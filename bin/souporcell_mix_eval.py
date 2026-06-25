#!/usr/bin/env python3
"""Score souporcell deconvolution against known sample-of-origin.

Built for the SOUPORCELL_ONLY validation experiments. The pipeline retags every cell barcode with a
``<sample>__`` prefix before merging samples into a mix, so for the artificial mixes the true donor of
each cell is recoverable directly from the barcode. This script reads each ``clusters.tsv``, recovers
the truth label from the prefix, and reports per-mix / per-K:

  * n cells, singlet / doublet / unassigned counts
  * confusion matrix  (true sample  x  souporcell cluster)  over singlets
  * Adjusted Rand Index (truth vs souporcell assignment) over singlets
  * best 1:1 cluster<->sample mapping accuracy
  * for solo (single-sample) negative controls: purity = fraction of singlets in the dominant cluster
    (high purity at k>=2 means souporcell did NOT spuriously split one individual)

Usage:
  souporcell_mix_eval.py --results <dir> [--outdir <dir>]
      globs <dir>/**/souporcell/<mix>/k<K>/clusters.tsv
  souporcell_mix_eval.py --clusters a/clusters.tsv b/clusters.tsv ...
"""
import argparse
import glob
import os
import re
import sys
from itertools import permutations

import numpy as np
import pandas as pd


def adjusted_rand_index(labels_true, labels_pred):
    """ARI without a sklearn dependency."""
    ct = pd.crosstab(pd.Series(labels_true), pd.Series(labels_pred)).to_numpy(dtype=float)
    n = ct.sum()
    if n < 2:
        return float("nan")
    a = ct.sum(axis=1)
    b = ct.sum(axis=0)
    comb = lambda x: x * (x - 1.0) / 2.0
    sum_ij = comb(ct).sum()
    sum_a = comb(a).sum()
    sum_b = comb(b).sum()
    expected = sum_a * sum_b / comb(n)
    maxi = 0.5 * (sum_a + sum_b)
    if maxi == expected:
        return 1.0
    return (sum_ij - expected) / (maxi - expected)


def best_mapping_accuracy(ct):
    """Max accuracy over all 1:1 cluster->sample assignments (rows=truth, cols=cluster)."""
    mat = ct.to_numpy(dtype=float)
    n_true, n_pred = mat.shape
    total = mat.sum()
    if total == 0:
        return float("nan")
    best = 0.0
    # iterate over assignments of clusters to (a subset of) true samples
    rows = range(n_true)
    for perm in permutations(range(n_pred), min(n_true, n_pred)):
        hit = sum(mat[r, perm[i]] for i, r in enumerate(rows) if i < len(perm))
        best = max(best, hit)
    return best / total


def parse_mix_k(path):
    m = re.search(r"souporcell/([^/]+)/k(\d+)/clusters\.tsv$", path)
    if m:
        return m.group(1), int(m.group(2))
    # fall back: <mix>/k<K>/clusters.tsv anywhere
    m = re.search(r"/([^/]+)/k(\d+)/clusters\.tsv$", path)
    return (m.group(1), int(m.group(2))) if m else (os.path.dirname(path), -1)


def evaluate(path):
    mix, k = parse_mix_k(path)
    df = pd.read_csv(path, sep="\t")
    df["truth"] = df["barcode"].str.split("__", n=1).str[0]
    n = len(df)
    status = df["status"].value_counts().to_dict()
    singlets = df[df["status"] == "singlet"].copy()
    singlets["assignment"] = singlets["assignment"].astype(str)

    row = {
        "mix": mix, "k": k, "n_cells": n,
        "n_singlet": int(status.get("singlet", 0)),
        "n_doublet": int(status.get("doublet", 0)),
        "n_unassigned": int(status.get("unassigned", 0)),
        "doublet_rate": round(status.get("doublet", 0) / n, 4) if n else float("nan"),
        "n_true_samples": singlets["truth"].nunique(),
    }
    ct = None
    if len(singlets):
        ct = pd.crosstab(singlets["truth"], singlets["assignment"])
        row["ari"] = round(adjusted_rand_index(singlets["truth"], singlets["assignment"]), 4)
        row["mapping_accuracy"] = round(best_mapping_accuracy(ct), 4)
        dom = singlets["assignment"].value_counts(normalize=True).iloc[0]
        row["dominant_cluster_purity"] = round(dom, 4)
    else:
        row["ari"] = row["mapping_accuracy"] = row["dominant_cluster_purity"] = float("nan")
    row["is_solo"] = row["n_true_samples"] == 1
    return row, ct


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", help="dir to glob <dir>/**/souporcell/<mix>/k<K>/clusters.tsv")
    ap.add_argument("--clusters", nargs="*", default=[], help="explicit clusters.tsv paths")
    ap.add_argument("--outdir", default=".", help="where to write summary + confusion matrices")
    args = ap.parse_args()

    paths = list(args.clusters)
    if args.results:
        paths += glob.glob(os.path.join(args.results, "**", "clusters.tsv"), recursive=True)
    paths = sorted(set(paths))
    if not paths:
        sys.exit("no clusters.tsv found")

    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    for p in paths:
        try:
            row, ct = evaluate(p)
        except Exception as e:  # noqa: BLE001
            print(f"WARN: failed on {p}: {e}", file=sys.stderr)
            continue
        rows.append(row)
        if ct is not None:
            ct.to_csv(os.path.join(args.outdir, f"confusion_{row['mix']}_k{row['k']}.csv"))
        print(f"== {row['mix']} k{row['k']} ==")
        if ct is not None:
            print(ct.to_string())
        # For true mixes (>=2 donors), success = high ARI of souporcell clusters vs donor identity.
        # For solo negative controls there is no true split: souporcell still returns K clusters
        # (it does NOT refuse to split — see ~50/50 purity on a pure sample), so purity/ARI here are
        # NOT a pass/fail signal. Real specificity needs the cluster_genotypes.vcf distinctness, which
        # is out of scope for this barcode-based scorer; we just report it.
        if row["is_solo"]:
            flag = "SOLO neg-ctrl (no true split; inspect cluster_genotypes.vcf for distinctness)"
        else:
            flag = "PASS" if row["ari"] >= 0.9 else "CHECK"
        print(f"  ARI={row['ari']} map_acc={row['mapping_accuracy']} "
              f"doublet_rate={row['doublet_rate']} solo={row['is_solo']} -> {flag}\n")

    summary = pd.DataFrame(rows).sort_values(["mix", "k"])
    out = os.path.join(args.outdir, "souporcell_mix_eval_summary.csv")
    summary.to_csv(out, index=False)
    print(f"wrote {out} ({len(summary)} rows)")


if __name__ == "__main__":
    main()

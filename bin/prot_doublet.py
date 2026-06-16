#!/usr/bin/env python3
"""
Protein/ADT doublet detection — ported from scripts/14_prot_doublet_detection.ipynb.

Marker-based doublet detection: cells co-expressing markers from incompatible lineages
(T: CD3/CD8, B: CD19/CD20, Monocyte: CD14) above a normalized-expression threshold are flagged
and removed. Falls back to keeping all cells if fewer than two lineage markers are present.
Checkpoint: prot_03_doublet_filtered.h5mu.

Usage: prot_doublet.py --in IN.h5mu --sample S [--threshold 2.5] --out OUT.h5mu
"""
import argparse
import os
import warnings

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
ad.settings.allow_write_nullable_strings = True


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--threshold", type=float, default=2.5)
    p.add_argument("--out", required=True)
    return p.parse_args()


def pick_markers(proteins):
    """One marker per lineage (T / B / Monocyte), preferring canonical names."""
    markers = []
    t = [p for p in proteins if "CD3" in p] or [p for p in proteins if "CD8" in p]
    b = [p for p in proteins if "CD19" in p or "CD20" in p]
    mono = [p for p in proteins if "CD14" in p]
    if t:
        markers.append(t[0])
    if b:
        markers.append(b[0])
    if mono:
        markers.append(mono[0])
    return markers


def main():
    args = parse_args()
    mdata = mu.read(args.inp)
    proteins = list(mdata["prot"].var_names)
    genes = pick_markers(proteins)
    print(f"Doublet markers: {genes}")

    n = mdata.shape[0]
    if len(genes) < 2:
        print("Fewer than 2 lineage markers; treating all cells as singlets")
        doublets = np.zeros(n, dtype=bool)
    else:
        expr = mdata["prot"][:, genes].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        temp = expr.T  # markers x cells
        thr = args.threshold
        pair_t_b = (temp[0] > thr) & (temp[1] > thr)
        doublets = pair_t_b
        if len(genes) >= 3:
            doublets = pair_t_b | ((temp[0] > thr) & (temp[2] > thr))

    mdata["prot"].obs["doublets_markers"] = pd.Series(
        doublets.astype(bool), index=mdata["prot"].obs_names
    ).astype(str)
    n_doublets = int(doublets.sum())
    mdata.update()

    mdata = mdata[mdata["prot"].obs["doublets_markers"] == "False"].copy()
    mdata.update()
    print(f"Protein doublets {args.sample}: removed {n_doublets} ({100*n_doublets/max(n,1):.1f}%); "
          f"{mdata.n_obs} singlets kept")

    mdata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

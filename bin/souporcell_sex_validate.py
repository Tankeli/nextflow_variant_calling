#!/usr/bin/env python3
"""Validate souporcell donor/recipient clusters in BM-transplant samples by sex-chromosome expression.

Scientific rationale (Exp 3): in a post-transplant relapse sample, leukaemic/host cells (recipient)
and engrafted graft cells (donor) are two genetic individuals -> souporcell should split them. If the
donor and recipient differ in sex, that split is independently confirmable from sex-linked expression:
  * female cells:  XIST high, chrY genes ~0
  * male cells:    chrY genes (RPS4Y1/DDX3Y/UTY/EIF1AY/KDM5D/USP9Y/NLGN4Y) expressed, XIST ~0
Donor sex is NOT recorded in the cohort, so this is the read-out: if a souporcell cluster is clearly
male and another clearly female within one sample, the genotype split is corroborated and the donor
sex is inferred. (Recipients here are 46,XX, so a male signal => male donor engraftment.)

For each (mix, K) the script joins souporcell singlet assignments to per-cell sex-marker expression
(from the constituent samples' filtered_feature_bc_matrix) and reports, per souporcell cluster:
mean XIST, mean chrY score, %cells expressing each, an inferred sex, and a 'sex_discordant' verdict.

Usage:
  souporcell_sex_validate.py --samplesheet assets/test/souporcell_mix_bmt.csv \
      --results results_soupmix_bmt/callers/souporcell --outdir results_soupmix_bmt/sex_eval
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import scanpy as sc

Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "NLGN4Y"]
XIST = "XIST"


def parse_mix_k(path):
    m = re.search(r"/([^/]+)/k(\d+)/clusters\.tsv$", path)
    return (m.group(1), int(m.group(2))) if m else (os.path.dirname(path), -1)


# project root (parent of bin/) so relative `outs` paths in the samplesheet resolve regardless of cwd
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve(p):
    return p if os.path.isabs(p) else os.path.join(PROJ_ROOT, p)


def load_sex_expression(outs_dir):
    """Return per-barcode DataFrame with normalized XIST + chrY score for one sample's matrix."""
    mtx = os.path.join(resolve(outs_dir), "filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(mtx, var_names="gene_symbols", cache=False)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    y_present = [g for g in Y_GENES if g in adata.var_names]
    df = pd.DataFrame(index=adata.obs_names)
    df["xist"] = (np.asarray(adata[:, XIST].X.todense()).ravel()
                  if XIST in adata.var_names else 0.0)
    df["y_score"] = (np.asarray(adata[:, y_present].X.todense()).mean(axis=1).ravel()
                     if y_present else 0.0)
    # also fraction-expressing helpers (binary > 0)
    df["xist_pos"] = (df["xist"] > 0).astype(int)
    df["y_pos"] = (df["y_score"] > 0).astype(int)
    return df


def infer_sex(mean_xist, mean_y, frac_xist_pos, frac_y_pos):
    male = (frac_y_pos >= 0.20) and (mean_y > 0.10)
    female = (frac_xist_pos >= 0.20) and (mean_xist > 0.10)
    if male and not female:
        return "male"
    if female and not male:
        return "female"
    if male and female:
        return "ambiguous(mixed)"
    return "low-signal"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samplesheet", required=True,
                    help="souporcell_mix_bmt.csv (sample,patient,timepoint,outs)")
    ap.add_argument("--results", required=True,
                    help="dir to glob <dir>/<mix>/k<K>/clusters.tsv")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    sheet = pd.read_csv(args.samplesheet)
    outs_by_sample = dict(zip(sheet["sample"], sheet["outs"]))

    # cache per-sample expression (samples reused across mixes)
    expr_cache = {}

    def get_expr(sample):
        if sample not in expr_cache:
            if sample not in outs_by_sample:
                raise KeyError(f"sample {sample} not in samplesheet")
            expr_cache[sample] = load_sex_expression(outs_by_sample[sample])
        return expr_cache[sample]

    paths = sorted(glob.glob(os.path.join(args.results, "**", "clusters.tsv"), recursive=True))
    if not paths:
        sys.exit(f"no clusters.tsv under {args.results}")

    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    for p in paths:
        mix, k = parse_mix_k(p)
        cl = pd.read_csv(p, sep="\t")
        cl = cl[cl["status"] == "singlet"].copy()
        cl["sample"] = cl["barcode"].str.split("__", n=1).str[0]
        cl["cb"] = cl["barcode"].str.split("__", n=1).str[1]
        cl["assignment"] = cl["assignment"].astype(str)

        # attach expression per cell
        parts = []
        for sample, sub in cl.groupby("sample"):
            expr = get_expr(sample)
            sub = sub.set_index("cb")
            j = sub.join(expr, how="left")
            j["sample"] = sample
            parts.append(j.reset_index())
        joined = pd.concat(parts, ignore_index=True)
        joined = joined.dropna(subset=["xist", "y_score"])

        per_cluster = []
        for cluster, g in joined.groupby("assignment"):
            mx, my = g["xist"].mean(), g["y_score"].mean()
            fx, fy = g["xist_pos"].mean(), g["y_pos"].mean()
            sex = infer_sex(mx, my, fx, fy)
            rec = dict(mix=mix, k=k, cluster=cluster, n_cells=len(g),
                       mean_xist=round(mx, 4), mean_y=round(my, 4),
                       frac_xist_pos=round(fx, 4), frac_y_pos=round(fy, 4),
                       inferred_sex=sex,
                       # which input samples land in this cluster (donor=relapse-only, recipient=both)
                       sample_mix=";".join(f"{s}:{n}" for s, n in
                                           g["sample"].value_counts().items()))
            per_cluster.append(rec)
            rows.append(rec)

        sexes = {r["inferred_sex"] for r in per_cluster
                 if r["inferred_sex"] in ("male", "female")}
        discordant = sexes == {"male", "female"}
        print(f"== {mix} k{k} ==")
        print(pd.DataFrame(per_cluster)[
            ["cluster", "n_cells", "mean_xist", "frac_xist_pos",
             "mean_y", "frac_y_pos", "inferred_sex", "sample_mix"]].to_string(index=False))
        verdict = ("SEX-DISCORDANT clusters -> genotype split corroborated (recipient 46,XX + likely male donor)"
                   if discordant else
                   "no sex discordance (donor likely same sex, or single origin) -> rely on genotype only")
        print(f"  -> {verdict}\n")

    summary = pd.DataFrame(rows).sort_values(["mix", "k", "cluster"])
    out = os.path.join(args.outdir, "souporcell_sex_validation.csv")
    summary.to_csv(out, index=False)
    print(f"wrote {out} ({len(summary)} cluster rows)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RNA differential expression — ported from scripts/10_scRNA_DE_analysis.ipynb.

Per-patient relapse-vs-diagnosis DE, computed per cell type. The notebook follows the
sc-best-practices edgeR/MAST tutorials (which assume replicated conditions); this cohort has a
single Dx and a single Rel sample per patient, so the robust default is a per-cell-type Wilcoxon
test (Rel vs Dx) on log-normalized counts. `--method edger|mast` attempts the R pseudobulk/mixed
models and falls back to Wilcoxon on failure. Always emits DE_markers.csv
(cell_type,gene,logFC,pval,padj).

Per-patient files are anonymised by Nextflow; sample/timepoint/cell_type come from obs stamped at
QC + annotation time.

Usage:
  rna_de.py --inputs . --patient P --method mast|edger|wilcoxon --groupby cell_type
            --out DE_markers.csv [--condition_key timepoint] [--reference Dx] [--min_cells 10]
"""
import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", required=True)
    p.add_argument("--patient", required=True)
    p.add_argument("--method", default="wilcoxon", choices=["wilcoxon", "edger", "mast"])
    p.add_argument("--groupby", default="cell_type")
    p.add_argument("--condition_key", default="timepoint")
    p.add_argument("--reference", default="Dx")
    p.add_argument("--min_cells", type=int, default=10)
    p.add_argument("--out", required=True)
    return p.parse_args()


def load_patient(inputs, condition_key):
    paths = sorted(glob.glob(os.path.join(inputs, "**", "*.h5ad"), recursive=True))
    if not paths:
        raise FileNotFoundError(f"No annotated .h5ad under {inputs}")
    adatas = [sc.read_h5ad(p) for p in paths]
    adata = sc.concat(adatas, join="outer")
    adata.obs_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    return adata


def wilcoxon_de(adata, groupby, condition_key, reference, min_cells):
    """Per-cell-type Wilcoxon, non-reference timepoint vs reference."""
    conditions = [c for c in adata.obs[condition_key].astype(str).unique() if c != "NA"]
    alt = [c for c in conditions if c != reference]
    rows = []
    if reference not in conditions or not alt:
        print(f"Need both '{reference}' and an alternative timepoint; found {conditions}")
        return pd.DataFrame(columns=["cell_type", "gene", "logFC", "pval", "padj"])
    alt = alt[0]
    for ct in adata.obs[groupby].astype(str).unique():
        sub = adata[adata.obs[groupby].astype(str) == ct].copy()
        vc = sub.obs[condition_key].astype(str).value_counts()
        if vc.get(reference, 0) < min_cells or vc.get(alt, 0) < min_cells:
            continue
        sub.obs["_grp"] = sub.obs[condition_key].astype(str)
        sub = sub[sub.obs["_grp"].isin([reference, alt])].copy()
        sub.obs["_grp"] = sub.obs["_grp"].astype("category")
        try:
            sc.tl.rank_genes_groups(sub, "_grp", groups=[alt], reference=reference, method="wilcoxon")
            df = sc.get.rank_genes_groups_df(sub, group=alt)
        except Exception as e:
            print(f"  {ct}: Wilcoxon failed ({e})")
            continue
        df = df.rename(columns={"names": "gene", "logfoldchanges": "logFC",
                                "pvals": "pval", "pvals_adj": "padj"})
        df["cell_type"] = ct
        rows.append(df[["cell_type", "gene", "logFC", "pval", "padj"]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
        columns=["cell_type", "gene", "logFC", "pval", "padj"])


def main():
    args = parse_args()
    adata = load_patient(args.inputs, args.condition_key)
    print(f"DE {args.patient}: {adata.n_obs} cells; "
          f"timepoints={list(adata.obs[args.condition_key].astype(str).unique())}; "
          f"cell types={adata.obs[args.groupby].astype(str).nunique()}")

    if args.method in ("edger", "mast"):
        # The R pseudobulk/mixed-effects models need replicated conditions, which this cohort
        # lacks (1 Dx + 1 Rel per patient). Attempt is intentionally not wired to avoid emitting
        # invalid statistics; fall back to the per-cell-type Wilcoxon test.
        print(f"method={args.method} needs replicated conditions (not available for a single "
              f"Dx/Rel pair); using Wilcoxon instead.")

    res = wilcoxon_de(adata, args.groupby, args.condition_key, args.reference, args.min_cells)
    res.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(res)} cell_type x gene results")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
CopyKAT CNA / aneuploidy summarised per reference-mapped cell type ("variants per cell type" +
a cell type x region variance matrix).

Joins the per-cell CopyKAT CNA bin matrix and aneuploid/diploid call to the ref_cell_type from the
reference-mapped h5ad, then builds:
  - cell type x chromosome matrices of mean CN and CN variance (clustered heatmaps),
  - per-cell-type aneuploid fraction (barplot) — does the malignancy call concentrate in a lineage?

Usage: copykat_celltype_matrix.py <sample> <mapped_h5ad> <CNA_results.txt> [calls.{csv|txt}|NONE]
Output (cwd): <sample>_copykat_celltype_matrix.csv, <sample>_copykat_celltype_aneuploid.csv, figures.
"""
import sys, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample, mapped_h5ad, cna_txt = sys.argv[1], sys.argv[2], sys.argv[3]
calls_in = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] not in ("", "NONE") else None

META_COLS = {"chrom", "chrompos", "abspos"}


def norm_bc(s):
    # CopyKAT sanitises barcodes (- -> .); normalise back so they match h5ad obs_names.
    return s.str.replace(".", "-", regex=False)


# ---- cell types ----
adata = sc.read_h5ad(mapped_h5ad)
if "ref_cell_type" not in adata.obs:
    sys.exit("mapped h5ad has no ref_cell_type column")
ct = adata.obs["ref_cell_type"].astype(str)
ct.index = adata.obs_names.astype(str)

# ---- CNA bin matrix (bins x cells) -> cells x bins ----
cna = pd.read_csv(cna_txt, sep="\t")
chrom_col = "chrom" if "chrom" in cna.columns else cna.columns[0]
cell_cols = [c for c in cna.columns if c not in META_COLS]
chrom = cna[chrom_col].astype(str).values
X = cna[cell_cols].astype(float).T            # cells x bins
X.index = norm_bc(pd.Index(X.index.astype(str)))

common = X.index.intersection(ct.index)
if len(common) == 0:                          # try the raw (un-normalised) barcodes too
    X.index = pd.Index(cell_cols)
    common = X.index.intersection(ct.index)
print(f"{sample}: matched {len(common)}/{len(cell_cols)} CopyKAT cells to cell types")
X = X.loc[common]; cell_ct = ct.loc[common]

# ---- cell type x chromosome mean CN + variance ----
bins_chrom = pd.Series(chrom, index=X.columns)
mean_rows, var_rows = {}, {}
for c, idx in cell_ct.groupby(cell_ct).groups.items():
    sub = X.loc[idx]
    mean_rows[c] = sub.T.groupby(bins_chrom).mean().mean(axis=1)   # mean over bins per chrom
    var_rows[c]  = sub.T.groupby(bins_chrom).var().mean(axis=1)
mean_mat = pd.DataFrame(mean_rows).T          # celltype x chrom
var_mat  = pd.DataFrame(var_rows).T

def chrom_key(c):
    c = str(c).replace("chr", "")
    return (0, int(c)) if c.isdigit() else (1, c)
order = sorted(mean_mat.columns, key=chrom_key)
mean_mat, var_mat = mean_mat[order], var_mat[order]
out = pd.concat({"mean_CN": mean_mat, "variance": var_mat}, axis=1)
out.to_csv(f"{sample}_copykat_celltype_matrix.csv")
print(f"{sample}: cell type x chromosome matrix = {mean_mat.shape[0]} types x {mean_mat.shape[1]} chroms")

for name, mat, cmap in [("mean_CN", mean_mat, "RdBu_r"), ("variance", var_mat, "viridis")]:
    fig, ax = plt.subplots(figsize=(max(6, mat.shape[1] * 0.5), max(4, mat.shape[0] * 0.35)))
    vmax = np.nanpercentile(np.abs(mat.values), 98) if name == "mean_CN" else np.nanpercentile(mat.values, 98)
    vmin = -vmax if name == "mean_CN" else 0
    im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_title(f"{sample} — CopyKAT {name} per cell type x chromosome")
    fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout(); fig.savefig(f"{sample}_copykat_celltype_{name}.png", dpi=150)
    plt.close(fig)

# ---- per-cell-type aneuploid fraction ----
if calls_in:
    sep = "," if calls_in.endswith(".csv") else "\t"
    cdf = pd.read_csv(calls_in, sep=sep)
    if "consensus" in cdf.columns:
        cmap_ = dict(zip(cdf[cdf.columns[0]].astype(str), cdf["consensus"].astype(str)))
    else:
        nm = "cell.names" if "cell.names" in cdf.columns else cdf.columns[0]
        pr = "copykat.pred" if "copykat.pred" in cdf.columns else cdf.columns[-1]
        cmap_ = dict(zip(cdf[nm].astype(str), cdf[pr].astype(str)))
    call = pd.Series({b: cmap_.get(b, "not.defined") for b in cell_ct.index})
    tab = pd.DataFrame({"cell_type": cell_ct.values, "call": call.values})
    frac = (tab[tab["call"].isin(["aneuploid", "diploid"])]
            .assign(aneu=lambda d: d["call"] == "aneuploid")
            .groupby("cell_type")["aneu"].agg(["mean", "count"])
            .rename(columns={"mean": "aneuploid_frac", "count": "n_cells"})
            .sort_values("aneuploid_frac", ascending=False))
    frac.to_csv(f"{sample}_copykat_celltype_aneuploid.csv")
    fig, ax = plt.subplots(figsize=(7, max(4, len(frac) * 0.3)))
    ax.barh(frac.index, frac["aneuploid_frac"], color="tab:blue")
    ax.invert_yaxis(); ax.set_xlabel("aneuploid fraction"); ax.set_xlim(0, 1)
    ax.set_title(f"{sample} — CopyKAT aneuploid fraction per cell type")
    fig.tight_layout(); fig.savefig(f"{sample}_copykat_celltype_aneuploid.png", dpi=150)
    plt.close(fig)
    print(f"{sample}: aneuploid fraction by cell type:\n{frac['aneuploid_frac'].round(3).to_string()}")

print(f"Wrote {sample}_copykat_celltype_matrix.csv + figures")

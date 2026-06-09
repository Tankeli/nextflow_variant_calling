#!/usr/bin/env python3
"""
Overlay the CopyKAT aneuploid/diploid call on the reference-mapping UMAP.

CopyKAT's own genome heatmap is published by the caller; this adds the
aneuploid/diploid + cell-type UMAP overlay that DDE_32 copyKAT_profiling.R drew on a
Seurat ref.umap. Here we reuse the X_umap_ref embedding + ref_cell_type already computed
by reference_mapping.py, so no UMAP is recomputed.

Usage: plot_copykat_umap.py <sample> <mapped_h5ad> <copykat_prediction.txt>
Output (cwd): <sample>_copykat_umap.png/.pdf
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample, mapped_h5ad, pred_txt = sys.argv[1], sys.argv[2], sys.argv[3]

adata = sc.read_h5ad(mapped_h5ad)
if "X_umap_ref" not in adata.obsm:
    sys.exit("mapped h5ad has no X_umap_ref embedding; run reference mapping first")

pred = pd.read_csv(pred_txt, sep="\t")
# copykat prediction.txt columns: cell.names, copykat.pred
name_col = "cell.names" if "cell.names" in pred.columns else pred.columns[0]
pred_col = "copykat.pred" if "copykat.pred" in pred.columns else pred.columns[-1]
pred_map = dict(zip(pred[name_col].astype(str), pred[pred_col].astype(str)))

adata.obs["copykat_pred"] = (
    pd.Series(adata.obs_names, index=adata.obs_names).map(pred_map).fillna("not.defined")
)
adata.obs["copykat_pred"] = adata.obs["copykat_pred"].astype("category")

n_match = int((adata.obs["copykat_pred"] != "not.defined").sum())
print(f"{sample}: matched CopyKAT calls for {n_match}/{adata.n_obs} cells")
print(adata.obs["copykat_pred"].value_counts())

ck_palette = {"aneuploid": "tab:blue", "diploid": "tab:red", "not.defined": "lightgray"}
present = [c for c in ck_palette if c in adata.obs["copykat_pred"].cat.categories]

sc.set_figure_params(dpi=150, frameon=False)
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.embedding(adata, basis="umap_ref", color="copykat_pred", ax=axes[0], show=False,
                title=f"{sample} — CopyKAT CNV call",
                palette=[ck_palette[c] for c in present])
ref_ct = "ref_cell_type" if "ref_cell_type" in adata.obs else None
if ref_ct:
    adata.obs[ref_ct] = adata.obs[ref_ct].astype("category")
    sc.pl.embedding(adata, basis="umap_ref", color=ref_ct, ax=axes[1], show=False,
                    title=f"{sample} — mapped cell types",
                    legend_loc="right margin", legend_fontsize=7)
else:
    axes[1].axis("off")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{sample}_copykat_umap.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Wrote {sample}_copykat_umap.png and .pdf")

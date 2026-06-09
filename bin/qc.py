#!/usr/bin/env python3
"""
Per-sample scanpy QC for CITE-seq (GEX + Antibody Capture).
Ported/parameterized from DDE_23 scripts/stage1_qc.py.

Reads a cellranger filtered matrix dir (GEX + Antibody Capture), computes QC metrics,
runs Scrublet doublet detection, applies filters, and writes a filtered AnnData plus a
per-cell metrics table (all cells, with a pass/fail flag).

Usage: qc.py <matrix_dir> <sample> <min_genes> <min_umi> <max_mito_pct> <expected_doublet_rate>
Outputs (cwd): <sample>_qc.h5ad  (filtered),  <sample>_qc_metrics.csv  (all cells),
               <sample>_qc_panel.png/.pdf  (diagnostic 3x3 panel over all cells)
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import scrublet as scr
from scipy.sparse import issparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matrix_dir, sample = sys.argv[1], sys.argv[2]
min_genes      = int(sys.argv[3])
min_umi        = int(sys.argv[4])
max_mito_pct   = float(sys.argv[5])
exp_doublet    = float(sys.argv[6])

adata = sc.read_10x_mtx(matrix_dir, gex_only=False)
adata.var_names_make_unique()
adata.obs["sample_id"] = sample
print(f"Loaded {sample}: {adata.n_obs} cells x {adata.n_vars} features")

gex_mask = adata.var["feature_types"] == "Gene Expression"
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

def colsum(mask):
    v = adata[:, mask].X.sum(axis=1)
    return np.asarray(v).flatten()

adata.obs["gex_counts"]     = colsum(gex_mask)
adata.obs["protein_counts"] = colsum(~gex_mask)
ng = (adata[:, gex_mask].X > 0).sum(axis=1)
adata.obs["n_genes"] = np.asarray(ng).flatten()

# Scrublet on GEX counts
try:
    scrub = scr.Scrublet(adata[:, gex_mask].X, expected_doublet_rate=exp_doublet)
    scores, predicted = scrub.scrub_doublets(min_counts=2, min_cells=3,
                                             min_gene_variability_pctl=85, n_prin_comps=30)
except Exception as e:  # scrublet can fail on tiny/odd inputs
    print(f"Scrublet failed ({e}); marking no doublets")
    scores = np.zeros(adata.n_obs); predicted = np.zeros(adata.n_obs, dtype=bool)
adata.obs["doublet_score"]     = scores
adata.obs["predicted_doublet"] = predicted

# Per-cell pass/fail (recorded for ALL cells, then filter)
adata.obs["pass_qc"] = (
    (adata.obs["n_genes"] >= min_genes)
    & (adata.obs["gex_counts"] >= min_umi)
    & (adata.obs["pct_counts_mt"] <= max_mito_pct)
    & (~adata.obs["predicted_doublet"].astype(bool))
)

metrics_cols = ["sample_id", "n_genes", "gex_counts", "protein_counts",
                "pct_counts_mt", "doublet_score", "predicted_doublet", "pass_qc"]
adata.obs[metrics_cols].to_csv(f"{sample}_qc_metrics.csv")


def plot_qc_panel(ad, sample_id):
    """Diagnostic 3x3 QC panel over ALL cells (pre-filter), with threshold lines.
    Ported from DDE_23 scripts/stage1_qc.py:plot_qc_metrics."""
    obs = ad.obs
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle(f"QC metrics — {sample_id}", fontsize=16, y=0.995)

    axes[0, 0].hist(obs["n_genes"], bins=50, edgecolor="black")
    axes[0, 0].axvline(min_genes, color="red", linestyle="--", label="threshold")
    axes[0, 0].set(xlabel="Genes per cell", ylabel="Cells")
    axes[0, 0].legend()

    axes[0, 1].hist(np.log10(obs["gex_counts"] + 1), bins=50, edgecolor="black")
    axes[0, 1].axvline(np.log10(min_umi + 1), color="red", linestyle="--")
    axes[0, 1].set(xlabel="log10(GEX UMI + 1)", ylabel="Cells")

    axes[0, 2].hist(np.log10(obs["protein_counts"] + 1), bins=50, edgecolor="black")
    axes[0, 2].set(xlabel="log10(Protein UMI + 1)", ylabel="Cells")

    axes[1, 0].hist(obs["pct_counts_mt"], bins=50, edgecolor="black")
    axes[1, 0].axvline(max_mito_pct, color="red", linestyle="--")
    axes[1, 0].set(xlabel="% Mitochondrial", ylabel="Cells")

    axes[1, 1].hist(obs["doublet_score"], bins=50, edgecolor="black")
    axes[1, 1].set(xlabel="Scrublet doublet score", ylabel="Cells")

    sccat = axes[1, 2].scatter(np.log10(obs["gex_counts"] + 1), obs["n_genes"],
                               c=obs["pct_counts_mt"], s=1, alpha=0.5, cmap="viridis")
    axes[1, 2].set(xlabel="log10(GEX UMI + 1)", ylabel="Genes per cell")
    plt.colorbar(sccat, ax=axes[1, 2], label="% MT")

    axes[2, 0].scatter(np.log10(obs["gex_counts"] + 1),
                       np.log10(obs["protein_counts"] + 1), s=1, alpha=0.5)
    axes[2, 0].set(xlabel="log10(GEX UMI + 1)", ylabel="log10(Protein UMI + 1)")

    pass_mask = obs["pass_qc"].astype(bool).values
    axes[2, 1].scatter(obs["pct_counts_mt"][~pass_mask], obs["n_genes"][~pass_mask],
                       s=1, alpha=0.4, color="lightgray", label="fail")
    axes[2, 1].scatter(obs["pct_counts_mt"][pass_mask], obs["n_genes"][pass_mask],
                       s=1, alpha=0.5, color="tab:blue", label="pass")
    axes[2, 1].set(xlabel="% Mitochondrial", ylabel="Genes per cell")
    axes[2, 1].legend(markerscale=6)

    n_pass = int(pass_mask.sum())
    summary = (
        f"Summary — {sample_id}\n\n"
        f"Total cells:        {ad.n_obs:,}\n"
        f"Pass QC:            {n_pass:,} ({100*n_pass/max(ad.n_obs,1):.1f}%)\n"
        f"Median genes:       {np.median(obs['n_genes']):.0f}\n"
        f"Median GEX UMI:     {np.median(obs['gex_counts']):.0f}\n"
        f"Median Protein UMI: {np.median(obs['protein_counts']):.0f}\n"
        f"Median %% MT:        {np.median(obs['pct_counts_mt']):.2f}%%\n"
        f"Predicted doublets: {int(obs['predicted_doublet'].sum())} "
        f"({100*obs['predicted_doublet'].sum()/max(ad.n_obs,1):.1f}%%)\n\n"
        f"Thresholds: min_genes={min_genes}, min_umi={min_umi},\n"
        f"            max_mt={max_mito_pct}%%"
    )
    axes[2, 2].axis("off")
    axes[2, 2].text(0.0, 0.5, summary, fontsize=11, family="monospace",
                    verticalalignment="center")

    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{sample_id}_qc_panel.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {sample_id}_qc_panel.png and .pdf")


plot_qc_panel(adata, sample)

n0 = adata.n_obs
adata_f = adata[adata.obs["pass_qc"].values, :].copy()
print(f"QC {sample}: kept {adata_f.n_obs}/{n0} cells "
      f"({100*(n0-adata_f.n_obs)/max(n0,1):.1f}% removed)")

adata_f.layers["counts"] = adata_f.X.copy()
adata_f.write_h5ad(f"{sample}_qc.h5ad")
print(f"Wrote {sample}_qc.h5ad and {sample}_qc_metrics.csv")

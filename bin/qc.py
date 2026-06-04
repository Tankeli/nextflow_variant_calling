#!/usr/bin/env python3
"""
Per-sample scanpy QC for CITE-seq (GEX + Antibody Capture).
Ported/parameterized from DDE_23 scripts/stage1_qc.py.

Reads a cellranger filtered matrix dir (GEX + Antibody Capture), computes QC metrics,
runs Scrublet doublet detection, applies filters, and writes a filtered AnnData plus a
per-cell metrics table (all cells, with a pass/fail flag).

Usage: qc.py <matrix_dir> <sample> <min_genes> <min_umi> <max_mito_pct> <expected_doublet_rate>
Outputs (cwd): <sample>_qc.h5ad  (filtered),  <sample>_qc_metrics.csv  (all cells)
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import scanpy as sc
import scrublet as scr
from scipy.sparse import issparse

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

n0 = adata.n_obs
adata_f = adata[adata.obs["pass_qc"].values, :].copy()
print(f"QC {sample}: kept {adata_f.n_obs}/{n0} cells "
      f"({100*(n0-adata_f.n_obs)/max(n0,1):.1f}% removed)")

adata_f.layers["counts"] = adata_f.X.copy()
adata_f.write_h5ad(f"{sample}_qc.h5ad")
print(f"Wrote {sample}_qc.h5ad and {sample}_qc_metrics.csv")

#!/usr/bin/env python3
"""
Per-cell leukaemic-stem-cell scoring (weighted pLSC6 + LSC17).

Weighted signature scores = sum_g (lognorm_expr[g] * weight[g]), exactly as DDE_32
LSC_scoring.R. Runs on the reference-mapped h5ad (lognorm GEX), so it sits in the annotation
branch and feeds the Phase-0 master table / headline pLSC6 figures.

Usage: lsc_scoring.py <sample> <mapped_h5ad>
Output (cwd): <sample>_lsc.csv  (barcode, sample_id, pLSC6_score, LSC17_score)
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

sample, mapped_h5ad = sys.argv[1], sys.argv[2]

# Weighted signatures (Ng et al. pLSC6; Ng et al. LSC17). Gene-symbol aliases listed so the
# score is robust to whichever symbol the atlas/cellranger reference used.
PLSC6 = {"DNMT3B": 0.189, "ADGRG1": 0.054, "CD34": 0.0171,
         "SOCS2": 0.141, "SPINK2": 0.109, "FAM30A": 0.0516}
LSC17 = {"DNMT3B": 0.0874, "ZBTB46": -0.0347, "NYNRIN": 0.00865, "ARHGAP22": -0.0138,
         "LAPTM4B": 0.00582, "MMRN1": 0.0258, "DPYSL3": 0.0284, "KIAA0125": 0.0196,
         "CDK6": -0.0704, "CPXM1": -0.0258, "SOCS2": 0.0271, "SMIM24": -0.0226,
         "EMP1": 0.0146, "NGFRAP1": 0.0465, "CD34": 0.0338, "AKR1C3": -0.0402,
         "GPR56": 0.0501}
ALIASES = {"ADGRG1": ["ADGRG1", "GPR56"], "GPR56": ["GPR56", "ADGRG1"],
           "FAM30A": ["FAM30A", "KIAA0125"], "KIAA0125": ["KIAA0125", "FAM30A"],
           "NGFRAP1": ["NGFRAP1", "BEX3"], "SMIM24": ["SMIM24", "C19orf77"]}

adata = sc.read_h5ad(mapped_h5ad)
var_index = {g: i for i, g in enumerate(adata.var_names)}
X = adata.X.tocsr() if issparse(adata.X) else np.asarray(adata.X)


def gene_col(gene):
    for alias in ALIASES.get(gene, [gene]):
        if alias in var_index:
            col = X[:, var_index[alias]]
            return np.asarray(col.todense()).ravel() if issparse(col) else np.asarray(col).ravel()
    return None


def weighted_score(weights, min_genes, label):
    score = np.zeros(adata.n_obs)
    found = 0
    for gene, w in weights.items():
        v = gene_col(gene)
        if v is not None:
            score += v * w
            found += 1
    print(f"  {label}: {found}/{len(weights)} genes available")
    if found < min_genes:
        print(f"  {label}: too few genes ({found}<{min_genes}); emitting NaN")
        return np.full(adata.n_obs, np.nan)
    return score


print(f"LSC scoring {sample}: {adata.n_obs} cells x {adata.n_vars} genes")
out = pd.DataFrame({
    "sample_id":    adata.obs["sample_id"].values if "sample_id" in adata.obs else sample,
    "pLSC6_score":  weighted_score(PLSC6, 3, "pLSC6"),
    "LSC17_score":  weighted_score(LSC17, 10, "LSC17"),
}, index=adata.obs_names)
out.index.name = "barcode"
out.to_csv(f"{sample}_lsc.csv")
print(f"Wrote {sample}_lsc.csv "
      f"(pLSC6 mean={np.nanmean(out['pLSC6_score']):.3f}, LSC17 mean={np.nanmean(out['LSC17_score']):.3f})")

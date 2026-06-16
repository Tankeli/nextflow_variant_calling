#!/usr/bin/env python3
"""
RNA normalization — ported from scripts/02_scRNA_normalisation.ipynb.

Adds normalization layers to the QC'd AnnData:
  - log1p_norm                 (shifted-log of total-count normalization)
  - scran_normalization        (scran size-factor normalization + log1p)  [primary]
  - analytic_pearson_residuals (analytic Pearson residuals)
size_factors are stored in .obs.

Usage: rna_normalize.py --in IN.h5ad --sample S --out OUT.h5ad
"""
import argparse
import logging
import os
import warnings

import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix, issparse

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def setup_rpy2():
    import anndata2ri
    import rpy2.rinterface_lib.callbacks as rcb
    from rpy2.rinterface_lib import conversion

    rcb.logger.setLevel(logging.ERROR)
    _orig = conversion._get_cdata

    def _patched(obj):
        if isinstance(obj, np.str_):
            obj = str(obj)
        return _orig(obj)

    conversion._get_cdata = _patched
    anndata2ri.activate()


def scran_size_factors(adata):
    """Preliminary clustering + scran computeSumFactors. Mirrors notebook cells 8-12."""
    import rpy2.robjects as ro

    adata_pp = adata.copy()
    sc.pp.normalize_total(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp, n_comps=15)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="groups")

    data_mat = adata_pp.X.T
    if issparse(data_mat):
        data_mat = data_mat.tocoo() if data_mat.nnz > 2**31 - 1 else data_mat.tocsc()
    ro.globalenv["data_mat"] = data_mat
    ro.globalenv["input_groups"] = adata_pp.obs["groups"]
    del adata_pp

    size_factors = ro.r(
        """
        suppressPackageStartupMessages({ library(scran); library(BiocParallel) })
        sizeFactors(
            computeSumFactors(
                SingleCellExperiment(list(counts=data_mat)),
                clusters = input_groups, min.mean = 0.1, BPPARAM = MulticoreParam()
            )
        )
        """
    )
    return np.asarray(size_factors)


def main():
    args = parse_args()
    adata = sc.read(args.inp)

    # Shifted-log of total-count normalization
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    # scran size-factor normalization (primary)
    setup_rpy2()
    adata.obs["size_factors"] = scran_size_factors(adata)
    scran = adata.X / adata.obs["size_factors"].values[:, None]
    scran_csr = csr_matrix(scran)
    scran_csr.data = np.log1p(scran_csr.data)
    adata.layers["scran_normalization"] = scran_csr

    # Analytic Pearson residuals
    analytic = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic["X"])

    adata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

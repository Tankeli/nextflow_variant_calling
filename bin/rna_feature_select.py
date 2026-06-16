#!/usr/bin/env python3
"""
RNA feature selection — ported from scripts/03_scRNA_feature_selection.ipynb.

Deviance-based highly variable gene selection (scry::devianceFeatureSelection over raw counts in
.X). Flags the top --n_top_genes as `highly_deviant` and stores `binomial_deviance` in .var; also
runs scanpy's standard HVG on the scran layer for comparison.

Usage: rna_feature_select.py --in IN.h5ad --sample S --n_top_genes 4000 --out OUT.h5ad
"""
import argparse
import logging
import os
import warnings

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--n_top_genes", type=int, default=4000)
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


def main():
    args = parse_args()
    setup_rpy2()
    import rpy2.robjects as ro

    adata = sc.read(args.inp)

    ro.globalenv["adata"] = adata
    ro.r('suppressPackageStartupMessages(library(scry))')
    ro.r('sce = devianceFeatureSelection(adata, assay="X")')
    binomial_deviance = np.asarray(ro.r("rowData(sce)$binomial_deviance")).T

    idx = binomial_deviance.argsort()[-args.n_top_genes:]
    mask = np.zeros(adata.var_names.shape, dtype=bool)
    mask[idx] = True
    adata.var["highly_deviant"] = mask
    adata.var["binomial_deviance"] = binomial_deviance

    # scanpy HVG on the scran layer (kept for comparison, as in the notebook)
    sc.pp.highly_variable_genes(adata, layer="scran_normalization")

    adata.write(args.out)
    print(f"Wrote {args.out}: {int(mask.sum())} highly deviant genes")


if __name__ == "__main__":
    main()

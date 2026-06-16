#!/usr/bin/env python3
"""
RNA QC — ported from scripts/01_scRNA_quality_control.ipynb.

Reads a Cell Ranger filtered 10x h5 (GEX; ADT ignored here), computes QC metrics, removes MAD
outliers, optionally corrects ambient RNA with SoupX (needs the raw/table-of-droplets h5),
filters low-count genes, and optionally scores doublets with scDblFinder. Doublets are annotated,
not removed (parity with the notebook).

Usage:
  rna_qc.py --filtered_h5 F --sample S --out OUT.h5ad --metrics M.csv
            [--raw_h5 R] [--soupx] [--scdblfinder]
            [--nmads_counts 5] [--nmads_mt 3] [--max_mito_pct 8] [--min_cells 20]
"""
import argparse
import logging
import os
import warnings

import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--filtered_h5", required=True)
    p.add_argument("--raw_h5", default=None, help="raw_feature_bc_matrix.h5 for SoupX")
    p.add_argument("--sample", required=True)
    p.add_argument("--patient", default=None)
    p.add_argument("--timepoint", default=None)
    p.add_argument("--batch", default=None)
    p.add_argument("--out", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--soupx", action="store_true")
    p.add_argument("--scdblfinder", action="store_true")
    p.add_argument("--nmads_counts", type=int, default=5)
    p.add_argument("--nmads_mt", type=int, default=3)
    p.add_argument("--max_mito_pct", type=float, default=8.0)
    p.add_argument("--min_cells", type=int, default=20)
    return p.parse_args()


def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )


def setup_rpy2():
    """Activate anndata2ri conversion + the numpy.str_ workaround used in the notebooks."""
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


def run_soupx(adata, raw_h5):
    """Ambient-RNA correction. Mirrors notebook cells 14-22."""
    import rpy2.robjects as ro

    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups")
    soupx_groups = adata_pp.obs["soupx_groups"]
    del adata_pp

    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    adata_raw = sc.read_10x_h5(filename=raw_h5)
    adata_raw.var_names_make_unique()
    data_tod = adata_raw.X.T
    del adata_raw

    ro.globalenv["data"] = data
    ro.globalenv["data_tod"] = data_tod
    ro.globalenv["genes"] = ro.StrVector(list(genes))
    ro.globalenv["cells"] = ro.StrVector(list(cells))
    ro.globalenv["soupx_groups"] = ro.StrVector(list(soupx_groups.astype(str)))

    out = ro.r(
        """
        suppressPackageStartupMessages(library(SoupX))
        rownames(data) = genes
        colnames(data) = cells
        data <- as(data, "sparseMatrix")
        data_tod <- as(data_tod, "sparseMatrix")
        sc = SoupChannel(data_tod, data, calcSoupProfile = FALSE)
        soupProf = data.frame(row.names = rownames(data),
                              est = rowSums(data)/sum(data), counts = rowSums(data))
        sc = setSoupProfile(sc, soupProf)
        sc = setClusters(sc, soupx_groups)
        sc = autoEstCont(sc, doPlot=FALSE)
        out = adjustCounts(sc, roundToInt = TRUE)
        out
        """
    )
    adata.layers["soupX_counts"] = out.T
    adata.X = adata.layers["soupX_counts"]
    return adata


def run_scdblfinder(adata):
    """Doublet scoring. Mirrors notebook cells 24-27."""
    import rpy2.robjects as ro

    ro.globalenv["data_mat"] = adata.X.T
    res = ro.r(
        """
        suppressPackageStartupMessages({
            library(scater); library(scDblFinder); library(BiocParallel)
        })
        set.seed(123)
        sce = scDblFinder(SingleCellExperiment(list(counts=data_mat)))
        list(score=sce$scDblFinder.score, class=as.character(sce$scDblFinder.class))
        """
    )
    adata.obs["scDblFinder_score"] = np.asarray(res.rx2("score"))
    adata.obs["scDblFinder_class"] = list(res.rx2("class"))
    return adata


def main():
    args = parse_args()

    # DDE_33 boundary: Cell Ranger here publishes the 10x mtx *directory*
    # (filtered_feature_bc_matrix/) rather than the combined filtered .h5. Read either — a
    # directory via read_10x_mtx (gex_only=False keeps the feature_types column so the CITE-seq
    # GEX subset below still applies), a file via read_10x_h5.
    if os.path.isdir(args.filtered_h5):
        adata = sc.read_10x_mtx(args.filtered_h5, gex_only=False)
    else:
        adata = sc.read_10x_h5(filename=args.filtered_h5)
    adata.var_names_make_unique()
    # Restrict to Gene Expression if the h5 carries multiple feature types (CITE-seq).
    if "feature_types" in adata.var.columns:
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    # Stamp sample metadata onto .obs so cohort/per-patient stages (integration, DE,
    # composition) stay self-describing after the per-file names are anonymized by Nextflow.
    def _clean(v):
        return None if v in (None, "", "null", "NA") else v

    adata.obs["sample_id"] = args.sample
    adata.obs["patient"] = _clean(args.patient) or "NA"
    adata.obs["timepoint"] = _clean(args.timepoint) or "NA"
    adata.obs["batch"] = _clean(args.batch) or args.sample  # fall back to per-sample batch

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[ABDEGMQZ]\d*(?!\w)")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", args.nmads_counts)
        | is_outlier(adata, "log1p_n_genes_by_counts", args.nmads_counts)
        | is_outlier(adata, "pct_counts_in_top_20_genes", args.nmads_counts)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", args.nmads_mt) | (
        adata.obs["pct_counts_mt"] > args.max_mito_pct
    )
    n0 = adata.n_obs
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"QC {args.sample}: {adata.n_obs}/{n0} cells kept after outlier filtering")

    adata.layers["counts"] = adata.X.copy()

    need_r = args.soupx or args.scdblfinder
    if need_r:
        setup_rpy2()

    if args.soupx:
        if not args.raw_h5:
            print("--soupx requested but no --raw_h5; skipping ambient correction")
        else:
            adata = run_soupx(adata, args.raw_h5)

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"QC {args.sample}: {adata.n_vars} genes after min_cells={args.min_cells}")

    if args.scdblfinder:
        adata = run_scdblfinder(adata)

    obs_cols = [c for c in adata.obs.columns if c.startswith(("total_", "n_genes", "pct_counts", "scDblFinder"))]
    adata.obs.assign(sample_id=args.sample)[["sample_id"] + obs_cols].to_csv(args.metrics)
    adata.write(args.out)
    print(f"Wrote {args.out} and {args.metrics}")


if __name__ == "__main__":
    main()

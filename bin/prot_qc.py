#!/usr/bin/env python3
"""
Protein/ADT QC — ported from scripts/12_prot_quality_control.ipynb.

Loads the combined Cell Ranger filtered h5 as a MuData (rna + prot modalities), computes QC
metrics on the protein modality, filters extreme-count cells (potential doublets) and MAD
outliers, and stamps sample metadata onto .obs so cohort batch correction stays self-describing.
Checkpoint: prot_01_quality_control.h5mu.

Usage:
  prot_qc.py --filtered_h5 F --sample S --out OUT.h5mu --metrics M.csv
             [--patient P --timepoint T --batch B]
             [--nmads 5] [--max_total_counts 100000]
"""
import argparse
import os
import warnings

import anndata as ad
import muon as mu
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
ad.settings.allow_write_nullable_strings = True


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--filtered_h5", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--patient", default=None)
    p.add_argument("--timepoint", default=None)
    p.add_argument("--batch", default=None)
    p.add_argument("--nmads", type=int, default=5)
    p.add_argument("--max_total_counts", type=int, default=100000)
    p.add_argument("--out", required=True)
    p.add_argument("--metrics", required=True)
    return p.parse_args()


def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )


def main():
    args = parse_args()

    def _clean(v):
        return None if v in (None, "", "null", "NA") else v

    # DDE_33 boundary: Cell Ranger here publishes the 10x mtx *directory*
    # (filtered_feature_bc_matrix/) rather than the combined filtered .h5. muon reads either —
    # a directory via read_10x_mtx, a file via read_10x_h5 — and splits GEX/ADT into rna/prot.
    if os.path.isdir(args.filtered_h5):
        mdata = mu.read_10x_mtx(args.filtered_h5)
    else:
        mdata = mu.read_10x_h5(args.filtered_h5)
    print(f"Loaded {args.sample}: modalities {list(mdata.mod.keys())}")
    if "prot" not in mdata.mod:
        raise ValueError(f"No 'prot' (Antibody Capture) modality in {args.filtered_h5}")

    for mod in mdata.mod:
        sc.pp.calculate_qc_metrics(mdata[mod], inplace=True, percent_top=None)

    # Stamp metadata on each modality + the MuData-level obs.
    for obj in [mdata["prot"], mdata["rna"] if "rna" in mdata.mod else None]:
        if obj is None:
            continue
        obj.obs["sample_id"] = args.sample
        obj.obs["patient"] = _clean(args.patient) or "NA"
        obj.obs["timepoint"] = _clean(args.timepoint) or "NA"
        obj.obs["batch"] = _clean(args.batch) or args.sample

    n0 = mdata.n_obs
    # Extreme total-count cells (likely multiplets).
    high = (mdata["prot"].obs.total_counts > args.max_total_counts).sum()
    if high > 0:
        sc.pp.filter_cells(mdata["prot"], max_counts=args.max_total_counts)
        mdata.update()

    # MAD outliers on protein library size / complexity.
    mdata["prot"].obs["outlier"] = is_outlier(mdata["prot"], "log1p_total_counts", args.nmads) | is_outlier(
        mdata["prot"], "log1p_n_genes_by_counts", args.nmads
    )
    mdata = mdata[~mdata["prot"].obs["outlier"]].copy()
    mdata.update()
    print(f"Protein QC {args.sample}: {mdata.n_obs}/{n0} cells kept "
          f"({high} extreme-count, {n0 - mdata.n_obs} total removed)")

    cols = [c for c in ["sample_id", "patient", "timepoint", "batch",
                        "total_counts", "n_genes_by_counts", "log1p_total_counts"]
            if c in mdata["prot"].obs.columns]
    mdata["prot"].obs[cols].to_csv(args.metrics)
    mdata.write(args.out)
    print(f"Wrote {args.out} and {args.metrics}")


if __name__ == "__main__":
    main()

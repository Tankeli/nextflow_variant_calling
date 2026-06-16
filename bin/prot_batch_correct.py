#!/usr/bin/env python3
"""
Protein/ADT batch correction — ported from scripts/16_prot_batch_correction.ipynb.

Cohort step: concatenates the per-sample dimred protein modalities, recomputes PCA on the merged
object, and runs Harmony across batches. Neighbors + UMAP are recomputed on the corrected
representation (X_pcahm when Harmony ran, else X_pca). With a single batch it passes through with
no correction. Checkpoint: prot_05_batch_corrected.h5mu.

Per-sample files are anonymised by Nextflow, so batch identity comes from the obs `batch` column
stamped at QC time.

Usage: prot_batch_correct.py --inputs . --batch_key batch [--n_pcs 30] --out OUT.h5mu
"""
import argparse
import glob
import os
import warnings

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
ad.settings.allow_write_nullable_strings = True


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", required=True, help="dir containing the per-sample dimred h5mus")
    p.add_argument("--batch_key", default="batch")
    p.add_argument("--n_pcs", type=int, default=30)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    paths = sorted(glob.glob(os.path.join(args.inputs, "**", "*.h5mu"), recursive=True))
    if not paths:
        raise FileNotFoundError(f"No .h5mu files under {args.inputs}")
    print(f"Batch-correcting {len(paths)} protein samples")

    prot_list = [mu.read(p)["prot"] for p in paths]
    adata = sc.concat(prot_list, join="outer")
    adata.obs_names_make_unique()
    if args.batch_key not in adata.obs:
        adata.obs[args.batch_key] = "batch1"
    n_batches = adata.obs[args.batch_key].nunique()
    print(f"Merged protein: {adata.n_obs} cells x {adata.n_vars} proteins; batches={n_batches}")

    n_pcs = min(args.n_pcs, adata.n_vars - 1)
    sc.pp.pca(adata, svd_solver="arpack")
    use_rep = "X_pca"
    if n_batches > 1:
        try:
            import harmonypy as hm
            meta = pd.DataFrame({args.batch_key: adata.obs[args.batch_key].astype(str).values})
            ho = hm.run_harmony(adata.obsm["X_pca"], meta, args.batch_key)
            adata.obsm["X_pcahm"] = ho.Z_corr.T
            use_rep = "X_pcahm"
            print("Harmony batch correction complete")
        except Exception as e:
            print(f"Harmony failed/unavailable ({e}); proceeding on uncorrected PCA")
    else:
        print("Single batch — no correction needed")

    sc.pp.neighbors(adata, n_pcs=min(n_pcs, adata.obsm[use_rep].shape[1]), use_rep=use_rep)
    sc.tl.umap(adata)

    mdata = mu.MuData({"prot": adata})
    mdata.write(args.out)
    print(f"Wrote {args.out} (use_rep={use_rep})")


if __name__ == "__main__":
    main()

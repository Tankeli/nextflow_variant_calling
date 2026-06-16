#!/usr/bin/env python3
"""
Protein/ADT dimensionality reduction — ported from scripts/15_prot_dimensionality_reduction.ipynb.

PCA (arpack) on the normalized protein modality, neighborhood graph, and UMAP.
Checkpoint: prot_04_dimred.h5mu.

Usage: prot_dimred.py --in IN.h5mu --sample S [--n_pcs 20] --out OUT.h5mu
"""
import argparse
import os
import warnings

import anndata as ad
import muon as mu
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
ad.settings.allow_write_nullable_strings = True


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--n_pcs", type=int, default=20)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    mdata = mu.read(args.inp)

    n_pcs = min(args.n_pcs, mdata["prot"].n_vars - 1)
    sc.pp.pca(mdata["prot"], svd_solver="arpack")
    sc.pp.neighbors(mdata["prot"], n_pcs=n_pcs)
    sc.tl.umap(mdata["prot"])
    print(f"Protein dimred {args.sample}: PCA + UMAP on {n_pcs} PCs")

    mdata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RNA dimensionality reduction — ported from scripts/04_scRNA_dimension_reduction.ipynb.

Uses the log1p_norm layer as .X, restricts PCA to highly-deviant genes, then PCA / t-SNE / UMAP.

Usage: rna_dimred.py --in IN.h5ad --sample S [--n_pcs 50] [--n_neighbors 15] --out OUT.h5ad
"""
import argparse
import os
import warnings

import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--n_pcs", type=int, default=50)
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    adata = sc.read(args.inp)

    adata.X = adata.layers["log1p_norm"]
    # use the deviance-selected genes as scanpy's HVG flag
    adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
    sc.tl.tsne(adata, use_rep="X_pca")
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(adata)

    adata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RNA pseudotime — ported from scripts/08_scRNA_pseudotime.ipynb.

Diffusion pseudotime (DPT) on a per-sample annotated object. Ensures PCA/neighbors/UMAP exist,
computes the diffusion map, picks a root cell from a diffusion component, and runs DPT.
Checkpoint: rna_08_pseudotime.h5ad.

Usage: rna_pseudotime.py --in IN.h5ad --sample S --out OUT.h5ad [--cluster_key cell_type]
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
    p.add_argument("--cluster_key", default="cell_type")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    adata = sc.read(args.inp)

    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata)
    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_pcs=30)
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    sc.tl.diffmap(adata)
    # Root = cell with the smallest value on diffusion component 3 (notebook strategy 1).
    root_ix = int(adata.obsm["X_diffmap"][:, 3].argmin())
    adata.uns["iroot"] = root_ix
    sc.tl.dpt(adata)

    ck = args.cluster_key if args.cluster_key in adata.obs else None
    print(f"Pseudotime {args.sample}: root cell {root_ix}; "
          f"dpt range [{adata.obs['dpt_pseudotime'].min():.3f}, {adata.obs['dpt_pseudotime'].max():.3f}]"
          + (f"; cluster_key={ck}" if ck else ""))

    adata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

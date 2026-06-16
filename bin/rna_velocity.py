#!/usr/bin/env python3
"""
RNA velocity — ported from scripts/09_scRNA_RNA_velocity.ipynb.

scVelo over a velocyto loom merged with the per-sample annotated object: filter+normalize,
moments, velocity (stochastic by default; dynamical optional), and the velocity graph.
Checkpoint: rna_09_velocity.h5ad.

Usage:
  rna_velocity.py --in ANNOTATED.h5ad --loom LOOM --sample S --out OUT.h5ad
                  [--mode stochastic|deterministic|dynamical] [--n_top_genes 2000]
"""
import argparse
import os
import warnings

import scanpy as sc
import scvelo as scv

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0
scv.settings.verbosity = 1


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--loom", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--mode", default="stochastic",
                   choices=["stochastic", "deterministic", "dynamical"])
    p.add_argument("--n_top_genes", type=int, default=2000)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    adata_velo = scv.read(args.loom, cache=True)
    adata_velo.var_names_make_unique()
    adata_anno = sc.read_h5ad(args.inp)

    adata = scv.utils.merge(adata_velo, adata_anno)
    print(f"Velocity {args.sample}: merged loom+annotated -> {adata.n_obs} cells x {adata.n_vars} genes")

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=args.n_top_genes)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.tl.umap(adata)

    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata)
    scv.tl.velocity(adata, mode=args.mode)
    scv.tl.velocity_graph(adata)
    print(f"Velocity {args.sample}: mode={args.mode}; graph computed")

    adata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

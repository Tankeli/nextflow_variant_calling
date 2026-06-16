#!/usr/bin/env python3
"""
RNA clustering — ported from scripts/05_scRNA_clustering.ipynb.

Leiden clustering at multiple resolutions. Input is the dimred object (PCA/neighbors present);
neighbors are recomputed for determinism, then leiden at each requested resolution into a
`leiden_res_<r>` column.

Usage: rna_cluster.py --in IN.h5ad --sample S --resolutions 0.25,0.5,1.0 --out OUT.h5ad
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
    p.add_argument("--resolutions", default="0.25,0.5,1.0")
    p.add_argument("--out", required=True)
    return p.parse_args()


def res_key(r):
    # 0.25 -> leiden_res_0_25 ; 1.0 -> leiden_res_1
    s = ("%g" % r).replace(".", "_")
    return f"leiden_res_{s}"


def main():
    args = parse_args()
    adata = sc.read(args.inp)

    if "neighbors" not in adata.uns:
        sc.pp.neighbors(adata, n_pcs=30)
    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    resolutions = [float(x) for x in args.resolutions.split(",") if x.strip()]
    for r in resolutions:
        # flavor='igraph' uses python-igraph's community_leiden (this env ships igraph, not
        # leidenalg); n_iterations=2/directed=False are scanpy's recommended igraph defaults.
        sc.tl.leiden(adata, key_added=res_key(r), resolution=r,
                     flavor="igraph", n_iterations=2, directed=False)
    print(f"Clustered {args.sample} at resolutions {resolutions}")

    adata.write(args.out)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

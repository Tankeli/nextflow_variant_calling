#!/usr/bin/env python3
"""
Protein/ADT annotation — ported from scripts/17_prot_annotation.ipynb.

Leiden clustering on the batch-corrected protein modality + Wilcoxon marker proteins per cluster,
then a marker-based automatic cell-type label per cluster (assigns the lineage whose canonical
markers are the cluster's top-ranked proteins; unmatched clusters keep "Unknown"). The notebook
left annotation manual; here it is automated so the stage emits concrete labels.
Checkpoint: prot_06_annotated.h5mu + prot_celltypes.csv.

Usage: prot_annotate.py --in IN.h5mu [--resolution 1.0] --out OUT.h5mu --celltypes CT.csv
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

LINEAGE = {
    "T cells": ["CD3", "CD4", "CD8"],
    "B cells": ["CD19", "CD20"],
    "NK cells": ["CD56", "CD16"],
    "Monocytes": ["CD14", "CD11b"],
    "Dendritic": ["CD11c"],
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--out", required=True)
    p.add_argument("--celltypes", required=True)
    return p.parse_args()


def label_cluster(top_proteins):
    """Map a cluster to a lineage by matching its top-ranked proteins to canonical markers."""
    for ct, markers in LINEAGE.items():
        for m in markers:
            if any(m in p for p in top_proteins):
                return ct
    return "Unknown"


def main():
    args = parse_args()
    mdata = mu.read(args.inp)
    adata = mdata["prot"].copy()

    sc.tl.leiden(adata, resolution=args.resolution, key_added="leiden")
    print(f"Protein clustering: {adata.obs['leiden'].nunique()} clusters at res={args.resolution}")

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")
    names = adata.uns["rank_genes_groups"]["names"]

    cluster_to_ct = {}
    for cl in adata.obs["leiden"].cat.categories:
        top = [names[cl][i] for i in range(min(5, len(names[cl])))]
        cluster_to_ct[cl] = label_cluster(top)
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_to_ct).astype(str)
    print("Cluster -> cell type:", cluster_to_ct)

    out_cols = [c for c in ["sample_id", "batch", "leiden", "cell_type"] if c in adata.obs.columns]
    adata.obs[out_cols].to_csv(args.celltypes)

    mdata.mod["prot"] = adata
    mdata.write(args.out)
    print(f"Wrote {args.out} and {args.celltypes}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
RNA annotation — ported from scripts/06_scRNA_annotation.ipynb (CellTypist portion).

Annotates cells with CellTypist using two immune models: coarse (Immune_All_High) and fine
(Immune_All_Low), both with majority voting over the clustering. The model named by
--celltypist_model also populates a primary `cell_type` column used downstream (DE/composition).

The notebook additionally maps onto a Zeng bone-marrow reference via scArches/scVI; that requires a
separate reference atlas and is intentionally left out of this stage (it belongs with integration /
a dedicated reference-mapping module).

Usage:
  rna_annotate.py --in IN.h5ad --sample S --celltypist_model Immune_All_Low.pkl
                  [--resolution 1.0] --out OUT.h5ad --celltypes CT.csv
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
    p.add_argument("--celltypist_model", default="Immune_All_Low.pkl")
    p.add_argument("--resolution", type=float, default=1.0)
    p.add_argument("--out", required=True)
    p.add_argument("--celltypes", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    import celltypist
    from celltypist import models

    adata = sc.read(args.inp)

    # CellTypist expects log1p of CP10k over raw counts, dense.
    adata_ct = adata.copy()
    adata_ct.X = adata.layers["counts"]
    sc.pp.normalize_total(adata_ct, target_sum=10**4)
    sc.pp.log1p(adata_ct)
    adata_ct.X = adata_ct.X.toarray()

    wanted = {"Immune_All_High.pkl", "Immune_All_Low.pkl", args.celltypist_model}
    models.download_models(force_update=False, model=sorted(wanted))

    model_high = models.Model.load(model="Immune_All_High.pkl")
    model_low = models.Model.load(model="Immune_All_Low.pkl")

    pred_high = celltypist.annotate(adata_ct, model=model_high, majority_voting=True).to_adata()
    adata.obs["celltypist_cell_label_coarse"] = pred_high.obs.loc[adata.obs.index, "majority_voting"]
    adata.obs["celltypist_conf_score_coarse"] = pred_high.obs.loc[adata.obs.index, "conf_score"]

    pred_low = celltypist.annotate(adata_ct, model=model_low, majority_voting=True).to_adata()
    adata.obs["celltypist_cell_label_fine"] = pred_low.obs.loc[adata.obs.index, "majority_voting"]
    adata.obs["celltypist_conf_score_fine"] = pred_low.obs.loc[adata.obs.index, "conf_score"]

    # Primary cell_type from the requested model (fine = Low, coarse = High).
    primary = "celltypist_cell_label_fine" if "Low" in args.celltypist_model else "celltypist_cell_label_coarse"
    adata.obs["cell_type"] = adata.obs[primary].astype(str)

    adata.obs[["cell_type", "celltypist_cell_label_coarse", "celltypist_cell_label_fine"]].to_csv(args.celltypes)
    adata.write(args.out)
    print(f"Wrote {args.out} and {args.celltypes}")


if __name__ == "__main__":
    main()

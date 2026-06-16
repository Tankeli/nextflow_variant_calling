#!/usr/bin/env python3
"""
RNA integration — ported from scripts/07_scRNA_integration.ipynb.

Concatenates all per-sample annotated objects and runs a configurable set of batch-integration
methods (scVI / scANVI / BBKNN / Seurat), benchmarks them with scib fast metrics, and writes the
best-scoring embedding as the primary checkpoint (rna_07_integrated.h5ad). Per-method objects are
also written alongside.

The per-sample objects are anonymised by Nextflow (all named rna_06_annotated.h5ad under
input*/), so we rely on the obs columns stamped at QC time (sample_id / batch / patient /
timepoint / cell_type) rather than file names.

Usage:
  rna_integrate.py --inputs . --methods scvi,scanvi,bbknn --batch_key batch
                   --out rna_07_integrated.h5ad --metrics integration_metrics.csv
                   [--label_key cell_type] [--n_top_genes 2000]
"""
import argparse
import glob
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs", required=True, help="dir containing the per-sample annotated h5ads")
    p.add_argument("--methods", default="scvi,scanvi,bbknn")
    p.add_argument("--batch_key", default="batch")
    p.add_argument("--label_key", default="cell_type")
    p.add_argument("--n_top_genes", type=int, default=2000)
    p.add_argument("--out", required=True)
    p.add_argument("--metrics", default="integration_metrics.csv")
    return p.parse_args()


def load_cohort(inputs, batch_key):
    paths = sorted(glob.glob(os.path.join(inputs, "**", "*.h5ad"), recursive=True))
    paths = [p for p in paths if os.path.basename(p) != os.path.basename(inputs)]
    if not paths:
        raise FileNotFoundError(f"No annotated .h5ad files under {inputs}")
    print(f"Integrating {len(paths)} samples")
    adatas = [sc.read_h5ad(p) for p in paths]
    keys = []
    for i, a in enumerate(adatas):
        keys.append(str(a.obs["sample_id"].iloc[0]) if "sample_id" in a.obs else f"sample_{i}")
    adata = sc.concat(adatas, join="outer", label="_concat_key", keys=keys)
    adata.obs_names_make_unique()
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    if batch_key not in adata.obs:
        adata.obs[batch_key] = adata.obs["_concat_key"].astype(str)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["logcounts"] = adata.X.copy()
    print(f"Combined: {adata.n_obs} cells x {adata.n_vars} genes; batches={adata.obs[batch_key].nunique()}")
    return adata


def prep_hvg(adata, batch_key, n_top_genes):
    sc.pp.filter_genes(adata, min_cells=1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="cell_ranger", batch_key=batch_key)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata[:, adata.var["highly_variable"]].copy()


def run_scvi(adata_hvg, batch_key):
    import scvi
    a = adata_hvg.copy()
    scvi.model.SCVI.setup_anndata(a, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(a)
    max_epochs = int(np.min([round((20000 / a.n_obs) * 400), 400]))
    model.train(max_epochs=max_epochs)
    a.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(a, use_rep="X_scVI")
    sc.tl.umap(a)
    return a, model, max_epochs


def run_scanvi(adata_scvi, model_scvi, label_key, max_epochs_scvi):
    import scvi
    model = scvi.model.SCANVI.from_scvi_model(
        model_scvi, labels_key=label_key, unlabeled_category="unlabelled"
    )
    max_epochs = int(np.min([10, np.max([2, round(max_epochs_scvi / 3.0)])]))
    model.train(max_epochs=max_epochs)
    a = adata_scvi.copy()
    a.obsm["X_scANVI"] = model.get_latent_representation()
    sc.pp.neighbors(a, use_rep="X_scANVI")
    sc.tl.umap(a)
    return a


def run_bbknn(adata_hvg, batch_key):
    import bbknn
    a = adata_hvg.copy()
    a.X = a.layers["logcounts"].copy()
    sc.pp.pca(a)
    nwb = 25 if a.n_obs > 100000 else 3
    bbknn.bbknn(a, batch_key=batch_key, neighbors_within_batch=nwb)
    sc.tl.umap(a)
    return a


def score(adata_unint, adata_int, batch_key, label_key, embed=None):
    import scib
    try:
        return scib.metrics.metrics_fast(adata_unint, adata_int, batch_key, label_key, embed=embed)
    except Exception as e:  # scib is brittle on small/odd inputs
        print(f"  scib metrics failed ({e})")
        return None


def main():
    args = parse_args()
    methods = [m.strip().lower() for m in args.methods.split(",") if m.strip()]

    adata = load_cohort(args.inputs, args.batch_key)
    label_key = args.label_key if args.label_key in adata.obs else "_concat_key"
    adata_hvg = prep_hvg(adata, args.batch_key, args.n_top_genes)

    results = {}        # name -> adata
    metrics_cols = {}   # name -> scib Series
    results["Unintegrated"] = adata_hvg
    m = score(adata_hvg, adata_hvg, args.batch_key, label_key)
    if m is not None:
        metrics_cols["Unintegrated"] = m

    model_scvi = max_epochs_scvi = adata_scvi = None
    if "scvi" in methods or "scanvi" in methods:
        try:
            adata_scvi, model_scvi, max_epochs_scvi = run_scvi(adata_hvg, args.batch_key)
            results["scVI"] = adata_scvi
            m = score(adata_hvg, adata_scvi, args.batch_key, label_key, embed="X_scVI")
            if m is not None:
                metrics_cols["scVI"] = m
        except Exception as e:
            print(f"scVI failed: {e}")

    if "scanvi" in methods and adata_scvi is not None and model_scvi is not None:
        try:
            adata_scanvi = run_scanvi(adata_scvi, model_scvi, label_key, max_epochs_scvi)
            results["scANVI"] = adata_scanvi
            m = score(adata_hvg, adata_scanvi, args.batch_key, label_key, embed="X_scANVI")
            if m is not None:
                metrics_cols["scANVI"] = m
        except Exception as e:
            print(f"scANVI failed: {e}")

    if "bbknn" in methods:
        try:
            adata_bbknn = run_bbknn(adata_hvg, args.batch_key)
            results["BBKNN"] = adata_bbknn
            m = score(adata_hvg, adata_bbknn, args.batch_key, label_key)
            if m is not None:
                metrics_cols["BBKNN"] = m
        except Exception as e:
            print(f"BBKNN failed: {e}")

    if "seurat" in methods:
        print("Seurat integration requires R/Seurat (anndata2ri); skipped in this script "
              "unless explicitly ported — use scVI/scANVI/BBKNN.")

    # Benchmark table + best method selection (mirrors the notebook's scaled overall score).
    best = "scVI" if "scVI" in results else next(iter(m for m in results if m != "Unintegrated"), "Unintegrated")
    if len(metrics_cols) >= 2:
        metrics = pd.concat(metrics_cols.values(), axis="columns").set_axis(list(metrics_cols), axis="columns")
        keep = [r for r in ["ASW_label", "ASW_label/batch", "PCR_batch",
                            "isolated_label_silhouette", "graph_conn"] if r in metrics.index]
        metrics = metrics.loc[keep, :].T
        scaled = (metrics - metrics.min()) / (metrics.max() - metrics.min())
        batch_m = [c for c in ["ASW_label/batch", "PCR_batch", "graph_conn"] if c in scaled]
        bio_m = [c for c in ["ASW_label", "isolated_label_silhouette"] if c in scaled]
        scaled["Batch"] = scaled[batch_m].mean(axis=1) if batch_m else 0.0
        scaled["Bio"] = scaled[bio_m].mean(axis=1) if bio_m else 0.0
        scaled["Overall"] = 0.4 * scaled["Batch"] + 0.6 * scaled["Bio"]
        scaled.to_csv(args.metrics)
        cand = scaled.drop(index="Unintegrated", errors="ignore")["Overall"]
        if len(cand) and not cand.isna().all():
            best = cand.idxmax()
        print(f"Integration benchmark written; best method = {best}")
    else:
        pd.DataFrame({"method": list(results), "note": "metrics unavailable"}).to_csv(args.metrics, index=False)

    # Persist per-method objects + the primary checkpoint (best embedding).
    name_to_file = {"scVI": "rna_07_integrated_scvi.h5ad", "scANVI": "rna_07_integrated_scanvi.h5ad",
                    "BBKNN": "rna_07_integrated_bbknn.h5ad", "Unintegrated": "rna_07_unintegrated.h5ad"}
    for name, a in results.items():
        if name in name_to_file:
            a.write(name_to_file[name])
    results[best].write(args.out)
    print(f"Wrote {args.out} (method={best}) and per-method objects")


if __name__ == "__main__":
    main()

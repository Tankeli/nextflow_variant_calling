#!/usr/bin/env python3
"""
Stage 2 — visualisation. Python port of DDE_31 scripts/2a.Visualisation_HM_PCA.R and
2b.Visualisation_UMAP.R (folded).

2a: clustered heatmaps (raw + batch-corrected) + PCA (condition/batch/replicate) + PCA loadings.
2b: sample-level UMAP coloured by condition/batch/cluster, with diffusion-pseudotime (scanpy DPT)
    standing in for the Seurat+slingshot path (consistent with the RNA branch's rna_pseudotime.py).
    Heavy single-cell tooling (Seurat, slingshot) is intentionally not reproduced.

Usage:
  prot_ms_viz.py --raw matrix_raw.tsv --corrected matrix_combat.tsv --design design_corrected.tsv \
      [--method both] [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import prot_ms_utils as U
import prot_ms_plotting as P


def _load_expr(path, cfg):
    full = U.read_expression_matrix(path, U.cfg_get(cfg, "input.id_column", 1))
    meta_idx = U.meta_columns(cfg, full.shape[1])
    sample_cols = [c for i, c in enumerate(full.columns) if i not in meta_idx]
    return U.coerce_numeric(full[sample_cols])


def _finite_var(mat):
    m = mat.replace([np.inf, -np.inf], np.nan).dropna(how="any")
    return m.loc[m.var(axis=1) > 0]


def _run_pca(mat, design, cfg, tag, cond_col, batch_col, rep_col, color_map, label, out):
    m = _finite_var(mat)
    if m.shape[0] < 3:
        return
    pca = PCA(n_components=min(10, m.shape[1] - 1, m.shape[0]))
    scores = pca.fit_transform(m.T.values)  # samples x PCs
    sc = pd.DataFrame(scores, index=m.columns,
                      columns=[f"PC{i+1}" for i in range(scores.shape[1])])
    sc["sample"] = sc.index
    for col in (cond_col, batch_col, rep_col):
        if col in design.columns:
            sc[col] = design.loc[sc.index, col].astype(str).values
    ve = pca.explained_variance_ratio_
    for pair in U.cfg_get(cfg, "visual.pca_pairs", [[1, 2], [3, 4]]):
        x, y = int(pair[0]), int(pair[1])
        if max(x, y) > scores.shape[1]:
            continue
        P.plot_pca(sc, ve, x, y, cond_col, color_map,
                   f"PCA (PC{x} vs PC{y}) - {label}", out(f"{tag}_PC{x}_PC{y}_condition.png"))
    for grp in (batch_col, rep_col):
        if grp in sc.columns:
            nm = "batch" if grp == batch_col else "donor"
            P.plot_pca(sc, ve, 1, 2, grp, None, f"PCA (PC1 vs PC2) - {label} by {grp}",
                       out(f"{tag}_PC1_PC2_{nm}.png"), show_legend=grp != batch_col)
    loadings = pd.DataFrame(pca.components_.T, index=m.index,
                            columns=[f"PC{i+1}" for i in range(scores.shape[1])])
    loadings.insert(0, "Gene", loadings.index)
    loadings.to_csv(out(f"PCA_loadings_{tag}_all.csv"), index=False)
    for i in range(1, min(6, scores.shape[1]) + 1):
        P.plot_pca_loadings(loadings, f"PC{i}", 20, f"Top/Bottom PC{i} Loadings",
                            out(f"{tag}_PC{i}_loadings.png"))


def _umap(mat, design, cfg, cond_col, batch_col, out):
    try:
        import scanpy as sc
        import anndata as ad
    except ImportError:
        print("[warn] scanpy/anndata unavailable; skipping UMAP")
        return
    m = _finite_var(mat)
    if m.shape[1] < 4:
        print("[warn] too few samples for UMAP; skipping")
        return
    adata = ad.AnnData(m.T.values, obs=design.loc[m.columns].copy())
    adata.var_names = m.index.astype(str)
    adata.obs_names = m.columns.astype(str)
    sc.pp.scale(adata, max_value=10)
    n_comps = min(int(U.cfg_get(cfg, "umap.dims_end", 10)), adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=max(2, n_comps))
    nn = max(2, min(15, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=max(2, n_comps))
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=float(U.cfg_get(cfg, "umap.resolution", 0.6)),
                 flavor="igraph", n_iterations=2, directed=False)
    try:
        adata.uns["iroot"] = 0
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
    except Exception as e:  # noqa: BLE001
        print(f"[warn] DPT failed: {e}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [k for k in (cond_col, batch_col, "leiden", "dpt_pseudotime") if k in adata.obs]
    sc.pl.umap(adata, color=keys, show=False, ncols=2)
    plt.savefig(out("UMAP_overview.png"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    adata.obs.to_csv(out("umap_sample_metadata.csv"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True)
    ap.add_argument("--corrected", default=None)
    ap.add_argument("--design", required=True)
    ap.add_argument("--method", default=None, help="override batch_correction.method")
    ap.add_argument("--config", default=None)
    ap.add_argument("--default_config", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)

    cfg = U.load_config(a.config, a.default_config)
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")
    batch_col = U.cfg_get(cfg, "sample_design.batch_column", "Batch")
    rep_col = U.cfg_get(cfg, "sample_design.replicate_column", "replicate")
    method = (a.method or U.cfg_get(cfg, "batch_correction.method", "both")).lower()

    design = U.build_design(U.load_design(a.design, cfg), cfg)
    raw = _load_expr(a.raw, cfg)
    common = [s for s in raw.columns if s in design.index]
    raw, design = raw[common], design.loc[common]

    cond_levels = U.get_condition_levels(design, cfg)
    palette = dict(zip(cond_levels, ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"][: len(cond_levels)]))
    user_colors = U.cfg_get(cfg, "visual.condition_colors", {}) or {}
    palette.update(user_colors)

    annot = pd.DataFrame({
        "Population": design[cond_col].astype(str),
        "Replicate": design[rep_col].astype(str) if rep_col in design else "NA",
        "Batch": design[batch_col].astype(str) if batch_col in design else "NA",
    }, index=design.index)
    annot_colors = {"Population": palette}

    # heatmaps
    try:
        P.clustermap(_finite_var(raw), annot, annot_colors, "Heatmap of non-corrected data",
                     out("HM_noncorrected.png"))
    except SystemExit as e:
        print(e)

    run_corrected = method != "none" and a.corrected and os.path.exists(a.corrected)
    corrected = _load_expr(a.corrected, cfg)[common] if run_corrected else raw
    tag = "limma" if method == "limma" else "combat"
    label = "limma-corrected" if method == "limma" else "batch-corrected"

    if run_corrected:
        try:
            P.clustermap(_finite_var(corrected), annot, annot_colors,
                         f"Heatmap of {label} data", out(f"HM_{tag}.png"))
        except SystemExit as e:
            print(e)
        _run_pca(corrected, design, cfg, tag, cond_col, batch_col, rep_col, palette, label, out)
        _umap(corrected, design, cfg, cond_col, batch_col, out)

    _run_pca(raw, design, cfg, "noncorrected", cond_col, batch_col, rep_col, palette,
             "non-corrected", out)
    print("2 visualisation complete.")


if __name__ == "__main__":
    main()

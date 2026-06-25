#!/usr/bin/env python3
"""
Stage 4 — detailed interpretation. Python port of DDE_31 scripts 4a-4f (folded).

  4a HM_clusters        heatmap of significant proteins + cutree row-cluster assignments
  4b PCA_sig            PCA on significant proteins + loadings
  4c DEGreports         k-means on per-condition means + condition-mean heatmap
  4d Pseudotime         lightweight stand-in: diffusion pseudotime (scanpy) + per-protein Spearman
                        association (the slingshot+GAM / tradeSeq R path is not reproduced)
  4e Boxplots           per-protein boxplots for stage4.boxplot_proteins (GO-term sets omitted)
  4f Offset_plots       top sample-missingness patterns by condition (needs --norm pre-filter matrix)

Usage:
  prot_ms_stage4.py --analysis matrix_combat.tsv --design design_corrected.tsv --de de_results.csv \
      [--norm NORM.tsv] [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import prot_ms_utils as U
import prot_ms_plotting as P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--norm", default=None)
    ap.add_argument("--config", default=None)
    ap.add_argument("--default_config", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)

    cfg = U.load_config(a.config, a.default_config)
    id_col = U.cfg_get(cfg, "input.id_column", 1)
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")

    full = U.read_expression_matrix(a.analysis, id_col)
    meta_idx = U.meta_columns(cfg, full.shape[1])
    sample_cols = [c for i, c in enumerate(full.columns) if i not in meta_idx]
    expr = U.coerce_numeric(full[sample_cols]).replace([np.inf, -np.inf], np.nan).dropna(how="any")
    design = U.build_design(U.load_design(a.design, cfg), cfg)
    common = [s for s in expr.columns if s in design.index]
    expr, design = expr[common], design.loc[common]

    de = pd.read_csv(a.de)
    if "protein" not in de.columns:
        de = de.rename(columns={de.columns[0]: "protein"})
    sig = [p for p in U.filter_significant_proteins(de, cfg) if p in expr.index]

    cond_levels = U.get_condition_levels(design, cfg)
    palette = dict(zip(cond_levels, ["#1f77b4", "#d62728", "#2ca02c", "#9467bd"][: len(cond_levels)]))
    logp = U.cfg_get(cfg, "stage4.logp_threshold", -0.1)

    # 4a — significant-protein heatmap + row clusters
    if len(sig) >= 2:
        mat = expr.loc[sig, design.index]
        annot = pd.DataFrame({"Cell_type": design[cond_col].astype(str)}, index=design.index)
        try:
            g = P.clustermap(mat, annot, {"Cell_type": palette},
                             f"Significant proteins (-log10p>={logp})",
                             out(f"HM_clusters_logp{logp}.png"))
            k = int(U.cfg_get(cfg, "stage4.hm_cutree_rows", 7))
            z = np.nan_to_num(stats.zscore(mat.values, axis=1))
            link = linkage(z, method="ward")
            clusters = fcluster(link, t=k, criterion="maxclust")
            pd.DataFrame({"Protein": mat.index, "Cluster": clusters}).to_csv(
                out(f"Protein_clusters_logp{logp}.csv"), index=False)
        except SystemExit as e:
            print(e)

    # 4b — PCA on significant proteins
    if len(sig) >= 3:
        m = expr.loc[sig, design.index]
        pca = PCA(n_components=min(10, m.shape[1] - 1, m.shape[0]))
        scores = pca.fit_transform(m.T.values)
        sc = pd.DataFrame(scores, index=m.columns,
                          columns=[f"PC{i+1}" for i in range(scores.shape[1])])
        sc["sample"] = sc.index
        sc[cond_col] = design.loc[sc.index, cond_col].astype(str).values
        for pair in U.cfg_get(cfg, "visual.pca_pairs", [[1, 2]]):
            x, y = int(pair[0]), int(pair[1])
            if max(x, y) <= scores.shape[1]:
                P.plot_pca(sc, pca.explained_variance_ratio_, x, y, cond_col, palette,
                           f"PCA significant proteins (PC{x} vs PC{y})",
                           out(f"PC{x}_PC{y}_logp_{logp}.png"))
        load = pd.DataFrame(pca.components_.T, index=m.index,
                            columns=[f"PC{i+1}" for i in range(scores.shape[1])])
        load.insert(0, "Gene", load.index)
        load.to_csv(out(f"PCA_loadings_logp{logp}.csv"), index=False)
        for i in range(1, min(6, scores.shape[1]) + 1):
            P.plot_pca_loadings(load, f"PC{i}", int(U.cfg_get(cfg, "stage4.pca_top_loadings", 20)),
                                f"Top/Bottom PC{i} Loadings", out(f"PC{i}_loadings_logp{logp}.png"))

    # 4c — k-means on per-condition means
    if len(sig) >= 5:
        means = pd.DataFrame({
            c: expr.loc[sig, [s for s in design.index[design[cond_col] == c]]].mean(axis=1)
            for c in cond_levels})
        kk = min(6, means.shape[0])
        km = KMeans(n_clusters=kk, n_init=25, random_state=0).fit(
            np.nan_to_num(stats.zscore(means.values, axis=1)))
        pd.DataFrame({"Protein": means.index, "Cluster": km.labels_}).to_csv(
            out("DEGreports_cluster_assignments.csv"), index=False)
        try:
            P.clustermap(means, None, None, "Mean expression of significant proteins by condition",
                         out("DEGreports_condition_means.png"))
        except SystemExit as e:
            print(e)

    # 4d — lightweight pseudotime association (DPT + Spearman)
    _pseudotime(expr, design, cfg, out)

    # 4e — per-protein boxplots
    genes = full["Genes"] if "Genes" in full.columns else pd.Series(full.index, index=full.index)
    gene_to_prot = {str(g): p for p, g in genes.items()}
    bdir = os.path.join(a.outdir, "boxplots")
    for name in (U.cfg_get(cfg, "stage4.boxplot_proteins", []) or []):
        pid = gene_to_prot.get(str(name))
        if pid is None or pid not in expr.index:
            continue
        os.makedirs(bdir, exist_ok=True)
        df = pd.DataFrame({"value": expr.loc[pid, design.index].values,
                           "condition": design[cond_col].astype(str).values})
        fig, ax = plt.subplots(figsize=(4.5, 4))
        sns.boxplot(data=df, x="condition", y="value", order=cond_levels, palette=palette, ax=ax)
        sns.stripplot(data=df, x="condition", y="value", order=cond_levels, color="black", ax=ax)
        ax.set_title(str(name))
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in str(name))
        fig.savefig(os.path.join(bdir, f"{safe}_boxplot.png"), dpi=200, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)

    # 4f — missingness offset patterns (needs pre-filter norm matrix)
    if a.norm and os.path.exists(a.norm):
        _offset_patterns(a.norm, design, cfg, out)

    print("4 stage-4 analysis complete.")


def _pseudotime(expr, design, cfg, out):
    try:
        import scanpy as sc
        import anndata as ad
    except ImportError:
        print("[warn] scanpy unavailable; skipping pseudotime")
        return
    m = expr.replace([np.inf, -np.inf], np.nan).dropna(how="any")
    m = m.loc[m.var(axis=1) > 0]
    if m.shape[1] < 4 or m.shape[0] < 10:
        return
    adata = ad.AnnData(m.T.values)
    adata.obs_names = m.columns.astype(str)
    adata.var_names = m.index.astype(str)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(10, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=max(2, min(15, adata.n_obs - 1)))
    try:
        adata.uns["iroot"] = 0
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
    except Exception as e:  # noqa: BLE001
        print(f"[warn] DPT failed: {e}")
        return
    pt = adata.obs["dpt_pseudotime"].values
    rows = []
    for p in m.index:
        rho, pval = stats.spearmanr(m.loc[p].values, pt)
        rows.append((p, rho, pval))
    res = pd.DataFrame(rows, columns=["protein", "spearman_rho", "pvalue"])
    res["padj"] = stats.false_discovery_control(np.nan_to_num(res["pvalue"], nan=1.0), method="bh")
    res.to_csv(out("pseudotime_association_results.csv"), index=False)


def _offset_patterns(norm_path, design, cfg, out):
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")
    norm = U.read_expression_matrix(norm_path, U.cfg_get(cfg, "input.id_column", 1))
    cols = [c for c in norm.columns if c in design.index]
    norm = U.coerce_numeric(norm[cols])
    thr = float(U.cfg_get(cfg, "qc.protein_id_threshold", 2500))
    keep = [c for c in norm.columns if norm[c].notna().sum() >= thr]
    norm = norm[keep]
    d = design.loc[keep]
    cond_levels = U.get_condition_levels(d, cfg)
    presence = norm.notna().astype(int)
    patterns = presence.apply(lambda r: "".join(map(str, r.values)), axis=1)
    top = patterns.value_counts().head(12)
    rows = []
    for pat in top.index:
        vals = np.array([int(x) for x in pat])
        by_cond = {c: vals[[i for i, s in enumerate(norm.columns) if d.loc[s, cond_col] == c]].mean()
                   for c in cond_levels}
        rows.append(by_cond)
    pdf = pd.DataFrame(rows)
    pdf["Pattern"] = [f"Pattern_{i+1}" for i in range(len(pdf))]
    pdf["Count"] = top.values
    pdf.to_csv(out("offset_missingness_patterns.csv"), index=False)


if __name__ == "__main__":
    main()

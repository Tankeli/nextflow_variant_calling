#!/usr/bin/env python3
"""
Plotting helpers for the bulk-proteomics branch. Python port of DDE_31 scripts/plotting.utils.R
(ggplot2 + pheatmap -> matplotlib + seaborn). Figures are functional equivalents, not pixel-identical.

Imported by bin/prot_ms_*.py via PYTHONPATH=$projectDir/bin.
"""
from __future__ import annotations

import re
import sys

import numpy as np
import pandas as pd

# scipy/seaborn draw hierarchical dendrograms recursively; with thousands of leaves the default
# limit (1000) overflows. Bulk-proteomics heatmaps can carry several thousand proteins.
sys.setrecursionlimit(100000)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def _safe_color(c, default="grey"):
    """Translate R colour names like 'grey28'/'gray70' (config-carried from the R pipeline) to a
    matplotlib-valid grayscale tuple; pass valid colours through unchanged."""
    if c is None:
        return default
    m = re.fullmatch(r"gr[ae]y(\d{1,3})", str(c))
    if m:
        v = min(100, int(m.group(1))) / 100.0
        return (v, v, v)
    return c


def _save(fig, path: str, **kw):
    fig.savefig(path, dpi=kw.get("dpi", 300), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_protein_ids_all_samples(expr: pd.DataFrame, threshold: float, out_path: str):
    counts = expr.notna().sum(axis=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.bar(counts.index.astype(str), counts.values, color="black")
    ax.axhline(threshold, ls=":", color="red")
    ax.set_ylabel("Protein IDs")
    ax.set_xticklabels(counts.index.astype(str), rotation=90, ha="center")
    _save(fig, out_path)


def plot_protein_completeness(expr: pd.DataFrame, out_path: str):
    completeness = expr.notna().sum(axis=1) / expr.shape[1] * 100
    bins = np.arange(0, 105, 5)
    labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins) - 1)]
    binned = pd.cut(completeness, bins=bins, include_lowest=True, labels=labels)
    counts = binned.value_counts().reindex(labels, fill_value=0)
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.bar(labels, counts.values, color="black")
    ax.set_title("Protein completeness distribution")
    ax.set_xlabel("% of samples detected in")
    ax.set_ylabel("Number of proteins")
    ax.set_xticklabels(labels, rotation=45, ha="right")
    _save(fig, out_path)


def plot_by_condition(df: pd.DataFrame, metric_col: str, condition_col: str,
                      batch_col: str | None, condition_order, title: str, y_label: str,
                      out_path: str):
    """Mean bar + SE error bars + jittered points, grouped by batch when given."""
    d = df.copy()
    if condition_order is not None:
        d[condition_col] = pd.Categorical(d[condition_col], categories=condition_order, ordered=True)
        d = d.sort_values(condition_col)
    fig, ax = plt.subplots(figsize=(10, 8))
    hue = batch_col if (batch_col and batch_col in d.columns) else None
    sns.barplot(data=d, x=condition_col, y=metric_col, hue=hue, errorbar="se",
                edgecolor="black", ax=ax, capsize=0.1)
    sns.stripplot(data=d, x=condition_col, y=metric_col, hue=hue, dodge=hue is not None,
                  color="black", alpha=0.7, ax=ax, legend=False)
    ax.set_title(title)
    ax.set_xlabel("Population")
    ax.set_ylabel(y_label)
    for lab in ax.get_xticklabels():
        lab.set_rotation(45)
        lab.set_ha("right")
    _save(fig, out_path)


def plot_quant(expr: pd.DataFrame, method: str, out_path: str):
    """Boxplot/density of per-feature intensities by sample."""
    long = expr.melt(var_name="sample", value_name="intensity").dropna()
    fig, ax = plt.subplots(figsize=(10, 8))
    if method == "box":
        sns.boxplot(data=long, x="sample", y="intensity", ax=ax)
        for lab in ax.get_xticklabels():
            lab.set_rotation(45)
            lab.set_ha("right")
        ax.set_ylabel("Feature intensity")
        ax.set_xlabel("")
    else:  # density
        for s, sub in long.groupby("sample"):
            sns.kdeplot(sub["intensity"], ax=ax, label=str(s), warn_singular=False)
        ax.set_xlabel("Feature intensity")
        ax.set_ylabel("Density")
        ax.legend(fontsize=6, ncol=2)
    _save(fig, out_path)


def clustermap(mat: pd.DataFrame, col_annot: pd.DataFrame | None, annot_colors: dict | None,
               title: str, out_path: str, row_scale: bool = True, show_rownames: bool = False,
               width: float = 10, height: float = 8, max_rows: int = 3000):
    """Row-scaled clustered heatmap (pheatmap analogue) via seaborn.clustermap. Caps to the
    `max_rows` most-variable rows so very large matrices stay renderable."""
    m = mat.replace([np.inf, -np.inf], np.nan).dropna(how="any")
    if m.shape[0] == 0:
        raise SystemExit(f"No finite rows available for heatmap: {title}")
    if m.shape[0] > max_rows:
        m = m.loc[m.var(axis=1).sort_values(ascending=False).index[:max_rows]]
    col_colors = None
    if col_annot is not None and annot_colors is not None:
        col_colors = pd.DataFrame(index=col_annot.index)
        for ann in col_annot.columns:
            cmap = annot_colors.get(ann, {})
            col_colors[ann] = col_annot[ann].map(cmap)
        col_colors = col_colors.reindex(m.columns)
    g = sns.clustermap(
        m, z_score=0 if row_scale else None, cmap="RdBu_r", center=0,
        col_colors=col_colors, method="average", metric="correlation",
        xticklabels=True, yticklabels=show_rownames, figsize=(width, height),
    )
    g.fig.suptitle(title)
    g.savefig(out_path, dpi=300)
    plt.close(g.fig)
    return g


def plot_pca(scores: pd.DataFrame, var_exp, pc_x: int, pc_y: int, color_col: str,
             color_map: dict | None, title: str, out_path: str, show_legend: bool = True):
    fig, ax = plt.subplots(figsize=(8, 6))
    groups = scores[color_col].astype(str)
    for grp, sub in scores.groupby(groups):
        col = (color_map or {}).get(grp)
        ax.scatter(sub[f"PC{pc_x}"], sub[f"PC{pc_y}"], label=grp, s=40,
                   color=col if col else None)
    for _, row in scores.iterrows():
        ax.annotate(str(row["sample"]), (row[f"PC{pc_x}"], row[f"PC{pc_y}"]), fontsize=6)
    ax.set_xlabel(f"PC{pc_x} ({var_exp[pc_x-1]*100:.1f}%)")
    ax.set_ylabel(f"PC{pc_y} ({var_exp[pc_y-1]*100:.1f}%)")
    ax.set_title(title)
    if show_legend:
        ax.legend(title=color_col, fontsize=8)
    _save(fig, out_path)


def plot_pca_loadings(loadings: pd.DataFrame, pc_name: str, n_top: int, title: str, out_path: str):
    ordered = loadings.sort_values(pc_name)
    top = pd.concat([ordered.head(n_top), ordered.tail(n_top)]).drop_duplicates(subset="Gene")
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(top["Gene"].astype(str), top[pc_name], color="#dc7717", edgecolor="black")
    ax.set_xlabel("Gene")
    ax.set_ylabel(f"{pc_name} Loading")
    ax.set_title(title)
    for lab in ax.get_xticklabels():
        lab.set_rotation(45)
        lab.set_ha("right")
        lab.set_fontsize(8)
    _save(fig, out_path)


def plot_volcano(plot_df: pd.DataFrame, logfc_col: str, comp: str, colors: dict,
                 logp_threshold: float, label_threshold: float, max_overlaps: int, out_path: str):
    fig, ax = plt.subplots(figsize=(8, 5))
    cmap = {"up": _safe_color(colors.get("up"), "#00798c"),
            "down": _safe_color(colors.get("down"), "#d1495b"),
            "ns": _safe_color(colors.get("ns"), "grey")}
    for reg, sub in plot_df.groupby("regulation"):
        ax.scatter(sub[logfc_col], sub["logpval"], s=8, alpha=0.8, color=cmap.get(reg, "grey"))
    lab_df = plot_df[(plot_df["logpval"] > logp_threshold) &
                     (plot_df[logfc_col].abs() > label_threshold)]
    for _, row in lab_df.head(max_overlaps).iterrows():
        ax.annotate(str(row["Genes"]), (row[logfc_col], row["logpval"]), fontsize=6)
    ax.set_xlabel(f"Log2 Fold Change {comp}")
    ax.set_ylabel(f"-log10(adj.P.val) {comp}")
    ax.set_title(comp)
    _save(fig, out_path)

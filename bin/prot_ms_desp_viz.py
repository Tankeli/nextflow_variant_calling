#!/usr/bin/env python3
"""
Stage 6 (viz) — DESP visualisations. Python port of DDE_31 6a/6c (top-protein heatmap + direction
summary; DE-ranked delta heatmap + cell-type contribution) and 6e (per-patient row-percentile
contribution heatmaps). The ggalluvial sankey is rendered as a diverging direction bar.

Reads the DESP outputs from prot_ms_desp_run.R (per-condition profiles + delta matrix [+ per_patient/])
plus the bulk matrix, DE table, proportions and design.

Usage:
  prot_ms_desp_viz.py --bulk matrix_limma.tsv --de de_results.csv --proportions props.tsv \
      --design design_corrected.tsv --desp_dir DESP_OUT [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import prot_ms_utils as U


def _row_percentile(mat: np.ndarray) -> np.ndarray:
    out = np.zeros_like(mat, dtype=float)
    for i in range(mat.shape[0]):
        x = mat[i].astype(float)
        ok = np.isfinite(x)
        if ok.sum() == 0:
            continue
        r = pd.Series(x[ok]).rank(method="average").values
        out[i, ok] = (r - 1) / max(1, ok.sum() - 1)
    return out


def _gene_label(de: pd.DataFrame) -> dict:
    if "Genes" not in de.columns:
        return {p: p for p in de["protein"]}
    out = {}
    for _, row in de.iterrows():
        g = str(row.get("Genes", "")).split(";")[0]
        out[row["protein"]] = g if g and g != "nan" else row["protein"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--proportions", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--desp_dir", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--default_config", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)

    cfg = U.load_config(a.config, a.default_config)
    id_col = U.cfg_get(cfg, "input.id_column", 1)
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")
    rep_col = U.cfg_get(cfg, "sample_design.replicate_column", "replicate")
    top_n = int(U.cfg_get(cfg, "desp.top_n_proteins", 30))
    top_ct = int(U.cfg_get(cfg, "desp.top_celltypes", 25))

    full = U.read_expression_matrix(a.bulk, id_col)
    meta_idx = U.meta_columns(cfg, full.shape[1])
    sample_cols = [c for i, c in enumerate(full.columns) if i not in meta_idx]
    expr = U.coerce_numeric(full[sample_cols])

    de = pd.read_csv(a.de)
    if "protein" not in de.columns:
        de = de.rename(columns={de.columns[0]: "protein"})
    fc_col = "logFC [Relapse vs Diagnosis]"
    adj_col = "adj.P.Val [Relapse vs Diagnosis]"
    if fc_col not in de.columns:
        print("[warn] no 'Relapse vs Diagnosis' contrast in DE; DESP viz limited")
        return
    de = de.sort_values([adj_col, fc_col], key=lambda s: s.abs() if s.name == fc_col else s)
    labels = _gene_label(de)

    design = U.build_design(U.load_design(a.design, cfg), cfg)
    prop_df = pd.read_csv(a.proportions, sep="\t")
    prop_mat = prop_df.set_index("cell_type")
    props = prop_mat.T  # samples x cell_type

    # ---- 6a: top-protein heatmap + direction bar ----
    de_ranked = de.copy()
    de_ranked["abs_fc"] = de_ranked[fc_col].abs()
    top_hm = de_ranked.head(min(25, len(de_ranked)))
    common = [s for s in design.index if s in expr.columns]
    hm_prot = [p for p in top_hm["protein"] if p in expr.index]
    if hm_prot:
        m = expr.loc[hm_prot, common]
        annot = pd.DataFrame({"Condition": design.loc[common, cond_col].astype(str)}, index=common)
        try:
            g = sns.clustermap(m.apply(lambda r: (r - r.mean()) / (r.std() or 1), axis=1).fillna(0),
                               cmap="RdBu_r", center=0,
                               col_colors=annot["Condition"].map(
                                   {"Diagnosis": "#1f77b4", "Relapse": "#d62728"}),
                               yticklabels=[labels.get(p, p) for p in m.index], figsize=(9, 8))
            g.savefig(out("desp_top_proteins_heatmap.png"), dpi=200)
            plt.close(g.fig)
        except Exception as e:  # noqa: BLE001
            print(f"[warn] top-protein heatmap failed: {e}")
    top_dir = de_ranked.head(min(15, len(de_ranked))).copy()
    top_dir["label"] = top_dir["protein"].map(labels)
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["#d62728" if v >= 0 else "#1f77b4" for v in top_dir[fc_col]]
    ax.barh(top_dir["label"][::-1], top_dir[fc_col][::-1], color=colors[::-1])
    ax.set_xlabel("logFC (Relapse vs Diagnosis)")
    ax.set_title("Top proteins changing between Diagnosis and Relapse")
    fig.savefig(out("desp_protein_direction.png"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    top_hm[["protein", fc_col, adj_col]].to_csv(out("desp_top_proteins.csv"), index=False)

    # ---- 6c: DE-ranked delta heatmap + contribution ----
    delta_path = os.path.join(a.desp_dir, "desp_delta_matrix.tsv")
    diag_path = os.path.join(a.desp_dir, "desp_diagnosis_cell_state_profiles.tsv")
    rel_path = os.path.join(a.desp_dir, "desp_relapse_cell_state_profiles.tsv")
    if os.path.exists(delta_path):
        delta = pd.read_csv(delta_path, sep="\t").set_index("feature")
        top_proteins = [p for p in de_ranked["protein"].head(top_n) if p in delta.index]
        if top_proteins:
            dtop = delta.loc[top_proteins]
            ct_score = dtop.abs().mean(axis=0).sort_values(ascending=False)
            dtop = dtop[ct_score.index[:top_ct]]
            dtop.index = [labels.get(p, p) for p in dtop.index]
            try:
                g = sns.clustermap(dtop.apply(lambda r: (r - r.mean()) / (r.std() or 1), axis=1).fillna(0),
                                   cmap="RdBu_r", center=0, figsize=(11, 9))
                g.fig.suptitle("DESP delta (Relapse - Diagnosis): top DE proteins x cell types")
                g.savefig(out("desp_delta_heatmap.png"), dpi=200)
                plt.close(g.fig)
            except Exception as e:  # noqa: BLE001
                print(f"[warn] delta heatmap failed: {e}")
            delta.loc[top_proteins].to_csv(out("desp_delta_matrix_top.csv"))

        # contribution = profile x mean cell-type proportion per condition
        if os.path.exists(diag_path) and os.path.exists(rel_path):
            _contribution(diag_path, rel_path, top_proteins[:12], props, design, cond_col,
                          labels, out)

    # ---- 6e: per-patient row-percentile contribution heatmaps ----
    for pdir in sorted(glob.glob(os.path.join(a.desp_dir, "per_patient", "*"))):
        if os.path.isdir(pdir):
            _per_patient(pdir, props, design, rep_col, cond_col, os.path.basename(pdir), a.outdir)

    print("6 DESP viz complete.")


def _contribution(diag_path, rel_path, proteins, props, design, cond_col, labels, out):
    if not proteins:
        return
    diag = pd.read_csv(diag_path, sep="\t").set_index("feature")
    rel = pd.read_csv(rel_path, sep="\t").set_index("feature")
    diag_ids = [s for s in design.index[design[cond_col] == "Diagnosis"] if s in props.index]
    rel_ids = [s for s in design.index[design[cond_col] == "Relapse"] if s in props.index]
    mean_d = props.loc[diag_ids].mean() if diag_ids else props.mean()
    mean_r = props.loc[rel_ids].mean() if rel_ids else props.mean()
    rows = []
    for cond, prof, mp in [("Diagnosis", diag, mean_d), ("Relapse", rel, mean_r)]:
        for p in proteins:
            if p not in prof.index:
                continue
            for ct in prof.columns:
                if ct in mp.index:
                    rows.append({"protein": labels.get(p, p), "condition": cond, "cell_type": ct,
                                 "contribution": float(prof.loc[p, ct]) * float(mp[ct])})
    cdf = pd.DataFrame(rows)
    if cdf.empty:
        return
    top_fill = (cdf.groupby("cell_type")["contribution"].apply(lambda x: x.abs().sum())
                .sort_values(ascending=False).head(10).index)
    cdf["cell_type_plot"] = np.where(cdf["cell_type"].isin(top_fill), cdf["cell_type"], "Other")
    cdf.to_csv(out("desp_celltype_contributions.csv"), index=False)
    agg = cdf.groupby(["condition", "protein", "cell_type_plot"])["contribution"].sum().reset_index()
    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    for ax, cond in zip(axes, ["Diagnosis", "Relapse"]):
        sub = agg[agg["condition"] == cond].pivot_table(index="protein", columns="cell_type_plot",
                                                         values="contribution", fill_value=0)
        sub.plot(kind="bar", stacked=True, ax=ax, legend=(cond == "Diagnosis"), width=0.85)
        ax.set_title(cond)
        ax.set_ylabel("Estimated contribution")
    fig.suptitle("Cell-type contribution breakdown for top proteins")
    fig.savefig(out("desp_celltype_contributions.png"), dpi=160, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)


def _per_patient(pdir, props, design, rep_col, cond_col, pid, outroot):
    odir = os.path.join(outroot, "per_patient", pid)
    os.makedirs(odir, exist_ok=True)
    for cond in ("diagnosis", "relapse"):
        path = os.path.join(pdir, f"cell_state_profiles_{cond}.tsv")
        if not os.path.exists(path):
            continue
        prof = pd.read_csv(path, sep="\t").set_index("feature")
        ids = [s for s in design.index[(design[rep_col].astype(str) == pid) &
               (design[cond_col].str.lower() == cond)] if s in props.index]
        if not ids:
            continue
        prop_vec = props.loc[ids].mean()
        cts = [c for c in prof.columns if c in prop_vec.index]
        if not cts:
            continue
        contrib = prof[cts].values * prop_vec[cts].values
        keep = np.nansum(np.abs(contrib), axis=1) > 0
        contrib = contrib[keep]
        if contrib.shape[0] == 0:
            continue
        pct = _row_percentile(contrib)
        fig, ax = plt.subplots(figsize=(min(16, 1 + 0.3 * len(cts)), 10))
        sns.heatmap(pct, cmap="coolwarm", cbar=True, xticklabels=cts, yticklabels=False, ax=ax)
        ax.set_title(f"{pid} all proteins by cell type (row percentile) - {cond}")
        fig.savefig(os.path.join(odir, f"all_proteins_by_celltype_percentile_{cond}.png"),
                    dpi=140, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        pd.DataFrame(contrib, index=prof.index[keep], columns=cts).to_csv(
            os.path.join(odir, f"contribution_{cond}.csv"))


if __name__ == "__main__":
    main()

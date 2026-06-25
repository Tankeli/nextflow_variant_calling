#!/usr/bin/env python3
"""
Stage 3 — differential expression + volcano. Python port of DDE_31 scripts/3a.DE_analysis.R and
3b.Volcano_plots.R (folded into one stage).

limma is reimplemented: OLS fit of the design (default ~0 + condition + replicate), all pairwise
condition contrasts, and empirical-Bayes variance moderation (squeezeVar / fitFDist / trigammaInverse,
faithful to limma) for moderated t-statistics and BH-adjusted p-values. Per-patient paired logFC
(Relapse - Diagnosis, mean difference) is emitted when a replicate column is present. Then writes
volcano PNGs per comparison + a per-patient logFC heatmap.

Usage:
  prot_ms_de.py --matrix filtered_log2_imputed.tsv --design design_reduced.tsv \
      [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import special, stats

import prot_ms_utils as U
import prot_ms_plotting as P


def trigamma_inverse(x: np.ndarray) -> np.ndarray:
    """limma::trigammaInverse — solve trigamma(y)=x by Newton iteration."""
    x = np.asarray(x, dtype=float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(np.abs(dif / y)) < 1e-8:
            break
    return y


def squeeze_var(s2: np.ndarray, df1: float):
    """limma::squeezeVar via fitFDist. Returns (s2_post, df_prior, var_prior)."""
    s2 = np.asarray(s2, dtype=float)
    ok = np.isfinite(s2) & (s2 > 0)
    x = s2[ok]
    n = x.size
    z = np.log(x)
    e = z - special.digamma(df1 / 2) + np.log(df1 / 2)
    emean = e.mean()
    evar = (e.var(ddof=1)) - special.polygamma(1, df1 / 2)
    if evar > 0:
        df2 = 2 * trigamma_inverse(np.array([evar]))[0]
        s20 = np.exp(emean + special.digamma(df2 / 2) - np.log(df2 / 2))
        s2_post = (df1 * s2 + df2 * s20) / (df1 + df2)
    else:
        df2 = np.inf
        s20 = np.exp(emean)
        s2_post = np.full_like(s2, s20)
    return s2_post, df2, s20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--default_config", default=None)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    out = lambda f: os.path.join(a.outdir, f)
    vol = os.path.join(a.outdir, "volcano")
    os.makedirs(vol, exist_ok=True)

    cfg = U.load_config(a.config, a.default_config)
    id_col = U.cfg_get(cfg, "input.id_column", 1)
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")
    rep_col = U.cfg_get(cfg, "sample_design.replicate_column", "replicate")

    full = U.read_expression_matrix(a.matrix, id_col)
    design = U.load_design(a.design, cfg)
    meta_idx = U.meta_columns(cfg, full.shape[1])
    genes = full["Genes"] if "Genes" in full.columns else None

    sample_cols = [c for c in full.columns if c in design["sample"].values]
    expr = U.coerce_numeric(full[sample_cols])
    common = [s for s in expr.columns if s in design.index]
    expr, design = expr[common], design.loc[common]

    # exclude patients
    excl = [str(p).rstrip("_") for p in (U.cfg_get(cfg, "de.exclude_patients", []) or [])]
    if rep_col in design.columns and excl:
        design = design[~design[rep_col].astype(str).str.rstrip("_").isin(excl)]
        expr = expr[design.index.tolist()]

    design[cond_col] = design[cond_col].astype(str)
    conds = sorted(design[cond_col].unique())
    if len(conds) < 2:
        raise SystemExit("Need at least two conditions for DE.")

    # design matrix: ~0 + condition (+ replicate when in formula and present)
    formula = U.cfg_get(cfg, "de.design_formula", "~ 0 + condition + replicate")
    parts = [cond_col]
    if rep_col in design.columns and "replicate" in formula:
        if design[rep_col].nunique() > 1:
            parts.append(rep_col)
    Xd = pd.get_dummies(design[parts].astype(str), columns=parts, drop_first=False)
    # drop one replicate dummy to avoid collinearity with condition block (treatment contrast)
    if rep_col in parts:
        rep_dummies = [c for c in Xd.columns if c.startswith(rep_col + "_")]
        Xd = Xd.drop(columns=rep_dummies[:1])
    Xd = Xd.astype(float)

    Y = expr.loc[:, design.index].values            # proteins x samples
    Xmat = Xd.values                                  # samples x p
    XtX_inv = np.linalg.pinv(Xmat.T @ Xmat)
    B = (XtX_inv @ Xmat.T @ Y.T)                       # p x proteins
    resid = Y.T - Xmat @ B
    df_resid = Xmat.shape[0] - np.linalg.matrix_rank(Xmat)
    sigma2 = (resid ** 2).sum(axis=0) / df_resid       # per protein
    s2_post, df_prior, _ = squeeze_var(sigma2, df_resid)
    df_total = df_resid + (0 if not np.isfinite(df_prior) else df_prior)

    cond_cols = [f"{cond_col}_{c}" for c in conds]
    col_index = {c: i for i, c in enumerate(Xd.columns)}
    results = pd.DataFrame({"protein": expr.index})

    for c1, c2 in combinations(conds, 2):
        contrast = np.zeros(Xmat.shape[1])
        contrast[col_index[f"{cond_col}_{c2}"]] = 1.0
        contrast[col_index[f"{cond_col}_{c1}"]] = -1.0
        logfc = contrast @ B                          # proteins
        v = float(contrast @ XtX_inv @ contrast)
        se = np.sqrt(s2_post * v)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = logfc / se
        p = 2 * stats.t.sf(np.abs(t), df=df_total if np.isfinite(df_total) else df_resid)
        padj = stats.false_discovery_control(np.nan_to_num(p, nan=1.0), method="bh")
        tag = f"{c2} vs {c1}"
        results[f"logFC [{tag}]"] = logfc
        results[f"P.Value [{tag}]"] = p
        results[f"adj.P.Val [{tag}]"] = padj

    if genes is not None:
        results = results.merge(
            pd.DataFrame({"protein": expr.index, "Genes": genes.reindex(expr.index).values}),
            on="protein", how="left")
    results.to_csv(out("de_results.csv"), index=False)

    # per-patient paired logFC (Relapse - Diagnosis), mean difference
    per_patient = None
    if rep_col in design.columns:
        pats = []
        cvals = design[cond_col].values
        pvals = design[rep_col].astype(str).values
        for pid in pd.unique(pvals):
            cs = set(cvals[pvals == pid])
            if {"Diagnosis", "Relapse"} <= cs:
                pats.append(pid)
        if pats:
            per_patient = pd.DataFrame({"protein": expr.index})
            for pid in pats:
                idx = design.index[pvals == pid]
                pd_ = design.loc[idx]
                dx = pd_.index[pd_[cond_col] == "Diagnosis"]
                rl = pd_.index[pd_[cond_col] == "Relapse"]
                per_patient[f"logFC [{pid} Relapse vs Diagnosis]"] = (
                    expr[rl].mean(axis=1).values - expr[dx].mean(axis=1).values)
            if genes is not None:
                per_patient = per_patient.merge(
                    pd.DataFrame({"protein": expr.index, "Genes": genes.reindex(expr.index).values}),
                    on="protein", how="left")
            per_patient.to_csv(out("de_results_per_patient.csv"), index=False)

    # ---- volcano figures ----
    vcfg = cfg.get("volcano", {})
    colors = {"up": vcfg.get("up_color"), "down": vcfg.get("down_color"), "ns": vcfg.get("ns_color")}
    logp_thr = float(vcfg.get("logp_threshold", 1.3))

    def make_volcano(df, tag):
        if "Genes" not in df.columns:
            df = df.copy(); df["Genes"] = df["protein"]
        df["Genes"] = df["Genes"].fillna(df["protein"]).replace("", np.nan).fillna(df["protein"])
        for fc in [c for c in df.columns if c.startswith("logFC [")]:
            comp = fc[len("logFC ["):-1]
            adj = f"adj.P.Val [{comp}]"
            if adj not in df.columns:
                continue
            d = df.copy()
            d["logpval"] = -np.log10(pd.to_numeric(d[adj], errors="coerce"))
            d["regulation"] = np.where((d["logpval"] > logp_thr) & (d[fc] > 0), "up",
                              np.where((d["logpval"] > logp_thr) & (d[fc] < 0), "down", "ns"))
            d = d[np.isfinite(d[fc]) & np.isfinite(d["logpval"])]
            P.plot_volcano(d, fc, comp, colors, logp_thr,
                           float(vcfg.get("logfc_label_threshold", 0)),
                           int(vcfg.get("max_overlaps", 20)),
                           os.path.join(vol, f"Volcano_{tag}_{comp.replace(' ', '_')}.png"))

    make_volcano(results, "all")
    if per_patient is not None:
        make_volcano(per_patient, "per_patient")
        # per-patient logFC heatmap (top 150 by max |logFC|)
        fc_cols = [c for c in per_patient.columns if c.startswith("logFC [")]
        mat = per_patient.set_index("protein")[fc_cols].apply(pd.to_numeric, errors="coerce")
        mat = mat.dropna(how="any")
        if mat.shape[0] >= 2:
            mat = mat.loc[mat.abs().max(axis=1).sort_values(ascending=False).index[:150]]
            mat.columns = [c[len("logFC ["):-1] for c in mat.columns]
            try:
                P.clustermap(mat, None, None, "Per-patient logFC (Relapse vs Diagnosis)",
                             os.path.join(vol, "Per_patient_logFC_heatmap.png"),
                             row_scale=False, width=8, height=10)
            except SystemExit:
                pass

    print(f"3 DE complete: {results.shape[0]} proteins, {len(conds)} conditions")


if __name__ == "__main__":
    main()

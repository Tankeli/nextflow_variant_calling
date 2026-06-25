#!/usr/bin/env python3
"""
Stage 1b — batch-effect correction. Python port of DDE_31 scripts/1b.Data_Batch_correction.R.

limma::removeBatchEffect is reimplemented directly (it is a linear model: fit X ~ [design | batch]
per protein, then subtract the fitted batch component). ComBat uses inmoose.pycombat_norm when
available. Optionally replaces condition labels with balanced random labels (modelling-artefact
control). Emits matrix_limma / matrix_combat / matrix_raw + design_corrected, preserving protein
metadata columns, matching the R outputs (incl. the "downstream-ready copy under combat path" rule).

Usage:
  prot_ms_batch.py --matrix filtered_log2_imputed.tsv --design design_reduced.tsv \
      [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import prot_ms_utils as U


def remove_batch_effect(X: pd.DataFrame, batch: pd.Series, design: np.ndarray) -> pd.DataFrame:
    """limma::removeBatchEffect. X: proteins x samples; batch aligned to X.columns;
    design: samples x p covariates to PRESERVE (incl. intercept)."""
    batch_mm = pd.get_dummies(batch.astype(str), drop_first=True).astype(float)
    batch_mm = batch_mm.reindex(X.columns)
    if batch_mm.shape[1] == 0:
        return X.copy()
    M = np.hstack([design, batch_mm.values])             # samples x (p + nbatch)
    coefs, *_ = np.linalg.lstsq(M, X.values.T, rcond=None)  # (p+nbatch) x proteins
    batch_beta = coefs[design.shape[1]:, :]              # nbatch x proteins
    corrected = X.values - (batch_mm.values @ batch_beta).T
    return pd.DataFrame(corrected, index=X.index, columns=X.columns)


def run_combat(X: pd.DataFrame, batch: pd.Series, mod: pd.DataFrame) -> pd.DataFrame:
    try:
        from inmoose.pycombat import pycombat_norm
    except ImportError:
        raise SystemExit("ComBat requested but inmoose is not installed (pip install inmoose).")
    covars = mod.reset_index(drop=True)
    res = pycombat_norm(X.values, batch.reindex(X.columns).astype(str).values, covar_mod=covars.values)
    return pd.DataFrame(res, index=X.index, columns=X.columns)


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

    cfg = U.load_config(a.config, a.default_config)
    id_col = U.cfg_get(cfg, "input.id_column", 1)
    cond_col = U.cfg_get(cfg, "sample_design.condition_column", "condition")
    batch_col = U.cfg_get(cfg, "sample_design.batch_column", "Batch")

    full = U.read_expression_matrix(a.matrix, id_col)
    design = U.build_design(U.load_design(a.design, cfg), cfg)

    meta_idx = U.meta_columns(cfg, full.shape[1])
    meta = full.iloc[:, meta_idx]
    sample_cols = [c for c in full.columns if c in design["sample"].values]
    expr = U.coerce_numeric(full[sample_cols])

    # align
    common = [s for s in expr.columns if s in design.index]
    expr = expr[common]
    design = design.loc[common]

    # confounded-condition filtering
    excluded = set(U.cfg_get(cfg, "batch_correction.excluded_conditions", []) or [])
    if U.cfg_get(cfg, "batch_correction.remove_single_batch_conditions", False) and batch_col in design:
        nb = design.groupby(cond_col)[batch_col].nunique()
        excluded |= set(nb.index[nb < 2])
    rep_col = U.cfg_get(cfg, "sample_design.replicate_column", "replicate")
    if U.cfg_get(cfg, "batch_correction.remove_single_replicate_conditions", False) and rep_col in design:
        nr = design.groupby(cond_col)[rep_col].nunique()
        excluded |= set(nr.index[nr < 2])
    design = design[~design[cond_col].isin(excluded)]
    expr = expr[design.index.tolist()]

    # optional randomized-label control
    if U.cfg_get(cfg, "batch_correction.randomize_condition_labels", False):
        labels = U.balanced_random_condition_labels(
            design[cond_col],
            U.cfg_get(cfg, "batch_correction.random_condition_labels", ["Condition 1", "Condition 2"]),
            int(U.cfg_get(cfg, "batch_correction.randomize_condition_seed", 42)))
        design["original_condition"] = design[cond_col]
        design[cond_col] = labels
        print("Randomized condition labels (original preserved in original_condition).")

    if batch_col not in design or cond_col not in design:
        raise SystemExit(f"Missing batch ({batch_col}) or condition ({cond_col}) column.")

    # design matrix to preserve (~ condition, with intercept) for removeBatchEffect / ComBat covars
    mod = pd.get_dummies(design[cond_col].astype(str), drop_first=True).astype(float)
    mod.insert(0, "(Intercept)", 1.0)

    def save(mat: pd.DataFrame, path: str):
        df = pd.concat([meta.loc[mat.index], mat], axis=1)
        df.to_csv(path, sep="\t", index=False, na_rep="")

    method = str(U.cfg_get(cfg, "batch_correction.method", "both")).lower()
    if method not in ("none", "limma", "combat", "both"):
        raise SystemExit(f"Invalid batch_correction.method: {method}")

    downstream = expr
    if method in ("limma", "both"):
        limma = remove_batch_effect(expr, design[batch_col], mod.values)
        save(limma, out("matrix_limma.tsv"))
        if method == "limma":
            downstream = limma
    if method in ("combat", "both"):
        combat = run_combat(expr, design[batch_col], mod)
        save(combat, out("matrix_combat.tsv"))
        downstream = combat

    save(expr, out("matrix_raw.tsv"))
    # keep later stages working when method is none/limma (they read matrix_combat)
    if method in ("none", "limma"):
        save(downstream, out("matrix_combat.tsv"))

    design.to_csv(out("design_corrected.tsv"), sep="\t", index=False)
    print(f"1b batch correction complete (method={method}).")


if __name__ == "__main__":
    main()

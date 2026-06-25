#!/usr/bin/env python3
"""
Stage 5 — machine learning. Python port of DDE_31 scripts 5a.Decision_tree.R + 5b.Random_forest.R.

Trains a decision tree (rpart -> sklearn DecisionTreeClassifier) and a random forest (ranger ->
sklearn RandomForestClassifier with impurity importance, class-balanced weights) to classify
condition from protein abundance. Stratified train/test split (caret::createDataPartition analogue).
Emits confusion matrices, RF importance table + plot, and the decision-tree text rules.

Sample sizes here are tiny (~12), so this is exploratory; guards bail out cleanly on degenerate splits.

Usage:
  prot_ms_ml.py --matrix matrix_combat.tsv --design design_corrected.tsv [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier, export_text, plot_tree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import prot_ms_utils as U


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
    seed = int(U.cfg_get(cfg, "ml.seed", 42))

    full = U.read_expression_matrix(a.matrix, id_col)
    meta_idx = U.meta_columns(cfg, full.shape[1])
    sample_cols = [c for i, c in enumerate(full.columns) if i not in meta_idx]
    expr = U.coerce_numeric(full[sample_cols]).replace([np.inf, -np.inf], np.nan).dropna(how="any")
    design = U.build_design(U.load_design(a.design, cfg), cfg)
    common = [s for s in expr.columns if s in design.index]
    expr, design = expr[common], design.loc[common]

    genes = full["Genes"] if "Genes" in full.columns else pd.Series(full.index, index=full.index)

    X = expr.T                                  # samples x proteins
    y = design.loc[X.index, cond_col].astype(str)
    if y.nunique() < 2:
        raise SystemExit("Need >=2 conditions for ML.")

    train_frac = float(U.cfg_get(cfg, "ml.train_fraction", 0.8))
    try:
        Xtr, Xte, ytr, yte = train_test_split(X, y, train_size=train_frac, stratify=y,
                                              random_state=seed)
    except ValueError:
        raise SystemExit("Cannot make a stratified split (too few samples per class).")
    if len(Xte) == 0:
        raise SystemExit("Empty test set; reduce ml.train_fraction or add samples.")

    # decision tree
    dt = DecisionTreeClassifier(ccp_alpha=float(U.cfg_get(cfg, "ml.dt_cp", 0.01)),
                                max_depth=int(U.cfg_get(cfg, "ml.dt_maxdepth", 10)),
                                random_state=seed).fit(Xtr, ytr)
    _confusion(yte, dt.predict(Xte), out("decision_tree_confusion_matrix.csv"))
    fig, ax = plt.subplots(figsize=(12, 9))
    plot_tree(dt, feature_names=[str(g) for g in genes.reindex(X.columns)], class_names=dt.classes_,
              filled=True, fontsize=7, ax=ax)
    fig.savefig(out("decision_tree.png"), dpi=120, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with open(out("decision_tree_rules.txt"), "w") as fh:
        fh.write(export_text(dt, feature_names=list(X.columns.astype(str))))

    # random forest
    rf = RandomForestClassifier(n_estimators=int(U.cfg_get(cfg, "ml.rf_num_trees", 1000)),
                                class_weight="balanced", random_state=seed,
                                n_jobs=-1).fit(Xtr, ytr)
    _confusion(yte, rf.predict(Xte), out("random_forest_confusion_matrix.csv"))
    # Impurity (Gini) importance, computed during training — fast and equivalent in spirit to ranger's
    # built-in importance. Post-hoc permutation importance is impractical here (~6500 features).
    imp = pd.DataFrame({
        "Protein": X.columns,
        "Gene": genes.reindex(X.columns).fillna(pd.Series(X.columns, index=X.columns)).values,
        "Importance": rf.feature_importances_,
    }).sort_values("Importance", ascending=False)
    imp.to_csv(out("random_forest_importance.csv"), index=False)

    top_n = min(int(U.cfg_get(cfg, "ml.top_features", 30)), len(imp))
    top = imp.head(top_n)
    fig, ax = plt.subplots(figsize=(7, 8))
    ax.barh(top["Gene"].astype(str)[::-1], top["Importance"][::-1], color="steelblue")
    ax.set_title(f"Top {top_n} proteins by RF importance")
    ax.set_xlabel("Permutation importance")
    fig.savefig(out("random_forest_top_features.png"), dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("5 ML complete.")


def _confusion(y_true, y_pred, path):
    ct = pd.crosstab(pd.Series(y_pred, name="Prediction"), pd.Series(list(y_true), name="Reference"))
    ct.reset_index().melt(id_vars="Prediction", var_name="Reference", value_name="Freq").to_csv(
        path, index=False)


if __name__ == "__main__":
    main()

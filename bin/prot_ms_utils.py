#!/usr/bin/env python3
"""
Shared helpers for the bulk-proteomics branch (PROTEOMICS subworkflow).

Python port of DDE_31 scripts/utils.R + scripts/imputation.utils.R. The R pipeline routed every
output through a config-driven results tree (cfg$outputs$*). Here paths are explicit per Nextflow
process (each stage stages named inputs and emits named outputs, publishDir handles layout); the
YAML config drives only *parameters* (thresholds, column names, methods), not path routing.

Imported by the bin/prot_ms_*.py stage scripts via `PYTHONPATH=$projectDir/bin`.
"""
from __future__ import annotations

import os
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml


# --------------------------------------------------------------------------- config

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge `override` into `base` (mirrors merge_nested_lists in utils.R).
    Named-mapping values merge key-by-key; everything else (incl. lists) is replaced."""
    if not isinstance(base, dict) or not isinstance(override, dict):
        return override
    out = dict(base)
    for k, v in override.items():
        if isinstance(out.get(k), dict) and isinstance(v, dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out


def load_config(config_path: str | None, default_path: str | None = None) -> dict:
    """Load the default proteomics config and overlay an optional experiment config.

    `default_path` defaults to assets/proteomics_default.yaml relative to this file's repo root.
    Either argument may be None / "" / "[]" (Nextflow's empty-path sentinel)."""
    if default_path in (None, "", "[]"):
        here = os.path.dirname(os.path.abspath(__file__))
        default_path = os.path.join(os.path.dirname(here), "assets", "proteomics_default.yaml")
    with open(default_path) as fh:
        cfg = yaml.safe_load(fh) or {}
    if config_path not in (None, "", "[]") and os.path.exists(config_path):
        with open(config_path) as fh:
            user = yaml.safe_load(fh) or {}
        cfg = _deep_merge(cfg, user)
    return cfg


def cfg_get(cfg: dict, dotted: str, default: Any = None) -> Any:
    """Fetch a nested config value by dotted path, e.g. cfg_get(cfg, 'qc.protein_id_threshold')."""
    node: Any = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return default
        node = node[part]
    return default if node is None else node


# --------------------------------------------------------------------------- IO

def read_expression_matrix(path: str, id_column: Any = 1) -> pd.DataFrame:
    """Read a tab-delimited expression matrix and index it by `id_column` (1-indexed int or name).
    Duplicate IDs are made unique (R make.unique style)."""
    mat = pd.read_csv(path, sep="\t", header=0, dtype=str)
    if isinstance(id_column, (int, np.integer)):
        id_col = mat.columns[int(id_column) - 1]
    else:
        id_col = id_column
    if id_col not in mat.columns:
        raise SystemExit(f"ID column not found in matrix: {id_column}")
    mat.index = make_unique(mat[id_col].astype(str).tolist())
    return mat


def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """Coerce a (string) expression block to float, empty strings -> NaN."""
    return df.replace({"": np.nan}).apply(pd.to_numeric, errors="coerce")


def make_unique(values: Iterable[str]) -> list[str]:
    """Replicate R make.unique: append .1, .2 ... to duplicate entries in order."""
    seen: dict[str, int] = {}
    out: list[str] = []
    for v in values:
        v = str(v)
        if v in seen:
            seen[v] += 1
            out.append(f"{v}.{seen[v]}")
        else:
            seen[v] = 0
            out.append(v)
    return out


def load_contaminants(path: str | None) -> list[str]:
    if path in (None, "", "[]") or not os.path.exists(path):
        return []
    with open(path) as fh:
        vals = [ln.strip() for ln in fh]
    return [v for v in vals if v]


def meta_columns(cfg: dict, n_cols: int) -> list[int]:
    """1-indexed metadata column positions clamped to the matrix width -> 0-indexed list."""
    raw = cfg_get(cfg, "input.metadata_columns", [1, 2, 3, 4]) or []
    cols = [int(c) for c in raw if 1 <= int(c) <= n_cols]
    return [c - 1 for c in cols]


# --------------------------------------------------------------------------- design

def load_design(path: str, cfg: dict) -> pd.DataFrame:
    design = pd.read_csv(path, sep="\t", header=0, dtype=str)
    sample_col = cfg_get(cfg, "sample_design.sample_column", "sample")
    if sample_col not in design.columns:
        sample_col = design.columns[0]
    design["sample"] = design[sample_col].astype(str)
    design.index = design["sample"]
    return design


def build_design(design: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Apply condition mapping / cleanup (parse_from_sample_id path is not used by F157)."""
    sd = cfg.get("sample_design", {})
    cond_col = sd.get("condition_column", "condition")

    if cond_col in design.columns:
        mapping = sd.get("condition_mapping") or {}
        vals = design[cond_col].astype(str)
        if mapping:
            vals = vals.map(lambda x: mapping.get(x, x))
        for pat in (sd.get("condition_cleanup_patterns") or []):
            vals = vals.str.replace(pat, "", regex=True)
        design[cond_col] = vals.str.strip()
    return design


def reconcile_sample_ids(design_samples, matrix_cols) -> dict:
    """Map design sample IDs to matrix column names, exact-match first, else by integer value.

    The Spectronaut matrix often stores numeric sample IDs without leading zeros (e.g. '109', '984')
    while the design file keeps them ('0109', '0984'); R's read.delim silently coerced both to
    integers, hiding the mismatch. We replicate that: a design '0109' matches a matrix column '109'.
    Returns {design_sample: matrix_col} only for samples present in the matrix."""
    matcols = [str(c) for c in matrix_cols]
    exact = set(matcols)
    by_int = {}
    for c in matcols:
        if c.lstrip("-").isdigit():
            by_int.setdefault(str(int(c)), c)
    out = {}
    for s in (str(x) for x in design_samples):
        if s in exact:
            out[s] = s
        elif s.lstrip("-").isdigit() and str(int(s)) in by_int:
            out[s] = by_int[str(int(s))]
    return out


def get_condition_levels(design: pd.DataFrame, cfg: dict) -> list[str]:
    cond_col = cfg_get(cfg, "sample_design.condition_column", "condition")
    configured = cfg_get(cfg, "sample_design.condition_order", []) or []
    observed = list(pd.unique(design[cond_col].astype(str)))
    if not configured:
        return observed
    return [c for c in configured if c in observed] + [c for c in observed if c not in configured]


def balanced_random_condition_labels(condition_values, label_levels=("Condition 1", "Condition 2"),
                                     seed: int = 42) -> np.ndarray:
    """Balanced random relabelling within each original condition (control for 1b)."""
    condition_values = np.asarray([str(x) for x in condition_values])
    uniq = pd.unique(condition_values[pd.notna(condition_values)])
    if len(uniq) != 2:
        raise SystemExit("Balanced randomization expects exactly two condition levels.")
    if len(label_levels) != 2:
        raise SystemExit("random_condition_labels must contain exactly two labels.")
    rng = np.random.default_rng(seed)
    out = np.empty(len(condition_values), dtype=object)
    for cond in uniq:
        idx = np.where(condition_values == cond)[0]
        idx = rng.permutation(idx)
        split = max(1, len(idx) // 2)
        out[idx[:split]] = label_levels[0]
        if split < len(idx):
            out[idx[split:]] = label_levels[1]
    return out


# --------------------------------------------------------------------------- filtering / imputation

def filter_proteins_by_presence(expr: pd.DataFrame, design: pd.DataFrame, cond_col: str,
                                cfg: dict) -> pd.Index:
    """Keep proteins detected in >= min_presence_per_condition fraction of samples in ANY one
    condition group (port of filter_proteins_by_presence)."""
    thr = float(cfg_get(cfg, "qc.min_presence_per_condition", 0.75))
    groups = pd.unique(design[cond_col].astype(str))
    keep = pd.Series(False, index=expr.index)
    notna = expr.notna()
    for g in groups:
        samples = [s for s in design.index[design[cond_col].astype(str) == g] if s in expr.columns]
        if not samples:
            continue
        frac = notna[samples].sum(axis=1) / len(samples)
        keep = keep | (frac >= thr)
    return expr.index[keep.values]


def impute_missing_values(mat: pd.DataFrame, cfg_qc: dict) -> pd.DataFrame:
    """Port of imputation.utils.R::impute_missing_values (none|zero|min|perseus|knn).
    'mle' falls back to perseus (imputeLCMD has no light Python equivalent)."""
    method = str(cfg_qc.get("imputation_method") or "").lower()
    if not method:
        method = "zero" if cfg_qc.get("impute_na") else "none"
    if method == "none":
        return mat
    out = mat.copy()
    if method == "zero":
        return out.fillna(0.0)
    if method == "min":
        mn = np.nanmin(out.values)
        if not np.isfinite(mn):
            raise SystemExit("Cannot apply min imputation: no observed values.")
        return out.fillna(mn)
    if method in ("perseus", "mle"):
        params = cfg_qc.get("imputation_params") or {}
        downshift = float(params.get("perseus_downshift", 1.8))
        width = float(params.get("perseus_width", 0.3))
        seed = params.get("seed")
        rng = np.random.default_rng(None if seed is None else int(seed))
        for col in out.columns:
            v = out[col].values.astype(float)
            miss = np.isnan(v)
            if not miss.any():
                continue
            obs = v[~miss]
            if obs.size == 0:
                v[miss] = 0.0
                out[col] = v
                continue
            mu, sd = float(np.mean(obs)), float(np.std(obs, ddof=1)) if obs.size > 1 else 0.0
            if not np.isfinite(sd) or sd == 0:
                sd = float(np.median(np.abs(obs - np.median(obs)))) or 1.0
            v[miss] = rng.normal(mu - downshift * sd, width * sd, size=int(miss.sum()))
            out[col] = v
        return out
    if method == "knn":
        try:
            from sklearn.impute import KNNImputer
        except ImportError:
            raise SystemExit("kNN imputation requested but scikit-learn is not installed.")
        k = int((cfg_qc.get("imputation_params") or {}).get("knn_k", 10))
        imp = KNNImputer(n_neighbors=k)
        # impute across features per sample (sklearn imputes column-wise -> transpose protein x sample)
        vals = imp.fit_transform(out.values)
        return pd.DataFrame(vals, index=out.index, columns=out.columns)
    raise SystemExit(f"Unknown imputation method: {method}")


# --------------------------------------------------------------------------- DE helpers

def filter_significant_proteins(de: pd.DataFrame, cfg: dict) -> list[str]:
    """Proteins significant in ANY comparison by the stage4 -log10p / |logFC| thresholds."""
    adj_cols = [c for c in de.columns if c.startswith("adj.P.Val [")]
    fc_cols = [c for c in de.columns if c.startswith("logFC [")]
    if not adj_cols or not fc_cols:
        raise SystemExit("DE results missing expected logFC/adj.P.Val columns.")
    logp_thr = float(cfg_get(cfg, "stage4.logp_threshold", -0.1))
    fc_thr = float(cfg_get(cfg, "stage4.logfc_threshold", -10))
    mask = pd.Series(False, index=de.index)
    for a, l in zip(sorted(adj_cols), sorted(fc_cols)):
        with np.errstate(divide="ignore"):
            logp = -np.log10(pd.to_numeric(de[a], errors="coerce"))
        mask = mask | ((logp >= logp_thr) & (pd.to_numeric(de[l], errors="coerce").abs() >= fc_thr))
    return list(pd.unique(de.loc[mask.values, "protein"]))

#!/usr/bin/env python3
"""
Stage 1a — QC & preprocessing for the bulk-proteomics branch. Python port of DDE_31
scripts/1a.QC_and_prep.R.

Pipeline: align matrices to design -> per-sample protein-ID filter -> contaminant removal ->
presence-per-condition filter -> log2 transform -> imputation. Emits the filtered/log2/imputed
matrix (protein metadata + per-sample values), the reduced design, QC figures and summary tables.

The MSnbase peptide->protein "robust summarization" workshop block in the R script is intentionally
dropped: the Spectronaut input here is already protein-level (Protein.Ids), so it is a no-op.

Usage:
  prot_ms_prep.py --nonnorm NONNORM.tsv --norm NORM.tsv --design DESIGN.tsv \
      [--contaminants C.txt] [--config params.yaml] --outdir .
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import prot_ms_utils as U
import prot_ms_plotting as P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nonnorm", required=True)
    ap.add_argument("--norm", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--contaminants", default=None)
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

    nonnorm_raw = U.read_expression_matrix(a.nonnorm, id_col)
    norm_raw = U.read_expression_matrix(a.norm, id_col)
    design = U.build_design(U.load_design(a.design, cfg), cfg)

    # Reconcile design sample IDs to the matrix column names (handles leading-zero / int forms),
    # then rename the design to the matrix form so every downstream stage uses one canonical naming.
    recon = U.reconcile_sample_ids(design["sample"], norm_raw.columns)
    if not recon:
        raise SystemExit("No overlapping samples between matrix columns and design.")
    design = design[design["sample"].isin(recon.keys())].copy()
    design["sample"] = design["sample"].map(recon)
    design.index = design["sample"]
    sample_cols = list(design["sample"])
    norm = U.coerce_numeric(norm_raw[sample_cols])
    nonnorm = U.coerce_numeric(nonnorm_raw[sample_cols])
    if U.cfg_get(cfg, "qc.replace_na_in_nonnorm_with_zero", True):
        nonnorm = nonnorm.fillna(0.0)

    cond_order = U.get_condition_levels(design, cfg)
    thr = float(U.cfg_get(cfg, "qc.protein_id_threshold", 2500))

    P.plot_protein_ids_all_samples(norm, thr, out("protein_ids_across_samples.png"))
    P.plot_protein_completeness(norm, out("protein_completeness.png"))

    # filter samples by protein-ID count
    keep_samples = [s for s in norm.columns if norm[s].notna().sum() >= thr]
    norm_red, nonnorm_red = norm[keep_samples], nonnorm[keep_samples]
    design_red = design.loc[keep_samples]

    # remove contaminants
    contaminants = set(U.load_contaminants(a.contaminants))
    norm_filt = norm_red[~norm_red.index.isin(contaminants)]
    nonnorm_filt = nonnorm_red[~nonnorm_red.index.isin(contaminants)]

    # presence-per-condition filter
    keep_prot = U.filter_proteins_by_presence(norm_filt, design_red, cond_col, cfg)
    filtered_norm = norm_filt.loc[keep_prot]
    nonnorm_final = nonnorm_filt.loc[keep_prot]

    # QC dataframe (abundance / protein ids by condition+batch)
    qc_df = pd.DataFrame({
        "sample": filtered_norm.columns,
        "abundance": nonnorm_final.sum(axis=0).values,
        "condition": design_red.loc[filtered_norm.columns, cond_col].values,
        "batch": design_red.loc[filtered_norm.columns, batch_col].values
            if batch_col in design_red.columns else "1",
    })
    qc_df["protein_ids"] = filtered_norm.notna().sum(axis=0).values
    P.plot_by_condition(qc_df, "protein_ids", "condition", None, cond_order,
                        "Protein IDs per population", "Total Protein IDs",
                        out("protein_ids_by_condition.png"))
    P.plot_by_condition(qc_df, "abundance", "condition", "batch", cond_order,
                        "Protein abundance per population (filtered)", "Protein abundance",
                        out("abundance_by_condition_batch_filtered.png"))

    # filtering summary
    def mstats(m):
        v = m.values
        tot = v.size
        miss = int(np.isnan(v.astype(float)).sum())
        return miss, (miss / tot if tot else np.nan), int((~np.isnan(v.astype(float))).all(axis=1).sum())

    steps = {"initial_aligned": norm, "after_sample_filter": norm_red,
             "after_contaminant_filter": norm_filt, "after_presence_filter": filtered_norm}
    rows = []
    for name, m in steps.items():
        miss, frac, complete = mstats(m)
        rows.append({"step": name, "n_proteins": m.shape[0], "n_samples": m.shape[1],
                     "missing_values_total": miss, "missing_fraction": round(frac, 4),
                     "peptides_no_missing": complete,
                     "pct_peptides_no_missing": round(100 * complete / m.shape[0], 2) if m.shape[0] else np.nan})
    summary = pd.DataFrame(rows)
    summary["pct_proteins_retained"] = round(100 * summary["n_proteins"] / summary["n_proteins"].iloc[0], 2)
    summary["pct_samples_retained"] = round(100 * summary["n_samples"] / summary["n_samples"].iloc[0], 2)
    summary.to_csv(out("filtering_summary.tsv"), sep="\t", index=False)

    missing_by_sample = pd.DataFrame({
        "sample": filtered_norm.columns,
        "missing_values": filtered_norm.isna().sum(axis=0).values,
        "missing_fraction": filtered_norm.isna().mean(axis=0).round(4).values,
        "condition": design_red.loc[filtered_norm.columns, cond_col].values,
    })
    if batch_col in design_red.columns:
        missing_by_sample["batch"] = design_red.loc[filtered_norm.columns, batch_col].values
    missing_by_sample.to_csv(out("missing_by_sample.tsv"), sep="\t", index=False)

    # log2 + imputation
    offset = float(U.cfg_get(cfg, "qc.log2_pseudocount", 1) or 0)
    log2 = np.log2(filtered_norm + offset)
    log2 = U.impute_missing_values(log2, cfg.get("qc", {}))

    P.plot_quant(log2, "box", out("quant_boxplot.png"))
    P.plot_quant(log2, "density", out("quant_density.png"))
    P.plot_quant(np.log2(nonnorm_final + offset), "density", out("quant_density_nonnorm.png"))

    # write filtered_log2 matrix (protein metadata + per-sample values), R column layout
    meta_idx = U.meta_columns(cfg, norm_raw.shape[1])
    meta = norm_raw.iloc[:, meta_idx].loc[log2.index]
    final = pd.concat([meta, log2], axis=1)
    final.to_csv(out("filtered_log2_imputed.tsv"), sep="\t", index=False, na_rep="")
    design_red.drop(columns=[c for c in ["sample"] if c in design_red.columns and c != design_red.index.name],
                    errors="ignore")
    design_red.to_csv(out("design_reduced.tsv"), sep="\t", index=False)

    print(f"1a prep complete: {final.shape[0]} proteins x {log2.shape[1]} samples")


if __name__ == "__main__":
    main()

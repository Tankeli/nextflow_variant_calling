#!/usr/bin/env python3
"""
RNA compositional analysis — ported from scripts/11_scRNA_composition_analysis.ipynb.

Cell-type composition shifts across the condition (timepoint) on the integrated cohort object,
using scCODA (pertpy). Falls back to a descriptive proportion-shift table (mean %Rel - mean %Dx
per cell type) if scCODA is unavailable or fails to converge on this small cohort. Always emits
composition_results.csv (cell_type,effect,hdi_low,hdi_high).

Usage:
  rna_composition.py --in INTEGRATED.h5ad --out composition_results.csv
                     [--cell_type_key cell_type] [--sample_key sample_id]
                     [--condition_key timepoint] [--fdr 0.2]
"""
import argparse
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"
sc.settings.verbosity = 0

COLUMNS = ["cell_type", "effect", "hdi_low", "hdi_high"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--cell_type_key", default="cell_type")
    p.add_argument("--sample_key", default="sample_id")
    p.add_argument("--condition_key", default="timepoint")
    p.add_argument("--fdr", type=float, default=0.2)
    return p.parse_args()


def proportion_fallback(adata, cell_type_key, sample_key, condition_key):
    """Descriptive: mean per-sample %composition difference (Rel - Dx) per cell type."""
    df = adata.obs[[sample_key, condition_key, cell_type_key]].astype(str)
    prop = (
        pd.crosstab([df[sample_key], df[condition_key]], df[cell_type_key], normalize="index") * 100
    )
    prop = prop.reset_index()
    means = prop.groupby(condition_key).mean(numeric_only=True)
    conds = list(means.index)
    rows = []
    if "Rel" in conds and "Dx" in conds:
        diff = means.loc["Rel"] - means.loc["Dx"]
    else:
        diff = means.iloc[-1] - means.iloc[0] if len(conds) >= 2 else means.iloc[0] * 0
    for ct, val in diff.items():
        rows.append({"cell_type": ct, "effect": float(val), "hdi_low": np.nan, "hdi_high": np.nan})
    return pd.DataFrame(rows, columns=COLUMNS)


def sccoda_effects(adata, cell_type_key, sample_key, condition_key, fdr):
    import pertpy as pt

    sccoda = pt.tl.Sccoda()
    data = sccoda.load(
        adata, type="cell_level", generate_sample_level=True,
        cell_type_identifier=cell_type_key, sample_identifier=sample_key,
        covariate_obs=[condition_key],
    )
    data = sccoda.prepare(data, modality_key="coda", formula=condition_key,
                          reference_cell_type="automatic")
    sccoda.run_nuts(data, modality_key="coda", rng_key=1234)
    sccoda.set_fdr(data, fdr)
    coda = data["coda"]
    eff_keys = [k for k in coda.varm.keys() if k.startswith("effect_df")]
    if not eff_keys:
        raise RuntimeError("no effect_df in scCODA output")
    eff = coda.varm[eff_keys[0]]
    eff = eff if isinstance(eff, pd.DataFrame) else pd.DataFrame(eff, index=coda.var_names)
    out = pd.DataFrame({"cell_type": list(eff.index)})
    fc_col = next((c for c in eff.columns if "log2" in str(c).lower() or "fold" in str(c).lower()), None)
    eff_col = next((c for c in eff.columns if "effect" in str(c).lower() or "final" in str(c).lower()), None)
    out["effect"] = eff[eff_col].values if eff_col else (eff[fc_col].values if fc_col else eff.iloc[:, 0].values)
    lo = next((c for c in eff.columns if "hdi" in str(c).lower() and "3" in str(c)), None)
    hi = next((c for c in eff.columns if "hdi" in str(c).lower() and "97" in str(c)), None)
    out["hdi_low"] = eff[lo].values if lo else np.nan
    out["hdi_high"] = eff[hi].values if hi else np.nan
    return out[COLUMNS]


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.inp)

    for key in (args.cell_type_key, args.condition_key):
        if key not in adata.obs:
            print(f"obs missing '{key}'; writing empty result")
            pd.DataFrame(columns=COLUMNS).to_csv(args.out, index=False)
            return

    try:
        res = sccoda_effects(adata, args.cell_type_key, args.sample_key, args.condition_key, args.fdr)
        print(f"scCODA effects computed for {len(res)} cell types")
    except Exception as e:
        print(f"scCODA unavailable/failed ({e}); writing descriptive proportion shifts")
        res = proportion_fallback(adata, args.cell_type_key, args.sample_key, args.condition_key)

    res.to_csv(args.out, index=False)
    print(f"Wrote {args.out}: {len(res)} rows")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Cohort-level QC summary across all samples.

Reads every per-sample <sample>_qc_metrics.csv and draws a one-page panel: cells per
sample (total vs pass-QC), median genes, median %MT, and doublet rate — the pipeline-level
"are the samples comparable?" view (in the spirit of DDE_32 preliminary_plots.R).

Usage: cohort_summary.py <qc_metrics.csv>...
Output (cwd): cohort_summary.png/.pdf, cohort_summary.csv
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

csvs = sys.argv[1:]
if not csvs:
    sys.exit("no qc_metrics.csv inputs given")

rows = []
for f in csvs:
    df = pd.read_csv(f)
    if "sample_id" not in df.columns:
        print(f"[warn] {f} has no sample_id; skipping")
        continue
    sid = str(df["sample_id"].iloc[0])
    passed = df["pass_qc"].astype(str).str.lower().isin(["true", "1"]) if "pass_qc" in df else pd.Series(True, index=df.index)
    pf = df[passed]
    rows.append({
        "sample_id":      sid,
        "n_cells_total":  len(df),
        "n_cells_pass":   int(passed.sum()),
        "median_genes":   float(np.median(pf["n_genes"])) if "n_genes" in pf and len(pf) else np.nan,
        "median_mt_pct":  float(np.median(pf["pct_counts_mt"])) if "pct_counts_mt" in pf and len(pf) else np.nan,
        "doublet_rate":   float(100 * df["predicted_doublet"].astype(str).str.lower().isin(["true", "1"]).mean()) if "predicted_doublet" in df else np.nan,
    })

summary = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
summary.to_csv("cohort_summary.csv", index=False)
print(summary.to_string(index=False))

x = np.arange(len(summary))
fig, axes = plt.subplots(2, 2, figsize=(max(10, 1.2 * len(summary)), 10))
fig.suptitle("Cohort QC summary", fontsize=16, y=0.995)

axes[0, 0].bar(x - 0.2, summary["n_cells_total"], width=0.4, label="total", color="lightgray")
axes[0, 0].bar(x + 0.2, summary["n_cells_pass"], width=0.4, label="pass QC", color="tab:blue")
axes[0, 0].set(ylabel="Cells", title="Cells per sample")
axes[0, 0].legend()

axes[0, 1].bar(x, summary["median_genes"], color="tab:green")
axes[0, 1].set(ylabel="Median genes (pass cells)", title="Gene complexity")

axes[1, 0].bar(x, summary["median_mt_pct"], color="tab:red")
axes[1, 0].set(ylabel="Median % MT (pass cells)", title="Mitochondrial fraction")

axes[1, 1].bar(x, summary["doublet_rate"], color="tab:purple")
axes[1, 1].set(ylabel="% predicted doublets", title="Doublet rate")

for ax in axes.flat:
    ax.set_xticks(x)
    ax.set_xticklabels(summary["sample_id"], rotation=45, ha="right", fontsize=9)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"cohort_summary.{ext}", bbox_inches="tight")
plt.close(fig)
print("Wrote cohort_summary.png/.pdf and cohort_summary.csv")

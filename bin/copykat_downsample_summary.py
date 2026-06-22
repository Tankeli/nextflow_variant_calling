#!/usr/bin/env python3
"""
CopyKAT cell-number confounder summary: is the per-sample ARI driven by N (cell count)?

For each sample, collects the robustness ARI / call-agreement / seed-switch-rate at 100% (full-N,
from the original controls sweep) and at the 25/50/75% downsamples, then plots each metric vs N.
If ARI is ~flat across N within a sample (and samples overlap at matched N), the cross-sample ARI
differences are NOT a cell-count artefact.

Usage: copykat_downsample_summary.py <ds_analysis_dir> <full_analysis_dir> <out_prefix> \
                                     <sample:fullN> [<sample:fullN> ...]
  e.g. ... results_controls/.../downsample/_analysis results_controls/robustness/_analysis \
           results_controls/.../downsample/_analysis/downsample PBM_2:7302 PBM_1:6248 PBMMC_3:5548
Reads <ds_analysis>/<sample>_frac<f>_copykat_{ari,stability,seed_switch}.csv and
      <full_analysis>/<sample>_copykat_{ari,stability,seed_switch}.csv.
"""
import sys, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ds_dir, full_dir, out_prefix = sys.argv[1], sys.argv[2], sys.argv[3]
samples = [a.split(":")[0] for a in sys.argv[4:]]
fullN   = {a.split(":")[0]: int(a.split(":")[1]) for a in sys.argv[4:]}
FRACS = [0.25, 0.50, 0.75]


def mean_ari(path):
    if not os.path.exists(path):
        return np.nan
    m = pd.read_csv(path, index_col=0).values
    iu = np.triu_indices(len(m), 1)
    return float(np.nanmean(m[iu]))


def mean_agree(path):
    return float(pd.read_csv(path)["agreement"].mean()) if os.path.exists(path) else np.nan


def mean_switch(path):
    return float(pd.read_csv(path)["cell_switch_rate"].mean()) if os.path.exists(path) else np.nan


rows = []
for s in samples:
    # 100% from the full-N controls sweep
    rows.append({"sample": s, "fraction": 1.0, "N": fullN[s],
                 "mean_ARI": mean_ari(f"{full_dir}/{s}_copykat_ari.csv"),
                 "mean_agreement": mean_agree(f"{full_dir}/{s}_copykat_stability.csv"),
                 "seed_switch": mean_switch(f"{full_dir}/{s}_copykat_seed_switch.csv")})
    for f in FRACS:
        lab = f"{s}_frac{f:.2f}"
        rows.append({"sample": s, "fraction": f, "N": int(round(f * fullN[s])),
                     "mean_ARI": mean_ari(f"{ds_dir}/{lab}_copykat_ari.csv"),
                     "mean_agreement": mean_agree(f"{ds_dir}/{lab}_copykat_stability.csv"),
                     "seed_switch": mean_switch(f"{ds_dir}/{lab}_copykat_seed_switch.csv")})

df = pd.DataFrame(rows).sort_values(["sample", "N"])
df.to_csv(f"{out_prefix}_summary_table.csv", index=False)
print(df.to_string(index=False))

fig, ax = plt.subplots(1, 3, figsize=(18, 5.5))
metrics = [("mean_ARI", "mean pairwise ARI (seed/param reproducibility)"),
           ("mean_agreement", "mean per-cell call agreement"),
           ("seed_switch", "seed switch-rate")]
for a, (col, title) in zip(ax, metrics):
    for s, g in df.groupby("sample"):
        g = g.sort_values("N")
        a.plot(g["N"], g[col], marker="o", label=s)
    a.set_xlabel("number of cells (N)"); a.set_title(title); a.grid(alpha=0.3)
    if col != "seed_switch":
        a.set_ylim(0, 1)
ax[0].legend(title="sample", fontsize=9)
fig.suptitle("CopyKAT robustness vs cell number — downsampling test (100/75/50/25%)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{out_prefix}_robustness_vs_N.png", dpi=140)
print(f"Wrote {out_prefix}_summary_table.csv and {out_prefix}_robustness_vs_N.png")

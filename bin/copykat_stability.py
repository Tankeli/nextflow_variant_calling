#!/usr/bin/env python3
"""
CopyKAT seed/parameter stability + classification boundary, over a COPYKAT_ROBUSTNESS sweep.

Reads every `*_copykat_prediction.txt` under <sweep_dir> (one per param x seed combo, named
`ks<>_win<>_ng<>_<dist>_norm<0|1>_seed<>_copykat_prediction.txt`), and quantifies how stable the
per-cell aneuploid/diploid call is:
  - per-cell consensus call + fraction-agreement + call entropy (over all combos),
  - seed switch-rate: how often a cell's call flips across seeds within one parameter set,
  - pairwise Adjusted Rand Index between combos (seed stability) and the aneuploid fraction as a
    function of KS.cut (where the diploid/aneuploid boundary lies).

Usage: copykat_stability.py <sample> <sweep_dir> [mapped_h5ad]
Output (cwd): <sample>_copykat_stability.csv  + figures.
"""
import sys, re, glob, os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

sample, sweep_dir = sys.argv[1], sys.argv[2]
mapped_h5ad = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "NONE") else None

COMBO_RE = re.compile(
    r"ks(?P<ks>[0-9.]+)_win(?P<win>\d+)_ng(?P<ng>\d+)_(?P<dist>[a-zA-Z]+)_"
    r"norm(?P<norm>[01])_seed(?P<seed>\d+)")

files = sorted(glob.glob(os.path.join(sweep_dir, "**", "*_copykat_prediction.txt"),
                         recursive=True))
if not files:
    sys.exit(f"No *_copykat_prediction.txt under {sweep_dir}")

# Build a cell x combo table of calls.
calls, meta = {}, {}
for f in files:
    m = COMBO_RE.search(os.path.basename(f))
    if not m:
        print(f"[skip] cannot parse combo from {f}")
        continue
    combo = m.group(0)
    meta[combo] = {k: m.group(k) for k in ("ks", "win", "ng", "dist", "norm", "seed")}
    df = pd.read_csv(f, sep="\t")
    name_col = "cell.names" if "cell.names" in df.columns else df.columns[0]
    pred_col = "copykat.pred" if "copykat.pred" in df.columns else df.columns[-1]
    calls[combo] = pd.Series(df[pred_col].astype(str).values,
                             index=df[name_col].astype(str).values)

mat = pd.DataFrame(calls)            # cells x combos
meta = pd.DataFrame(meta).T
meta[["ks", "win", "ng", "seed"]] = meta[["ks", "win", "ng", "seed"]].apply(pd.to_numeric)
print(f"{sample}: {mat.shape[0]} cells x {mat.shape[1]} combos")

# ---- per-cell consensus + stability (over all combos) ----
def cell_stats(row):
    vals = row[row.notna() & (row != "not.defined")]
    n = len(vals)
    if n == 0:
        return pd.Series({"consensus": "not.defined", "agreement": np.nan,
                          "entropy": np.nan, "n_aneuploid": 0, "n_calls": 0})
    vc = vals.value_counts()
    p = vc / n
    ent = float(-(p * np.log2(p)).sum())
    return pd.Series({"consensus": vc.idxmax(), "agreement": float(vc.iloc[0] / n),
                      "entropy": ent, "n_aneuploid": int((vals == "aneuploid").sum()),
                      "n_calls": int(n)})

stab = mat.apply(cell_stats, axis=1)
stab.index.name = "barcode"
stab.to_csv(f"{sample}_copykat_stability.csv")
unstable = float((stab["agreement"] < 1.0).mean())
print(f"{sample}: mean per-cell agreement={stab['agreement'].mean():.3f}; "
      f"{unstable*100:.1f}% of cells flip across combos")

# ---- seed switch-rate within each non-seed parameter set ----
def to_binary(s):  # aneuploid=1, diploid=0, else NaN
    return s.map({"aneuploid": 1.0, "diploid": 0.0})

param_key = meta.assign(pk=meta[["ks", "win", "ng", "dist", "norm"]].astype(str).agg("_".join, axis=1))
switch_rows = []
for pk, grp in param_key.groupby("pk"):
    cols = grp.index.tolist()
    if len(cols) < 2:
        continue
    b = mat[cols].apply(to_binary)
    # a cell "switches" if it has both 0 and 1 calls across seeds of this param set
    switched = ((b.max(axis=1) == 1) & (b.min(axis=1) == 0))
    switch_rows.append({"param_set": pk, "n_seeds": len(cols),
                        "cell_switch_rate": float(switched.mean())})
switch = pd.DataFrame(switch_rows)
if not switch.empty:
    switch.to_csv(f"{sample}_copykat_seed_switch.csv", index=False)
    print(f"{sample}: per-param seed switch-rate "
          f"(mean over param sets)={switch['cell_switch_rate'].mean():.3f}")

# ---- classification boundary: aneuploid fraction vs KS.cut ----
frac = []
for combo in mat.columns:
    b = to_binary(mat[combo]).dropna()
    frac.append({**meta.loc[combo].to_dict(),
                 "aneuploid_frac": float(b.mean()) if len(b) else np.nan})
frac = pd.DataFrame(frac)
frac.to_csv(f"{sample}_copykat_aneuploid_fraction.csv", index=False)

fig, ax = plt.subplots(figsize=(7, 5))
for (win, ng, dist, norm), g in frac.groupby(["win", "ng", "dist", "norm"]):
    g = g.sort_values("ks")
    agg = g.groupby("ks")["aneuploid_frac"].agg(["mean", "std"]).reset_index()
    ax.errorbar(agg["ks"], agg["mean"], yerr=agg["std"].fillna(0), marker="o",
                capsize=3, label=f"win{win}/ng{ng}/{dist}/norm{norm}")
ax.set_xlabel("KS.cut"); ax.set_ylabel("aneuploid fraction")
ax.set_title(f"{sample} — classification boundary (mean ± SD over seeds)")
ax.set_ylim(-0.02, 1.02); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{sample}_copykat_boundary_curve.png", dpi=150)
plt.close(fig)

# ---- pairwise ARI across combos (seed/param stability) ----
combos = list(mat.columns)
n = len(combos)
ari = np.eye(n)
for i in range(n):
    for j in range(i + 1, n):
        a, b = to_binary(mat[combos[i]]), to_binary(mat[combos[j]])
        ok = a.notna() & b.notna()
        ari[i, j] = ari[j, i] = adjusted_rand_score(a[ok], b[ok]) if ok.sum() > 1 else np.nan
ari_df = pd.DataFrame(ari, index=combos, columns=combos)
ari_df.to_csv(f"{sample}_copykat_ari.csv")
print(f"{sample}: mean pairwise ARI={np.nanmean(ari[np.triu_indices(n, 1)]):.3f}")

fig, ax = plt.subplots(figsize=(max(6, n * 0.35), max(5, n * 0.3)))
im = ax.imshow(ari, vmin=0, vmax=1, cmap="viridis")
ax.set_xticks(range(n)); ax.set_yticks(range(n))
ax.set_xticklabels(combos, rotation=90, fontsize=5)
ax.set_yticklabels(combos, fontsize=5)
ax.set_title(f"{sample} — pairwise ARI of CopyKAT calls")
fig.colorbar(im, ax=ax, fraction=0.046)
fig.tight_layout(); fig.savefig(f"{sample}_copykat_ari_heatmap.png", dpi=150)
plt.close(fig)

# ---- optional per-cell stability overlay on the reference-map UMAP ----
if mapped_h5ad:
    try:
        import scanpy as sc
        adata = sc.read_h5ad(mapped_h5ad)
        if "X_umap_ref" in adata.obsm:
            common = adata.obs_names.intersection(stab.index)
            adata = adata[common].copy()
            adata.obs["copykat_consensus"] = stab.loc[common, "consensus"].astype("category")
            adata.obs["copykat_agreement"] = stab.loc[common, "agreement"].astype(float)
            sc.set_figure_params(dpi=150, frameon=False)
            fig, axes = plt.subplots(1, 2, figsize=(16, 7))
            sc.pl.embedding(adata, basis="umap_ref", color="copykat_consensus", ax=axes[0],
                            show=False, title=f"{sample} — consensus call")
            sc.pl.embedding(adata, basis="umap_ref", color="copykat_agreement", ax=axes[1],
                            show=False, cmap="RdYlGn", title=f"{sample} — call agreement")
            fig.tight_layout(); fig.savefig(f"{sample}_copykat_stability_umap.png", dpi=150)
            plt.close(fig)
        else:
            print("mapped h5ad lacks X_umap_ref; skipping UMAP overlay")
    except Exception as e:
        print(f"UMAP overlay skipped ({e})")

print(f"Wrote {sample}_copykat_stability.csv + boundary/ARI figures")

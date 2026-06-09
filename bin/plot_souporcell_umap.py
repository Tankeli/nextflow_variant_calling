#!/usr/bin/env python3
"""
Per-patient souporcell clone overlay in the paediatric reference-map space.

Souporcell is run jointly per patient on a CB-retagged merged BAM, so clone IDs are
comparable across timepoints. This figure projects those clones onto the SHARED frozen
reference UMAP (X_umap_ref, attached by reference_mapping.py) and shows the clone x
timepoint composition. Because the reference frame is shared, no UMAP is recomputed and
the same analysis runs whether the patient has a Dx+Rel pair or a single standalone sample.

Usage:
  plot_souporcell_umap.py <patient> <clusters.tsv> <K> <samples_csv> <timepoints_csv> <mapped_h5ad>...
Output (cwd): <patient>_souporcell_umap.png/.pdf, <patient>_souporcell_composition.png/.pdf
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

patient, clusters_tsv, K = sys.argv[1], sys.argv[2], sys.argv[3]
samples    = sys.argv[4].split(",")
timepoints = sys.argv[5].split(",")
h5ads      = sys.argv[6:]
tp_of      = dict(zip(samples, timepoints))

# ---- souporcell clones (retagged "<sample>__<barcode>") ----
clones = pd.read_csv(clusters_tsv, sep="\t")
clones = clones[clones["status"] == "singlet"].copy()
clones["clone"] = "S" + clones["assignment"].astype(str)
clone_of = dict(zip(clones["barcode"].astype(str), clones["clone"]))
print(f"{patient}: {len(clone_of)} singlet souporcell assignments (k={K})")

# ---- load the patient's mapped samples (carry the shared reference UMAP) ----
parts = []
for f in h5ads:
    a = sc.read_h5ad(f)
    if "X_umap_ref" not in a.obsm:
        sys.exit(f"{f} has no X_umap_ref; reference mapping must run before this plot")
    sid = str(a.obs["sample_id"].iloc[0]) if "sample_id" in a.obs else None
    if sid not in tp_of:  # fall back to matching the filename
        sid = next((s for s in samples if s in f), sid)
    a.obs["sample_id"] = sid
    a.obs["timepoint"] = tp_of.get(sid, "NA")
    a.obs["soup_key"]  = [f"{sid}__{bc}" for bc in a.obs_names]
    # keep only what the figure needs so concat is cheap and var-agnostic
    b = ad.AnnData(obs=a.obs[["sample_id", "timepoint", "soup_key"]].copy())
    b.obsm["X_umap_ref"] = np.asarray(a.obsm["X_umap_ref"])
    parts.append(b)

adata = ad.concat(parts, index_unique=None)
adata.obs["clone"] = adata.obs["soup_key"].map(clone_of).fillna("unassigned")
n_assigned = int((adata.obs["clone"] != "unassigned").sum())
print(f"{patient}: {n_assigned}/{adata.n_obs} cells received a clone label")
for col in ("clone", "timepoint", "sample_id"):
    adata.obs[col] = adata.obs[col].astype("category")

# ---- figure 1: clone + timepoint, on the shared reference-map UMAP ----
sc.set_figure_params(dpi=150, frameon=False)
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.embedding(adata, basis="umap_ref", color="clone", ax=axes[0], show=False,
                title=f"{patient} — souporcell clones (k={K})",
                legend_loc="on data", legend_fontsize=8)
sc.pl.embedding(adata, basis="umap_ref", color="timepoint", ax=axes[1], show=False,
                title=f"{patient} — timepoint", palette="Set1")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{patient}_souporcell_umap.{ext}", bbox_inches="tight")
plt.close(fig)

# ---- figure 2: clone x timepoint composition (stacked %) ----
ct = (pd.crosstab(adata.obs["clone"], adata.obs["timepoint"], normalize="columns") * 100)
ct = ct.reindex(sorted(ct.index, key=lambda c: (c == "unassigned", c)))
fig, ax = plt.subplots(figsize=(7, 6))
bottom = np.zeros(ct.shape[1])
cmap = plt.get_cmap("tab20")
for i, clone in enumerate(ct.index):
    ax.bar(ct.columns, ct.loc[clone], bottom=bottom, label=clone, color=cmap(i % 20))
    bottom += ct.loc[clone].values
ax.set(ylabel="% of cells", title=f"{patient} — clone composition by timepoint")
ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title="clone")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{patient}_souporcell_composition.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Wrote {patient}_souporcell_umap.* and {patient}_souporcell_composition.*")

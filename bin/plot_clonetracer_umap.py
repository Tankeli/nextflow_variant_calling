#!/usr/bin/env python3
"""
Per-patient CloneTracer clone overlay in the paediatric reference-map space.

CloneTracer infers clones jointly per patient, so clone IDs are comparable across timepoints.
This figure projects the per-cell clone assignment (argmax of the clone posterior) and its
posterior confidence onto the SHARED frozen reference UMAP (X_umap_ref, attached by
reference_mapping.py) and shows the clone x timepoint composition. The reference frame is shared,
so nothing is recomputed and the same analysis runs for a Dx+Rel pair or a standalone sample.

Usage:
  plot_clonetracer_umap.py <patient> <clone_assignments.csv> <samples_csv> <timepoints_csv> <mapped_h5ad>...
Output (cwd): <patient>_clonetracer_umap.png/.pdf, <patient>_clonetracer_composition.png/.pdf
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

patient, assign_csv = sys.argv[1], sys.argv[2]
samples    = sys.argv[3].split(",")
timepoints = sys.argv[4].split(",")
h5ads      = sys.argv[5:]
tp_of      = dict(zip(samples, timepoints))

# ---- CloneTracer per-cell assignments (barcode = "<sample>__<barcode>") ----
assign = pd.read_csv(assign_csv)
if assign.empty or "barcode" not in assign.columns:
    # trivial / single-clone result (e.g. healthy control): emit placeholder figures and exit 0
    for suffix in ("umap", "composition"):
        for ext in ("png", "pdf"):
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.text(0.5, 0.5, f"{patient}: no CloneTracer clones", ha="center", va="center")
            ax.axis("off")
            fig.savefig(f"{patient}_clonetracer_{suffix}.{ext}", bbox_inches="tight")
            plt.close(fig)
    print(f"{patient}: empty clone_assignments.csv — wrote placeholder figures")
    sys.exit(0)

assign["clone"] = "C" + assign["clone"].astype(str)
clone_of = dict(zip(assign["barcode"].astype(str), assign["clone"]))
prob_of  = dict(zip(assign["barcode"].astype(str), assign["max_prob"].astype(float)))
print(f"{patient}: {len(clone_of)} CloneTracer clone assignments")

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
    a.obs["ct_key"]    = [f"{sid}__{bc}" for bc in a.obs_names]
    b = ad.AnnData(obs=a.obs[["sample_id", "timepoint", "ct_key"]].copy())
    b.obsm["X_umap_ref"] = np.asarray(a.obsm["X_umap_ref"])
    parts.append(b)

adata = ad.concat(parts, index_unique=None)
adata.obs["clone"]    = adata.obs["ct_key"].map(clone_of).fillna("unassigned")
adata.obs["max_prob"] = adata.obs["ct_key"].map(prob_of).astype(float)
n_assigned = int((adata.obs["clone"] != "unassigned").sum())
print(f"{patient}: {n_assigned}/{adata.n_obs} cells received a clone label")
for col in ("clone", "timepoint", "sample_id"):
    adata.obs[col] = adata.obs[col].astype("category")

# ---- figure 1: clone + posterior confidence + timepoint, on the shared reference-map UMAP ----
sc.set_figure_params(dpi=150, frameon=False)
fig, axes = plt.subplots(1, 3, figsize=(26, 7))
sc.pl.embedding(adata, basis="umap_ref", color="clone", ax=axes[0], show=False,
                title=f"{patient} — CloneTracer clones", legend_loc="on data", legend_fontsize=8)
sc.pl.embedding(adata, basis="umap_ref", color="max_prob", ax=axes[1], show=False,
                title=f"{patient} — clone posterior", color_map="viridis")
sc.pl.embedding(adata, basis="umap_ref", color="timepoint", ax=axes[2], show=False,
                title=f"{patient} — timepoint", palette="Set1")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{patient}_clonetracer_umap.{ext}", bbox_inches="tight")
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
    fig.savefig(f"{patient}_clonetracer_composition.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Wrote {patient}_clonetracer_umap.* and {patient}_clonetracer_composition.*")

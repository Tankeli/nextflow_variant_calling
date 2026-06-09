#!/usr/bin/env python3
"""
Per-sample reference mapping via scanpy.tl.ingest.
Ported/parameterized from DDE_23 scripts/stage3_reference_mapping_individual.py.

Maps a QC'd sample onto a reference atlas (.h5ad with a cell-type obs column), transferring
cell-type labels and a mapping-confidence score. Atlas is configurable (defaults to the DDE_32
paediatric BM atlas; Zeng BoneMarrowMap or others can be swapped in via --refmap_atlas).

Usage: reference_mapping.py <sample_h5ad> <atlas_h5ad> <sample> <celltype_key> <confidence_threshold> <n_pcs>
Outputs (cwd): <sample>_mapped.h5ad,  <sample>_celltypes.csv,  <sample>_mapping_umap.png/.pdf
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample_h5ad, atlas_h5ad, sample = sys.argv[1], sys.argv[2], sys.argv[3]
celltype_key = sys.argv[4]
conf_thresh  = float(sys.argv[5])
n_pcs        = int(sys.argv[6])
# Optional frozen reference UMAP (cell, UMAP_1, UMAP_2[, celltype]). When supplied, every
# sample is projected into THIS embedding, so X_umap_ref is one shared "paediatric
# reference-map space" across all samples — paired or standalone. Defaults to NONE.
umap_tsv = sys.argv[7] if len(sys.argv) > 7 and sys.argv[7] not in ("", "[]") else "NONE"

# ---- reference ----
ref = sc.read_h5ad(atlas_h5ad)
if celltype_key not in ref.obs:
    for alt in ("cell_type", "celltype", "CellType", "annotation"):
        if alt in ref.obs:
            print(f"celltype key '{celltype_key}' absent; using '{alt}'")
            celltype_key = alt
            break
    else:
        sys.exit(f"No cell-type column in atlas obs: {list(ref.obs.columns)[:20]}")

# Attach the frozen reference-map embedding (the shared coordinate frame) when given.
if umap_tsv not in ("NONE", "none"):
    emb = pd.read_csv(umap_tsv, sep="\t", index_col=0)
    coords = emb[["UMAP_1", "UMAP_2"]].reindex(ref.obs_names)
    covered = int(coords.notna().all(axis=1).sum())
    if covered >= 0.5 * ref.n_obs:
        keep = coords.notna().all(axis=1).values
        if keep.sum() < ref.n_obs:
            print(f"[warn] {ref.n_obs - int(keep.sum())} atlas cells absent from frozen UMAP; dropping")
            ref = ref[keep].copy()
            coords = coords.loc[ref.obs_names]
        ref.obsm["X_umap"] = coords.to_numpy()
        print(f"Loaded frozen reference UMAP for {ref.n_obs} atlas cells (shared map space)")
    else:
        print(f"[warn] frozen UMAP covers only {covered}/{ref.n_obs} atlas cells; ignoring it")

# ---- sample (GEX only, shared genes) ----
adata = sc.read_h5ad(sample_h5ad)
if "feature_types" in adata.var:
    adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
adata.X = adata.layers["counts"].copy() if "counts" in adata.layers else adata.X
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# The atlas X is RAW counts — normalise it the SAME way as the query before PCA/ingest.
# Without this the query (lognorm) is projected through a raw-count PCA, so neighbours and
# transferred labels are meaningless and cells land in the wrong place on the reference UMAP.
if ref.X.max() > 50:   # heuristic: still raw counts, not already log-normalised
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    print("Normalised reference (raw counts -> lognorm) to match the query")

shared = ref.var_names.intersection(adata.var_names)
print(f"Shared genes: {len(shared)}")
sub = adata[:, shared].copy()
ref_sub = ref[:, shared].copy()              # column subset preserves obsm['X_umap']
sc.pp.pca(ref_sub, n_comps=min(50, len(shared) - 1))
sc.pp.neighbors(ref_sub, n_pcs=n_pcs)
# sc.tl.ingest needs a fitted UMAP model (ref_sub.uns['umap']['params']) to transform the query;
# only sc.tl.umap populates it. So always run it — but when a frozen embedding is attached,
# restore those coordinates afterwards so every sample lands in the same shared frame.
frozen_umap = ref_sub.obsm["X_umap"].copy() if "X_umap" in ref_sub.obsm else None
sc.tl.umap(ref_sub, random_state=123)
if frozen_umap is not None:
    ref_sub.obsm["X_umap"] = frozen_umap
    print("Projecting onto frozen reference UMAP (shared map space)")
else:
    print("[warn] no frozen reference UMAP; per-run embedding (NOT cross-sample stable)")

sc.tl.ingest(sub, ref_sub, obs=celltype_key)
adata.obs["ref_cell_type"] = sub.obs[celltype_key].values

# Mapping confidence from UMAP-space distance to nearest reference centroid (0-1, higher=better).
try:
    cent = pd.DataFrame(ref_sub.obsm["X_umap"], index=ref_sub.obs_names)
    cent["ct"] = ref_sub.obs[celltype_key].values
    centroids = cent.groupby("ct")[[0, 1]].mean()
    emb = sub.obsm["X_umap"]
    d = np.array([np.linalg.norm(emb[i] - centroids.loc[adata.obs["ref_cell_type"].iloc[i]].values)
                  for i in range(adata.n_obs)])
    conf = 1.0 / (1.0 + d)
except Exception as e:
    print(f"Confidence calc failed ({e}); setting NaN")
    conf = np.full(adata.n_obs, np.nan)

adata.obs["mapping_confidence"] = conf
adata.obs["poorly_mapped"] = adata.obs["mapping_confidence"] < conf_thresh

out_cols = ["sample_id", "ref_cell_type", "mapping_confidence", "poorly_mapped"]
out_cols = [c for c in out_cols if c in adata.obs]
adata.obs[out_cols].to_csv(f"{sample}_celltypes.csv")

# Attach the reference-projection embedding so it travels with the checkpoint.
try:
    adata.obsm["X_umap_ref"] = np.asarray(sub.obsm["X_umap"])
except Exception as e:
    print(f"Could not attach reference UMAP ({e})")

adata.write_h5ad(f"{sample}_mapped.h5ad")
print(f"Wrote {sample}_mapped.h5ad and {sample}_celltypes.csv "
      f"({adata.obs['ref_cell_type'].nunique()} cell types)")


# ---- diagnostic mapping UMAP (reference-projection space) ----
# Cells coloured by transferred cell type, mapping confidence, and poorly-mapped flag.
# Ported from DDE_23 stage3_reference_mapping_individual.py.
def plot_mapping_umap(ad, sample_id):
    if "X_umap_ref" not in ad.obsm:
        print("No reference UMAP embedding; skipping mapping figure")
        return
    sc.set_figure_params(dpi=150, frameon=False)
    ad.obs["ref_cell_type"] = ad.obs["ref_cell_type"].astype("category")
    ad.obs["poorly_mapped"] = ad.obs["poorly_mapped"].astype(str).astype("category")
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    sc.pl.embedding(ad, basis="umap_ref", color="ref_cell_type", ax=axes[0],
                    show=False, title=f"{sample_id} — mapped cell types",
                    legend_loc="right margin", legend_fontsize=7)
    sc.pl.embedding(ad, basis="umap_ref", color="mapping_confidence", ax=axes[1],
                    show=False, title=f"{sample_id} — mapping confidence",
                    cmap="viridis")
    sc.pl.embedding(ad, basis="umap_ref", color="poorly_mapped", ax=axes[2],
                    show=False, title=f"{sample_id} — poorly mapped (<{conf_thresh})",
                    palette=["lightgray", "red"])
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{sample_id}_mapping_umap.{ext}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {sample_id}_mapping_umap.png and .pdf")


plot_mapping_umap(adata, sample)

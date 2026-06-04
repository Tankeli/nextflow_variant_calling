#!/usr/bin/env python3
"""
Per-sample reference mapping via scanpy.tl.ingest.
Ported/parameterized from DDE_23 scripts/stage3_reference_mapping_individual.py.

Maps a QC'd sample onto a reference atlas (.h5ad with a cell-type obs column), transferring
cell-type labels and a mapping-confidence score. Atlas is configurable (defaults to the DDE_32
paediatric BM atlas; Zeng BoneMarrowMap or others can be swapped in via --refmap_atlas).

Usage: reference_mapping.py <sample_h5ad> <atlas_h5ad> <sample> <celltype_key> <confidence_threshold> <n_pcs>
Outputs (cwd): <sample>_mapped.h5ad,  <sample>_celltypes.csv
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc

sample_h5ad, atlas_h5ad, sample = sys.argv[1], sys.argv[2], sys.argv[3]
celltype_key = sys.argv[4]
conf_thresh  = float(sys.argv[5])
n_pcs        = int(sys.argv[6])

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

if "X_pca" not in ref.obsm or "X_umap" not in ref.obsm:
    print("Preprocessing reference (normalize/log1p/HVG/scale/PCA/neighbors/UMAP)...")
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    sc.pp.highly_variable_genes(ref, n_top_genes=2000)
    sc.pp.scale(ref, max_value=10)
    sc.tl.pca(ref, n_comps=50, svd_solver="arpack")
    sc.pp.neighbors(ref, n_pcs=n_pcs)
    sc.tl.umap(ref)

# ---- sample (GEX only, shared genes) ----
adata = sc.read_h5ad(sample_h5ad)
if "feature_types" in adata.var:
    adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
adata.X = adata.layers["counts"].copy() if "counts" in adata.layers else adata.X
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

shared = ref.var_names.intersection(adata.var_names)
print(f"Shared genes: {len(shared)}")
sub = adata[:, shared].copy()
ref_sub = ref[:, shared].copy()
sc.pp.pca(ref_sub, n_comps=min(50, len(shared) - 1))
sc.pp.neighbors(ref_sub, n_pcs=n_pcs)
sc.tl.umap(ref_sub)

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
adata.write_h5ad(f"{sample}_mapped.h5ad")
print(f"Wrote {sample}_mapped.h5ad and {sample}_celltypes.csv "
      f"({adata.obs['ref_cell_type'].nunique()} cell types)")

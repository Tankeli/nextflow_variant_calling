#!/usr/bin/env python3
"""
Cross-reference CopyKAT aneuploidy driver genes against (a) the reference-mapping anchor genes and
(b) cell-type signatures, to test whether CopyKAT is keying on the same expression axes used to
assign cell type / lineage (which would confound the malignancy call).

Gene sets:
  - drivers       : top driver genes from copykat_drivers.py (by |Δ| aneuploid vs diploid).
  - anchor genes  : top-|loading| genes over the first n_pcs PCs of the atlas PCA, recomputed with
    the SAME normalisation as bin/reference_mapping.py (normalize_total -> log1p -> pca). For the
    scanpy.tl.ingest mapping there is no explicit HVG step, so these loadings literally define the
    projection space the query is mapped through.
  - signatures    : per-cell-type markers (sc.tl.rank_genes_groups on the atlas) + the pLSC6/LSC17
    weighted-score gene sets (mirrors bin/lsc_scoring.py).
Overlap is reported as Jaccard + hypergeometric enrichment p-value over the shared gene universe.

Usage: copykat_crossref.py <sample> <drivers_csv> <atlas_h5ad> \
                           [n_top_drivers=200] [n_pcs=30] [n_top_anchor=200] \
                           [n_markers=50] [celltype_key=auto] [marker_method=wilcoxon]
Output (cwd): <sample>_copykat_crossref_overlap.csv, <sample>_copykat_crossref_celltype.csv, figures.
"""
import sys, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample, drivers_csv, atlas_h5ad = sys.argv[1], sys.argv[2], sys.argv[3]
n_top_drivers = int(sys.argv[4]) if len(sys.argv) > 4 else 200
n_pcs         = int(sys.argv[5]) if len(sys.argv) > 5 else 30
n_top_anchor  = int(sys.argv[6]) if len(sys.argv) > 6 else 200
n_markers     = int(sys.argv[7]) if len(sys.argv) > 7 else 50
celltype_key  = sys.argv[8] if len(sys.argv) > 8 else "auto"
marker_method = sys.argv[9] if len(sys.argv) > 9 else "wilcoxon"   # 't-test' is much faster

# pLSC6 / LSC17 gene symbols (mirrors bin/lsc_scoring.py — kept inline to avoid import side effects).
PLSC6 = {"DNMT3B", "ADGRG1", "GPR56", "CD34", "SOCS2", "SPINK2", "FAM30A", "KIAA0125"}
LSC17 = {"DNMT3B", "ZBTB46", "NYNRIN", "ARHGAP22", "LAPTM4B", "MMRN1", "DPYSL3", "KIAA0125",
         "CDK6", "CPXM1", "SOCS2", "SMIM24", "EMP1", "NGFRAP1", "BEX3", "CD34", "AKR1C3",
         "GPR56", "ADGRG1"}

# ---- driver genes ----
dr = pd.read_csv(drivers_csv)
gene_col = "gene" if "gene" in dr.columns else dr.columns[0]
if "abs_diff" in dr.columns:
    dr = dr.sort_values("abs_diff", ascending=False)
drivers = list(dict.fromkeys(dr[gene_col].astype(str)))[:n_top_drivers]
drivers = {g for g in drivers if g and g.lower() != "nan"}
print(f"{sample}: {len(drivers)} driver genes (top {n_top_drivers})")

# ---- atlas-derived gene sets (anchor genes + cell-type markers) ----
# These depend only on the atlas + (n_pcs, n_top_anchor, n_markers), NOT on the sample, and the
# atlas is 1-2 GB with a slow Wilcoxon marker step. Cache to JSON in cwd so a multi-sample driver
# (jobs/run_copykat_robustness.sh) recomputes once, not per sample.
import json, hashlib, os
cache_key = hashlib.md5(f"{os.path.realpath(atlas_h5ad)}|{n_pcs}|{n_top_anchor}|{n_markers}|"
                        f"{celltype_key}|{marker_method}".encode()).hexdigest()[:10]
cache_f = f".atlas_genesets_{cache_key}.json"

if os.path.exists(cache_f):
    with open(cache_f) as fh:
        cached = json.load(fh)
    anchor    = set(cached["anchor"])
    sig_by_ct = {k: set(v) for k, v in cached["sig_by_ct"].items()}
    universe  = set(cached["universe"])
    print(f"{sample}: loaded atlas gene sets from cache {cache_f}")
else:
    ref = sc.read_h5ad(atlas_h5ad)
    if celltype_key == "auto" or celltype_key not in ref.obs:
        for alt in ("cell_type", "celltype", "CellType", "annotation", "ref_cell_type"):
            if alt in ref.obs:
                celltype_key = alt; break
    if ref.X.max() > 50:
        sc.pp.normalize_total(ref, target_sum=1e4); sc.pp.log1p(ref)
    sc.pp.pca(ref, n_comps=min(n_pcs, ref.n_vars - 1, 50))

    # anchor genes = top |loading| symbols across the leading PCs
    loadings = np.abs(ref.varm["PCs"][:, :min(n_pcs, ref.varm["PCs"].shape[1])])
    anchor_rank = pd.Series(loadings.max(axis=1), index=ref.var_names).sort_values(ascending=False)
    anchor = set(anchor_rank.head(n_top_anchor).index.astype(str))
    print(f"{sample}: {len(anchor)} anchor genes (top {n_top_anchor} over {min(n_pcs, loadings.shape[1])} PCs)")

    sig_by_ct = {}
    try:
        sc.tl.rank_genes_groups(ref, celltype_key, method=marker_method, n_genes=n_markers)
        names = pd.DataFrame(ref.uns["rank_genes_groups"]["names"])
        for ct in names.columns:
            sig_by_ct[str(ct)] = set(names[ct].astype(str).head(n_markers))
    except Exception as e:
        print(f"rank_genes_groups failed ({e}); cell-type marker overlap skipped")
    universe = set(ref.var_names.astype(str))
    with open(cache_f, "w") as fh:
        json.dump({"anchor": sorted(anchor), "universe": sorted(universe),
                   "sig_by_ct": {k: sorted(v) for k, v in sig_by_ct.items()}}, fh)
    print(f"{sample}: cached atlas gene sets -> {cache_f}")

sig_all = set().union(*sig_by_ct.values()) if sig_by_ct else set()
drivers_u = drivers & universe


def overlap(a, b):
    a, b = (a & universe), (b & universe)
    inter = a & b
    M, n, N, k = len(universe), len(a), len(b), len(inter)
    p = hypergeom.sf(k - 1, M, n, N) if (M and n and N) else np.nan
    j = k / len(a | b) if (a | b) else np.nan
    return {"set_size_a": len(a), "set_size_b": len(b), "overlap": k,
            "jaccard": j, "hypergeom_p": p, "overlap_genes": ";".join(sorted(inter))}

rows = [
    {"comparison": "drivers_vs_anchor",  **overlap(drivers_u, anchor)},
    {"comparison": "drivers_vs_markers", **overlap(drivers_u, sig_all)},
    {"comparison": "drivers_vs_pLSC6",   **overlap(drivers_u, PLSC6)},
    {"comparison": "drivers_vs_LSC17",   **overlap(drivers_u, LSC17)},
    {"comparison": "anchor_vs_markers",  **overlap(anchor, sig_all)},
]
res = pd.DataFrame(rows)
res.to_csv(f"{sample}_copykat_crossref_overlap.csv", index=False)
print(res[["comparison", "set_size_a", "set_size_b", "overlap", "jaccard", "hypergeom_p"]].to_string(index=False))

# per-cell-type driver overlap (which lineages do the drivers resemble?)
ct_rows = []
for ct, genes in sig_by_ct.items():
    o = overlap(drivers_u, genes)
    ct_rows.append({"cell_type": ct, **o})
if ct_rows:
    ct_df = pd.DataFrame(ct_rows).sort_values("hypergeom_p")
    ct_df.to_csv(f"{sample}_copykat_crossref_celltype.csv", index=False)

# ---- figures ----
fig, ax = plt.subplots(figsize=(7, 4))
plot = res.copy()
plot["-log10p"] = -np.log10(plot["hypergeom_p"].clip(lower=1e-300))
ax.bar(plot["comparison"], plot["-log10p"], color="tab:purple")
ax.axhline(-np.log10(0.05), ls="--", c="grey", label="p=0.05")
ax.set_ylabel("-log10 hypergeometric p"); ax.set_title(f"{sample} — driver-gene overlap enrichment")
ax.tick_params(axis="x", rotation=30); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{sample}_copykat_crossref_enrichment.png", dpi=150)
plt.close(fig)

if ct_rows:
    cd = pd.DataFrame(ct_rows).sort_values("overlap", ascending=True).tail(20)
    fig, ax = plt.subplots(figsize=(7, max(4, len(cd) * 0.3)))
    ax.barh(cd["cell_type"], cd["overlap"], color="teal")
    ax.set_xlabel("driver genes shared"); ax.set_title(f"{sample} — drivers vs cell-type markers")
    fig.tight_layout(); fig.savefig(f"{sample}_copykat_crossref_celltype.png", dpi=150)
    plt.close(fig)

print(f"Wrote {sample}_copykat_crossref_overlap.csv + figures")

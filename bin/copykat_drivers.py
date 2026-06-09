#!/usr/bin/env python3
"""
Which genes / genomic regions drive the CopyKAT aneuploidy call?

Splits cells into the aneuploid vs diploid groups (from a CopyKAT prediction.txt or a
copykat_stability.csv consensus), then for each genomic feature ranks it by how strongly it
separates the two groups: abs_diff = |mean CN(aneuploid) - mean CN(diploid)|, alongside the
across-cell variance. Region drivers come from the CopyKAT CNA bin matrix; gene drivers (consumed
by copykat_crossref.py) come from the gene-by-cell raw matrix when available.

Usage: copykat_drivers.py <sample> <CNA_results.txt> [gene_by_cell.txt|NONE] [calls.{csv|txt}|NONE]
Output (cwd): <sample>_copykat_drivers.csv (gene-level if available else region-level),
              <sample>_copykat_region_drivers.csv, <sample>_copykat_driver_track.png
"""
import sys, os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample   = sys.argv[1]
cna_txt  = sys.argv[2]
gbc_txt  = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "NONE") else None
calls_in = sys.argv[4] if len(sys.argv) > 4 and sys.argv[4] not in ("", "NONE") else None

META_COLS = {"chrom", "chrompos", "abspos", "chromosome_name", "start_position",
             "end_position", "ensembl_gene_id", "hgnc_symbol", "gene", "symbol", "band",
             "ensembl_transcript_id", "entrezgene"}
SYMBOL_COLS = ["hgnc_symbol", "gene", "symbol"]


def load_calls(path):
    """Map barcode -> 'aneuploid'/'diploid' from a stability csv or a copykat prediction.txt."""
    if path is None:
        return None
    sep = "," if path.endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep)
    if "consensus" in df.columns:                       # stability csv
        idx = df.columns[0]
        return dict(zip(df[idx].astype(str), df["consensus"].astype(str)))
    name = "cell.names" if "cell.names" in df.columns else df.columns[0]
    pred = "copykat.pred" if "copykat.pred" in df.columns else df.columns[-1]
    return dict(zip(df[name].astype(str), df[pred].astype(str)))


def split_cols(df):
    meta = [c for c in df.columns if c in META_COLS]
    cells = [c for c in df.columns if c not in META_COLS]
    return meta, cells


def lookup_call(call_map, cell):
    # CopyKAT's CNA matrix sanitises barcodes (- -> .); prediction.txt keeps '-'. Try both.
    c = str(cell)
    return call_map.get(c) or call_map.get(c.replace(".", "-"), "not.defined")


def rank_features(mat, call_map, label_cols):
    """mat: features x cells DataFrame (cell columns). Returns ranked driver table."""
    cells = mat.columns
    grp = pd.Series({c: lookup_call(call_map, c) for c in cells})
    aneu = cells[grp.values == "aneuploid"]
    dip  = cells[grp.values == "diploid"]
    X = mat.astype(float)
    out = label_cols.copy()
    out["mean_aneuploid"] = X[aneu].mean(axis=1) if len(aneu) else np.nan
    out["mean_diploid"]   = X[dip].mean(axis=1) if len(dip) else np.nan
    out["abs_diff"]       = (out["mean_aneuploid"] - out["mean_diploid"]).abs()
    out["variance"]       = X.var(axis=1)
    out["n_aneuploid"], out["n_diploid"] = len(aneu), len(dip)
    rank_by = "abs_diff" if (len(aneu) and len(dip)) else "variance"
    return out.sort_values(rank_by, ascending=False).reset_index(drop=True), len(aneu), len(dip)


calls = load_calls(calls_in)
if calls is None:
    sys.exit("No calls file given — pass a copykat_stability.csv or prediction.txt as arg 4")

# ---- region drivers (CNA bin matrix) ----
cna = pd.read_csv(cna_txt, sep="\t")
meta_c, cell_c = split_cols(cna)
if not cell_c:
    sys.exit(f"No cell columns parsed from {cna_txt}")
region_labels = pd.DataFrame({
    "chrom":    cna[meta_c[0]] if meta_c else np.arange(len(cna)),
    "chrompos": cna["chrompos"] if "chrompos" in cna else np.nan,
    "abspos":   cna["abspos"] if "abspos" in cna else np.arange(len(cna)),
})
region_drivers, na, nd = rank_features(cna[cell_c], calls, region_labels)
region_drivers.to_csv(f"{sample}_copykat_region_drivers.csv", index=False)
print(f"{sample}: {len(cna)} CNA bins, {na} aneuploid / {nd} diploid cells; "
      f"top region |Δ|={region_drivers['abs_diff'].iloc[0]:.3f}")

# ---- gene drivers (gene-by-cell), the set crossref consumes ----
gene_drivers = None
if gbc_txt and os.path.exists(gbc_txt):
    gbc = pd.read_csv(gbc_txt, sep="\t")
    sym = next((c for c in SYMBOL_COLS if c in gbc.columns), None)
    meta_g, cell_g = split_cols(gbc)
    if sym is None:                       # symbols may be the row index
        gbc = gbc.reset_index().rename(columns={"index": "gene"})
        sym = "gene"
        meta_g, cell_g = split_cols(gbc)
    labels = pd.DataFrame({"gene": gbc[sym].astype(str),
                           "chrom": gbc["chromosome_name"] if "chromosome_name" in gbc else np.nan,
                           "start": gbc["start_position"] if "start_position" in gbc else np.nan})
    gene_drivers, _, _ = rank_features(gbc[cell_g], calls, labels)
    gene_drivers.to_csv(f"{sample}_copykat_drivers.csv", index=False)
    print(f"{sample}: {len(gbc)} genes; top driver = "
          f"{gene_drivers['gene'].iloc[0]} (|Δ|={gene_drivers['abs_diff'].iloc[0]:.3f})")
else:
    # No gene-by-cell — emit region drivers as the primary driver table so crossref can still run
    # on whatever symbol-like column exists (degraded; warns).
    region_drivers.to_csv(f"{sample}_copykat_drivers.csv", index=False)
    print(f"{sample}: no gene-by-cell file; wrote region drivers as the primary table")

# ---- driver track along the genome (region |Δ|) ----
rd = region_drivers.dropna(subset=["abs_diff"]).copy()
if "abspos" in rd and rd["abspos"].notna().any():
    rd = rd.sort_values("abspos")
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.plot(rd["abspos"].values, rd["abs_diff"].values, lw=0.6)
    ax.set_xlabel("genome position (abspos)")
    ax.set_ylabel("|mean CN aneuploid − diploid|")
    ax.set_title(f"{sample} — CopyKAT aneuploidy driver track")
    fig.tight_layout(); fig.savefig(f"{sample}_copykat_driver_track.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {sample}_copykat_driver_track.png")

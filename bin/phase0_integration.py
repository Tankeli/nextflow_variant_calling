#!/usr/bin/env python3
"""
Phase-0 per-patient master table.

Joins every modality into one patient-stable, per-cell table so Dx and Rel cells sit on the
same axes. Reduced port of DDE_32 joint_phase0_analysis.R — the mgatk mtDNA clone axis is
dropped (mgatk is excluded from DDE_33); clones come from joint Numbat + joint souporcell.

Per cell it records: Numbat-joint clone (N<k>), souporcell clone (S<k>), CopyKAT call, mapped
cell type + confidence, pLSC6 / LSC17, and timepoint. A malignant flag is derived exactly as
the headline script (>=2 of: CopyKAT aneuploid, blast-like cell type, pLSC6 z>1 vs normal,
has a Numbat clone).

Usage:
  phase0_integration.py <patient> <samples_csv> <timepoints_csv> <numbat_out_dir> \\
      <souporcell_clusters.tsv> <celltypes_csvs_csv> <copykat_txts_csv> <lsc_csvs_csv>
Output (cwd): <patient>_cells.tsv, <patient>_clone_QC.png/.pdf
"""
import sys
import warnings
warnings.filterwarnings("ignore")

import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

patient        = sys.argv[1]
samples        = sys.argv[2].split(",")
timepoints     = sys.argv[3].split(",")
numbat_dir     = sys.argv[4]
soup_clusters  = sys.argv[5]
celltype_files = [f for f in sys.argv[6].split(",") if f]
copykat_files  = [f for f in sys.argv[7].split(",") if f]
lsc_files      = [f for f in sys.argv[8].split(",") if f]
tp_of          = dict(zip(samples, timepoints))
# longest sample names first so prefix-stripping is unambiguous (names contain "_")
samples_by_len = sorted(samples, key=len, reverse=True)

NORMAL_CT = re.compile(r"NK cell|T cell|Treg|B cell|Pre.B|Plasma|Mast cell|pDC|DC|Stromal", re.I)
BLAST_CT  = re.compile(r"HSC|MPP|HSPC|GMP|CMP|MEP|LMPP|Myelocyte|Metamyelocyte|preNeutrophil|"
                       r"preMonocyte|Erythroblast|Megakaryocyte|blast|progenitor", re.I)
ENRICHED_LSC_Z = 1.0


def split_prefix(cell, sep):
    """Return (sample, barcode) by stripping a known '<sample><sep>' prefix."""
    for s in samples_by_len:
        if cell.startswith(s + sep):
            return s, cell[len(s) + len(sep):]
    return None, cell


def sample_of_file(path):
    base = path.split("/")[-1]
    for s in samples_by_len:
        if base.startswith(s + "_"):
            return s
    return None


# ---- per-sample phenotype layers (plain barcodes) ----
def load_per_sample(files, value_cols, rename):
    frames = []
    for f in files:
        sid = sample_of_file(f)
        df = pd.read_csv(f)
        bc_col = "barcode" if "barcode" in df.columns else df.columns[0]
        df = df.rename(columns={bc_col: "barcode"})
        keep = ["barcode"] + [c for c in value_cols if c in df.columns]
        df = df[keep].rename(columns=rename)
        df["sample"] = sid
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["barcode", "sample"])


celltypes = load_per_sample(celltype_files, ["ref_cell_type", "mapping_confidence"],
                            {"ref_cell_type": "celltype"})
lsc       = load_per_sample(lsc_files, ["pLSC6_score", "LSC17_score"],
                            {"pLSC6_score": "pLSC6", "LSC17_score": "LSC17"})

# copykat prediction.txt is tab-separated with cell.names / copykat.pred
ck_frames = []
for f in copykat_files:
    sid = sample_of_file(f)
    df = pd.read_csv(f, sep="\t")
    name_col = "cell.names" if "cell.names" in df.columns else df.columns[0]
    pred_col = "copykat.pred" if "copykat.pred" in df.columns else df.columns[-1]
    ck_frames.append(pd.DataFrame({"barcode": df[name_col].astype(str),
                                   "copykat_pred": df[pred_col].astype(str), "sample": sid}))
copykat = pd.concat(ck_frames, ignore_index=True) if ck_frames else pd.DataFrame(columns=["barcode", "copykat_pred", "sample"])

# ---- joint Numbat clones (cell = "<sample>_<barcode>") ----
cps = sorted(glob.glob(f"{numbat_dir}/clone_post_*.tsv"),
             key=lambda p: int(re.sub(r".*clone_post_(\d+)\.tsv", r"\1", p)))
numbat_rows = []
if cps:
    nb = pd.read_csv(cps[-1], sep="\t")
    if "cell" in nb.columns and "clone_opt" in nb.columns:
        for cell, clone in zip(nb["cell"].astype(str), nb["clone_opt"]):
            s, bc = split_prefix(cell, "_")
            numbat_rows.append({"sample": s, "barcode": bc, "numbat_clone_joint": f"N{clone}"})
numbat = pd.DataFrame(numbat_rows) if numbat_rows else pd.DataFrame(columns=["sample", "barcode", "numbat_clone_joint"])
print(f"{patient}: {len(numbat)} Numbat-clone cells")

# ---- joint souporcell clones (barcode = "<sample>__<barcode>") ----
soup_rows = []
if soup_clusters and glob.glob(soup_clusters):
    sp = pd.read_csv(soup_clusters, sep="\t")
    sp = sp[sp["status"] == "singlet"]
    for bc, a in zip(sp["barcode"].astype(str), sp["assignment"]):
        s, core = split_prefix(bc, "__")
        soup_rows.append({"sample": s, "barcode": core, "souporcell_clone": f"S{a}"})
soup = pd.DataFrame(soup_rows) if soup_rows else pd.DataFrame(columns=["sample", "barcode", "souporcell_clone"])
print(f"{patient}: {len(soup)} souporcell-clone singlets")

# ---- merge on (sample, barcode); celltypes is the cell spine ----
cells = celltypes
for df in (lsc, copykat, numbat, soup):
    if len(df):
        cells = cells.merge(df, on=["sample", "barcode"], how="outer")
cells["patient"]   = patient
cells["timepoint"] = cells["sample"].map(tp_of)

# ---- malignant flag (mirror headline prepare_cells) ----
cells["is_normal_ct"] = cells["celltype"].fillna("").str.contains(NORMAL_CT)
cells["is_blast_ct"]  = cells["celltype"].fillna("").str.contains(BLAST_CT)
ref = cells[cells["is_normal_ct"] & cells["pLSC6"].notna()]["pLSC6"]
ref_mean, ref_sd = (ref.mean(), ref.std()) if len(ref) > 5 else (cells["pLSC6"].mean(), cells["pLSC6"].std())
cells["pLSC6_z"] = (cells["pLSC6"] - ref_mean) / ref_sd if ref_sd and ref_sd > 0 else np.nan
mscore = (
    (cells["copykat_pred"].fillna("").str.lower() == "aneuploid").astype(int)
    + cells["is_blast_ct"].astype(int)
    + (cells["pLSC6_z"] > ENRICHED_LSC_Z).fillna(False).astype(int)
    + cells["numbat_clone_joint"].notna().astype(int)
)
cells["malignant"] = mscore >= 2

cols = ["barcode", "sample", "patient", "timepoint", "celltype", "mapping_confidence",
        "copykat_pred", "numbat_clone_joint", "souporcell_clone", "pLSC6", "LSC17",
        "pLSC6_z", "malignant"]
cells = cells[[c for c in cols if c in cells.columns]]
cells.to_csv(f"{patient}_cells.tsv", sep="\t", index=False)
print(f"{patient}: wrote master table — {len(cells)} cells, {int(cells['malignant'].sum())} malignant")

# ---- clone QC sanity panel ----
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
fig.suptitle(f"{patient} — Phase-0 clone QC", fontsize=15)
for ax, col, title in ((axes[0], "numbat_clone_joint", "Numbat clone × timepoint"),
                       (axes[1], "souporcell_clone", "Souporcell clone × timepoint")):
    if col in cells and cells[col].notna().any():
        ct = pd.crosstab(cells[col], cells["timepoint"])
        ct.plot(kind="bar", stacked=True, ax=ax)
        ax.set(title=title, ylabel="cells")
    else:
        ax.text(0.5, 0.5, f"no {col}", ha="center"); ax.axis("off")
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{patient}_clone_QC.{ext}", bbox_inches="tight")
plt.close(fig)
print(f"Wrote {patient}_clone_QC.png/.pdf")

#!/usr/bin/env python3
"""
Numbat specificity / false-positive characterisation on healthy controls.

The counterpart to `copykat_stability.py` / analysis 02, but for Numbat: rather than asking how
*reproducible* a per-cell call is, it asks how *specific* Numbat is — does the joint clone caller
correctly stay silent on true-normal samples, and when it does emit clones on healthy cells, how
confident (LLR) and how lineage-confounded are those calls versus a real tumour?

For every `numbat_joint/<sample>/numbat_out/` under each results dir it:
  - classifies the sample SILENT (no `segs_consensus_*.tsv` -> Numbat found no CNV after LLR
    filtering, the "No CNV remains" log path) vs CALLED;
  - for CALLED samples reads the latest `segs_consensus_*.tsv` (per-CNV LLR + state + length) and
    `clone_post_*.tsv` (per-cell clone_opt + p_cnv) -> n_clones, clone sizes, fraction of cells with
    p_cnv>0.5, CNV LLR distribution, total CNV span (Mb);
  - (optional) overlays `reference_mapping/<sample>/<sample>_celltypes.csv` to test whether the
    aneuploid (non-reference) cells concentrate in erythroid lineage, the known BM/PBMC confound that
    drives CopyKAT's over-call (analysis 02).

Usage:
  numbat_specificity.py --out-dir DIR --results healthy=results_controls \
      --results tumour=results_patients [--celltype-subdir reference_mapping]

Outputs (in --out-dir): numbat_specificity_summary.csv + figures.
"""
import argparse, glob, os, re, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Cell-type substrings counted as erythroid (the CopyKAT/Numbat BM confound lineage).
ERYTHROID = re.compile(r"eryth|erythro|HBG|GYPA|MEP|megakaryocyte.?eryth|proeryth", re.I)
MIN_LLR_USED = 3.0      # the relaxed min_LLR this cohort was called at (numbat_min_llr)


def latest(d, pat):
    fs = sorted(glob.glob(os.path.join(d, pat)),
                key=lambda f: int(re.search(r"_(\d+)\.tsv", f).group(1))
                if re.search(r"_(\d+)\.tsv", f) else 0)
    return fs[-1] if fs else None


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--results", action="append", required=True,
                    help="label=results_dir (repeatable); label is the group, e.g. healthy=results_controls")
    ap.add_argument("--celltype-subdir", default="reference_mapping",
                    help="subdir under each results dir holding <sample>/<sample>_celltypes.csv")
    return ap.parse_args()


def celltype_overlay(res_dir, ct_subdir, sample, clone_post):
    """Return (erythroid_frac_of_aneuploid, top_celltype) for cells assigned to a non-reference
    clone, or (None, None) if no celltype file. Reference clone is the most populous clone_opt."""
    ct = glob.glob(os.path.join(res_dir, ct_subdir, sample, f"{sample}_celltypes.csv"))
    if not ct:
        return None, None
    cells = pd.read_csv(ct[0], index_col=0)
    col = "ref_cell_type" if "ref_cell_type" in cells.columns else cells.columns[-1]
    cp = clone_post.copy()
    cp["bc"] = cp["cell"].str.split("_").str[-1]            # strip <sample>_ prefix -> raw barcode
    cells = cells.rename_axis("bc").reset_index()
    cells["bc"] = cells["bc"].astype(str)
    ref_clone = cp["clone_opt"].value_counts().idxmax()
    aneu = cp[(cp["clone_opt"] != ref_clone) & (cp["p_cnv"] > 0.5)]
    merged = aneu.merge(cells[["bc", col]], on="bc", how="left").dropna(subset=[col])
    if merged.empty:
        return None, None
    ery = merged[col].str.contains(ERYTHROID).mean()
    top = merged[col].value_counts().idxmax()
    return round(float(ery), 3), top


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sources = [r.split("=", 1) for r in args.results]

    rows, seg_llr = [], {}   # seg_llr[sample] = array of per-CNV LLR
    for group, res_dir in sources:
        for nd in sorted(glob.glob(os.path.join(res_dir, "numbat_joint", "*", "numbat_out"))):
            sample = os.path.basename(os.path.dirname(nd))
            segf = latest(nd, "segs_consensus_*.tsv")
            cpf = latest(nd, "clone_post_*.tsv")
            row = dict(sample=sample, group=group, status="silent", n_cells=0, n_clones=0,
                       frac_aneuploid=0.0, n_cnv_segs=0, llr_max=np.nan, llr_median=np.nan,
                       cnv_span_mb=np.nan, cnv_states="", erythroid_frac=np.nan, top_aneu_celltype="")
            if segf and cpf:
                row["status"] = "called"
                seg = pd.read_csv(segf, sep="\t")
                cp = pd.read_csv(cpf, sep="\t")
                row["n_cells"] = len(cp)
                # reference clone = most populous; non-reference clones are the aneuploid calls
                ref_clone = cp["clone_opt"].value_counts().idxmax()
                row["n_clones"] = int(cp["clone_opt"].nunique())
                row["frac_aneuploid"] = round(
                    float(((cp["clone_opt"] != ref_clone) & (cp["p_cnv"] > 0.5)).mean()), 3)
                cnv = seg[seg["cnv_state_post"].astype(str).str.lower() != "neu"]
                row["n_cnv_segs"] = int(len(cnv))
                if len(cnv):
                    row["llr_max"] = round(float(cnv["LLR"].max()), 1)
                    row["llr_median"] = round(float(cnv["LLR"].median()), 1)
                    row["cnv_states"] = ",".join(sorted(cnv["cnv_state_post"].astype(str).unique()))
                    if "seg_length" in cnv:
                        row["cnv_span_mb"] = round(float(cnv["seg_length"].sum()) / 1e6, 1)
                    seg_llr[sample] = cnv["LLR"].to_numpy()
                ery, top = celltype_overlay(res_dir, args.celltype_subdir, sample, cp)
                row["erythroid_frac"], row["top_aneu_celltype"] = ery, top
            rows.append(row)

    summary = pd.DataFrame(rows).sort_values(["group", "status", "sample"])
    out_csv = os.path.join(args.out_dir, "numbat_specificity_summary.csv")
    summary.to_csv(out_csv, index=False)
    pd.set_option("display.width", 220, "display.max_colwidth", 40)
    print(summary.to_string(index=False))

    n_silent = (summary["status"] == "silent").sum()
    n_total = len(summary)
    print(f"\nSILENT (no CNV / correctly negative): {n_silent}/{n_total}")

    _plot_status(summary, args.out_dir)
    _plot_llr(summary, seg_llr, args.out_dir)
    if summary["erythroid_frac"].notna().any():
        _plot_erythroid(summary, args.out_dir)
    print(f"\nWrote {out_csv} + figures to {args.out_dir}")


def _palette(groups):
    cols = {"healthy": "#2c7fb8", "tumour": "#d7301f"}
    return [cols.get(g, "#888888") for g in groups]


def _plot_status(summary, out_dir):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    s = summary.sort_values(["group", "n_clones", "sample"])
    y = np.arange(len(s))
    ax.barh(y, s["n_clones"], color=_palette(s["group"]))
    for i, (_, r) in enumerate(s.iterrows()):
        lab = "silent (no CNV)" if r["status"] == "silent" else \
              f"{r['n_clones']} clones · LLR med {r['llr_median']:.0f}"
        ax.text(max(r["n_clones"], 0.05) + 0.05, i, lab, va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels(s["sample"], fontsize=8)
    ax.set_xlabel("Numbat clones called (0 = correctly silent)")
    ax.set_title("Numbat specificity: clones called per sample\n"
                 "blue = healthy control, red = tumour")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "numbat_specificity_status.png"), dpi=150)
    plt.close(fig)


def _plot_llr(summary, seg_llr, out_dir):
    called = summary[summary["status"] == "called"].sort_values(["group", "llr_median"])
    if called.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, (_, r) in enumerate(called.iterrows()):
        vals = seg_llr.get(r["sample"], np.array([]))
        col = _palette([r["group"]])[0]
        if len(vals):
            jit = np.random.RandomState(0).uniform(-0.12, 0.12, len(vals))
            ax.scatter(np.clip(vals, 1, None), i + jit, s=18, alpha=0.6, color=col)
        ax.scatter(max(np.median(vals), 1) if len(vals) else 1, i, marker="|", s=400,
                   color="black", zorder=5)
    ax.set_ylim(-0.6, len(called) - 0.4)
    ax.axvline(MIN_LLR_USED, ls="--", color="grey", lw=1)
    ax.text(MIN_LLR_USED, -0.55, f" min_LLR={MIN_LLR_USED:g} (called at)",
            fontsize=8, color="grey", va="bottom", ha="left")
    ax.axvline(20, ls=":", color="green", lw=1)
    ax.text(20, -0.55, " LLR=20 (suggested)", fontsize=8, color="green", va="bottom", ha="left")
    ax.set_xscale("log")
    ax.set_yticks(range(len(called)))
    ax.set_yticklabels(called["sample"], fontsize=8)
    ax.set_xlabel("per-CNV-segment LLR (log scale); black | = median")
    ax.set_title("CNV confidence (LLR) of called clones: healthy false-positives sit near\n"
                 "the min_LLR floor; tumour CNVs are orders of magnitude higher", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "numbat_specificity_llr.png"), dpi=150)
    plt.close(fig)


def _plot_erythroid(summary, out_dir):
    s = summary[summary["erythroid_frac"].notna()].sort_values("erythroid_frac")
    if s.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    y = np.arange(len(s))
    ax.barh(y, s["erythroid_frac"], color=_palette(s["group"]))
    ax.set_yticks(y); ax.set_yticklabels(s["sample"], fontsize=8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("erythroid fraction of Numbat-aneuploid cells")
    ax.set_title("Lineage confound check: are 'aneuploid' cells erythroid?\n"
                 "(parallels the CopyKAT erythroid driver finding, analysis 02)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "numbat_specificity_erythroid.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

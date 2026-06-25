#!/usr/bin/env python3
"""
Numbat reproducibility across the seed x min_LLR sweep (bin/numbat_sweep.R outputs).

The Numbat counterpart to `copykat_stability.py`: rather than a single-run specificity snapshot
(analysis 03), it quantifies how *reproducible* the clone call is when only the RNG seed changes, and
how the call responds to the min_LLR threshold. For each sample it reads every combo dir
`numbat_robustness/<label>/seed<>_llr<>_ent<>/` and computes:

  SEED reproducibility (at each fixed min_LLR):
    - called/silent agreement across seeds (reproducibly silent = the good negative-control case);
    - n_clones and aneuploid fraction per seed (range across seeds);
    - mean pairwise Adjusted Rand Index of per-cell clone_opt across seeds (1.0 = identical
      partition every seed) — the headline "does the clone survive re-seeding" number.

  min_LLR response (at fixed seed): n_clones / aneuploid frac / n_CNV_segs / median LLR vs min_LLR.

Usage:
  numbat_reproducibility.py --out-dir DIR --results healthy=results_controls \
      --results tumour=results_patients
Outputs: numbat_reproducibility_{seed,llr}.csv + figures.
"""
import argparse, glob, os, re, warnings, itertools
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import adjusted_rand_score

COMBO_RE = re.compile(r"seed(?P<seed>\d+)_llr(?P<llr>[0-9.]+)_ent(?P<ent>[0-9.]+)")


def latest(d, pat):
    fs = sorted(glob.glob(os.path.join(d, pat)),
                key=lambda f: int(re.search(r"_(\d+)\.tsv", f).group(1))
                if re.search(r"_(\d+)\.tsv", f) else 0)
    return fs[-1] if fs else None


def read_combo(combo_dir):
    """Return dict of per-combo metrics + the per-cell clone_opt Series (or None if silent)."""
    segf = latest(combo_dir, "segs_consensus_*.tsv")
    cpf = latest(combo_dir, "clone_post_*.tsv")
    status_f = os.path.join(combo_dir, "_sweep_status.txt")
    status = open(status_f).read().strip() if os.path.exists(status_f) else "unknown"
    out = dict(status=status, called=False, n_clones=0, frac_aneuploid=0.0,
               n_cnv_segs=0, llr_median=np.nan, clones=None)
    if segf and cpf:
        cp = pd.read_csv(cpf, sep="\t")
        seg = pd.read_csv(segf, sep="\t")
        ref = cp["clone_opt"].value_counts().idxmax()
        out.update(called=True, n_clones=int(cp["clone_opt"].nunique()),
                   frac_aneuploid=round(float(((cp["clone_opt"] != ref) &
                                               (cp["p_cnv"] > 0.5)).mean()), 3))
        cnv = seg[seg["cnv_state_post"].astype(str).str.lower() != "neu"]
        out["n_cnv_segs"] = int(len(cnv))
        if len(cnv):
            out["llr_median"] = round(float(cnv["LLR"].median()), 1)
        out["clones"] = pd.Series(cp["clone_opt"].values, index=cp["cell"].astype(str))
    return out


def mean_pairwise_ari(clone_series_list):
    """Mean ARI of clone_opt partitions over common cells, across all seed pairs."""
    series = [s for s in clone_series_list if s is not None]
    if len(series) < 2:
        return np.nan
    aris = []
    for a, b in itertools.combinations(series, 2):
        common = a.index.intersection(b.index)
        if len(common) < 50:
            continue
        aris.append(adjusted_rand_score(a.loc[common].values, b.loc[common].values))
    return round(float(np.mean(aris)), 3) if aris else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--results", action="append", required=True, help="label=results_dir")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Gather every combo: sample -> (seed, llr) -> metrics
    combos = {}   # (group, sample) -> {(seed, llr): metrics}
    for src in args.results:
        group, res_dir = src.split("=", 1)
        for cdir in sorted(glob.glob(os.path.join(res_dir, "numbat_robustness", "*", "*"))):
            m = COMBO_RE.search(os.path.basename(cdir))
            if not m:
                continue
            sample = os.path.basename(os.path.dirname(cdir))
            seed, llr = int(m["seed"]), float(m["llr"])
            combos.setdefault((group, sample), {})[(seed, llr)] = read_combo(cdir)

    seed_rows, llr_rows = [], []
    for (group, sample), cmap in sorted(combos.items()):
        llrs = sorted({k[1] for k in cmap})
        seeds = sorted({k[0] for k in cmap})
        # --- SEED reproducibility at each min_LLR ---
        for llr in llrs:
            entries = [cmap[(s, llr)] for s in seeds if (s, llr) in cmap]
            n_called = sum(e["called"] for e in entries)
            ari = mean_pairwise_ari([e["clones"] for e in entries])
            ncl = [e["n_clones"] for e in entries]
            fa = [e["frac_aneuploid"] for e in entries]
            seed_rows.append(dict(
                group=group, sample=sample, min_LLR=llr, n_seeds=len(entries),
                n_called=n_called,
                reproducibly_silent=(n_called == 0 and len(entries) > 0),
                n_clones_min=min(ncl), n_clones_max=max(ncl),
                frac_aneuploid_min=min(fa), frac_aneuploid_max=max(fa),
                seed_ari=ari))
        # --- min_LLR response at each seed ---
        for seed in seeds:
            for llr in llrs:
                if (seed, llr) not in cmap:
                    continue
                e = cmap[(seed, llr)]
                llr_rows.append(dict(group=group, sample=sample, seed=seed, min_LLR=llr,
                                     called=e["called"], n_clones=e["n_clones"],
                                     frac_aneuploid=e["frac_aneuploid"],
                                     n_cnv_segs=e["n_cnv_segs"], llr_median=e["llr_median"]))

    seed_df = pd.DataFrame(seed_rows).sort_values(["group", "sample", "min_LLR"])
    llr_df = pd.DataFrame(llr_rows).sort_values(["group", "sample", "seed", "min_LLR"])
    seed_df.to_csv(os.path.join(args.out_dir, "numbat_reproducibility_seed.csv"), index=False)
    llr_df.to_csv(os.path.join(args.out_dir, "numbat_reproducibility_llr.csv"), index=False)
    pd.set_option("display.width", 220, "display.max_colwidth", 30)
    print("=== SEED reproducibility (per sample x min_LLR) ===")
    print(seed_df.to_string(index=False))

    _plot_seed_ari(seed_df, args.out_dir)
    _plot_llr_response(llr_df, args.out_dir)
    print(f"\nWrote CSVs + figures to {args.out_dir}")


def _plot_seed_ari(seed_df, out_dir):
    # headline: seed ARI at production min_LLR=3 for samples that ever called
    d = seed_df[(seed_df["min_LLR"] == 3) & (seed_df["n_called"] > 0)].copy()
    if d.empty:
        return
    d = d.sort_values("seed_ari")
    fig, ax = plt.subplots(figsize=(8, 4.5))
    cols = ["#d7301f" if g == "tumour" else "#2c7fb8" for g in d["group"]]
    y = np.arange(len(d))
    ax.barh(y, d["seed_ari"].fillna(0), color=cols)
    for i, (_, r) in enumerate(d.iterrows()):
        txt = f"clones {r['n_clones_min']}–{r['n_clones_max']} over {int(r['n_called'])}/{int(r['n_seeds'])} seeds"
        ax.text(0.02, i, txt, va="center", fontsize=8)
    ax.axvline(1.0, ls=":", color="grey", lw=1)
    ax.set_xlim(0, 1.05)
    ax.set_yticks(y); ax.set_yticklabels(d["sample"], fontsize=8)
    ax.set_xlabel("mean pairwise seed ARI at min_LLR=3 (1.0 = identical clones every seed)")
    ax.set_title("Numbat seed reproducibility: does the clone survive re-seeding?\n"
                 "blue = healthy control, red = tumour", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "numbat_reproducibility_seed_ari.png"), dpi=150)
    plt.close(fig)


def _plot_llr_response(llr_df, out_dir):
    samples = llr_df.groupby("sample")["n_clones"].max()
    active = samples[samples > 0].index.tolist()
    if not active:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for sample in active:
        d = llr_df[llr_df["sample"] == sample].groupby("min_LLR").agg(
            n_clones=("n_clones", "mean"), frac=("frac_aneuploid", "mean")).reset_index()
        axes[0].plot(d["min_LLR"], d["n_clones"], marker="o", label=sample)
        axes[1].plot(d["min_LLR"], d["frac"], marker="o", label=sample)
    axes[0].set_xlabel("min_LLR"); axes[0].set_ylabel("mean n_clones (over seeds)")
    axes[1].set_xlabel("min_LLR"); axes[1].set_ylabel("mean aneuploid fraction")
    axes[0].set_title("Clones vs min_LLR"); axes[1].set_title("Aneuploid fraction vs min_LLR")
    axes[1].legend(fontsize=7, loc="best")
    fig.suptitle("min_LLR threshold response (active samples)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "numbat_reproducibility_llr_response.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()

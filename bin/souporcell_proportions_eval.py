#!/usr/bin/env python3
"""Score the souporcell proportion-titration sweep against ground truth.

For each (pair, minority%) run, recovers the true origin from the <sample>__ barcode prefix and
reports, over singlets:
  * ARI (souporcell vs true origin) — overall deconvolution accuracy
  * minority recall / precision — of the souporcell cluster best matching the minority sample, how
    many true-minority cells it captures (recall) and how pure it is (precision). This is the key
    read-out: can souporcell still pull out a rare population?
Then plots ARI and minority-recall vs minority fraction, one line per pair.

Usage:
  souporcell_proportions_eval.py --manifest results_soupmix_prop/manifest.csv \
      --results results_soupmix_prop --outdir <figdir> --csvout results_soupmix_prop/eval
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def adjusted_rand_index(a, b):
    ct = pd.crosstab(pd.Series(a), pd.Series(b)).to_numpy(float)
    n = ct.sum()
    if n < 2:
        return float("nan")
    comb = lambda x: (x * (x - 1) / 2).sum()
    si, sa, sb = comb(ct), comb(ct.sum(1)), comb(ct.sum(0))
    exp = sa * sb / (n * (n - 1) / 2)
    mx = 0.5 * (sa + sb)
    return 1.0 if mx == exp else (si - exp) / (mx - exp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--csvout", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    os.makedirs(a.csvout, exist_ok=True)

    man = pd.read_csv(a.manifest)
    rows = []
    for _, m in man.iterrows():
        cl_path = os.path.join(a.results, m["pair"], f"r{int(m['pct'])}", "clusters.tsv")
        if not os.path.exists(cl_path):
            print(f"MISSING {cl_path}"); continue
        cl = pd.read_csv(cl_path, sep="\t")
        cl["truth"] = cl["barcode"].str.split("__", n=1).str[0]
        s = cl[cl["status"] == "singlet"].copy()
        s["assignment"] = s["assignment"].astype(str)
        if not len(s):
            continue
        ari = adjusted_rand_index(s["truth"], s["assignment"])
        # minority cluster = souporcell cluster with most minority cells
        mins = m["minority_sample"]
        is_min = s["truth"] == mins
        n_min_true = int(is_min.sum())
        recall = precision = float("nan")
        if n_min_true:
            ct = pd.crosstab(s["assignment"], is_min)
            if True in ct.columns:
                min_clu = ct[True].idxmax()
                in_clu = s["assignment"] == min_clu
                recall = (in_clu & is_min).sum() / n_min_true
                precision = (in_clu & is_min).sum() / max(int(in_clu.sum()), 1)
        rows.append(dict(pair=m["pair"], pct=int(m["pct"]),
                         minority_sample=mins, n_cells=len(s),
                         n_minority_true=n_min_true,
                         ari=round(ari, 4),
                         minority_recall=round(float(recall), 4),
                         minority_precision=round(float(precision), 4)))
        print(f"{m['pair']} r{int(m['pct'])}%: ARI={ari:.3f} recall={recall:.3f} prec={precision:.3f} "
              f"(n_min={n_min_true})")

    df = pd.DataFrame(rows).sort_values(["pair", "pct"])
    out = os.path.join(a.csvout, "souporcell_proportions_summary.csv")
    df.to_csv(out, index=False)
    print(f"wrote {out}")

    # ---- figure: ARI + minority recall vs minority % ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for pair, g in df.groupby("pair"):
        axes[0].plot(g["pct"], g["ari"], "o-", label=pair)
        axes[1].plot(g["pct"], g["minority_recall"], "o-", label=pair)
    for ax, ttl, yl in [(axes[0], "Deconvolution accuracy vs mixing ratio", "ARI vs true origin"),
                        (axes[1], "Rare-population recovery vs mixing ratio", "minority recall")]:
        ax.set_xscale("log"); ax.set_xticks([1, 5, 10, 25, 50]); ax.set_xticklabels([1, 5, 10, 25, 50])
        ax.set_xlabel("minority sample (% of cells)"); ax.set_ylabel(yl)
        ax.set_ylim(-0.02, 1.02); ax.set_title(ttl); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    fig_out = os.path.join(a.outdir, "stats_proportions.png")
    fig.savefig(fig_out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {fig_out}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()

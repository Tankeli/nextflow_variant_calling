#!/usr/bin/env python3
"""Figures for the souporcell deconvolution-validation report.

For each mix it builds ONE shared UMAP embedding from the pooled cells and shows:
  * "before"  — coloured by TRUE sample-of-origin (the ground truth, from the <sample>__ barcode)
  * "after"   — coloured by souporcell cluster assignment, one panel per K (singlets only; doublet /
                unassigned greyed). Panel titles carry the ARI vs the true origin.
Plus cohort-level statistics figures:
  * ARI per mix at the correct K (n clusters == n true donors), controls vs patients
  * per-mix confusion-matrix heatmaps (true origin x souporcell cluster, row-normalised)

Usage:
  souporcell_mix_figures.py --samplesheet assets/test/souporcell_mix_controls.csv \
      --results results_soupmix_controls/callers/souporcell --label controls --outdir <figdir>
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 0
TAB = plt.get_cmap("tab10").colors


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


def load_mix_adata(samples_outs):
    """Pool one mix's samples into a single AnnData with obs['origin'] + souporcell-style barcodes."""
    ads = []
    for sample, outs in samples_outs:
        a = sc.read_10x_mtx(os.path.join(outs, "filtered_feature_bc_matrix"),
                            var_names="gene_symbols", cache=False)
        a.var_names_make_unique()
        a.obs["origin"] = sample
        a.obs_names = [f"{sample}__{bc}" for bc in a.obs_names]
        ads.append(a)
    adata = sc.concat(ads, join="inner") if len(ads) > 1 else ads[0]
    adata.obs_names_make_unique()
    return adata


def embed(adata):
    sc.pp.filter_genes(adata, min_cells=3)
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    a = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=30)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30)
    sc.tl.umap(a)
    adata.obsm["X_umap"] = a.obsm["X_umap"]
    return adata


def load_clusters(results, mix, k):
    p = os.path.join(results, mix, f"k{k}", "clusters.tsv")
    if not os.path.exists(p):
        return None
    cl = pd.read_csv(p, sep="\t")
    cl["label"] = np.where(cl["status"] == "singlet",
                           "C" + cl["assignment"].astype(str), cl["status"])
    return cl.set_index("barcode")


def scatter(ax, xy, labels, title, palette=None, order=None):
    cats = order or sorted(pd.unique(labels), key=str)
    pal = palette or {c: TAB[i % 10] for i, c in enumerate(cats)}
    grey = {"doublet": "0.6", "unassigned": "0.85"}
    for c in cats:
        m = labels == c
        ax.scatter(xy[m, 0], xy[m, 1], s=3, lw=0,
                   c=[grey.get(c, pal.get(c, "0.5"))], label=f"{c} (n={int(m.sum())})")
    ax.set_title(title, fontsize=10)
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(markerscale=3, fontsize=7, loc="best", frameon=False)


def fig_mix(mix, samples_outs, results, ks, outdir):
    adata = embed(load_mix_adata(samples_outs))
    xy = adata.obsm["X_umap"]
    origin = adata.obs["origin"].astype(str).values
    o_cats = sorted(pd.unique(origin))
    o_pal = {c: TAB[i % 10] for i, c in enumerate(o_cats)}

    panels = [("before — true origin", origin, o_pal, o_cats)]
    for k in ks:
        cl = load_clusters(results, mix, k)
        if cl is None:
            continue
        lab = cl.reindex(adata.obs_names)["label"].fillna("unassigned").values
        sing = ~pd.Series(lab).isin(["doublet", "unassigned"]).values
        ari = adjusted_rand_index(origin[sing], lab[sing]) if sing.sum() else float("nan")
        panels.append((f"after — souporcell k={k}  (ARI={ari:.3f})", lab, None, None))

    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(4.4 * n, 4.4))
    if n == 1:
        axes = [axes]
    for ax, (title, lab, pal, order) in zip(axes, panels):
        scatter(ax, xy, np.asarray(lab), title, palette=pal, order=order)
    fig.suptitle(mix, fontsize=12, y=1.02)
    fig.tight_layout()
    out = os.path.join(outdir, f"umap_{mix}.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--results", required=True, help="<outdir>/callers/souporcell")
    ap.add_argument("--label", required=True, help="controls|patients|bmt (used in stats fig names)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ks", default="2,3")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ks = [int(x) for x in args.ks.split(",")]
    sheet = pd.read_csv(args.samplesheet)

    for mix, g in sheet.groupby("patient", sort=False):
        print(f"== {mix} ==")
        samples_outs = list(zip(g["sample"], g["outs"]))
        try:
            fig_mix(mix, samples_outs, args.results, ks, args.outdir)
        except Exception as e:  # noqa: BLE001
            print(f"  WARN {mix}: {e}")

    # ---- stats: ARI at correct K + confusion heatmaps (read eval summary if present) ----
    summ = os.path.join(os.path.dirname(args.results.rstrip("/")), "..", "eval",
                        "souporcell_mix_eval_summary.csv")
    summ = os.path.normpath(summ)
    if os.path.exists(summ):
        df = pd.read_csv(summ)
        corr = df[df["k"] == df["n_true_samples"]].copy()
        corr = corr[~corr["is_solo"]].sort_values("ari")
        if len(corr):
            fig, ax = plt.subplots(figsize=(6, 0.5 * len(corr) + 1.5))
            ax.barh(corr["mix"] + " (k" + corr["k"].astype(str) + ")", corr["ari"], color="#3b7dd8")
            ax.set_xlim(0, 1.02); ax.set_xlabel("Adjusted Rand Index vs true origin")
            ax.set_title(f"Deconvolution accuracy at correct K — {args.label}")
            for y, v in enumerate(corr["ari"]):
                ax.text(v - 0.02, y, f"{v:.3f}", va="center", ha="right", color="white", fontsize=8)
            fig.tight_layout()
            out = os.path.join(args.outdir, f"stats_ari_{args.label}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
            print(f"  wrote {out}")

        # confusion heatmaps at correct K
        eval_dir = os.path.dirname(summ)
        for _, r in corr.iterrows():
            cpath = os.path.join(eval_dir, f"confusion_{r['mix']}_k{int(r['k'])}.csv")
            if not os.path.exists(cpath):
                continue
            cm = pd.read_csv(cpath, index_col=0)
            norm = cm.div(cm.sum(1), axis=0)
            fig, ax = plt.subplots(figsize=(0.7 * cm.shape[1] + 2, 0.5 * cm.shape[0] + 1.8))
            im = ax.imshow(norm.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
            ax.set_xticks(range(cm.shape[1])); ax.set_xticklabels(cm.columns)
            ax.set_yticks(range(cm.shape[0])); ax.set_yticklabels(cm.index)
            ax.set_xlabel("souporcell cluster"); ax.set_ylabel("true origin")
            ax.set_title(f"{r['mix']} (k{int(r['k'])})", fontsize=10)
            for i in range(cm.shape[0]):
                for j in range(cm.shape[1]):
                    ax.text(j, i, f"{int(cm.values[i, j])}", ha="center", va="center",
                            color="white" if norm.values[i, j] > 0.5 else "black", fontsize=8)
            fig.colorbar(im, ax=ax, fraction=0.046, label="row fraction")
            fig.tight_layout()
            out = os.path.join(args.outdir, f"confusion_{r['mix']}.png")
            fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
            print(f"  wrote {out}")


if __name__ == "__main__":
    main()

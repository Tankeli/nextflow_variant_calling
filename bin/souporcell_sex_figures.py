#!/usr/bin/env python3
"""Figures for Exp 3 (transplant donor/recipient sex validation).

Two kinds:
  * summary  — per souporcell cluster, mean XIST vs mean chrY score, one marker per sample, sized by
               n cells, coloured by inferred sex. One glance shows which samples are male/female and
               where a single sample carries BOTH (sex-discordant => donor/recipient of different sex).
               Reads results_soupmix_bmt/sex_eval/souporcell_sex_validation.csv (no heavy compute).
  * umap     — per sample, shared embedding coloured by: souporcell cluster (k2) | XIST | chrY score
               [| timepoint if the run pooled Dx+Rel]. Needs the Cell Ranger matrices (sbatch).

Usage:
  souporcell_sex_figures.py summary --sexcsv <...> --outdir <figdir>
  souporcell_sex_figures.py umap --samplesheet assets/test/souporcell_mix_bmt.csv \
      --results results_soupmix_bmt/callers/souporcell --outdir <figdir> [--k 2]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "NLGN4Y"]
SEX_COLOR = {"male": "#2c7fb8", "female": "#de2d26", "ambiguous(mixed)": "#756bb1",
             "low-signal": "0.6"}


def resolve(p):
    return p if os.path.isabs(p) else os.path.join(PROJ_ROOT, p)


def fig_summary(sexcsv, outdir):
    df = pd.read_csv(sexcsv)
    # one headline K per sample: smallest K available (k2) keeps it readable
    df = df[df["k"] == df.groupby("mix")["k"].transform("min")]
    fig, ax = plt.subplots(figsize=(8, 6))
    for _, r in df.iterrows():
        ax.scatter(r["mean_xist"], r["mean_y"], s=20 + r["n_cells"] / 15,
                   c=SEX_COLOR.get(r["inferred_sex"], "0.5"), alpha=0.8, edgecolor="k", lw=0.4)
        ax.annotate(f"{r['mix'].replace('_rel_solo','').replace('_dx_solo','').replace('_dxrel','+DxRel')}"
                    f"·C{r['cluster']}",
                    (r["mean_xist"], r["mean_y"]), fontsize=7,
                    xytext=(4, 3), textcoords="offset points")
    ax.set_xlabel("mean XIST (log-norm)  →  female")
    ax.set_ylabel("mean chrY score (log-norm)  →  male")
    ax.set_title("Souporcell clusters by sex-chromosome expression (BMT samples, k2)\n"
                 "point size ∝ n cells; a single sample with both a red & blue point = sex-discordant")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=s)
               for s, c in SEX_COLOR.items()]
    ax.legend(handles=handles, title="inferred sex", fontsize=8, loc="best", frameon=False)
    fig.tight_layout()
    out = os.path.join(outdir, "sex_summary_bmt.png")
    fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}")


def fig_umaps(samplesheet, results, outdir, k):
    import scanpy as sc
    sc.settings.verbosity = 0
    sheet = pd.read_csv(samplesheet)
    for mix, g in sheet.groupby("patient", sort=False):
        print(f"== {mix} ==")
        ads = []
        for _, row in g.iterrows():
            a = sc.read_10x_mtx(os.path.join(resolve(row["outs"]), "filtered_feature_bc_matrix"),
                                var_names="gene_symbols", cache=False)
            a.var_names_make_unique()
            a.obs["timepoint"] = row["timepoint"]
            a.obs_names = [f"{row['sample']}__{bc}" for bc in a.obs_names]
            ads.append(a)
        adata = sc.concat(ads, join="inner") if len(ads) > 1 else ads[0]
        adata.obs_names_make_unique()
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
        adata.obs["XIST"] = (np.asarray(adata[:, "XIST"].X.todense()).ravel()
                             if "XIST" in adata.var_names else 0.0)
        yg = [x for x in Y_GENES if x in adata.var_names]
        adata.obs["chrY"] = np.asarray(adata[:, yg].X.todense()).mean(1).ravel() if yg else 0.0
        hv = adata.copy()
        sc.pp.highly_variable_genes(hv, n_top_genes=2000)
        hv = hv[:, hv.var.highly_variable].copy()
        sc.pp.scale(hv, max_value=10); sc.tl.pca(hv, n_comps=30)
        sc.pp.neighbors(hv, n_neighbors=15, n_pcs=30); sc.tl.umap(hv)
        xy = hv.obsm["X_umap"]

        cl = pd.read_csv(os.path.join(resolve(results), mix, f"k{k}", "clusters.tsv"), sep="\t")
        cl["label"] = np.where(cl["status"] == "singlet", "C" + cl["assignment"].astype(str),
                               cl["status"])
        lab = cl.set_index("barcode").reindex(adata.obs_names)["label"].fillna("unassigned").values

        multi = g["sample"].nunique() > 1
        ncol = 4 if multi else 3
        fig, axes = plt.subplots(1, ncol, figsize=(4.3 * ncol, 4.3))
        # souporcell cluster
        cats = sorted(pd.unique(lab), key=str)
        pal = {c: plt.get_cmap("tab10").colors[i % 10] for i, c in enumerate(cats)}
        for c in cats:
            m = lab == c
            col = {"doublet": "0.6", "unassigned": "0.85"}.get(c, pal[c])
            axes[0].scatter(xy[m, 0], xy[m, 1], s=3, lw=0, c=[col], label=f"{c} ({int(m.sum())})")
        axes[0].set_title(f"souporcell k={k}"); axes[0].legend(markerscale=3, fontsize=7, frameon=False)
        # XIST + chrY continuous
        for ax, key, cmap in [(axes[1], "XIST", "Reds"), (axes[2], "chrY", "Blues")]:
            scat = ax.scatter(xy[:, 0], xy[:, 1], s=3, lw=0, c=adata.obs[key], cmap=cmap)
            ax.set_title(key + (" (→female)" if key == "XIST" else " (→male)"))
            fig.colorbar(scat, ax=ax, fraction=0.046)
        if multi:
            tp = adata.obs["timepoint"].astype(str).values
            for i, t in enumerate(sorted(pd.unique(tp))):
                m = tp == t
                axes[3].scatter(xy[m, 0], xy[m, 1], s=3, lw=0,
                                c=[plt.get_cmap("Set1").colors[i]], label=f"{t} ({int(m.sum())})")
            axes[3].set_title("timepoint"); axes[3].legend(markerscale=3, fontsize=7, frameon=False)
        for ax in axes:
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(mix, fontsize=12, y=1.02); fig.tight_layout()
        out = os.path.join(outdir, f"umap_sex_{mix}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
        print(f"  wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("summary"); s.add_argument("--sexcsv", required=True); s.add_argument("--outdir", required=True)
    u = sub.add_parser("umap")
    u.add_argument("--samplesheet", required=True); u.add_argument("--results", required=True)
    u.add_argument("--outdir", required=True); u.add_argument("--k", type=int, default=2)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    if args.cmd == "summary":
        fig_summary(args.sexcsv, args.outdir)
    else:
        fig_umaps(args.samplesheet, args.results, args.outdir, args.k)


if __name__ == "__main__":
    main()

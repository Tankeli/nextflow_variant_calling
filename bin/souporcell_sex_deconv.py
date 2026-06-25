#!/usr/bin/env python3
"""Positive control for the sex-validation method: a known male x female diagnosis mix.

A male and a female single-origin sample are pooled. This gives TWO independent ground truths:
  * barcode origin  (<sample>__ prefix) — which individual a cell really came from
  * sex             (each individual is one known sex)
We then deconvolute the pool THREE ways and check they all agree:
  1. souporcell      — genotype-based clusters (its normal mode)
  2. sex-expression  — classify each cell male/female from XIST vs chrY ALONE (no genotype)
  3. true origin     — the barcode prefix
If (1)≈(2)≈(3) then the sex-expression read-out is a valid orthogonal confirmation of souporcell's
genotype split — which is what Exp 3 relies on for transplant donor/recipient validation.

Reports per mix: souporcell ARI vs origin, sex-call accuracy vs origin, sex-call vs souporcell
agreement; and a figure (XIST-vs-chrY coloured by origin; UMAP by origin / souporcell / sex-call).

Usage:
  souporcell_sex_deconv.py --samplesheet assets/test/souporcell_mix_sex.csv \
      --results results_soupmix_sex/callers/souporcell \
      --sexmap Sample_8178:M,Sample_2977:F,Sample_2395:M,Sample_2958:F \
      --outdir <figdir> --csvout <dir> [--k 2]
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.verbosity = 0
PROJ_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "EIF1AY", "KDM5D", "USP9Y", "NLGN4Y"]
ORIG_COL = plt.get_cmap("tab10").colors
SEXC = {"male": "#2c7fb8", "female": "#de2d26", "ambiguous": "0.6"}


def resolve(p):
    return p if os.path.isabs(p) else os.path.join(PROJ_ROOT, p)


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


def classify_sex(xist, ychr, xist_thr=0.5, y_thr=0.1):
    """Per-cell sex from expression ALONE."""
    out = np.full(len(xist), "ambiguous", dtype=object)
    out[(ychr > y_thr) & (xist <= xist_thr)] = "male"
    out[(xist > xist_thr) & (ychr <= y_thr)] = "female"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--sexmap", required=True, help="Sample:M,Sample:F,...")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--csvout", required=True)
    ap.add_argument("--k", type=int, default=2)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.csvout, exist_ok=True)

    sexmap = dict(kv.split(":") for kv in args.sexmap.split(","))
    sexmap = {k: ("male" if v.upper().startswith("M") else "female") for k, v in sexmap.items()}
    sheet = pd.read_csv(args.samplesheet)
    rows = []

    for mix, g in sheet.groupby("patient", sort=False):
        print(f"== {mix} ==")
        ads = []
        for _, r in g.iterrows():
            a = sc.read_10x_mtx(os.path.join(resolve(r["outs"]), "filtered_feature_bc_matrix"),
                                var_names="gene_symbols", cache=False)
            a.var_names_make_unique()
            a.obs["origin"] = r["sample"]
            a.obs_names = [f"{r['sample']}__{bc}" for bc in a.obs_names]
            ads.append(a)
        adata = sc.concat(ads, join="inner")
        adata.obs_names_make_unique()
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
        adata.obs["XIST"] = (np.asarray(adata[:, "XIST"].X.todense()).ravel()
                             if "XIST" in adata.var_names else 0.0)
        yg = [x for x in Y_GENES if x in adata.var_names]
        adata.obs["chrY"] = np.asarray(adata[:, yg].X.todense()).mean(1).ravel() if yg else 0.0
        adata.obs["sexcall"] = classify_sex(adata.obs["XIST"].values, adata.obs["chrY"].values)
        adata.obs["true_sex"] = adata.obs["origin"].map(sexmap).astype(str)

        # embedding
        hv = adata.copy()
        sc.pp.highly_variable_genes(hv, n_top_genes=2000)
        hv = hv[:, hv.var.highly_variable].copy()
        sc.pp.scale(hv, max_value=10); sc.tl.pca(hv, n_comps=30)
        sc.pp.neighbors(hv, n_neighbors=15, n_pcs=30); sc.tl.umap(hv)
        xy = hv.obsm["X_umap"]

        # souporcell labels
        cl = pd.read_csv(os.path.join(resolve(args.results), mix, f"k{args.k}", "clusters.tsv"), sep="\t")
        cl["label"] = np.where(cl["status"] == "singlet", "C" + cl["assignment"].astype(str), cl["status"])
        adata.obs["souporcell"] = cl.set_index("barcode").reindex(adata.obs_names)["label"].fillna("unassigned").values

        # ---- metrics on confidently-classified singlets ----
        o = adata.obs
        sing = ~o["souporcell"].isin(["doublet", "unassigned"])
        sx_conf = o["sexcall"] != "ambiguous"
        # souporcell vs true origin
        ari_soup = adjusted_rand_index(o.loc[sing, "origin"], o.loc[sing, "souporcell"])
        # sex-call vs true sex (== true origin here): accuracy among confidently-sexed cells
        sex_acc = (o.loc[sx_conf, "sexcall"].values == o.loc[sx_conf, "true_sex"].values).mean()
        sex_cover = sx_conf.mean()
        # sex-call vs souporcell agreement
        both = sing & sx_conf
        ari_sex_soup = adjusted_rand_index(o.loc[both, "sexcall"], o.loc[both, "souporcell"])
        rows.append(dict(mix=mix, k=args.k, n_cells=len(o),
                         samples=";".join(f"{s}({sexmap[s][0].upper()})" for s in g["sample"]),
                         souporcell_ARI_vs_origin=round(ari_soup, 4),
                         sexcall_accuracy_vs_origin=round(float(sex_acc), 4),
                         sexcall_coverage=round(float(sex_cover), 4),
                         sexcall_vs_souporcell_ARI=round(ari_sex_soup, 4)))
        print(f"  souporcell ARI vs origin = {ari_soup:.4f} | sex-call acc vs origin = {sex_acc:.4f} "
              f"(coverage {sex_cover:.2f}) | sex-call vs souporcell ARI = {ari_sex_soup:.4f}")

        # ---- figure ----
        fig, axes = plt.subplots(1, 4, figsize=(17.5, 4.4))
        # XIST vs chrY coloured by true origin
        for i, s in enumerate(sorted(o["origin"].unique())):
            m = o["origin"] == s
            axes[0].scatter(o.loc[m, "XIST"], o.loc[m, "chrY"], s=4, lw=0,
                            c=[ORIG_COL[i]], label=f"{s} ({sexmap[s]})", alpha=0.5)
        axes[0].set_xlabel("XIST (→female)"); axes[0].set_ylabel("chrY score (→male)")
        axes[0].set_title("sex genes by true origin"); axes[0].legend(fontsize=7, markerscale=2, frameon=False)
        # UMAPs
        for ax, key, title in [(axes[1], "origin", "true origin"),
                               (axes[2], "souporcell", f"souporcell k={args.k}"),
                               (axes[3], "sexcall", "sex-expression call")]:
            cats = sorted(o[key].unique(), key=str)
            for i, c in enumerate(cats):
                m = (o[key] == c).values
                col = SEXC.get(c, {"doublet": "0.6", "unassigned": "0.85"}.get(c, ORIG_COL[i % 10]))
                ax.scatter(xy[m, 0], xy[m, 1], s=4, lw=0, c=[col], label=f"{c} ({int(m.sum())})")
            ax.set_title(title); ax.legend(fontsize=7, markerscale=2, frameon=False)
            ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"{mix}  —  souporcell ARI={ari_soup:.3f} · sex-call acc={sex_acc:.3f}", y=1.03)
        fig.tight_layout()
        out = os.path.join(args.outdir, f"sexdeconv_{mix}.png")
        fig.savefig(out, dpi=150, bbox_inches="tight"); plt.close(fig)
        print(f"  wrote {out}")

    summ = pd.DataFrame(rows)
    out = os.path.join(args.csvout, "souporcell_sex_deconv_summary.csv")
    summ.to_csv(out, index=False)
    print(f"wrote {out}")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()

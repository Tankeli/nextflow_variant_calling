#!/usr/bin/env python3
"""
Headline Dx->Rel clonal-tracing figures (Python/plotly reimplementation of DDE_32
headline_clonal_tracing.R).

  FIG 1 — Numbat-clone Sankey: each clone tracked Dx->Rel, node height = cell count, node
          colour = mean pLSC6 z (vs normal), stratum label = clone + CNV signature. Differing
          Dx/Rel clone sizes are balanced with "(lost)"/"(gained)" nodes so heights are true.
  FIG 2 — pLSC6-quartile fingerprint Sankey: Q1-Q4 thresholds from Dx malignant cells applied to
          Rel too; flows = malignant cells moving between quartiles Dx->Rel.

Both emit absolute-count and normalised-% variants. Single-timepoint patients (no Dx+Rel pair)
get a per-quartile / per-clone composition bar instead of a flow.

Usage: headline_sankey.py <patient> <cells.tsv> <numbat_out_dir>
Output (cwd): <patient>_fig1_numbat_sankey[_pct].{pdf,png}, <patient>_fig2_pLSC6_fingerprint[_pct].{pdf,png}
"""
import sys
import glob
import re
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import plotly.graph_objects as go

patient, cells_tsv, numbat_dir = sys.argv[1], sys.argv[2], sys.argv[3]
cells = pd.read_csv(cells_tsv, sep="\t")

# Map Dx/Rel labels -> ordered display names.
TP_ORDER = ["Diagnosis", "Relapse"]
tp_map = {"Dx": "Diagnosis", "Rel": "Relapse", "Diagnosis": "Diagnosis", "Relapse": "Relapse"}
cells["tp"] = cells["timepoint"].map(tp_map).fillna(cells["timepoint"])
present_tps = [t for t in TP_ORDER if t in set(cells["tp"])]
paired = len(present_tps) == 2

MAGMA = [(-1.5, "#3b0f70"), (0.5, "#8c2981"), (2.0, "#de4968"), (4.0, "#fe9f6d")]


def zcolor(z):
    if pd.isna(z):
        return "#cccccc"
    z = max(-1.5, min(4.0, z))
    for (lo, c1), (hi, c2) in zip(MAGMA[:-1], MAGMA[1:]):
        if lo <= z <= hi:
            return c1 if (z - lo) < (hi - lo) / 2 else c2
    return MAGMA[-1][1]


def save(fig, stem):
    for ext in ("pdf", "png"):
        try:
            fig.write_image(f"{stem}.{ext}", width=1100, height=750, scale=2)
        except Exception as e:
            print(f"[warn] could not write {stem}.{ext} ({e})")


# ---------------------------------------------------------------------------
# Numbat per-clone CNV signature (mean event posterior > 0.5 -> +/=/- tokens)
# ---------------------------------------------------------------------------
def numbat_cnv_labels():
    def latest(prefix):
        fs = glob.glob(f"{numbat_dir}/{prefix}_*.tsv") + glob.glob(f"{numbat_dir}/{prefix}_*.tsv.gz")
        if not fs:
            return None
        return max(fs, key=lambda p: int(re.sub(rf".*{prefix}_(\d+)\.tsv.*", r"\1", p)))
    geno_f, cp_f, segs_f = latest("geno"), latest("clone_post"), latest("segs_consensus")
    if not geno_f or not cp_f:
        return {}
    geno = pd.read_csv(geno_f, sep="\t")
    cp = pd.read_csv(cp_f, sep="\t")
    if "cell" not in geno or "clone_opt" not in cp:
        return {}
    sym = {}
    if segs_f:
        segs = pd.read_csv(segs_f, sep="\t")
        if {"seg", "cnv_state"}.issubset(segs.columns):
            for _, r in segs.iterrows():
                sym[str(r["seg"])] = {"amp": "+", "bamp": "+", "del": "−", "bdel": "−",
                                      "loh": "="}.get(str(r["cnv_state"]), "?")
    ev_cols = [c for c in geno.columns if c != "cell"]
    cp = cp[["cell", "clone_opt"]].copy()
    cp["clone"] = "N" + cp["clone_opt"].astype(str)
    j = geno.merge(cp[["cell", "clone"]], on="cell", how="inner")
    labels = {}
    for clone, g in j.groupby("clone"):
        means = g[ev_cols].mean()
        toks = [f"{sym.get(seg, '?')}{seg}" for seg, p in means.items() if p > 0.5]
        labels[clone] = " ".join(toks) if toks else "no CNV"
    return labels


# ---------------------------------------------------------------------------
# FIG 1 — Numbat clone tracking
# ---------------------------------------------------------------------------
def fig1(normalize):
    malig = cells[cells["malignant"] & cells["numbat_clone_joint"].notna()]
    if len(malig) < 30:
        print(f"{patient}: <30 malignant Numbat cells; skipping fig1")
        return
    cnv = numbat_cnv_labels()
    summ = (malig.groupby(["numbat_clone_joint", "tp"])
            .agg(n=("barcode", "size"), z=("pLSC6_z", "mean")).reset_index())
    clones = summ.groupby("numbat_clone_joint")["n"].sum().sort_values(ascending=False).index.tolist()

    def cnt(clone, tp):
        r = summ[(summ.numbat_clone_joint == clone) & (summ.tp == tp)]
        return int(r["n"].iloc[0]) if len(r) else 0

    def meanz(clone):
        r = summ[summ.numbat_clone_joint == clone]
        return float(np.nansum(r["n"] * r["z"]) / max(r["n"].sum(), 1)) if len(r) else np.nan

    if not paired:
        # single timepoint -> composition bar
        tp = present_tps[0]
        import plotly.express as px
        df = summ[summ.tp == tp]
        fig = px.bar(df, x="numbat_clone_joint", y="n", color="z",
                     color_continuous_scale="magma",
                     title=f"{patient} — Numbat clone composition ({tp})")
        save(fig, f"{patient}_fig1_numbat_sankey" + ("_pct" if normalize else ""))
        return

    labels, colors, node_index = [], [], {}

    def node(name, clone_for_color=None):
        if name not in node_index:
            node_index[name] = len(labels)
            lab = name
            if clone_for_color is not None:
                sig = cnv.get(clone_for_color, "")
                lab = f"{clone_for_color}<br>{sig}" if sig else clone_for_color
            labels.append(lab)
            colors.append(zcolor(meanz(clone_for_color)) if clone_for_color else "#dddddd")
        return node_index[name]

    src, tgt, val = [], [], []
    for clone in clones:
        nd, nr = cnt(clone, "Diagnosis"), cnt(clone, "Relapse")
        dx = node(f"Dx:{clone}", clone)
        rl = node(f"Rel:{clone}", clone)
        cont = min(nd, nr)
        if cont > 0:
            src.append(dx); tgt.append(rl); val.append(cont)
        if nd > nr:  # shrank
            lost = node("(lost)")
            src.append(dx); tgt.append(lost); val.append(nd - nr)
        elif nr > nd:  # grew
            gained = node("(gained)")
            src.append(gained); tgt.append(rl); val.append(nr - nd)

    if normalize:
        tot = {"Diagnosis": sum(cnt(c, "Diagnosis") for c in clones),
               "Relapse": sum(cnt(c, "Relapse") for c in clones)}
        val = [v / max(tot["Diagnosis"], tot["Relapse"], 1) * 100 for v in val]

    fig = go.Figure(go.Sankey(
        node=dict(label=labels, color=colors, pad=18, thickness=18,
                  line=dict(color="black", width=0.5)),
        link=dict(source=src, target=tgt, value=val, color="rgba(120,120,120,0.35)")))
    fig.update_layout(
        title=f"{patient} — Joint Numbat clone evolution{' (normalised %)' if normalize else ''}"
              f"<br><sub>Malignant cells (n={len(malig)}); node colour = mean pLSC6 z; "
              f"label = clone + CNV (+/=/−)</sub>",
        font_size=11)
    save(fig, f"{patient}_fig1_numbat_sankey" + ("_pct" if normalize else ""))


# ---------------------------------------------------------------------------
# FIG 2 — pLSC6 quartile fingerprint
# ---------------------------------------------------------------------------
def fig2(normalize):
    malig = cells[cells["malignant"] & cells["pLSC6"].notna()].copy()
    dx = malig[malig.tp == "Diagnosis"]
    if len(dx) < 20:
        print(f"{patient}: <20 Dx malignant cells with pLSC6; skipping fig2")
        return
    qs = dx["pLSC6"].quantile([0.25, 0.5, 0.75]).values
    def to_q(v):
        return "Q1" if v <= qs[0] else "Q2" if v <= qs[1] else "Q3" if v <= qs[2] else "Q4"
    malig["Q"] = malig["pLSC6"].map(to_q)
    Qs = ["Q1", "Q2", "Q3", "Q4"]

    if not paired:
        import plotly.express as px
        tp = present_tps[0]
        comp = malig[malig.tp == tp]["Q"].value_counts().reindex(Qs).fillna(0).reset_index()
        comp.columns = ["Q", "n"]
        fig = px.bar(comp, x="Q", y="n", color="Q",
                     title=f"{patient} — pLSC6 quartile composition ({tp})")
        save(fig, f"{patient}_fig2_pLSC6_fingerprint" + ("_pct" if normalize else ""))
        return

    # Sankey: Dx quartile -> Rel quartile transitions, matched per clone where possible.
    # Without per-cell tracking across timepoints, aggregate the clone-level Q distributions:
    # a clone contributes flow Dx_Qi -> Rel_Qj proportional to its Dx and Rel quartile shares.
    labels = [f"Dx {q}" for q in Qs] + [f"Rel {q}" for q in Qs]
    idx = {l: i for i, l in enumerate(labels)}
    src, tgt, val = [], [], []
    for clone, g in malig.groupby(malig["numbat_clone_joint"].fillna(malig.get("souporcell_clone", "NA"))):
        gd = g[g.tp == "Diagnosis"]["Q"].value_counts()
        gr = g[g.tp == "Relapse"]["Q"].value_counts()
        nd, nr = gd.sum(), gr.sum()
        if nd == 0 or nr == 0:
            continue
        for qi in Qs:
            for qj in Qs:
                v = (gd.get(qi, 0) / nd) * gr.get(qj, 0)  # expected cells flowing Qi->Qj
                if v > 0:
                    src.append(idx[f"Dx {qi}"]); tgt.append(idx[f"Rel {qj}"]); val.append(v)
    if not val:
        print(f"{patient}: no paired clone quartile flows; skipping fig2")
        return
    if normalize:
        s = sum(val)
        val = [v / s * 100 for v in val]
    qcol = {"Q1": "#2c7bb6", "Q2": "#abd9e9", "Q3": "#fdae61", "Q4": "#d7191c"}
    node_colors = [qcol[q] for q in Qs] * 2
    fig = go.Figure(go.Sankey(
        node=dict(label=labels, color=node_colors, pad=18, thickness=18,
                  line=dict(color="black", width=0.5)),
        link=dict(source=src, target=tgt, value=val, color="rgba(150,150,150,0.35)")))
    fig.update_layout(
        title=f"{patient} — pLSC6-quartile fingerprint Dx→Rel{' (normalised %)' if normalize else ''}"
              f"<br><sub>Q1–Q4 thresholds from Dx malignant cells; Q4 = most stem-like</sub>",
        font_size=11)
    save(fig, f"{patient}_fig2_pLSC6_fingerprint" + ("_pct" if normalize else ""))


for norm in (False, True):
    fig1(norm)
    fig2(norm)
print(f"{patient}: headline figures done (paired={paired})")

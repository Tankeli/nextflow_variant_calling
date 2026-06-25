#!/usr/bin/env python3
"""
Cohort result visualisations for the CloneTracer report. Per patient, a Dx->Rel clonal-dynamics
alluvial: clone fractions at each timepoint as stacked bars with connecting ribbons (same clone
colour), so clonal expansion/contraction across diagnosis->relapse is visible at a glance. Also a
compact cohort summary (clones per patient, leukaemia fraction Dx vs Rel).

Usage:
  clonetracer_report_viz.py <outdir> <patient>:<assignments_csv>:<dx_sample>:<rel_sample> [...]
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

outdir = sys.argv[1]
os.makedirs(outdir, exist_ok=True)
specs = []
for a in sys.argv[2:]:
    pat, csv, dx, rel = a.split(":")
    specs.append((pat, csv, dx, rel))

CMAP = plt.get_cmap("tab10")


def clone_key(c):
    s = str(c)
    return int(s) if s.isdigit() else s


def fractions(df, sample, clones):
    sub = df[df["sample"] == sample]
    n = len(sub)
    return ({c: (sub["clone"] == c).sum() / n if n else 0.0 for c in clones}, n)


def alluvial(ax, df, pat, dx, rel):
    clones = sorted(df["clone"].unique(), key=clone_key)
    fdx, ndx = fractions(df, dx, clones)
    frel, nrel = fractions(df, rel, clones)
    xL, xR, w = 0.0, 1.0, 0.18

    def cum(fr):
        y, pos = 0.0, {}
        for c in clones:
            pos[c] = (y, y + fr[c]); y += fr[c]
        return pos

    pL, pR = cum(fdx), cum(frel)
    for i, c in enumerate(clones):
        col = CMAP(i % 10)
        ax.fill_between([xL - w, xL], pL[c][0], pL[c][1], color=col, lw=0)
        ax.fill_between([xR, xR + w], pR[c][0], pR[c][1], color=col, lw=0)
        ax.add_patch(Polygon([(xL, pL[c][0]), (xL, pL[c][1]), (xR, pR[c][1]), (xR, pR[c][0])],
                             closed=True, color=col, alpha=0.32, lw=0))
        # clone label at the larger of the two bands
        if max(fdx[c], frel[c]) > 0.03:
            if fdx[c] >= frel[c]:
                ax.text(xL - w - 0.04, sum(pL[c]) / 2, f"C{c}", ha="right", va="center", fontsize=8, color=col, fontweight="bold")
            else:
                ax.text(xR + w + 0.04, sum(pR[c]) / 2, f"C{c}", ha="left", va="center", fontsize=8, color=col, fontweight="bold")
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(1.02, -0.02)
    ax.set_xticks([xL - w / 2, xR + w / 2])
    ax.set_xticklabels([f"{dx}\nDx  (n={ndx})", f"{rel}\nRel  (n={nrel})"], fontsize=9)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["0", "25", "50", "75", "100%"], fontsize=8)
    ax.set_ylabel("clone fraction")
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_title(pat, fontsize=11, fontweight="bold")
    # leukaemia (non-healthy = clone != 0) fraction per timepoint
    leuk = lambda fr: 1 - fr.get(min(clones, key=clone_key), 0.0)
    return leuk(fdx), leuk(frel)


# ---- per-patient alluvials (one combined figure) ----
dfs = {pat: pd.read_csv(csv) for pat, csv, _, _ in specs}
fig, axes = plt.subplots(1, len(specs), figsize=(4.6 * len(specs), 5.2), squeeze=False)
leuk_rows = []
for ax, (pat, csv, dx, rel) in zip(axes[0], specs):
    ldx, lrel = alluvial(ax, dfs[pat], pat, dx, rel)
    leuk_rows.append((pat, ldx, lrel))
fig.suptitle("CloneTracer — clonal dynamics across diagnosis → relapse", fontsize=13, y=1.02)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(outdir, f"cohort_clone_dynamics.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
print("wrote cohort_clone_dynamics.png/.pdf")

# ---- per-patient single alluvial (for embedding individually) ----
for pat, csv, dx, rel in specs:
    f, a = plt.subplots(figsize=(4.6, 5.2))
    alluvial(a, dfs[pat], pat, dx, rel)
    for ext in ("png", "pdf"):
        f.savefig(os.path.join(outdir, f"{pat}_clone_dynamics.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(f)
    print(f"wrote {pat}_clone_dynamics.png/.pdf")

# ---- cohort leukaemia-fraction summary ----
f, a = plt.subplots(figsize=(1.6 + 1.3 * len(specs), 4))
x = np.arange(len(leuk_rows)); bw = 0.36
a.bar(x - bw / 2, [r[1] for r in leuk_rows], bw, label="Dx", color="#5b8def")
a.bar(x + bw / 2, [r[2] for r in leuk_rows], bw, label="Rel", color="#e05c5c")
a.set_xticks(x); a.set_xticklabels([r[0] for r in leuk_rows])
a.set_ylabel("non-healthy clone fraction"); a.set_ylim(0, 1)
a.set_title("Malignant fraction Dx vs Rel"); a.legend(fontsize=8)
for s in ("top", "right"):
    a.spines[s].set_visible(False)
for ext in ("png", "pdf"):
    f.savefig(os.path.join(outdir, f"cohort_malignant_fraction.{ext}"), dpi=150, bbox_inches="tight")
plt.close(f)
print("wrote cohort_malignant_fraction.png/.pdf")
print("done")

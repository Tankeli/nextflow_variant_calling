#!/usr/bin/env python3
"""
Downstream CloneTracer figures from a <name>_out.pickle, ported to Python (matplotlib/seaborn/
networkx) from the veltenlab CloneTracer R vignette funct_clonal_analysis.R
(https://github.com/veltenlab/CloneTracer/tree/master/clonal_inference/vignettes).

Generates, from the pickle alone:
  1. <name>_elbo.png/pdf       - ELBO across candidate trees in the final heuristic round
                                 (lower = higher evidence); selected trees highlighted.
  2. <name>_trees.png/pdf      - the selected clonal hierarchy/hierarchies (nodes = clones labelled
                                 by the mutations gained; root = Healthy).
  3. <name>_heatmap.png/pdf    - single cells x mutations VAF = M/(M+N), cells ordered by assigned
                                 clone, with clone + leukaemia-probability side bars.

Pickle schema (dict): trees, tree_indices, children, parents, ELBO, clonal_prob, mutations_tree,
mutations_matrix, M, N, cell_barcode. Robust to the degenerate single-clone (healthy) case.

Usage: clonetracer_figures.py <out.pickle> <name> [outdir]
"""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import seaborn as sns

pkl_path = sys.argv[1]
name     = sys.argv[2]
outdir   = sys.argv[3] if len(sys.argv) > 3 else "."
os.makedirs(outdir, exist_ok=True)

with open(pkl_path, "rb") as f:
    P = pickle.load(f)

mut_tree = list(P["mutations_tree"])           # column order of P['trees']
trees    = P["trees"]                           # list of selected tree matrices (n_clones x n_muts)
tidx     = list(P["tree_indices"])              # indices into children/parents/potential_trees


def mito_rename(s):
    # mt_7366_C_A -> mt:7366C>A ; leave nuclear SNV / CNV names as-is
    if s.startswith("mt_") or s.startswith("mt."):
        parts = s.replace(".", "_").split("_")
        if len(parts) >= 4:
            return f"mt:{parts[1]}{parts[2]}>{parts[3]}"
    return s

mut_lab = [mito_rename(m) for m in mut_tree]


def tree_attr(sel_pos):
    """Port of get_tree_attr: returns (edges, node_label, node_order) for selected tree sel_pos."""
    tmat = np.asarray(trees[sel_pos], dtype=int)          # n_nodes x n_muts
    ti   = tidx[sel_pos]
    parents = [list(p) for p in P["parents"][ti]]          # per node: list of ALL ancestors
    n = tmat.shape[0]
    # order nodes by depth (number of ancestors)
    order = sorted(range(n), key=lambda i: len(parents[i]))
    label = {}
    edges = []
    placed = []                                            # node ids in placement order
    for rank, node in enumerate(order):
        if rank == 0:
            label[node] = "Healthy"
        else:
            # immediate parent = deepest already-placed ancestor
            anc = [a for a in parents[node] if a in placed]
            par = max(anc, key=lambda a: len(parents[a])) if anc else order[0]
            gained = [mut_lab[j] for j in range(tmat.shape[1]) if tmat[node, j] - tmat[par, j] == 1]
            label[node] = "\n".join(gained) if gained else f"clone {node}"
            edges.append((par, node))
        placed.append(node)
    return edges, label, order


def fig_trees():
    fig, axes = plt.subplots(1, len(trees), figsize=(5.5 * len(trees), 4.8), squeeze=False)
    for k in range(len(trees)):
        ax = axes[0][k]
        edges, label, order = tree_attr(k)
        # assign depth (y) and spread siblings (x)
        depth = {order[0]: 0}
        for par, node in edges:
            depth[node] = depth.get(par, 0) + 1
        bylevel = {}
        for nd, d in depth.items():
            bylevel.setdefault(d, []).append(nd)
        # spread siblings widely enough that multi-line node labels don't collide
        maxsibs = max(len(nds) for nds in bylevel.values())
        xgap = 2.6
        pos = {}
        for d, nds in bylevel.items():
            for i, nd in enumerate(sorted(nds)):
                pos[nd] = ((i - (len(nds) - 1) / 2.0) * xgap, -d)
        for par, node in edges:
            x0, y0 = pos[par]; x1, y1 = pos[node]
            ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                         mutation_scale=14, color="#555555", lw=1.4,
                                         shrinkA=26, shrinkB=26, zorder=1))
        for nd, (x, y) in pos.items():
            healthy = label[nd] == "Healthy"
            ax.scatter([x], [y], s=2200, c=("#cfd8dc" if healthy else "#e57373"),
                       edgecolors="#37474f", linewidths=1.3, zorder=2)
            ax.text(x, y, label[nd], ha="center", va="center", fontsize=7,
                    fontweight="bold", zorder=3)
        ax.set_title(f"tree {tidx[k]}" + (" (selected)" if len(trees) == 1 else ""), fontsize=10)
        xs = [p[0] for p in pos.values()]
        ax.set_xlim(min(xs) - xgap, max(xs) + xgap)
        ax.set_ylim(min(p[1] for p in pos.values()) - 0.8, 0.8)
        ax.axis("off")
    fig.suptitle(f"{name} — CloneTracer clonal hierarchy", fontsize=12, y=1.02)
    _save(fig, "trees")


def fig_elbo():
    elbo = P.get("ELBO", {})
    if not elbo:
        return
    # last-iteration ELBO per candidate tree (lower = better)
    last = {int(k): float(np.asarray(v)[-1]) for k, v in elbo.items() if len(np.asarray(v))}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
    ks = sorted(last)
    cols = ["#e57373" if k in tidx else "#90a4ae" for k in ks]
    a1.bar([str(k) for k in ks], [last[k] for k in ks], color=cols)
    a1.set_xlabel("candidate tree"); a1.set_ylabel("ELBO (last iter, lower = better)")
    a1.set_title("Model evidence per tree")
    # ELBO trajectories for selected trees
    for k in ks:
        v = np.asarray(elbo[k]); first = min(100, len(v) // 2)
        style = dict(lw=2.0) if k in tidx else dict(lw=0.8, alpha=0.4)
        a2.plot(range(first, len(v)), v[first:], label=f"tree {k}" + (" *" if k in tidx else ""), **style)
    a2.set_xlabel("SVI iteration"); a2.set_ylabel("ELBO")
    a2.set_title("ELBO trajectory"); a2.legend(fontsize=7, ncol=2)
    fig.suptitle(f"{name} — CloneTracer ELBO (selected trees highlighted)", y=1.02)
    _save(fig, "elbo")


def fig_heatmap():
    sel = 1 if 1 in P["clonal_prob"] else sorted(P["clonal_prob"])[0]
    cp = np.asarray(P["clonal_prob"][sel])                 # cells x clones
    clone_n = cp.argmax(1)                                 # 0 = healthy
    leuk = 1 - cp[:, 0]
    # VAF over the tree mutations; map M/N (mutations_matrix order) -> mutations_tree order
    mm = list(P["mutations_matrix"])
    col = [mm.index(m) for m in mut_tree]
    M = np.asarray(P["M"])[:, col]; N = np.asarray(P["N"])[:, col]
    with np.errstate(invalid="ignore", divide="ignore"):
        vaf = np.where((M + N) > 0, M / (M + N), np.nan)
    order = np.argsort(clone_n, kind="stable")
    vaf, clone_n, leuk = vaf[order], clone_n[order], leuk[order]

    fig, axes = plt.subplots(1, 3, figsize=(9, 6.5), gridspec_kw={"width_ratios": [10, 0.5, 0.5]})
    sns.heatmap(vaf, ax=axes[0], cmap="RdBu_r", vmin=0, vmax=1, cbar_kws={"label": "VAF", "shrink": .4},
                xticklabels=mut_lab, yticklabels=False)
    axes[0].set_xticklabels(mut_lab, rotation=45, ha="right", fontsize=8)
    axes[0].set_ylabel(f"{vaf.shape[0]} cells (ordered by clone)"); axes[0].set_title("Mutation VAF")
    ncl = cp.shape[1]
    sns.heatmap(clone_n.reshape(-1, 1), ax=axes[1], cmap="tab10", vmin=0, vmax=max(9, ncl - 1),
                cbar=False, xticklabels=["clone"], yticklabels=False)
    sns.heatmap(leuk.reshape(-1, 1), ax=axes[2], cmap="magma", vmin=0, vmax=1,
                cbar_kws={"label": "leukaemia prob", "shrink": .4}, xticklabels=["leuk."], yticklabels=False)
    fig.suptitle(f"{name} — single-cell VAF + CloneTracer clone assignment (tree {tidx[0]})", y=1.01)
    _save(fig, "heatmap")


def _save(fig, suffix):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(outdir, f"{name}_{suffix}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {name}_{suffix}.png/.pdf")


fig_trees()
fig_elbo()
fig_heatmap()
print("done")

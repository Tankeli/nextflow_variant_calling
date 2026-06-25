#!/usr/bin/env python3
"""
Build a CloneTracer input JSON for one patient from this pipeline's caller outputs.

CloneTracer (veltenlab) is a Bayesian *integration* model: it does not discover variants,
it consumes per-cell mutant (M) and reference (N) read-count matrices over a curated set of
mutations and infers a clonal hierarchy. Because this cohort has standard CITE-seq only (no
Optimized 10x targeted libraries), we synthesise M/N from the standard 10x outputs:

  * CNV         (mut_type 0): Numbat consensus segments -> sum of UMIs over the affected
                              region (M) vs total UMIs per cell (N); r_cnv = gain/loss.
  * nuclear SNV (mut_type 1): souporcell alt.mtx / ref.mtx at clonally-informative sites.
  * mtDNA  SNV  (mut_type 2): cellsnp-lite / mgatk AD (M) and DP-AD (N) over chrM variants.

All cells are keyed in the joint `<sample>__<barcode>` namespace (matching souporcell-prep), and
the JSON is built over the *union* of barcodes, zero-filled where a mutation type was not measured
in a cell (CloneTracer treats M=N=0 as "no information"). At least one source must yield a mutation.

Output JSON keys (see CloneTracer/clonal_inference/README.md):
  M, N, mut_type, mut_names, r_cnv, cell_barcode, class_assign, class_names
"""

import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread


def log(msg):
    print(f"[clonetracer_build_json] {msg}", file=sys.stderr)


def opt(path):
    """Treat 'NA'/'' / missing paths as absent."""
    if path is None:
        return None
    path = str(path).strip()
    if path in ("", "NA", "null", "None"):
        return None
    return path if os.path.exists(path) else None


# ----------------------------------------------------------------------------- barcodes / matrices

def read_matrix_barcodes(mtx_dir):
    """Return the list of cell barcodes from a 10x filtered matrix dir."""
    for name in ("barcodes.tsv.gz", "barcodes.tsv"):
        p = os.path.join(mtx_dir, name)
        if os.path.exists(p):
            opener = gzip.open if p.endswith(".gz") else open
            with opener(p, "rt") as fh:
                return [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    raise FileNotFoundError(f"no barcodes.tsv[.gz] in {mtx_dir}")


# ----------------------------------------------------------------------------------------- CNV (Numbat)

def parse_gtf_genes(gtf_path):
    """Map gene symbol / Ensembl id -> (chrom, start, end) from a GTF 'gene' feature."""
    genes = {}
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = f[0].replace("chr", "")
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            attrs = f[8]
            for key in ("gene_name", "gene_id"):
                marker = key + ' "'
                i = attrs.find(marker)
                if i >= 0:
                    j = attrs.find('"', i + len(marker))
                    val = attrs[i + len(marker):j]
                    genes.setdefault(val, (chrom, start, end))
                    if key == "gene_id":
                        genes.setdefault(val.split(".")[0], (chrom, start, end))
    return genes


GAIN_STATES = {"amp", "bamp", "gain", "1", "2", "amp_2"}
LOSS_STATES = {"del", "bdel", "loss", "-1"}


def cnv_mutations(numbat_dir, sample_matrices, gtf_path, cell_index):
    """
    Build CNV M/N columns from Numbat consensus segments + per-sample expression.

    Returns (names, mut_type, r_cnv, M_cols, N_cols) aligned to `cell_index`.
    """
    import scanpy as sc

    segs_path = None
    for fn in os.listdir(numbat_dir):
        if fn.startswith("segs_consensus") and fn.endswith(".tsv"):
            segs_path = os.path.join(numbat_dir, fn)
            break
    if segs_path is None:
        log("CNV: no segs_consensus_*.tsv in Numbat dir; skipping CNV mutations")
        return [], [], [], [], []

    segs = pd.read_csv(segs_path, sep="\t")
    state_col = "cnv_state_post" if "cnv_state_post" in segs.columns else "cnv_state"
    if state_col not in segs.columns:
        log(f"CNV: no cnv_state column in {segs_path}; skipping CNV mutations")
        return [], [], [], [], []

    def classify(s):
        s = str(s).lower()
        if any(g in s for g in ("amp", "gain")):
            return 1.5
        if any(l in s for l in ("del", "loss")):
            return 0.5
        return None

    segs = segs.copy()
    segs["r_cnv"] = segs[state_col].map(classify)
    segs = segs[segs["r_cnv"].notna()]
    if segs.empty:
        log("CNV: no gain/loss segments after filtering; skipping CNV mutations")
        return [], [], [], [], []

    genes = parse_gtf_genes(gtf_path) if opt(gtf_path) else {}
    if not genes:
        log("CNV: no usable GTF gene annotation; cannot map genes to segments, skipping CNV")
        return [], [], [], [], []

    # Per-cell total UMIs and per-segment region UMIs, accumulated across samples.
    n_cells = len(cell_index)
    total_umi = np.zeros(n_cells)
    seg_keys = [
        (str(r["CHROM"]).replace("chr", ""), int(r["seg_start"]), int(r["seg_end"]),
         f"CNV_chr{r['CHROM']}_{r[state_col]}", float(r["r_cnv"]))
        for _, r in segs.iterrows()
        if {"CHROM", "seg_start", "seg_end"}.issubset(segs.columns)
    ]
    if not seg_keys:
        log("CNV: segs_consensus lacks CHROM/seg_start/seg_end; skipping CNV mutations")
        return [], [], [], [], []
    region_umi = {k[3]: np.zeros(n_cells) for k in seg_keys}

    for sample, mtx_dir in sample_matrices:
        adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
        bc = [f"{sample}__{b}" for b in adata.obs_names]
        rows = np.array([cell_index.get(b, -1) for b in bc])
        keep = rows >= 0
        if not keep.any():
            continue
        rows = rows[keep]
        X = adata.X[keep]
        total_umi[rows] += np.asarray(X.sum(axis=1)).ravel()
        # gene -> column once per sample
        gpos = {g: i for i, g in enumerate(adata.var_names)}
        for chrom, start, end, name, _ in seg_keys:
            cols = [gpos[g] for g, (gc, gs, ge) in genes.items()
                    if g in gpos and gc == chrom and gs >= start and ge <= end]
            if cols:
                region_umi[name][rows] += np.asarray(X[:, cols].sum(axis=1)).ravel()

    names, mut_type, r_cnv, M_cols, N_cols = [], [], [], [], []
    for chrom, start, end, name, r in seg_keys:
        names.append(name)
        mut_type.append(0)
        r_cnv.append(r)
        M_cols.append(region_umi[name])
        N_cols.append(total_umi.copy())
    log(f"CNV: added {len(names)} segment(s): {names}")
    return names, mut_type, r_cnv, M_cols, N_cols


# ------------------------------------------------------------------------------- nuclear SNV (souporcell)

def souporcell_mutations(soup_dir, cell_index, max_snvs):
    """Nuclear-SNV M/N from souporcell alt.mtx / ref.mtx at informative sites."""
    alt_p = os.path.join(soup_dir, "alt.mtx")
    ref_p = os.path.join(soup_dir, "ref.mtx")
    clu_p = os.path.join(soup_dir, "clusters.tsv")
    vcf_p = os.path.join(soup_dir, "cluster_genotypes.vcf")
    if not (os.path.exists(alt_p) and os.path.exists(ref_p) and os.path.exists(clu_p)):
        log("SNV: missing alt.mtx/ref.mtx/clusters.tsv; skipping nuclear SNVs")
        return [], [], [], [], []

    alt = mmread(alt_p).tocsr()   # variants x cells
    ref = mmread(ref_p).tocsr()
    clusters = pd.read_csv(clu_p, sep="\t")
    bc_col = clusters.columns[0]
    cell_bc = clusters[bc_col].astype(str).tolist()   # already <sample>__<bc>
    if alt.shape[1] != len(cell_bc):
        log(f"SNV: alt.mtx cols ({alt.shape[1]}) != clusters rows ({len(cell_bc)}); skipping")
        return [], [], [], [], []

    n_var = alt.shape[0]
    # Variant identity + informativeness from cluster_genotypes.vcf (GT differs across clusters).
    var_names, informative = [], np.zeros(n_var, dtype=bool)
    if os.path.exists(vcf_p):
        vlines = [ln for ln in open(vcf_p) if not ln.startswith("#")]
        for i, ln in enumerate(vlines[:n_var]):
            f = ln.rstrip("\n").split("\t")
            var_names.append(f"SNV_{f[0]}_{f[1]}_{f[3]}_{f[4]}")
            gts = {s.split(":")[0] for s in f[9:]}
            informative[i] = len(gts) > 1
    while len(var_names) < n_var:
        var_names.append(f"SNV_var{len(var_names)}")

    idx = np.where(informative)[0]
    if idx.size == 0:
        log("SNV: no GT-differential sites in VCF; falling back to coverage-ranked sites")
        idx = np.arange(n_var)
    # rank by total coverage, keep the most-covered informative sites
    cov = np.asarray((alt + ref).sum(axis=1)).ravel()
    idx = idx[np.argsort(-cov[idx])][:max_snvs]
    if idx.size == 0:
        return [], [], [], [], []

    rows = np.array([cell_index[b] for b in cell_bc if b in cell_index])
    src = np.array([k for k, b in enumerate(cell_bc) if b in cell_index])
    n_cells = len(cell_index)

    names, mut_type, r_cnv, M_cols, N_cols = [], [], [], [], []
    for v in idx:
        M = np.zeros(n_cells); N = np.zeros(n_cells)
        M[rows] = np.asarray(alt[v].todense()).ravel()[src]
        N[rows] = np.asarray(ref[v].todense()).ravel()[src]
        names.append(var_names[v]); mut_type.append(1); r_cnv.append(0.0)
        M_cols.append(M); N_cols.append(N)
    log(f"SNV: added {len(names)} nuclear SNV site(s) (cap {max_snvs})")
    return names, mut_type, r_cnv, M_cols, N_cols


# ------------------------------------------------------------------------------------- mtDNA (cellsnp/mgatk)

def mtdna_mutations(mt_dirs, cell_index, min_cells, max_sites):
    """mtDNA-SNV M/N from cellsnp-lite output dirs (AD = M, DP-AD = N)."""
    names, mut_type, r_cnv, M_cols, N_cols = [], [], [], [], []
    n_cells = len(cell_index)
    # accumulate per-site columns across samples; key sites by chrom:pos:ref:alt
    site_M, site_N, site_cells = {}, {}, {}

    for sample, d in mt_dirs:
        ad_p = os.path.join(d, "cellSNP.tag.AD.mtx")
        dp_p = os.path.join(d, "cellSNP.tag.DP.mtx")
        vcf_p = os.path.join(d, "cellSNP.base.vcf")
        bc_p = os.path.join(d, "cellSNP.samples.tsv")
        if not (os.path.exists(ad_p) and os.path.exists(dp_p) and os.path.exists(bc_p)):
            log(f"mtDNA: incomplete cellSNP output in {d}; skipping this sample")
            continue
        AD = mmread(ad_p).tocsr()   # sites x cells
        DP = mmread(dp_p).tocsr()
        bc = [f"{sample}__{ln.strip()}" for ln in open(bc_p) if ln.strip()]
        site_ids = []
        if os.path.exists(vcf_p):
            for ln in open(vcf_p):
                if ln.startswith("#"):
                    continue
                f = ln.rstrip("\n").split("\t")
                site_ids.append(f"mt_{f[1]}_{f[3]}_{f[4]}")
        while len(site_ids) < AD.shape[0]:
            site_ids.append(f"mt_site{len(site_ids)}")

        rows = np.array([cell_index[b] for b in bc if b in cell_index])
        src = np.array([k for k, b in enumerate(bc) if b in cell_index])
        if rows.size == 0:
            continue
        for s in range(AD.shape[0]):
            sid = site_ids[s]
            if sid not in site_M:
                site_M[sid] = np.zeros(n_cells); site_N[sid] = np.zeros(n_cells)
                site_cells[sid] = 0
            ad = np.asarray(AD[s].todense()).ravel()[src]
            dp = np.asarray(DP[s].todense()).ravel()[src]
            site_M[sid][rows] += ad
            site_N[sid][rows] += (dp - ad)
            site_cells[sid] += int((dp > 0).sum())

    # keep informative, well-covered sites
    sites = [s for s in site_M if site_cells[s] >= min_cells]
    sites = sorted(sites, key=lambda s: -site_cells[s])[:max_sites]
    for s in sites:
        names.append(s); mut_type.append(2); r_cnv.append(0.0)
        M_cols.append(site_M[s]); N_cols.append(site_N[s])
    log(f"mtDNA: added {len(names)} site(s) (>= {min_cells} cells covered, cap {max_sites})")
    return names, mut_type, r_cnv, M_cols, N_cols


# ------------------------------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patient", required=True)
    ap.add_argument("--samples", required=True, help="comma-separated sample ids, ordered")
    ap.add_argument("--timepoints", required=True, help="comma-separated Dx/Rel, same order")
    ap.add_argument("--matrices", required=True, help="comma-separated 10x matrix dirs, same order")
    ap.add_argument("--numbat-dir", default="NA")
    ap.add_argument("--souporcell-dir", default="NA")
    ap.add_argument("--mtdna-dirs", default="NA", help="comma-separated cellsnp dirs, same sample order")
    ap.add_argument("--gtf", default="NA")
    ap.add_argument("--max-snvs", type=int, default=50)
    ap.add_argument("--mtdna-min-cells", type=int, default=10)
    ap.add_argument("--mtdna-max-sites", type=int, default=50)
    ap.add_argument("--max-total-muts", type=int, default=6,
                    help="hard cap on total mutations fed to the model (priority CNV>SNV>mtDNA). "
                         "CloneTracer's tree search explodes super-exponentially in mutation count, "
                         "so keep this small (~3-6); 0 disables the cap.")
    ap.add_argument("--pseudobulk", action="store_true",
                    help="synthesise per-class bulk_M/bulk_N (column sums per timepoint) so the "
                         "multi-sample model (-s) can run without external bulk exome/karyotype data")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    samples = args.samples.split(",")
    timepoints = args.timepoints.split(",")
    matrices = args.matrices.split(",")
    assert len(samples) == len(matrices) == len(timepoints), "samples/timepoints/matrices length mismatch"
    sample_matrices = list(zip(samples, matrices))

    # ---- cell universe (union of all sample barcodes, joint namespace) ----
    all_bc = []
    sample_of = {}
    for sample, mtx_dir in sample_matrices:
        for b in read_matrix_barcodes(mtx_dir):
            tag = f"{sample}__{b}"
            if tag not in sample_of:
                sample_of[tag] = sample
                all_bc.append(tag)
    cell_index = {b: i for i, b in enumerate(all_bc)}
    log(f"{args.patient}: {len(all_bc)} cells across {len(samples)} sample(s)")

    names, mut_type, r_cnv, M_cols, N_cols = [], [], [], [], []

    if opt(args.numbat_dir):
        try:
            r = cnv_mutations(opt(args.numbat_dir), sample_matrices, args.gtf, cell_index)
            for dst, s in zip((names, mut_type, r_cnv, M_cols, N_cols), r):
                dst.extend(s)
        except Exception as e:  # never let one source sink the run
            log(f"CNV: skipped due to error: {e}")

    if opt(args.souporcell_dir):
        try:
            r = souporcell_mutations(opt(args.souporcell_dir), cell_index, args.max_snvs)
            for dst, s in zip((names, mut_type, r_cnv, M_cols, N_cols), r):
                dst.extend(s)
        except Exception as e:
            log(f"SNV: skipped due to error: {e}")

    mt_raw = (args.mtdna_dirs or "").strip()
    mt_dirs = [d for d in mt_raw.split(",") if opt(d)] if mt_raw not in ("", "NA", "null", "None") else []
    if mt_dirs:
        mt_pairs = list(zip(samples, mt_dirs))
        try:
            r = mtdna_mutations(mt_pairs, cell_index, args.mtdna_min_cells, args.mtdna_max_sites)
            for dst, s in zip((names, mut_type, r_cnv, M_cols, N_cols), r):
                dst.extend(s)
        except Exception as e:
            log(f"mtDNA: skipped due to error: {e}")

    if not names:
        log("ERROR: no mutations from any source (CNV/SNV/mtDNA). Cannot build CloneTracer input.")
        sys.exit(2)

    # ---- hard total-mutation cap (priority CNV > SNV > mtDNA) -------------------------------------
    # CloneTracer's infer_hierarchy re-fits SVI to every candidate tree and the candidate count grows
    # super-exponentially with mutation count, so an unbounded set never finishes. Sources are already
    # best-first within type; a stable sort by mut_type keeps all CNVs, then fills with SNVs, then mtDNA.
    if args.max_total_muts and len(names) > args.max_total_muts:
        keep = sorted(range(len(names)), key=lambda i: mut_type[i])[:args.max_total_muts]
        keep_set = set(keep)
        dropped = len(names) - len(keep)
        names    = [names[i]    for i in range(len(names)) if i in keep_set]
        mut_type = [mut_type[i] for i in range(len(mut_type)) if i in keep_set]
        r_cnv    = [r_cnv[i]    for i in range(len(r_cnv)) if i in keep_set]
        M_cols   = [M_cols[i]   for i in range(len(M_cols)) if i in keep_set]
        N_cols   = [N_cols[i]   for i in range(len(N_cols)) if i in keep_set]
        log(f"total-cap: kept {len(names)}/{len(names)+dropped} mutations "
            f"(cap {args.max_total_muts}, priority CNV>SNV>mtDNA): {names}")

    # cells x mutations
    M = np.vstack(M_cols).T
    N = np.vstack(N_cols).T

    # drop cells with zero information across all mutations (M+N == 0 everywhere)
    informative = (M + N).sum(axis=1) > 0
    M, N = M[informative], N[informative]
    kept_bc = [b for b, k in zip(all_bc, informative) if k]
    tp_map = {s: t for s, t in zip(samples, timepoints)}
    class_names = list(dict.fromkeys(timepoints))   # preserve Dx,Rel order, unique
    class_idx = {c: i for i, c in enumerate(class_names)}
    class_assign = [class_idx[tp_map[sample_of[b]]] for b in kept_bc]

    out = {
        "M": M.astype(int).tolist(),
        "N": N.astype(int).tolist(),
        "mut_type": [int(x) for x in mut_type],
        "mut_names": names,
        "r_cnv": [float(x) for x in r_cnv],
        "cell_barcode": kept_bc,
    }
    # only add multi-sample class info when >1 timepoint is present
    if len(class_names) > 1:
        out["class_assign"] = [int(x) for x in class_assign]
        out["class_names"] = class_names
        # Optional pseudobulk: per-class column sums of M / N. CloneTracer's multi-sample model
        # (-s) needs bulk_M/bulk_N (a list of per-class vectors); we have no exome/karyotype here,
        # so synthesise it from the single-cell data. Without this, the model must run pooled (no -s).
        if args.pseudobulk:
            ca = np.asarray(class_assign)
            bulk_M, bulk_N = [], []
            for c in range(len(class_names)):
                sel = ca == c
                bulk_M.append(M[sel].sum(axis=0).astype(int).tolist())
                bulk_N.append(N[sel].sum(axis=0).astype(int).tolist())
            out["bulk_M"] = bulk_M
            out["bulk_N"] = bulk_N
            log(f"pseudobulk: added per-class bulk_M/bulk_N for {class_names}")

    with open(args.output, "w") as fh:
        json.dump(out, fh)
    log(f"wrote {args.output}: {M.shape[0]} cells x {M.shape[1]} mutations "
        f"(CNV={mut_type.count(0)}, SNV={mut_type.count(1)}, mtDNA={mut_type.count(2)}); "
        f"classes={class_names if len(class_names) > 1 else 'single-sample'}")


if __name__ == "__main__":
    main()

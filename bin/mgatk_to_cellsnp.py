#!/usr/bin/env python3
"""
Convert mgatk per-base count output into cellSNP-style matrices so the rest of the
CloneTracer mtDNA path is method-agnostic. Opt-in path (default is cellsnp-lite).

Reads an mgatk `final/` directory containing the sparse per-base files
  <sample>.{A,C,G,T}.txt[.gz]   (rows: "pos,cell_idx,count_fw,count_rev")
  <sample>_refAllele.txt        (pos, ref base)   [or chrM_refAllele.txt]
  <sample>.barcodes / .depthTable / .qc ...
and writes, into the output dir:
  cellSNP.tag.AD.mtx  (alt depth, sites x cells)
  cellSNP.tag.DP.mtx  (total depth, sites x cells)
  cellSNP.base.vcf    (one line per called site: CHROM POS . REF ALT ...)
  cellSNP.samples.tsv (cell barcodes, one per line)

A site is called per position by taking the most common non-reference base across cells
(simple majority-alt). This mirrors what cellsnp-lite emits closely enough for the
CloneTracer M/N builder (mtdna_mutations() in clonetracer_build_json.py).

NOTE: mgatk output layout varies across versions; this handles the common `tenx` layout.
If the expected files are absent it exits non-zero so the caller can fall back to cellsnp.
"""
import glob
import gzip
import os
import sys

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

BASES = ["A", "C", "G", "T"]


def _open(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def find(final_dir, suffix):
    hits = glob.glob(os.path.join(final_dir, f"*{suffix}")) + \
           glob.glob(os.path.join(final_dir, f"*{suffix}.gz"))
    return hits[0] if hits else None


def main():
    final_dir, out_dir = sys.argv[1], sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    # barcodes (cell order)
    bc_file = find(final_dir, ".barcodes.txt") or find(final_dir, "barcodes.tsv")
    if not bc_file:
        sys.exit("mgatk_to_cellsnp: no barcodes file found")
    with _open(bc_file) as fh:
        barcodes = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    ncells = len(barcodes)

    # ref allele per position
    ref_file = find(final_dir, "_refAllele.txt")
    refbase = {}
    if ref_file:
        with _open(ref_file) as fh:
            for ln in fh:
                p = ln.rstrip("\n").split("\t")
                if len(p) >= 2:
                    try:
                        refbase[int(p[0])] = p[1].upper()
                    except ValueError:
                        pass

    # per-base count matrices: counts[base][pos] = np.array(ncells)
    counts = {b: {} for b in BASES}
    maxpos = 0
    for b in BASES:
        f = find(final_dir, f".{b}.txt")
        if not f:
            sys.exit(f"mgatk_to_cellsnp: missing per-base file for {b}")
        with _open(f) as fh:
            for ln in fh:
                parts = ln.rstrip("\n").split(",")
                if len(parts) < 3:
                    continue
                pos = int(parts[0]); cell = int(parts[1]) - 1
                cnt = sum(int(float(x)) for x in parts[2:])
                if cell < 0 or cell >= ncells:
                    continue
                counts[b].setdefault(pos, np.zeros(ncells, dtype=np.int32))[cell] += cnt
                maxpos = max(maxpos, pos)

    positions = sorted({p for b in BASES for p in counts[b]})
    AD_rows, DP_rows, vcf = [], [], []
    for pos in positions:
        per_base = {b: counts[b].get(pos, np.zeros(ncells, dtype=np.int32)) for b in BASES}
        dp = sum(per_base.values())
        ref = refbase.get(pos)
        alt_candidates = [b for b in BASES if b != ref]
        # choose alt = non-ref base with the highest total count
        alt = max(alt_candidates, key=lambda b: int(per_base[b].sum())) if alt_candidates else "N"
        ad = per_base.get(alt, np.zeros(ncells, dtype=np.int32))
        if int(ad.sum()) == 0:
            continue
        AD_rows.append(ad)
        DP_rows.append(dp)
        vcf.append((pos, ref or "N", alt))

    if not AD_rows:
        sys.exit("mgatk_to_cellsnp: no alt-bearing sites found")

    AD = csr_matrix(np.vstack(AD_rows))
    DP = csr_matrix(np.vstack(DP_rows))
    mmwrite(os.path.join(out_dir, "cellSNP.tag.AD.mtx"), AD)
    mmwrite(os.path.join(out_dir, "cellSNP.tag.DP.mtx"), DP)
    with open(os.path.join(out_dir, "cellSNP.samples.tsv"), "w") as fh:
        fh.write("\n".join(barcodes) + "\n")
    with open(os.path.join(out_dir, "cellSNP.base.vcf"), "w") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref, alt in vcf:
            fh.write(f"chrM\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n")
    print(f"mgatk_to_cellsnp: wrote {len(vcf)} sites x {ncells} cells to {out_dir}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

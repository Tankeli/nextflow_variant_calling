#!/usr/bin/env python3
"""
Downsample a 10x filtered_feature_bc_matrix to a fraction of its CELLS (columns), keeping all genes.

Used by the CopyKAT cell-number confounder test: to check whether the per-sample ARI differences in
the robustness sweep are driven by cell count, we hold a sample fixed and vary only N. One fixed
cell subset per (sample, fraction) — selected with `seed` — so the downstream CopyKAT seed replicates
all see the SAME cells and the ARI we measure reflects N, not a re-draw of which cells.

Usage: downsample_10x.py <matrix_dir> <fraction> <seed> <out_dir>
Writes <out_dir>/{matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz} (cellranger layout).
"""
import sys, os, gzip, subprocess
import numpy as np
from scipy.io import mmread, mmwrite

matrix_dir, fraction, seed, out_dir = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
os.makedirs(out_dir, exist_ok=True)

mtx = mmread(os.path.join(matrix_dir, "matrix.mtx.gz")).tocsc()      # genes x cells
with gzip.open(os.path.join(matrix_dir, "barcodes.tsv.gz"), "rt") as fh:
    barcodes = [l.rstrip("\n") for l in fh]
n = len(barcodes)
k = max(1, round(fraction * n))
rng = np.random.default_rng(seed)
idx = np.sort(rng.choice(n, size=k, replace=False))

sub = mtx[:, idx]
sub_bc = [barcodes[i] for i in idx]
print(f"downsample {matrix_dir}: {n} -> {k} cells ({fraction:.0%}, seed={seed})")

# matrix.mtx.gz — integer counts, genes x cells (readMM-compatible)
tmp = os.path.join(out_dir, "matrix.mtx")
mmwrite(tmp, sub, field="integer")
subprocess.run(["gzip", "-f", tmp], check=True)

with gzip.open(os.path.join(out_dir, "barcodes.tsv.gz"), "wt") as fh:
    fh.write("\n".join(sub_bc) + "\n")

# features unchanged (all genes kept) — copy the original gz verbatim
import shutil
shutil.copyfile(os.path.join(matrix_dir, "features.tsv.gz"),
                os.path.join(out_dir, "features.tsv.gz"))
print(f"wrote {out_dir} ({sub.shape[0]} genes x {sub.shape[1]} cells)")

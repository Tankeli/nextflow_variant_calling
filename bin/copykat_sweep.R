#!/usr/bin/env Rscript
# CopyKAT robustness sweep: one parameterised + seeded CopyKAT run.
# Generalises bin/copykat.R (which runs all defaults, no seed) so the COPYKAT_ROBUSTNESS
# subworkflow can fan out over a parameter x seed grid and we can measure where the
# aneuploid/diploid boundary lies and how stable the per-cell call is across random seeds.
#
# copykat() defines its confident-diploid baseline with a stochastic clustering step, so the
# call is NOT deterministic — set.seed() before the call is what makes a replicate reproducible
# and lets seed-to-seed flips be measured.
#
# Usage: copykat_sweep.R <matrix_dir> <comboid> <ncores> <ks_cut> <win_size> <ngene_chr> \
#                        <distance> <seed> <norm_barcodes|NONE>
# Output (cwd, named by comboid): <comboid>_copykat_prediction.txt, <comboid>_copykat_CNA_results.txt,
#         <comboid>_copykat_*  (heatmaps, clustering rds, gene-by-cell raw)

suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

args        <- commandArgs(trailingOnly = TRUE)
matrix_dir  <- args[[1]]
comboid     <- args[[2]]                 # ks<>_win<>_ng<>_<dist>_norm<0|1>_seed<> ; also sam.name
ncores      <- as.integer(args[[3]])
ks_cut      <- as.numeric(args[[4]])
win_size    <- as.integer(args[[5]])
ngene_chr   <- as.integer(args[[6]])
distance    <- args[[7]]
seed        <- as.integer(args[[8]])
norm_arg    <- if (length(args) >= 9) args[[9]] else "NONE"

# ---- load matrix (identical to bin/copykat.R) ----
mtx  <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
bc   <- readLines(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")))
feat <- read.delim(gzfile(file.path(matrix_dir, "features.tsv.gz")),
                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gene_idx <- which(feat$V3 == "Gene Expression")
m <- as.matrix(mtx[gene_idx, , drop = FALSE])
rownames(m) <- make.unique(feat$V2[gene_idx])
colnames(m) <- bc

# ---- optional known-normal baseline (norm.cell.names) ----
norm_cells <- character(0)
if (!is.null(norm_arg) && !(norm_arg %in% c("NONE", "none", "[]", ""))) {
  if (file.exists(norm_arg)) {
    norm_cells <- readLines(norm_arg)
    norm_cells <- intersect(norm_cells, colnames(m))   # only barcodes present in this matrix
  }
}

cat(sprintf(paste0("CopyKAT sweep: %s — %d genes x %d cells (ncores=%d) | ",
                   "KS.cut=%s win.size=%d ngene.chr=%d distance=%s seed=%d norm.cells=%d\n"),
            comboid, nrow(m), ncol(m), ncores,
            format(ks_cut), win_size, ngene_chr, distance, seed, length(norm_cells)))

# Reproducibility: copykat()'s baseline clustering samples cells stochastically.
set.seed(seed)

copykat(rawmat       = m,
        id.type      = "symbol",
        sam.name     = comboid,
        KS.cut       = ks_cut,
        win.size     = win_size,
        ngene.chr    = ngene_chr,
        distance     = distance,
        norm.cell.names = norm_cells,
        n.cores      = ncores)

cat("CopyKAT sweep done:", comboid, "\n")

#!/usr/bin/env Rscript
# CopyKAT aneuploid/diploid gate, run directly on a 10X filtered matrix (Gene Expression).
# Leaner than DDE_32 copyKAT_profiling.R: that read downstream Seurat objects (for UMAP/celltype
# overlays, which are downstream analysis). Here we only need the per-cell malignancy call.
#
# Usage: copykat.R <matrix_dir> <sample> <ncores>
# Output (written to cwd by copykat): <sample>_copykat_prediction.txt, <sample>_copykat_CNA_results.txt,
#         <sample>_copykat_*heatmap*, <sample>_copykat_clustering_results.rds

suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

args       <- commandArgs(trailingOnly = TRUE)
matrix_dir <- args[[1]]
sample     <- args[[2]]
ncores     <- as.integer(args[[3]])

mtx  <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
bc   <- readLines(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")))
feat <- read.delim(gzfile(file.path(matrix_dir, "features.tsv.gz")),
                   header = FALSE, sep = "\t", stringsAsFactors = FALSE)

gene_idx <- which(feat$V3 == "Gene Expression")
m <- as.matrix(mtx[gene_idx, , drop = FALSE])
rownames(m) <- make.unique(feat$V2[gene_idx])
colnames(m) <- bc

cat(sprintf("CopyKAT: %s — %d genes x %d cells (ncores=%d)\n",
            sample, nrow(m), ncol(m), ncores))

copykat(rawmat = m, id.type = "symbol", sam.name = sample, n.cores = ncores)

cat("CopyKAT done:", sample, "\n")

#!/usr/bin/env Rscript
# Numbat reproducibility SWEEP runner — one (sample x seed x threshold) combo.
# The Numbat counterpart to bin/copykat_sweep.R: re-runs run_numbat() over a seed x min_LLR grid
# REUSING the existing per-sample pileup (allele_counts.tsv.gz) + Cell Ranger matrices, so only the
# cheap clone-calling step is repeated (no re-pileup, no re-Cell Ranger). Adds set.seed() — the
# production run_numbat.R sets no seed, so this is what makes the seed axis meaningful.
#
# Usage:
#   numbat_sweep.R <label> <out_dir> <samples_csv> <allele_counts_csv> <matrix_dirs_csv> \
#                  <ncores> <max_entropy> <min_LLR> <genome> <seed>
# (identical positional contract to run_numbat.R, plus a trailing <seed>.)

suppressPackageStartupMessages({
  library(numbat)
  library(Matrix)
  library(dplyr)
})

seurat_available <- requireNamespace("Seurat", quietly = TRUE)
if (seurat_available) library(Seurat)

args        <- commandArgs(trailingOnly = TRUE)
label       <- args[[1]]
out_dir     <- args[[2]]
samples     <- strsplit(args[[3]], ",")[[1]]
allele_files<- strsplit(args[[4]], ",")[[1]]
matrix_dirs <- strsplit(args[[5]], ",")[[1]]
ncores      <- as.integer(args[[6]])
max_entropy <- as.numeric(args[[7]])
min_LLR     <- as.numeric(args[[8]])
genome      <- if (length(args) >= 9) args[[9]] else "hg38"
seed        <- if (length(args) >= 10) as.integer(args[[10]]) else 1L

stopifnot(length(samples) == length(allele_files),
          length(samples) == length(matrix_dirs))

read_10x_gene_expression <- function(matrix_dir) {
  mtx <- readMM(file.path(matrix_dir, "matrix.mtx.gz"))
  bc  <- readLines(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")))
  feat <- read.delim(gzfile(file.path(matrix_dir, "features.tsv.gz")),
                     header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gene_idx <- which(feat$V3 == "Gene Expression")
  m <- as(mtx[gene_idx, , drop = FALSE], "dgCMatrix")
  rownames(m) <- make.unique(feat$V2[gene_idx])
  colnames(m) <- bc
  m
}

cat("=== Numbat sweep:", label, "| seed", seed,
    "| min_LLR", min_LLR, "| max_entropy", max_entropy, "===\n")
cat("Start:", format(Sys.time()), "\n")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)   # the whole point of the sweep — make the RNG path deterministic per combo

# Concatenate per-sample allele counts (cell prefixed with "<sample>_").
df_allele_list <- list()
for (i in seq_along(samples)) {
  sname <- samples[i]
  d <- read.delim(allele_files[i])
  d$cell   <- paste0(sname, "_", d$cell)
  d$sample <- sname
  df_allele_list[[sname]] <- d
}
df_allele <- do.call(rbind, df_allele_list)
cat(sprintf("Joint df_allele: %d rows, %d cells\n",
            nrow(df_allele), length(unique(df_allele$cell))))

# Joint count matrix.
count_mats <- list()
for (i in seq_along(samples)) {
  sample_name <- samples[i]
  matrix_dir  <- matrix_dirs[i]
  cm <- if (seurat_available) {
    x <- Read10X(matrix_dir)
    if (is.list(x)) x[["Gene Expression"]] else x
  } else {
    read_10x_gene_expression(matrix_dir)
  }
  cm <- as(cm, "dgCMatrix")
  colnames(cm) <- paste0(sample_name, "_", colnames(cm))
  count_mats[[sample_name]] <- cm
}
common_genes <- Reduce(intersect, lapply(count_mats, rownames))
count_mat <- do.call(cbind, lapply(count_mats, function(m) m[common_genes, , drop = FALSE]))
common_cells <- intersect(colnames(count_mat), unique(df_allele$cell))
cat(sprintf("Joint matrix: %d cells x %d genes; %d with allele data\n",
            ncol(count_mat), nrow(count_mat), length(common_cells)))
if (length(common_cells) < 100) stop("Too few common cells (<100) — aborting")
count_mat <- count_mat[, common_cells]

# Numbat's tree search (run_numbat) is itself parallel + stochastic; re-seed immediately before it.
set.seed(seed)
ok <- tryCatch({
  run_numbat(
    count_mat, ref_hca, df_allele,
    genome = genome, t = 1e-5, ncores = ncores, plot = TRUE,
    out_dir = out_dir, max_iter = 2,
    max_entropy = max_entropy, min_LLR = min_LLR
  )
  TRUE
}, error = function(e) { cat("run_numbat error:", conditionMessage(e), "\n"); FALSE })

# A combo that finds no CNV is a valid SILENT result, not a failure — mark it so the array task
# exits 0 and the downstream parser can distinguish silent-but-ran from crashed.
writeLines(if (ok) "ok" else "error", file.path(out_dir, "_sweep_status.txt"))
cat("=== done:", label, "seed", seed, "min_LLR", min_LLR,
    "->", if (ok) "ok" else "error", "| End:", format(Sys.time()), "===\n")

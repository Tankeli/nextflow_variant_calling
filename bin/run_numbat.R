#!/usr/bin/env Rscript
# Joint per-patient Numbat CNV/clone calling.
# Ported from DDE_32 scripts/joint_numbat_analysis.R; parameterized for one patient
# (the patient->sample cohort now comes from the samplesheet, not a hardcoded list).
#
# Usage:
#   run_numbat.R <patient> <out_dir> <samples_csv> <allele_counts_csv> <matrix_dirs_csv> \
#                <ncores> <max_entropy> <min_LLR> <genome>
# where *_csv are comma-separated and aligned by sample order, e.g.
#   samples_csv        = Sample_2395,Sample_3001
#   allele_counts_csv  = Sample_2395_allele_counts.tsv.gz,Sample_3001_allele_counts.tsv.gz
#   matrix_dirs_csv    = Sample_2395_filtered_feature_bc_matrix,Sample_3001_filtered_feature_bc_matrix
#
# Pooling Dx+Rel into one run_numbat() call gives cross-timepoint comparable clone IDs
# and more cells for the HMM/pseudobulk filters (per-sample runs failed at default thresholds).

suppressPackageStartupMessages({
  library(numbat)
  library(Matrix)
  library(dplyr)
})

seurat_available <- requireNamespace("Seurat", quietly = TRUE)
if (seurat_available) library(Seurat)

args        <- commandArgs(trailingOnly = TRUE)
patient     <- args[[1]]
out_dir     <- args[[2]]
samples     <- strsplit(args[[3]], ",")[[1]]
allele_files<- strsplit(args[[4]], ",")[[1]]
matrix_dirs <- strsplit(args[[5]], ",")[[1]]
ncores      <- as.integer(args[[6]])
max_entropy <- as.numeric(args[[7]])
min_LLR     <- as.numeric(args[[8]])
genome      <- if (length(args) >= 9) args[[9]] else "hg38"

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

cat("=== Joint Numbat:", patient, "===\n")
cat("Start:", format(Sys.time()), "\n")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Concatenate per-sample allele counts, prefixing `cell` with "<sample>_" so cells
# are unique across samples and match the joint count matrix's column names.
cat("Loading + concatenating per-sample allele counts...\n")
df_allele_list <- list()
for (i in seq_along(samples)) {
  sname <- samples[i]
  d <- read.delim(allele_files[i])
  d$cell   <- paste0(sname, "_", d$cell)
  d$sample <- sname
  df_allele_list[[sname]] <- d
  cat(sprintf("  %s: %d rows, %d cells\n", sname, nrow(d), length(unique(d$cell))))
}
df_allele <- do.call(rbind, df_allele_list)
cat(sprintf("Joint df_allele: %d rows, %d cells\n",
            nrow(df_allele), length(unique(df_allele$cell))))

# Joint count matrix: load each sample's 10X matrix, prefix barcodes with "<sample>_".
cat("Loading + concatenating 10X matrices...\n")
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
  cat(sprintf("  %s: %d cells x %d genes\n", sample_name, ncol(cm), nrow(cm)))
}

# Common gene set across samples (numbat needs a single feature space).
common_genes <- Reduce(intersect, lapply(count_mats, rownames))
count_mat <- do.call(cbind, lapply(count_mats, function(m) m[common_genes, , drop = FALSE]))
cat(sprintf("Joint count matrix: %d cells x %d genes (common)\n",
            ncol(count_mat), nrow(count_mat)))

# Restrict to cells present in df_allele.
common_cells <- intersect(colnames(count_mat), unique(df_allele$cell))
cat(sprintf("Cells with allele data: %d\n", length(common_cells)))
if (length(common_cells) < 100) stop("Too few common cells (<100) — aborting")
count_mat <- count_mat[, common_cells]

cat(sprintf("Running joint Numbat (max_entropy=%.2f, min_LLR=%g, genome=%s)...\n",
            max_entropy, min_LLR, genome))
run_numbat(
  count_mat,
  ref_hca,
  df_allele,
  genome      = genome,
  t           = 1e-5,
  ncores      = ncores,
  plot        = TRUE,
  out_dir     = out_dir,
  max_iter    = 2,
  max_entropy = max_entropy,
  min_LLR     = min_LLR
)

cat("=== Joint Numbat complete:", patient, "===\n")
cat("End:", format(Sys.time()), "\n")

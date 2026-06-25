#!/usr/bin/env Rscript
# Stage 6 (run) — DESP cell-state demixing. The one step kept in R: it wraps the external DESP
# package (Cell Reports Methods 2024; demixes bulk into per-cell-state profiles given known cell-type
# proportions). Port of the run_desp() logic in DDE_31 6c/6d + the per-patient deep-dive.
#
# Given a batch-corrected bulk matrix (features x samples) and a cell-type proportions matrix
# (cell_type x sample), runs DESP::DESP per condition (Diagnosis / Relapse), writes per-condition
# cell-state profiles + a Relapse-Diagnosis delta matrix + run_summary, and (optionally) the same
# per patient (replicate) that has both timepoints in the overlap.
#
# Usage:
#   prot_ms_desp_run.R --bulk matrix_limma.tsv --proportions props.tsv --design design_corrected.tsv \
#       --id_col Protein.Ids --condition_col condition --replicate_col replicate \
#       [--per_patient TRUE] --outdir .

suppressPackageStartupMessages({ library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bulk"),
  make_option("--proportions"),
  make_option("--design"),
  make_option("--id_col", default = "Protein.Ids"),
  make_option("--condition_col", default = "condition"),
  make_option("--replicate_col", default = "replicate"),
  make_option("--per_patient", default = "TRUE"),
  make_option("--outdir", default = ".")
)))

if (!requireNamespace("DESP", quietly = TRUE)) {
  stop("Package DESP is required (keep-DESP-as-R-step). Install it in the proteomics R env.")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
out <- function(f) file.path(opt$outdir, f)

bulk_df <- read.delim(opt$bulk, check.names = FALSE)
prop_df <- read.delim(opt$proportions, check.names = FALSE)
design  <- read.delim(opt$design, check.names = FALSE, colClasses = "character")
# Numeric sample IDs (e.g. 2977) must stay character, else they index columns BY POSITION below.
design$sample <- as.character(design$sample)

id_col <- if (opt$id_col %in% colnames(bulk_df)) opt$id_col else colnames(bulk_df)[1]
sample_cols <- intersect(design$sample, colnames(bulk_df))
bulk <- as.matrix(bulk_df[, sample_cols, drop = FALSE])
rownames(bulk) <- make.unique(as.character(bulk_df[[id_col]]))
storage.mode(bulk) <- "numeric"

if (!"cell_type" %in% colnames(prop_df)) stop("Proportions file needs a cell_type column.")
prop_mat <- as.matrix(prop_df[, setdiff(colnames(prop_df), "cell_type"), drop = FALSE])
rownames(prop_mat) <- as.character(prop_df$cell_type)
storage.mode(prop_mat) <- "numeric"
props <- t(prop_mat)   # samples x cell_types

common <- Reduce(intersect, list(colnames(bulk), rownames(props), design$sample))
if (length(common) < 2) stop("Too few overlapping samples across bulk, proportions, design.")
bulk <- bulk[, common, drop = FALSE]
props <- props[common, , drop = FALSE]
design <- design[match(common, design$sample), , drop = FALSE]

run_desp <- function(sample_ids) {
  b <- bulk[, sample_ids, drop = FALSE]
  p <- props[sample_ids, , drop = FALSE]
  # DESP cannot take NA: row-median impute (mirrors 6d non-corrected handling).
  for (i in seq_len(nrow(b))) {
    if (anyNA(b[i, ])) {
      fill <- median(b[i, ], na.rm = TRUE)
      if (!is.finite(fill)) fill <- 0
      b[i, is.na(b[i, ])] <- fill
    }
  }
  fit <- DESP::DESP(bulk = b, proportions = p)
  if (nrow(fit) == nrow(b)) rownames(fit) <- rownames(b)
  fit
}

write_profiles <- function(mat, path) {
  df <- data.frame(feature = rownames(mat), mat, check.names = FALSE)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

cond_col <- opt$condition_col
conditions <- intersect(c("Diagnosis", "Relapse"), unique(as.character(design[[cond_col]])))
if (length(conditions) < 1) stop("No Diagnosis/Relapse samples found.")

profiles <- list()
for (cond in conditions) {
  ids <- design$sample[as.character(design[[cond_col]]) == cond]
  if (length(ids) < 1) next
  message("Running DESP for ", cond, " (", length(ids), " samples)")
  fit <- run_desp(ids)
  write_profiles(fit, out(paste0("desp_", tolower(cond), "_cell_state_profiles.tsv")))
  profiles[[cond]] <- fit
}

# Relapse - Diagnosis delta on shared features/cell types.
if (all(c("Diagnosis", "Relapse") %in% names(profiles))) {
  d <- profiles[["Diagnosis"]]; r <- profiles[["Relapse"]]
  sf <- intersect(rownames(d), rownames(r)); sc <- intersect(colnames(d), colnames(r))
  delta <- r[sf, sc, drop = FALSE] - d[sf, sc, drop = FALSE]
  write_profiles(delta, out("desp_delta_matrix.tsv"))
}

# Per-patient (replicate) DESP for patients with both timepoints in the overlap.
rep_col <- opt$replicate_col
if (toupper(opt$per_patient) %in% c("TRUE", "T", "1") && rep_col %in% colnames(design)) {
  pats <- unique(as.character(design[[rep_col]]))
  for (pid in pats) {
    pd_idx <- which(as.character(design[[rep_col]]) == pid)
    pconds <- unique(as.character(design[[cond_col]][pd_idx]))
    if (!all(c("Diagnosis", "Relapse") %in% pconds)) next
    pdir <- file.path(opt$outdir, "per_patient", pid)
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    for (cond in c("Diagnosis", "Relapse")) {
      ids <- design$sample[pd_idx][as.character(design[[cond_col]][pd_idx]) == cond]
      fit <- run_desp(ids)
      write_profiles(fit, file.path(pdir, paste0("cell_state_profiles_", tolower(cond), ".tsv")))
    }
    writeLines(c(paste0("patient=", pid),
                 paste0("conditions=", paste(pconds, collapse = ",")),
                 "status=completed"),
               file.path(pdir, "run_summary.txt"))
  }
}

writeLines(c(paste0("status=completed"),
             paste0("n_samples_overlap=", length(common)),
             paste0("samples_overlap=", paste(common, collapse = ",")),
             paste0("conditions=", paste(conditions, collapse = ",")),
             paste0("n_features=", nrow(bulk)),
             paste0("n_celltypes=", ncol(props))),
           out("desp_run_summary.txt"))
message("DESP run complete.")

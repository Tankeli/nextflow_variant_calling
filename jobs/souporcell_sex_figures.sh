#!/bin/bash
#SBATCH --job-name=dde33_soup_sexfig
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/souporcell_sex_figures_%j.log
#SBATCH --error=logs/souporcell_sex_figures_%j.err
#
# Exp 3 figures: per-sample shared-embedding UMAPs coloured by souporcell cluster | XIST | chrY
# [| timepoint], for the BMT donor/recipient samples.
set -euo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FIGDIR="${PROJ}/docs/reports/souporcell_deconvolution/figures"
PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
cd "${PROJ}"; mkdir -p "${FIGDIR}" logs
module load Miniconda3/23.5.2-0 2>/dev/null || true

"${PY}" bin/souporcell_sex_figures.py umap \
    --samplesheet assets/test/souporcell_mix_bmt.csv \
    --results results_soupmix_bmt/callers/souporcell \
    --k 2 --outdir "${FIGDIR}"

echo "[$(date)] BMT sex UMAPs written to ${FIGDIR}"

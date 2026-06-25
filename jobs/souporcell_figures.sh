#!/bin/bash
#SBATCH --job-name=dde33_soup_fig
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/souporcell_figures_%j.log
#SBATCH --error=logs/souporcell_figures_%j.err
#
# Figures for the souporcell deconvolution report: per-mix shared-embedding UMAPs (true origin vs
# souporcell cluster per K) + ARI/confusion statistics. Reads the published clusters + Cell Ranger
# matrices; writes PNGs into the Obsidian report's figures/ folder.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FIGDIR="${PROJ}/docs/reports/souporcell_deconvolution/figures"
PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
cd "${PROJ}"
mkdir -p "${FIGDIR}" logs

module load Miniconda3/23.5.2-0 2>/dev/null || true

"${PY}" bin/souporcell_mix_figures.py \
    --samplesheet assets/test/souporcell_mix_controls.csv \
    --results results_soupmix_controls/callers/souporcell \
    --label controls --outdir "${FIGDIR}"

"${PY}" bin/souporcell_mix_figures.py \
    --samplesheet assets/test/souporcell_mix_patients.csv \
    --results results_soupmix_patients/callers/souporcell \
    --label patients --outdir "${FIGDIR}"

echo "[$(date)] figures written to ${FIGDIR}"
ls -la "${FIGDIR}"

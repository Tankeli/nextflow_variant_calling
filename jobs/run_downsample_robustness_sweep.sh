#!/bin/bash
#SBATCH --job-name=dde33_ds_cksweep
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-180%24
#SBATCH --output=logs/ds_cksweep_%A_%a.log
#SBATCH --error=logs/ds_cksweep_%A_%a.err
#
# CopyKAT cell-number confounder test: is the per-sample ARI difference driven by N?
# For 3 high-cell healthy controls (PBM_2, PBM_1, PBMMC_3) downsampled to 75/50/25% of cells
# (fixed subset per fraction), re-run the SAME 20-combo grid (KS.cut × seed) so ARI at each N is
# directly comparable to the full-N (100%) sweep already in results_controls/copykat_robustness/.
# Standalone (no Nextflow); short 4h walltime (downsampled = faster) to clear maintenance YOR796.
set -euo pipefail
PROJECT=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJECT}"
export PATH="${PROJECT}/bin:${PATH}"
module load Miniconda3/23.5.2-0
source activate snv

MANIFEST=jobs/patients_downsample_sweep_manifest.tsv
read -r SAMPLE FRAC COMBO KS WIN NG DIST SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "no manifest line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

BASE="${PROJECT}/results_controls/copykat_robustness_downsample/${SAMPLE}/frac${FRAC}"
MTX="${BASE}/matrix"
OUT="${BASE}/sweep/${COMBO}"
mkdir -p "${OUT}"; cd "${OUT}"

echo "[$(date)] ${SAMPLE} frac${FRAC} ${COMBO} (ncores=${SLURM_CPUS_PER_TASK})"
copykat_sweep.R "${MTX}" "${COMBO}" "${SLURM_CPUS_PER_TASK}" "${KS}" "${WIN}" "${NG}" "${DIST}" "${SEED}" NONE
echo "[$(date)] done ${SAMPLE} frac${FRAC} ${COMBO}"

#!/bin/bash
#SBATCH --job-name=dde33_pat_cksweep
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=6
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=1-80%24
#SBATCH --output=logs/pat_cksweep_%A_%a.log
#SBATCH --error=logs/pat_cksweep_%A_%a.err
#
# Standalone CopyKAT robustness sweep for the 2 DDE_32 prototype patients (4 samples, Dx+Rel).
# One array task per (sample × combo) from jobs/patients_sweep_manifest.tsv. Runs copykat_sweep.R
# DIRECTLY on the published CellRanger matrices — does NOT use Nextflow, so it neither locks the
# active work_patients run nor overwrites the production CopyKAT calls in results_patients/copykat/.
# Short walltime (6h) so tasks backfill before maintenance YOR796 (2026-06-11 09:00).
set -euo pipefail
PROJECT=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJECT}"
export PATH="${PROJECT}/bin:${PATH}"
module load Miniconda3/23.5.2-0
source activate snv   # copykat R package (same env as the COPYKAT process on the viking profile)

MANIFEST=jobs/patients_sweep_manifest.tsv
read -r SAMPLE COMBO KS WIN NG DIST SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
[ -n "${SAMPLE:-}" ] || { echo "no manifest line ${SLURM_ARRAY_TASK_ID}"; exit 1; }

MTX="${PROJECT}/results_patients/cellranger/${SAMPLE}/outs/filtered_feature_bc_matrix"
OUT="${PROJECT}/results_patients/copykat_robustness/${SAMPLE}/sweep/${COMBO}"
mkdir -p "${OUT}"; cd "${OUT}"

echo "[$(date)] ${SAMPLE} ${COMBO} (ncores=${SLURM_CPUS_PER_TASK})"
copykat_sweep.R "${MTX}" "${COMBO}" "${SLURM_CPUS_PER_TASK}" "${KS}" "${WIN}" "${NG}" "${DIST}" "${SEED}" NONE
echo "[$(date)] done ${SAMPLE} ${COMBO}"

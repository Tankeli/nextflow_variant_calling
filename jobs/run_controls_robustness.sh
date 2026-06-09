#!/bin/bash
#SBATCH --job-name=dde33_ck_sweep
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=40:00:00
#SBATCH --output=logs/copykat_sweep_orchestrator_%j.log
#SBATCH --error=logs/copykat_sweep_orchestrator_%j.err
#
# CopyKAT robustness SWEEP orchestrator (the Nextflow half of the hybrid track) for the controls.
# Reuses the cached CellRanger from -work-dir work (-resume); the other callers are disabled so the
# only new tasks are COPYKAT_SWEEP (parameter x seed grid). Downstream analysis is the separate
# jobs/run_copykat_robustness.sh.
# Usage: sbatch jobs/run_controls_robustness.sh [samplesheet.csv]
set -euo pipefail
cd /mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

SHEET="${1:-assets/controls_samplesheet_all9.csv}"
echo "CopyKAT sweep orchestrator start $(date) — samplesheet: ${SHEET}"

nextflow run . \
    -profile viking \
    -params-file params-controls.yaml \
    --input "${SHEET}" \
    --run_copykat_robustness \
    --run_numbat false --run_souporcell false --run_clonetracer false \
    --run_qc false --run_reference_mapping false --run_integration false \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir work \
    -resume

echo "CopyKAT sweep orchestrator finished $(date)"

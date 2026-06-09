#!/bin/bash
#SBATCH --job-name=dde33_controls
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/controls_orchestrator_%j.log
#SBATCH --error=logs/controls_orchestrator_%j.err
#
# Nextflow orchestrator for the Caron PBMMC control run (-profile viking).
# Lightweight head job; the slurm executor dispatches each pipeline process as its own job.
# Usage: sbatch jobs/run_controls.sh [samplesheet.csv]
set -euo pipefail

cd /mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

SHEET="${1:-assets/controls_samplesheet.csv}"
echo "Orchestrator start $(date) — samplesheet: ${SHEET}"

# Drop nf-co2footprint for this run: numbat.sif lacks `ps`, and the co2footprint plugin forces
# Nextflow's ps-based metric collection. nf-schema (validation) + nf-prov (BCO) only.
nextflow run . \
    -profile viking \
    -params-file params-controls.yaml \
    --input "${SHEET}" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir work \
    -resume

echo "Orchestrator finished $(date)"

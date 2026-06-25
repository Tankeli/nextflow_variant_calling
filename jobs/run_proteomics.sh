#!/bin/bash
#SBATCH --job-name=dde33_proteomics
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/proteomics_orchestrator_%j.log
#SBATCH --error=logs/proteomics_orchestrator_%j.err
#
# Nextflow orchestrator for the F157 bulk-proteomics branch (-entry PROTEOMICS, -profile viking).
# Lightweight head job; the slurm executor dispatches each PROT_MS_* process as its own job
# (Python stages -> conda env aml_proteomics, DESP -> conda env desp_r).
# Usage: sbatch jobs/run_proteomics.sh
set -euo pipefail

cd /mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

echo "Orchestrator start $(date) — F157 proteomics (-entry PROTEOMICS)"

# nf-schema (validation skipped for entry workflows anyway) + nf-prov (BCO); drop nf-co2footprint
# to avoid its ps-based metric collection, matching the other viking runs.
nextflow run . \
    -entry PROTEOMICS \
    -profile viking \
    -params-file params-proteomics.yaml \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir work_proteomics \
    --outdir results_proteomics \
    -resume

echo "Orchestrator finished $(date)"

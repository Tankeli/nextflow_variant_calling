#!/bin/bash
#SBATCH --job-name=dde33_ct
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/controls_ct_orchestrator_%j.log
#SBATCH --error=logs/controls_ct_orchestrator_%j.err
#
# Live CloneTracer test on the Caron PBMMC controls (-profile viking), resuming the cached
# CellRanger / Numbat / souporcell / reference-mapping work. Only the new mtDNA pileup +
# CloneTracer build + model (+ figure) run fresh.
#
# Model runs on GPU: conf/viking.config routes CLONETRACER to the `clonetracer_gpu` conda env
# (torch 1.13.1+cu117 / pyro 1.8.4) on the a40 `gpu` partition when --clonetracer_gpu is set.
#
# Usage: sbatch jobs/run_controls_clonetracer.sh [samplesheet.csv]
set -euo pipefail

cd /mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

SHEET="${1:-assets/controls_samplesheet_batch1.csv}"
GTF=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
echo "Orchestrator start $(date) — samplesheet: ${SHEET}"

nextflow run . \
    -profile viking \
    -params-file params-controls.yaml \
    --input "${SHEET}" \
    --run_clonetracer \
    --clonetracer_gpu \
    --clonetracer_gtf "${GTF}" \
    --clonetracer_mtdna_chrom chrM \
    -c conf/maint_cap.config \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir work \
    -resume

echo "Orchestrator finished $(date)"

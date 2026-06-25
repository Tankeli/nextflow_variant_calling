#!/bin/bash
#SBATCH --job-name=dde33_pat_ct
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/patients_ct_orchestrator_%j.log
#SBATCH --error=logs/patients_ct_orchestrator_%j.err
#
# Live CloneTracer run on the prototype AML cohort (Patient_1 + Patient_2) through DDE_33.
# Resumes the cached CellRanger / Numbat / souporcell / reference-mapping work; only the new
# mtDNA pileup + CloneTracer build + model + figures run fresh.
#
# Model runs on CPU (the `nodes` partition): conf/viking.config routes CLONETRACER to the pinned
# `clonetracer_gpu` conda env (pyro 1.8.4 / torch 1.13.1), which runs correctly+quickly on CPU, with
# OMP/MKL threads capped to the allocation (avoids the OpenMP deadlock). Mutations are hard-capped
# (clonetracer_max_total_muts=6) so the tree search cannot explode. No --clonetracer_gpu => no GPU
# queue. Isolated launchDir + work_patients so it can run alongside other orchestrators.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
# MUST reuse the original patient launchDir: -resume reads the session cache from <launchDir>/.nextflow,
# so launching from a fresh dir re-runs Cell Ranger + everything. run_patients.sh used .launch_patients.
LAUNCH="${PROJ}/.launch_patients"
mkdir -p "${LAUNCH}" "${PROJ}/logs"

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

GTF=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz

cd "${LAUNCH}"
echo "Patients CloneTracer orchestrator start $(date) — launchDir ${LAUNCH}"

# nf-co2footprint dropped (numbat.sif lacks `ps`); nf-schema + nf-prov only.
nextflow run "${PROJ}" \
    -profile viking \
    -params-file "${PROJ}/params-patients.yaml" \
    --run_clonetracer \
    --clonetracer_gtf "${GTF}" \
    --clonetracer_mtdna_chrom chrM \
    -c "${PROJ}/conf/maint_cap.config" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir "${PROJ}/work_patients" \
    -resume

echo "Patients CloneTracer orchestrator finished $(date)"

#!/bin/bash
#SBATCH --job-name=dde33_downstream
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/downstream_patients_orchestrator_%j.log
#SBATCH --error=logs/downstream_patients_orchestrator_%j.err
#
# Nextflow orchestrator for the DOWNSTREAM-only RNA best-practices stack (ported from DDE_27),
# run OFF the already-published Cell Ranger matrices (no Cell Ranger re-run, no variant callers).
# Stages: RNA core + cohort integration + advanced (pseudotime/DE/composition) + protein/ADT.
# RNA processes run from DDE_27's built .sif (conf/viking.config). Isolated launch dir + work dir
# so it does NOT share work/.nextflow with the patients/controls caller runs.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
LAUNCH="${PROJ}/.launch_downstream"   # isolated launchDir (own .nextflow + nextflow.log)
mkdir -p "${LAUNCH}" "${PROJ}/logs"

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

cd "${LAUNCH}"
echo "Downstream orchestrator start $(date) — launchDir ${LAUNCH}"

# nf-schema only (paramsSummaryLog); no prov/co2 for this analysis-only run.
nextflow run "${PROJ}" \
    -entry DOWNSTREAM \
    -profile viking \
    -params-file "${PROJ}/params-downstream-patients.yaml" \
    -plugins nf-schema@2.7.2 \
    -work-dir "${PROJ}/work_downstream" \
    -resume

echo "Downstream orchestrator finished $(date)"

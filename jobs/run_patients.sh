#!/bin/bash
#SBATCH --job-name=dde33_patients
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/patients_orchestrator_%j.log
#SBATCH --error=logs/patients_orchestrator_%j.err
#
# Nextflow orchestrator for the DDE_32 prototype AML cohort (Patient_1 + Patient_2) through DDE_33.
# Runs the FULL pipeline incl. Phase-2 integration + headline Sankeys (run_integration=true).
#
# Isolated from the controls run: launches from its own dir with -work-dir work_patients so it does
# NOT share work/ or .nextflow with jobs/run_controls.sh (both can run concurrently).
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
LAUNCH="${PROJ}/.launch_patients"   # isolated launchDir (own .nextflow + nextflow.log)
mkdir -p "${LAUNCH}" "${PROJ}/logs"

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

cd "${LAUNCH}"
echo "Patients orchestrator start $(date) — launchDir ${LAUNCH}"

# Drop nf-co2footprint (numbat.sif lacks `ps`); nf-schema + nf-prov only.
nextflow run "${PROJ}" \
    -profile viking \
    -params-file "${PROJ}/params-patients.yaml" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir "${PROJ}/work_patients" \
    -resume

echo "Patients orchestrator finished $(date)"

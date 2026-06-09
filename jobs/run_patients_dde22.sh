#!/bin/bash
#SBATCH --job-name=dde33_patients_dde22
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=logs/patients_dde22_orchestrator_%j.log
#SBATCH --error=logs/patients_dde22_orchestrator_%j.err
#
# Nextflow orchestrator for the extended Vivobank AML cohort (DDE_22 samples) through DDE_33.
# Runs the FULL pipeline incl. Phase-2 integration + headline Sankeys (run_integration=true).
#
# Isolated from the other runs: launches from its own dir with -work-dir work_patients_dde22 so it
# does NOT share work/ or .nextflow with jobs/run_patients.sh or jobs/run_controls.sh.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
LAUNCH="${PROJ}/.launch_patients_dde22"   # isolated launchDir (own .nextflow + nextflow.log)
mkdir -p "${LAUNCH}" "${PROJ}/logs"

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home

cd "${LAUNCH}"
echo "Patients(DDE_22) orchestrator start $(date) — launchDir ${LAUNCH}"

# Drop nf-co2footprint (numbat.sif lacks `ps`); nf-schema + nf-prov only.
nextflow run "${PROJ}" \
    -profile viking \
    -params-file "${PROJ}/params-patients-dde22.yaml" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -work-dir "${PROJ}/work_patients_dde22" \
    -resume

echo "Patients(DDE_22) orchestrator finished $(date)"

#!/bin/bash
#SBATCH --job-name=dde33_soup_mix
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/souporcell_mixing_orchestrator_%j.log
#SBATCH --error=logs/souporcell_mixing_orchestrator_%j.err
#
# Souporcell deconvolution-validation experiments (Exp 1 + Exp 2).
#   Exp 1: artificial mixes of healthy controls (within- and cross-study) -> assets/test/souporcell_mix_controls.csv
#   Exp 2: artificial mixes of patient diagnosis samples (single-origin)   -> assets/test/souporcell_mix_patients.csv
# Runs the standalone -entry SOUPORCELL_ONLY off ALREADY-PUBLISHED Cell Ranger outs (no Cell Ranger
# re-run, no other callers). Each mix is a fake `patient` grouping several samples; the <sample>__
# barcode prefix preserves ground truth so bin/souporcell_mix_eval.py can score the result.
#
# This is the lightweight orchestrator; souporcell itself is dispatched per-mix as its own SLURM job
# by the viking profile. Keep this off the login node (it is an sbatch job).
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FASTA=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home
# Compute nodes have no internet: point at the pre-pulled apptainer image cache so no runtime pull
# is attempted (samtools image pre-pulled on the login node into this dir).
export NXF_SINGULARITY_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
export NXF_APPTAINER_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
mkdir -p "${PROJ}/logs"

run_mix () {
    local tag=$1 sheet=$2
    local launch="${PROJ}/.launch_soupmix_${tag}"
    mkdir -p "${launch}"
    cd "${launch}"
    echo "[$(date)] souporcell mixing '${tag}' start — launchDir ${launch}"
    nextflow run "${PROJ}" \
        -entry SOUPORCELL_ONLY \
        -profile viking \
        --souporcell_samplesheet "${PROJ}/${sheet}" \
        --souporcell_fasta "${FASTA}" \
        --souporcell_k '2,3' \
        --outdir "${PROJ}/results_soupmix_${tag}" \
        -work-dir "${PROJ}/work_soupmix_${tag}" \
        -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
        -resume
    echo "[$(date)] souporcell mixing '${tag}' finished"
}

run_mix controls assets/test/souporcell_mix_controls.csv
run_mix patients assets/test/souporcell_mix_patients.csv

# Score both experiments (env with pandas/numpy).
EVAL_PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
for tag in controls patients; do
    "${EVAL_PY}" "${PROJ}/bin/souporcell_mix_eval.py" \
        --results "${PROJ}/results_soupmix_${tag}/callers/souporcell" \
        --outdir  "${PROJ}/results_soupmix_${tag}/eval"
done

echo "[$(date)] all souporcell mixing experiments + eval done"

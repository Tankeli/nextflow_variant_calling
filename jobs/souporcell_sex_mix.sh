#!/bin/bash
#SBATCH --job-name=dde33_soup_sexmix
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/souporcell_sex_mix_orchestrator_%j.log
#SBATCH --error=logs/souporcell_sex_mix_orchestrator_%j.err
#
# Positive control for the sex-validation method: artificially mix a MALE and a FEMALE diagnosis
# sample, deconvolute, and check souporcell (genotype) ≈ sex-expression call ≈ true barcode origin.
#   MIX_sex_M8178_F2977  Sample_8178 (AML107, M) + Sample_2977 (AML152, F)
#   MIX_sex_M2395_F2958  Sample_2395 (AML057, M) + Sample_2958 (AML155, F)
set -euo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FASTA=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa
SHEET=assets/test/souporcell_mix_sex.csv
FIGDIR="${PROJ}/docs/reports/souporcell_deconvolution/figures"
PY=/users/hbp534/.conda/envs/aml_scrna/bin/python

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home
export NXF_SINGULARITY_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
export NXF_APPTAINER_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
mkdir -p "${PROJ}/logs" "${FIGDIR}"

LAUNCH="${PROJ}/.launch_soupmix_sex"; mkdir -p "${LAUNCH}"; cd "${LAUNCH}"
echo "[$(date)] sex-mix souporcell start"
nextflow run "${PROJ}" \
    -entry SOUPORCELL_ONLY -profile viking \
    --souporcell_samplesheet "${PROJ}/${SHEET}" \
    --souporcell_fasta "${FASTA}" \
    --souporcell_k '2,3' \
    --outdir "${PROJ}/results_soupmix_sex" \
    -work-dir "${PROJ}/work_soupmix_sex" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 -resume
echo "[$(date)] sex-mix souporcell finished"

module load Miniconda3/23.5.2-0 2>/dev/null || true
# barcode ground-truth ARI
"${PY}" "${PROJ}/bin/souporcell_mix_eval.py" \
    --results "${PROJ}/results_soupmix_sex/callers/souporcell" \
    --outdir  "${PROJ}/results_soupmix_sex/eval"
# 3-way triangulation: souporcell vs sex-expression vs true origin (+ figures)
"${PY}" "${PROJ}/bin/souporcell_sex_deconv.py" \
    --samplesheet "${PROJ}/${SHEET}" \
    --results "${PROJ}/results_soupmix_sex/callers/souporcell" \
    --sexmap Sample_8178:M,Sample_2977:F,Sample_2395:M,Sample_2958:F \
    --outdir "${FIGDIR}" --csvout "${PROJ}/results_soupmix_sex/sex_deconv" --k 2

echo "[$(date)] sex-mix done"

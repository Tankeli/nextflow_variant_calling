#!/bin/bash
#SBATCH --job-name=dde33_soup_bmt
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=logs/souporcell_bmt_orchestrator_%j.log
#SBATCH --error=logs/souporcell_bmt_orchestrator_%j.err
#
# Souporcell deconvolution-validation Exp 3: real biological donor/recipient mixtures in
# BM-transplant samples (no artificial mixing). -entry SOUPORCELL_ONLY off published Cell Ranger outs.
#   AML066_rel_solo  Sample_1386            female recipient (46,XX), relapse alone  -> expect donor+recipient split
#   AML066_dxrel     Sample_3652 + 1386     Dx+Rel combined: recipient spans both, donor relapse-only
#   AML163_dx_solo   Sample_1187            female, Allograft, Dx-only (pre-transplant) -> single-origin negative ctrl
#   AML079_rel_solo  Sample_1255            male recipient, relapse (sex-discordance comparison)
#   AML107_rel_solo  Sample_1894            male recipient, relapse (sex-discordance comparison)
# After clustering, validate the split by sex-chromosome expression (XIST vs chrY genes): a male+female
# cluster pair within one sample corroborates the genotype split and infers donor sex (unrecorded).
# Lightweight orchestrator; souporcell dispatches per-mix as its own SLURM job via the viking profile.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FASTA=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa
SHEET=assets/test/souporcell_mix_bmt.csv

module load Java/17.0.6
export PATH="/users/hbp534/.local/bin:$PATH"
export NXF_HOME=/mnt/scratch/users/hbp534/nextflow_home
# Compute nodes have no internet: use the pre-pulled apptainer image cache (no runtime pull).
export NXF_SINGULARITY_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
export NXF_APPTAINER_CACHEDIR=/mnt/scratch/users/hbp534/apptainer/nf-apptainer
mkdir -p "${PROJ}/logs"

LAUNCH="${PROJ}/.launch_soupmix_bmt"
mkdir -p "${LAUNCH}"
cd "${LAUNCH}"
echo "[$(date)] souporcell BMT (Exp 3) start — launchDir ${LAUNCH}"
nextflow run "${PROJ}" \
    -entry SOUPORCELL_ONLY \
    -profile viking \
    --souporcell_samplesheet "${PROJ}/${SHEET}" \
    --souporcell_fasta "${FASTA}" \
    --souporcell_k '2,3' \
    --outdir "${PROJ}/results_soupmix_bmt" \
    -work-dir "${PROJ}/work_soupmix_bmt" \
    -plugins nf-schema@2.7.2,nf-prov@1.7.0 \
    -resume
echo "[$(date)] souporcell BMT clustering finished"

EVAL_PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
# generic ARI/doublet scoring (low ARI expected for true single-individual solo, higher where two genotypes)
"${EVAL_PY}" "${PROJ}/bin/souporcell_mix_eval.py" \
    --results "${PROJ}/results_soupmix_bmt/callers/souporcell" \
    --outdir  "${PROJ}/results_soupmix_bmt/eval"
# sex-chromosome corroboration of the donor/recipient split
"${EVAL_PY}" "${PROJ}/bin/souporcell_sex_validate.py" \
    --samplesheet "${PROJ}/${SHEET}" \
    --results "${PROJ}/results_soupmix_bmt/callers/souporcell" \
    --outdir  "${PROJ}/results_soupmix_bmt/sex_eval"

echo "[$(date)] souporcell BMT Exp 3 + validation done"

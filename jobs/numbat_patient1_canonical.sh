#!/bin/bash
#SBATCH --job-name=dde33_p1_numbat
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --output=logs/numbat_p1_canonical_%j.log
#SBATCH --error=logs/numbat_p1_canonical_%j.err
#
# One-off: produce Patient_1's canonical NUMBAT_RUN checkpoint directly (the main-pipeline -resume
# cache-missed CellRanger after the interrupted 06-09 run, so the orchestrator re-ran CellRanger
# needlessly). Invokes the SAME pipeline script (bin/run_numbat.R) the NUMBAT_RUN module calls, at the
# new operating point (min_LLR=5, max_entropy=0.8) with the bumped resources (400 G / 8 cores) the
# 15k-cell joint needs — reusing the existing published pileup + Cell Ranger matrices (no CellRanger).
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"; mkdir -p logs
NUMBAT_SIF=/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/singularity/numbat.sif
module load Apptainer/latest

PILE="$PROJ/results_patients/numbat_joint/Patient_1/Patient_1_pileup"
CR="$PROJ/results_patients/cellranger"
OUT="$PROJ/results_patients/numbat_joint/Patient_1/numbat_out"
mkdir -p "$OUT"

apptainer exec -B /mnt/scratch/users/hbp534 "$NUMBAT_SIF" \
    Rscript "$PROJ/bin/run_numbat.R" \
        Patient_1 \
        "$OUT" \
        Sample_2395,Sample_3001 \
        "$PILE/Sample_2395_allele_counts.tsv.gz,$PILE/Sample_3001_allele_counts.tsv.gz" \
        "$CR/Sample_2395/outs/filtered_feature_bc_matrix,$CR/Sample_3001/outs/filtered_feature_bc_matrix" \
        8 0.8 5 hg38

echo "[$(date)] Patient_1 numbat done -> $(ls "$OUT"/segs_consensus_*.tsv 2>/dev/null | tail -1)"

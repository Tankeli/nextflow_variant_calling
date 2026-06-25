#!/bin/bash
#SBATCH --job-name=pat_mtdna
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/pat_mtdna_%A_%a.log
#SBATCH --error=logs/pat_mtdna_%A_%a.err
# Standalone mtDNA pileup (cellsnp-lite via numbat.sif) on the published patient BAMs — mirrors the
# MTDNA_PILEUP module. Feeds the CloneTracer mtDNA axis. Bypasses the pipeline (full -resume would
# re-run cellranger/numbat/souporcell from scratch due to a task-hash change).
set -euo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
SIF=/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/singularity/numbat.sif
cd "$PROJ"
module load Apptainer/latest

SAMPLES=(Sample_2395 Sample_3001 Sample_2977 Sample_0109)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
BAM="$PROJ/results_patients/cellranger/$S/outs/possorted_genome_bam.bam"
MTX="$PROJ/results_patients/cellranger/$S/outs/filtered_feature_bc_matrix"
OUT="$PROJ/results_patients/clonetracer/mtdna/${S}_mtdna"
echo "[$S] start $(date +%T) bam=$BAM"

WORK=$(mktemp -d)
zcat "$MTX/barcodes.tsv.gz" > "$WORK/barcodes.txt"
apptainer exec -B "$PROJ" -B "$WORK" "$SIF" cellsnp-lite \
    -s "$BAM" -b "$WORK/barcodes.txt" -O "$OUT" -p "${SLURM_CPUS_PER_TASK:-4}" \
    --chrom chrM --minMAF 0.01 --minCOUNT 20 --cellTAG CB --UMItag UB --genotype
echo "[$S] done $(date +%T); outputs:"; ls -la "$OUT" 2>/dev/null | head

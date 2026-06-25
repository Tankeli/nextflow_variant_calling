#!/bin/bash
#SBATCH --job-name=dde33_nb_p1
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=8
#SBATCH --mem=400G
#SBATCH --time=12:00:00
#SBATCH --array=46-54%5
#SBATCH --output=logs/numbat_sweep_%A_%a.log
#SBATCH --error=logs/numbat_sweep_%A_%a.err
#
# Numbat reproducibility SWEEP — standalone SLURM array (one task per sample x seed x min_LLR combo).
# The Numbat counterpart to the CopyKAT robustness sweep, but deliberately NOT via Nextflow: it
# consumes the already-published per-sample pileups (allele_counts.tsv.gz) + Cell Ranger matrices
# directly, so it re-runs ONLY run_numbat() and cannot re-draw / overwrite the production calls the
# way the CopyKAT NF sweep did (analysis 02 caveat). Invokes numbat.sif via apptainer directly,
# which also sidesteps the numbat.sif procps issue that only bites Nextflow's metric launcher.
#
# Prereq: python3 bin/numbat_sweep_manifest.py  (writes assets/numbat_sweep_manifest.tsv)
# Usage:  sbatch jobs/numbat_sweep.sh   (edit --array=1-N%C to match the manifest row count)
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"
mkdir -p logs
MANIFEST="$PROJ/assets/numbat_sweep_manifest.tsv"
NUMBAT_SIF=/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/singularity/numbat.sif

module load Apptainer/latest

# Pull this task's row (skip header; column 1 == idx).
read -r idx label out_dir samples allele matrices ncores max_entropy min_LLR genome seed < <(
    awk -F'\t' -v id="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==id {print; exit}' "$MANIFEST")

if [ -z "${label:-}" ]; then
    echo "No manifest row for task $SLURM_ARRAY_TASK_ID" >&2; exit 1
fi

echo "[$(date)] task $SLURM_ARRAY_TASK_ID: $label seed=$seed min_LLR=$min_LLR ent=$max_entropy"
echo "  out_dir: $out_dir"
mkdir -p "$out_dir"

apptainer exec -B /mnt/scratch/users/hbp534 "$NUMBAT_SIF" \
    Rscript "$PROJ/bin/numbat_sweep.R" \
        "$label" "$out_dir" "$samples" "$allele" "$matrices" \
        "8" "$max_entropy" "$min_LLR" "$genome" "$seed"

echo "[$(date)] task $SLURM_ARRAY_TASK_ID done -> $(cat "$out_dir/_sweep_status.txt" 2>/dev/null || echo no-status)"

#!/bin/bash
#SBATCH --job-name=dde33_ds_analysis
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ds_analysis_%j.log
#SBATCH --error=logs/ds_analysis_%j.err
#
# Downsample confounder analysis: per-(sample,fraction) stability, then the ARI-vs-N summary.
# Run after the downsample sweep array completes. Reuses the existing full-N (100%) stability from
# results_controls/copykat_robustness/_analysis as the anchor point.
set -euo pipefail
PROJECT=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJECT}"
export PATH="${PROJECT}/bin:${PATH}"
module load Miniconda3/23.5.2-0
source activate aml_scrna

DS=results_controls/copykat_robustness_downsample
OUT=${DS}/_analysis
mkdir -p "${OUT}"

# 1. per-(sample,fraction) stability over the downsampled sweep
for s in PBM_2 PBM_1 PBMMC_3; do
  for f in 0.75 0.50 0.25; do
    SW=${DS}/${s}/frac${f}/sweep
    if [ -d "${SW}" ]; then
      ( cd "${OUT}" && copykat_stability.py "${s}_frac${f}" "${PROJECT}/${SW}" )
    else
      echo "[skip] no sweep dir ${SW}"
    fi
  done
done

# 2. ARI-vs-N summary (+ figure), reusing the full-N points
copykat_downsample_summary.py \
    "${PROJECT}/${OUT}" \
    "${PROJECT}/results_controls/copykat_robustness/_analysis" \
    "${PROJECT}/${OUT}/downsample" \
    PBM_2:7302 PBM_1:6248 PBMMC_3:5548

echo "Done — outputs in ${OUT}"

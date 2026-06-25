#!/bin/bash
#SBATCH --job-name=dde33_soup_prop_orch
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/souporcell_prop_orch_%j.log
#SBATCH --error=logs/souporcell_prop_orch_%j.err
#
# Proportion-titration experiment: does an unbalanced mix reduce souporcell accuracy?
# Reuses the cached merged BAMs from the 50-50 mixes (no re-prep), and for each pair re-runs
# souporcell at minority fractions 1/5/10/25/50 % with TOTAL cells held constant, so any accuracy
# drop is purely the proportion effect. Submits a job array (pair x ratio) then an eval job.
#   ctrl_PBMMC23   minority PBMMC_2  / majority PBMMC_3   (healthy controls)
#   pat_2395_2977  minority 2395     / majority 2977      (patient diagnosis)
set -euo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
OUTBASE="${PROJ}/results_soupmix_prop"
PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
RATIOS="1 5 10 25 50"
mkdir -p "${OUTBASE}" "${PROJ}/logs"

# pick the valid cached merged BAM for a mix: newest task dir with .exitcode==0 holding the BAM+bai+barcodes
find_cached () {  # $1 workdir  $2 mixname
    local wd="$1" mix="$2" best="" bestt=0
    while IFS= read -r bam; do
        local d; d=$(dirname "${bam}")
        [ -f "${d}/.exitcode" ] && [ "$(cat "${d}/.exitcode")" = "0" ] || continue
        [ -f "${bam}.bai" ] && [ -f "${d}/${mix}.barcodes.tsv" ] || continue
        local t; t=$(stat -c %Y "${bam}")
        if [ "${t}" -gt "${bestt}" ]; then bestt="${t}"; best="${bam}"; fi
    done < <(find "${wd}" -name "${mix}.merged.sorted.bam" 2>/dev/null)
    echo "${best}"
}

stage_pair () {  # $1 pair  $2 workdir  $3 mixname  $4 minority  $5 majority  $6 total
    local pair="$1" wd="$2" mix="$3" minority="$4" majority="$5" total="$6"
    local bam; bam=$(find_cached "${wd}" "${mix}")
    [ -n "${bam}" ] || { echo "ERROR: no valid cached merged BAM for ${mix} under ${wd}" >&2; exit 1; }
    local d; d=$(dirname "${bam}")
    mkdir -p "${OUTBASE}/${pair}"
    ln -sf "${bam}"        "${OUTBASE}/${pair}/merged.sorted.bam"
    ln -sf "${bam}.bai"    "${OUTBASE}/${pair}/merged.sorted.bam.bai"
    ln -sf "${d}/${mix}.barcodes.tsv" "${OUTBASE}/${pair}/full.barcodes.tsv"
    echo "  ${pair}: using ${bam}"
    for pct in ${RATIOS}; do
        echo "${pair},${pct},${OUTBASE}/${pair}/merged.sorted.bam,${OUTBASE}/${pair}/full.barcodes.tsv,${minority},${majority},${total}" >> "${OUTBASE}/manifest.csv"
    done
}

echo "pair,pct,merged,bclist,minority,majority,total" > "${OUTBASE}/manifest.csv"
# minority = the smaller sample so all ratios up to 50% are reachable at constant total per pair
stage_pair ctrl_PBMMC23  "${PROJ}/work_soupmix_controls" MIX_within_PBMMC23 PBMMC_2     PBMMC_3     4000
stage_pair pat_2395_2977 "${PROJ}/work_soupmix_patients" MIX_dx_2395_2977   Sample_2977 Sample_2395 2500

N=$(( $(wc -l < "${OUTBASE}/manifest.csv") - 1 ))
echo "manifest has ${N} runs"
ARRAY_ID=$(sbatch --parsable --array=0-$((N-1)) "${PROJ}/jobs/souporcell_proportions_run.sh")
echo "submitted run array ${ARRAY_ID}"
sbatch --dependency=afterok:${ARRAY_ID} --account=biol-stem-2022 --partition=nodes \
    --job-name=dde33_soup_prop_eval --cpus-per-task=2 --mem=8G --time=00:30:00 \
    --output=logs/souporcell_prop_eval_%j.log --error=logs/souporcell_prop_eval_%j.err \
    --wrap "${PY} ${PROJ}/bin/souporcell_proportions_eval.py \
        --manifest ${OUTBASE}/manifest.csv --results ${OUTBASE} \
        --outdir ${PROJ}/docs/reports/souporcell_deconvolution/figures \
        --csvout ${OUTBASE}/eval"
echo "submitted eval (afterok:${ARRAY_ID})"

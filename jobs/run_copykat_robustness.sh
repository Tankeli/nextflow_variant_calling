#!/bin/bash
#SBATCH --job-name=dde33_ck_robust
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/copykat_robustness_%j.log
#SBATCH --error=logs/copykat_robustness_%j.err
#
# Downstream CopyKAT robustness analysis (the standalone half of the hybrid track).
# Runs AFTER the COPYKAT_ROBUSTNESS sweep has published its combos. Operates on already-published
# checkpoints only, so it is cheap to re-run / iterate:
#   - stability + boundary  <- the sweep combos      (results_*/robustness/<s>/sweep)
#   - drivers / crossref / celltype  <- the production CopyKAT call + reference mapping + atlas
#
# Usage: sbatch jobs/run_copykat_robustness.sh [RESULTS_DIR] [ATLAS_H5AD] [SAMPLE ...]
#   RESULTS_DIR  default results_controls
#   ATLAS_H5AD   default = atlas from params-controls.yaml
#   SAMPLE ...   default = every sample dir under <RESULTS_DIR>/callers/copykat
set -euo pipefail
PROJECT=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJECT}"
# bin/ is only on PATH inside Nextflow tasks; add it here so the copykat_*.py scripts are callable.
export PATH="${PROJECT}/bin:${PATH}"

module load Miniconda3/23.5.2-0
source activate aml_scrna

RESULTS="${1:-results_controls}"
ATLAS="${2:-/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/references/paediatric_bm_reference/bone_marrow_atlas.h5ad}"
shift || true; shift || true
# Atlas-derived gene sets are computed once and cached (CKROB_* tune the crossref). Wilcoxon is the
# default marker method; export CKROB_MARKER_METHOD=t-test for a much faster (coarser) run.
MARKER_METHOD="${CKROB_MARKER_METHOD:-wilcoxon}"

if [ "$#" -gt 0 ]; then SAMPLES=("$@"); else
    mapfile -t SAMPLES < <(ls -1 "${RESULTS}/callers/copykat" 2>/dev/null | grep -v '^figures$' || true)
fi
[ "${#SAMPLES[@]}" -gt 0 ] || { echo "No samples found under ${RESULTS}/callers/copykat"; exit 1; }

OUT="${RESULTS}/robustness/_analysis"
mkdir -p "${OUT}"; cd "${OUT}"
RES="../../.."   # back to project root from results_*/robustness/_analysis

echo "CopyKAT robustness analysis — results=${RESULTS} atlas=${ATLAS} samples=${SAMPLES[*]}"

for S in "${SAMPLES[@]}"; do
    echo "=== ${S} ==="
    CK="${RES}/${RESULTS}/callers/copykat/${S}"
    SWEEP="${RES}/${RESULTS}/robustness/${S}/sweep"
    MAPPED="${RES}/${RESULTS}/annotation/reference_mapping/${S}/${S}_mapped.h5ad"
    CNA="${CK}/${S}_copykat_CNA_results.txt"
    GBC="${CK}/${S}_copykat_CNA_raw_results_gene_by_cell.txt"
    PRED="${CK}/${S}_copykat_prediction.txt"
    [ -f "${MAPPED}" ] || MAPPED=NONE
    [ -f "${GBC}" ]    || GBC=NONE

    # 1) seed/param stability + classification boundary (over the sweep)
    if [ -d "${SWEEP}" ]; then
        copykat_stability.py "${S}" "${SWEEP}" "${MAPPED}"
    else
        echo "  [skip] no sweep dir ${SWEEP} — run the COPYKAT_ROBUSTNESS sweep first"
    fi
    CONS="${S}_copykat_stability.csv"; [ -f "${CONS}" ] || CONS="${PRED}"

    # 2) driver genes/regions for the production call
    if [ -f "${CNA}" ]; then
        copykat_drivers.py "${S}" "${CNA}" "${GBC}" "${CONS}"
        # 3) cross-reference drivers vs anchor genes + cell-type signatures
        [ -f "${S}_copykat_drivers.csv" ] && copykat_crossref.py "${S}" "${S}_copykat_drivers.csv" "${ATLAS}" 200 30 200 50 auto "${MARKER_METHOD}"
        # 4) per-cell-type variance matrix + aneuploid fraction
        [ "${MAPPED}" != "NONE" ] && copykat_celltype_matrix.py "${S}" "${MAPPED}" "${CNA}" "${CONS}"
    else
        echo "  [skip] no CNA results ${CNA}"
    fi
done

echo "Done — outputs in ${RESULTS}/robustness/_analysis"

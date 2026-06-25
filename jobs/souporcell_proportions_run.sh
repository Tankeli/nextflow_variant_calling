#!/bin/bash
#SBATCH --job-name=dde33_soup_prop
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/souporcell_prop_%A_%a.log
#SBATCH --error=logs/souporcell_prop_%A_%a.err
#
# One souporcell run per (pair, minority%) of the proportion titration. Reads its row from the
# manifest by SLURM_ARRAY_TASK_ID, subsamples the combined barcode list to the target ratio (total
# cells held constant), and runs souporcell -k2 on the cached merged BAM.
set -euo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
FASTA=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa
SIF=/mnt/scratch/users/hbp534/DDE_24_CITE_seq_reference_map/singularity/souporcell_release.sif
MANIFEST="${PROJ}/results_soupmix_prop/manifest.csv"
PY=/users/hbp534/.conda/envs/aml_scrna/bin/python
module load Apptainer/latest

# manifest line 1 is the header; array task 0 -> data row 2
row=$(( SLURM_ARRAY_TASK_ID + 2 ))
IFS=, read -r pair pct merged bclist minority majority total < <(sed -n "${row}p" "${MANIFEST}")
echo "[$(date)] task ${SLURM_ARRAY_TASK_ID}: pair=${pair} pct=${pct} minority=${minority} total=${total}"

OUT="${PROJ}/results_soupmix_prop/${pair}/r${pct}"
mkdir -p "${OUT}"
BC="${OUT}/barcodes.txt"
"${PY}" "${PROJ}/bin/souporcell_subsample_barcodes.py" \
    --bclist "${bclist}" --minority-sample "${minority}" \
    --total "${total}" --minority-frac "$(${PY} -c "print(${pct}/100)")" --seed 0 --out "${BC}"

[ -f "${FASTA}.fai" ] || apptainer exec --bind /mnt/scratch "${SIF}" samtools faidx "${FASTA}"
apptainer exec --bind /mnt/scratch "${SIF}" souporcell_pipeline.py \
    -i "${merged}" -b "${BC}" -f "${FASTA}" -t "${SLURM_CPUS_PER_TASK}" \
    -o "${OUT}" -k 2 --no_umi true --skip_remap True --ignore True

echo "[$(date)] done ${pair} r${pct}"

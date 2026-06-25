#!/bin/bash
#SBATCH --job-name=ct_refig
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/ct_refig_%j.log
#SBATCH --error=logs/ct_refig_%j.err
set -uo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"; module load Miniconda3/23.5.2-0; source activate aml_scrna
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
for P in Patient_1 Patient_2; do
  python -u bin/clonetracer_figures.py results_patients/clonetracer/$P/${P}_out.pickle $P results_patients/clonetracer/$P/figures
done
echo "DONE"

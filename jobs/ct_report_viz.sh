#!/bin/bash
#SBATCH --job-name=ct_rviz
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=logs/ct_rviz_%j.log
#SBATCH --error=logs/ct_rviz_%j.err
set -uo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"; module load Miniconda3/23.5.2-0; source activate aml_scrna
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
python -u bin/clonetracer_report_viz.py results_patients/clonetracer/report_viz \
  Patient_2:results_patients/clonetracer/Patient_2/Patient_2_clone_assignments.csv:Sample_2977:Sample_0109 \
  Patient_1:results_patients/clonetracer/Patient_1/Patient_1_clone_assignments.csv:Sample_2395:Sample_3001
echo "DONE"; ls -la results_patients/clonetracer/report_viz/

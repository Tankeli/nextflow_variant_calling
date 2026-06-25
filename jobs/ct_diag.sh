#!/bin/bash
#SBATCH --job-name=ct_diag
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ct_diag_%j.log
#SBATCH --error=logs/ct_diag_%j.err
set -uo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"
module load Miniconda3/23.5.2-0
source activate clonetracer_gpu
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "host=$(hostname) start=$(date)"
python -u bin/ct_diag.py validate_ct/HD_BM_3_sub8.json 60 validate_ct
echo "exit=$? end=$(date)"

#!/bin/bash
#SBATCH --job-name=ct_real
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ct_real_%j.log
#SBATCH --error=logs/ct_real_%j.err
set -uo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "$PROJ"
module load Miniconda3/23.5.2-0
source activate clonetracer_gpu
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export PYTHONPATH="$PROJ/bin:${PYTHONPATH:-}"

echo "===== sub3 (3 muts), DEFAULT num_iter (300, init=200) $(date +%T) ====="
SECONDS=0
python -u bin/run_clonetracer_noredir.py -i validate_ct/HD_BM_3_sub3.json -n HD_BM_3_sub3 -o validate_ct
echo "sub3 exit=$? wall=${SECONDS}s $(date +%T)"
python -u bin/clonetracer_assignments.py validate_ct/HD_BM_3_sub3_out.pickle validate_ct/HD_BM_3_sub3_clone_assignments.csv && echo "assignments OK"
echo "=== outputs ==="; ls -la validate_ct/HD_BM_3_sub3*
echo "=== assignments head ==="; head validate_ct/HD_BM_3_sub3_clone_assignments.csv
echo "DONE $(date +%T)"

#!/bin/bash
#SBATCH --job-name=ct_sub8
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/ct_sub8_%j.log
#SBATCH --error=logs/ct_sub8_%j.err
#
# Run the full CloneTracer model on an 8-mutation subset of the real HD_BM_3 JSON (1645 cells,
# 4 SNV + 4 mtDNA, no CNVs -> 300 SVI iters). Earlier 100-mut/8-mut runs hung at ~2s CPU: that was
# an OpenMP oversubscription deadlock (OMP_NUM_THREADS unset -> torch spawns ~96 threads onto 4
# allocated cores). Capping OMP/MKL threads to the allocation fixes it (the stepwise diagnostic with
# OMP_NUM_THREADS=4 ran clean). CPU partition + fast pinned env (pyro 1.8.4 / torch 1.13.1), no -g.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJ}"
mkdir -p logs validate_ct

module load Miniconda3/23.5.2-0
source activate clonetracer_gpu
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
export PYTHONPATH="${PROJ}/bin:${PYTHONPATH:-}"
# Match thread pools to the SLURM allocation to avoid the OpenMP oversubscription deadlock.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

echo "Start $(date)  host=$(hostname)"
python -c "import torch; print('torch', torch.__version__, 'cuda?', torch.cuda.is_available())"

JSON=validate_ct/HD_BM_3_sub8.json
SECONDS=0
python -u bin/run_clonetracer.py -i "${JSON}" -n HD_BM_3_sub8 -o validate_ct
echo "model wall: ${SECONDS}s"
bin/clonetracer_assignments.py validate_ct/HD_BM_3_sub8_out.pickle validate_ct/HD_BM_3_sub8_clone_assignments.csv

echo "=== model stdout (CloneTracer writes its own log) ==="
tail -30 logs/HD_BM_3_sub8/HD_BM_3_sub8_stdout.txt 2>/dev/null
echo "=== outputs ==="
ls -la validate_ct/
echo "=== assignments head ==="
head validate_ct/HD_BM_3_sub8_clone_assignments.csv
echo "Finished $(date)"

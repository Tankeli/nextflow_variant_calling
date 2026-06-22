#!/bin/bash
#SBATCH --job-name=ct_sub8
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --output=logs/ct_sub8_%j.log
#SBATCH --error=logs/ct_sub8_%j.err
#
# Diagnosis test: run the CloneTracer model on an 8-mutation subset of the real HD_BM_3 JSON
# (1645 cells, 4 SNV + 4 mtDNA, no CNVs -> 300 SVI iters). The full 100-mut JSON timed out at 2h
# because infer_hierarchy does a combinatorial per-mutation tree search; this confirms a small
# mutation set runs quickly. CPU partition (all a40 GPUs were fully allocated, ~16h queue) using
# the fast pinned env (pyro 1.8.4 / torch 1.13.1) without -g. Unbuffered so progress is visible.
set -euo pipefail

PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
cd "${PROJ}"
mkdir -p logs validate_ct

module load Miniconda3/23.5.2-0
source activate clonetracer_gpu
export PYTHONNOUSERSITE=1
export PYTHONUNBUFFERED=1
export PYTHONPATH="${PROJ}/bin:${PYTHONPATH:-}"

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

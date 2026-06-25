#!/bin/bash
#SBATCH --job-name=pat_ct_run
#SBATCH --account=biol-stem-2022
#SBATCH --partition=nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/pat_ct_run_%j.log
#SBATCH --error=logs/pat_ct_run_%j.err
# Standalone CloneTracer build + model + figures for both patients, off the existing results_patients
# outputs (full-pipeline -resume is unusable: task hashes changed since the cached run). Mirrors
# CLONETRACER_WF + PLOT_CLONETRACER. Cap 4 total muts (Patient_2 has CNVs -> 500 iters; the heuristic
# tree search explodes, so keep small). Patient_1 = mtDNA-only (no numbat_out / souporcell).
set -uo pipefail
PROJ=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling
GTF=/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz
CAP=4
cd "$PROJ"
mtx(){ echo "$PROJ/results_patients/cellranger/$1/outs/filtered_feature_bc_matrix"; }
mt(){  echo "$PROJ/results_patients/clonetracer/mtdna/$1_mtdna"; }
module load Miniconda3/23.5.2-0

run_patient(){
  local PAT="$1" DX="$2" REL="$3" EXTRA="$4"
  local OUT="$PROJ/results_patients/clonetracer/$PAT"; mkdir -p "$OUT"
  echo "============ $PAT ($DX,$REL) $(date +%T) ============"

  # ---- 1. build JSON (aml_scrna) ----
  source activate aml_scrna
  export PYTHONNOUSERSITE=1
  python -u "$PROJ/bin/clonetracer_build_json.py" \
    --patient "$PAT" --samples "$DX,$REL" --timepoints Dx,Rel \
    --matrices "$(mtx $DX),$(mtx $REL)" \
    --mtdna-dirs "$(mt $DX),$(mt $REL)" \
    --gtf "$GTF" $EXTRA \
    --max-snvs "$CAP" --mtdna-max-sites "$CAP" --max-total-muts "$CAP" \
    --output "$OUT/$PAT.json" || { echo "$PAT: build FAILED"; return 1; }
  conda deactivate

  # ---- 2. model (clonetracer_gpu env, CPU; OMP capped to avoid deadlock) ----
  source activate clonetracer_gpu
  export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1
  export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
  export PYTHONPATH="$PROJ/bin:${PYTHONPATH:-}"
  python -u "$PROJ/bin/run_clonetracer_noredir.py" -i "$OUT/$PAT.json" -n "$PAT" -o "$OUT" \
    || { echo "$PAT: model FAILED"; conda deactivate; return 1; }
  python -u "$PROJ/bin/clonetracer_assignments.py" "$OUT/${PAT}_out.pickle" "$OUT/${PAT}_clone_assignments.csv"
  conda deactivate

  # ---- 3. figures (aml_scrna): trees/ELBO/heatmap + UMAP overlay ----
  source activate aml_scrna
  export PYTHONNOUSERSITE=1
  mkdir -p "$OUT/figures"
  python -u "$PROJ/bin/clonetracer_figures.py" "$OUT/${PAT}_out.pickle" "$PAT" "$OUT/figures"
  ( cd "$OUT/figures" && python -u "$PROJ/bin/plot_clonetracer_umap.py" "$PAT" \
      "$OUT/${PAT}_clone_assignments.csv" "$DX,$REL" Dx,Rel \
      "$PROJ/results_patients/reference_mapping/$DX/${DX}_mapped.h5ad" \
      "$PROJ/results_patients/reference_mapping/$REL/${REL}_mapped.h5ad" )
  conda deactivate
  echo "$PAT DONE $(date +%T); figures:"; ls "$OUT/figures" 2>/dev/null
}

# Patient_2: Numbat CNV + souporcell SNV + mtDNA (richest)
run_patient Patient_2 Sample_2977 Sample_0109 \
  "--numbat-dir $PROJ/results_patients/numbat_joint/Patient_2/numbat_out --souporcell-dir $PROJ/results_patients/souporcell/Patient_2/k3"

# Patient_1: mtDNA-only (no numbat_out / souporcell available)
run_patient Patient_1 Sample_2395 Sample_3001 ""

echo "ALL DONE $(date +%T)"

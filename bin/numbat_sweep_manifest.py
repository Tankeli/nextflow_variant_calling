#!/usr/bin/env python3
"""
Build the (sample x seed x min_LLR) manifest the Numbat sweep SLURM array consumes.

Screen-then-replicate design: the 5 samples that called clones in the production run get the full
seed x min_LLR grid (where the reproducibility variance + compute cost actually live); the 6 silent
controls get a seeds-only negative-control check at the production min_LLR (raising min_LLR keeps a
no-CNV sample trivially silent, so the ladder there is wasted). See lab-book session for the DoE
rationale (factor roles + power-for-instability) behind not running a full factorial.

One TSV row per combo with the exact run_numbat inputs (reusing existing published pileups +
Cell Ranger matrices). Columns:
  idx label out_dir samples_csv allele_csv matrix_dirs_csv ncores max_entropy min_LLR genome seed
"""
import os, sys, itertools

PROJ = "/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling"
SEEDS = [1, 2, 3]
MIN_LLR = [3, 5, 10]      # full ladder, ACTIVE samples only
SILENT_MIN_LLR = [3]      # negative-control: production floor only
MAX_ENTROPY = [0.8]
GENOME = "hg38"
NCORES = 8

# sample -> (results_root, [member_sample, ...]); single-member = control "joint" of one sample.
CONTROLS = ["HD_BM_1", "HD_BM_2", "HD_BM_3", "HD_BM_4",
            "PBM_1", "PBM_2", "PBMMC_1", "PBMMC_2", "PBMMC_3"]
SAMPLES = {s: ("results_controls", [s]) for s in CONTROLS}
SAMPLES["Patient_1"] = ("results_patients", ["Sample_2395", "Sample_3001"])
SAMPLES["Patient_2"] = ("results_patients", ["Sample_2977", "Sample_0109"])

# Samples that called clones at production thresholds (analysis 03) -> full grid. Everything else is
# silent -> seeds-only. Patient_1 never finished NUMBAT_RUN but its pileup exists and it is a tumour,
# so it is treated as ACTIVE (the sweep also yields its first clone calls).
ACTIVE = {"PBM_2", "PBMMC_2", "PBMMC_3", "Patient_1", "Patient_2"}


def allele_path(root, label, member):
    return f"{PROJ}/{root}/numbat_joint/{label}/{label}_pileup/{member}_allele_counts.tsv.gz"


def matrix_path(root, member):
    return f"{PROJ}/{root}/cellranger/{member}/outs/filtered_feature_bc_matrix"


def main():
    out = os.path.join(PROJ, "assets", "numbat_sweep_manifest.tsv")
    rows, missing = [], []
    for label, (root, members) in SAMPLES.items():
        alleles = [allele_path(root, label, m) for m in members]
        matrices = [matrix_path(root, m) for m in members]
        for f in alleles + matrices:
            if not os.path.exists(f):
                missing.append(f)
        llr_grid = MIN_LLR if label in ACTIVE else SILENT_MIN_LLR
        for seed, llr, ent in itertools.product(SEEDS, llr_grid, MAX_ENTROPY):
            combo = f"seed{seed}_llr{llr}_ent{ent}"
            out_dir = f"{PROJ}/{root}/numbat_robustness/{label}/{combo}"
            rows.append([label, out_dir, ",".join(members), ",".join(alleles),
                         ",".join(matrices), NCORES, ent, llr, GENOME, seed])

    if missing:
        print("MISSING inputs (fix before submitting):", file=sys.stderr)
        for m in sorted(set(missing)):
            print("  ", m, file=sys.stderr)
        sys.exit(1)

    header = ["idx", "label", "out_dir", "samples_csv", "allele_csv", "matrix_dirs_csv",
              "ncores", "max_entropy", "min_LLR", "genome", "seed"]
    with open(out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for i, r in enumerate(rows, 1):
            fh.write("\t".join([str(i)] + [str(x) for x in r]) + "\n")
    n_active = len(ACTIVE)
    n_silent = len(SAMPLES) - n_active
    print(f"Wrote {out}: {len(rows)} combos "
          f"({n_active} active x {len(SEEDS)*len(MIN_LLR)} + "
          f"{n_silent} silent x {len(SEEDS)*len(SILENT_MIN_LLR)})")


if __name__ == "__main__":
    main()

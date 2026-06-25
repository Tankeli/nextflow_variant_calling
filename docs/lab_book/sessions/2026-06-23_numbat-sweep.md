---
date: 2026-06-23
project: DDE_33
type: session
status: done
tags: [numbat, sweep, reproducibility, seed, min_llr, doe, slurm, oom, apptainer, robustness]
related: ["[[04_numbat_reproducibility]]", "[[03_numbat_specificity]]", "[[02_copykat_robustness]]"]
---

# 2026-06-23 — Numbat reproducibility sweep (seed × min_LLR)

> **Objective:** the "repeat the CopyKAT robustness on Numbat" the 2026-06-12 meeting flagged — sweep
> seed × min_LLR over controls + tumour to test whether Numbat clones survive re-seeding, off the
> existing pileups (no re-pileup, no Nextflow). Resolves the two open items in [[03_numbat_specificity]].

## Work done

### 1. Design — DoE instead of full factorial

Started at a full factorial (3 seeds × 3 min_LLR × 11 samples = 99). Reduced it by **factor role**:
seed = replication/nuisance, min_LLR = treatment curve, sample = block → don't full-cross. **Power for
detecting instability** with R seeds at flip-rate f: `R ≥ ln(β)/ln(1−f)` → R=3 gives ~80 % power at
f≥0.5 (gross), ~8 needed for f≈0.2. Combined with the fact that **6/11 samples are silent** (analysis
03), adopted a **screen-then-replicate** allocation: full 3×3 on the 5 clone-callers, seeds-only at
min_LLR=3 on the 6 silent controls → **63 combos** (−36 % vs 99) with no loss where the variance lives.

### 2. Standalone sweep harness (not Nextflow — by design)

- `bin/numbat_sweep.R` — `run_numbat()` + `set.seed()` (production sets none) + a `_sweep_status.txt`
  marker (silent-but-ran vs crashed); reuses published `allele_counts.tsv.gz` + Cell Ranger matrices.
- `bin/numbat_sweep_manifest.py` — emits `assets/numbat_sweep_manifest.tsv` (63 rows; ACTIVE vs silent).
- `jobs/numbat_sweep.sh` — SLURM array, `apptainer exec numbat.sif Rscript …` (direct apptainer
  sidesteps the numbat.sif procps issue that only bites Nextflow's metric launcher).
- Smoke-tested **one** fast element first (`35140582_1`, HD_BM_1, COMPLETED 3:47, silent → status ok)
  before releasing the array — a plumbing bug would otherwise have wasted ~63 jobs.

### 3. Running it — two resource bugs, both fixed

- **Array `35140733` (1-63%20, 24G):** 15 silent controls COMPLETED; PBMMC_2 finished (flagged for
  exceeding the *soft* 24G request); **39 combos OUT_OF_MEMORY**. Per-task MaxRSS showed the
  tree-building joints peak at **100–130 G** at ncores=8 (numbat forks workers → cores multiply RAM).
  The `nodes` partition has **503 G/node**, so 24 G was needlessly tight.
- **Hi-mem rerun `35176483` (200G, ncores=4):** PBM_1/PBM_2/PBMMC_3/Patient_2 recovered (peaked ~94–96 G
  at ncores=4); **Patient_1 (15,098 cells) TIMEOUT** at the 6 h wall — too slow at 4 cores.
- **Patient_1:** discovered the reproducibility parser already read its outputs — `run_numbat` had
  written clone_post + segs (57/57/57 segs across seeds) **before** the timeout kill; only the trailing
  `_sweep_status.txt` was missing. So a queued `35215594` (8c/400G/12h) Patient_1 rerun was **redundant
  → cancelled**; backfilled the status markers. **63/63 combos valid.**

### 4. Analysis — `bin/numbat_reproducibility.py`

Mean pairwise ARI of per-cell `clone_opt` across seeds + clone/aneuploid-fraction range + min_LLR
response. Outputs `results_controls/numbat_robustness/_analysis/` (+ figures committed to
`docs/lab_book/assets/04_numbat_reproducibility/`).

## Results / outcome

- **Specificity reproducible:** 6/6 silent controls silent in 3/3 seeds.
- **Seed-instability is confined to the worst over-call:** PBMMC_2 (median LLR 8) ARI **0.47**, clones
  flip 4↔5; **every other call ARI 1.00** — including the other healthy FPs and both tumours.
- **min_LLR response:** PBM_2 noise gone at ≥5; PBMMC_2 needs 10 to clear (but 10 also kills the real
  Patient_2); PBMMC_3 + Patient_1 persist at all thresholds.
- **PBMMC_3 = genuine CNV** (persists at llr10, ARI 1.0, LLR 463) — analysis-03 open question resolved.
- **Patient_1 now has clone calls** (4 clones, ARI 1.0) — tumour comparator n=1 → n=2.
- **Operating point:** `min_LLR=5` + a **seed-ARI gate** (run ≥3 seeds, flag ARI < ~0.9). Two orthogonal
  filters (LLR confidence + seed reproducibility) separate signal from noise where neither alone can.

Distilled into [[04_numbat_reproducibility]].

## Decisions & rationale

- **Standalone array, not a Nextflow sweep subworkflow.** Reuses published pileups directly → re-runs
  only `run_numbat()`, cannot re-draw production calls (the CopyKAT NF-sweep failure mode), and full
  control of the grid. Cost: not wired into `main.nf` (acceptable for an analysis track).
- **Screen-then-replicate over full factorial.** DoE + power + the 6 silent samples → 63 not 99 combos,
  no loss on the samples that carry the variance/cost.
- **ncores as the memory lever.** numbat forks → fewer cores = lower peak RAM; dropped 8→4 for the
  hi-mem rerun. Patient_1 alone needed cores back (8) for speed within walltime.

## Issues / blockers

- Soft `--mem` on this cluster: tasks exceed the request and are flagged OOM but may still finish if the
  node has RAM — non-deterministic. Mitigation: request true peak (≥128 G for tree-builders, ≥200 G for
  big joints) + modest concurrency, don't rely on the cap.
- Patient_1 outputs came from a timed-out (not clean) run — verified complete, but re-run clean if it
  becomes load-bearing.

## Next steps

- [ ] Adopt `min_LLR=5` + seed-ARI gate in the Numbat output contract (`run_numbat.R` / `numbat_run.nf`);
      surface mean clone ARI as a per-patient QC number.
- [ ] Inspect PBMMC_3 segments (`segs_consensus`, `bulk_clones`) — characterise the genuine CNV.
- [ ] (Optional) escalate to ~8 seeds only on borderline-ARI samples (sequential design).
- [ ] Push `results_controls/numbat_robustness/` to Longship.

## Environment

- git: DDE_33 `main`@2521d03 + uncommitted (`bin/numbat_sweep*.{R,py}`, `bin/numbat_reproducibility.py`,
  `jobs/numbat_sweep*.sh`, lab book).
- numbat.sif `/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/singularity/numbat.sif` via
  Apptainer/latest; analysis in conda `aml_scrna`. SLURM `nodes` partition (503 G/node), biol-stem-2022.
- Jobs: `35140582` (smoke), `35140733` (sweep), `35176483` (hi-mem rerun), `35215594` (cancelled).

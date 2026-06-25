---
title: Numbat reproducibility — seed × min_LLR sweep on controls + tumour
project: DDE_33
type: analysis
status: active
updated: 2026-06-23
tags: [numbat, cnv, reproducibility, seed-stability, min_llr, sweep, robustness, healthy-baseline, ari, operating-point]
related: ["[[03_numbat_specificity]]", "[[02_copykat_robustness]]", "[[2026-06-23_numbat-sweep]]", "[[2026-06-09_patient-cohort-runs]]"]
---

# 04 — Numbat reproducibility (seed × min_LLR sweep)

## Question / goal

The reproducibility counterpart to the specificity snapshot ([[03_numbat_specificity]]) and the direct
Numbat analogue of the CopyKAT robustness sweep ([[02_copykat_robustness]]) — the "repeat the above
on Numbat" the 2026-06-12 meeting anticipated. Analysis 03 showed Numbat is *specific* (silent on 6/9
healthy) and that its false-positives grade by CNV LLR, but off **single runs**: it could not say
whether a clone is **reproducible**. Here: does the clone call survive **re-seeding** (Numbat's tree
search is stochastic and the production `run_numbat.R` sets no seed), and how does it respond to the
**min_LLR** threshold? Resolves analysis 03's two open items — the PBMMC_3 "real or artefact?" flag
and the missing seed/threshold sweep.

## Data & provenance

Standalone seed × min_LLR sweep (`bin/numbat_sweep.R`, `set.seed()` + grid), **reusing the existing
published pileups + Cell Ranger matrices** so only `run_numbat()` re-runs — it cannot re-draw the
production calls the way the CopyKAT NF sweep did (analysis 02 caveat). **Screen-then-replicate**
design (DoE rationale + power analysis in [[2026-06-23_numbat-sweep]]): the 5 samples that called
clones in production get the full grid; the 6 silent controls get a seeds-only negative control.

| Class | Samples | grid | combos |
|---|---|---|---|
| active (clone-callers) | PBM_2, PBMMC_2, PBMMC_3, Patient_1, Patient_2 | seed{1,2,3} × min_LLR{3,5,10} | 5×9 = 45 |
| silent controls | HD_BM_1-4, PBM_1, PBMMC_1 | seed{1,2,3} × min_LLR{3} | 6×3 = 18 |

**63/63 combos completed** (SLURM array `35140733` + hi-mem rerun `35176483`; numbat.sif via apptainer).
Reproducibility parsed by `bin/numbat_reproducibility.py` (mean pairwise Adjusted Rand Index of
per-cell `clone_opt` across seeds; clone counts + aneuploid fraction range across seeds).

> **New positive comparator:** Patient_1 (Sample_2395 + Sample_3001, 15,098 cells) — its `NUMBAT_RUN`
> never finished in the main pipeline ([[2026-06-09_patient-cohort-runs]]); the sweep produced its
> **first clone calls**, raising the tumour class from n=1 (the analysis-03 limitation) to n=2.

## Results

`results_controls/numbat_robustness/_analysis/numbat_reproducibility_{seed,llr}.csv`.

### Seed reproducibility at production min_LLR=3 (the headline)

![Numbat seed reproducibility](../assets/04_numbat_reproducibility/numbat_reproducibility_seed_ari.png)

| Sample | group | seeds called | n_clones (min–max) | aneuploid frac (range) | **seed ARI** |
|---|---|---|---|---|---|
| HD_BM_1-4, PBM_1, PBMMC_1 | healthy | 0/3 | 0 | 0 | reproducibly **silent** |
| **PBMMC_2** | healthy | 3/3 | **4–5** | 0.096–0.181 | **0.47** |
| PBM_2 | healthy | 3/3 | 2–2 | 0.195 | 1.00 |
| PBMMC_3 | healthy | 3/3 | 3–3 | 0.189 | 1.00 |
| Patient_1 | tumour | 3/3 | 4–4 | 0.371 | 1.00 |
| Patient_2 | tumour | 3/3 | 2–2 | 0.374 | 1.00 |

Two clean findings:
1. **Specificity is reproducible.** All 6 truly-silent controls are silent in **3/3 seeds** — Numbat's
   correct negatives don't flicker into existence on a re-seed (contrast CopyKAT, where 18–98 % of
   cells flipped class across the param×seed sweep, [[02_copykat_robustness]]).
2. **Only the worst over-call is seed-unstable.** **PBMMC_2** — the sample that emitted 4–8 "clones" at
   a median CNV LLR of just 8 (analysis 03) — is the **sole non-reproducible** sample: ARI **0.47**,
   with the clone count itself flipping **4↔5** between seeds. Its low-LLR segments re-partition
   differently every seed: this *is* noise, now proven by re-seeding. Every other call — including the
   other two healthy false-positives (PBM_2, PBMMC_3) and both tumours — is **ARI 1.00** (bit-identical
   partitions across seeds).

### min_LLR threshold response

![min_LLR response](../assets/04_numbat_reproducibility/numbat_reproducibility_llr_response.png)

Mean n_clones over seeds, by min_LLR (3 → 5 → 10):

| Sample | llr3 | llr5 | llr10 | reading |
|---|---|---|---|---|
| PBM_2 (healthy, LLR 4.6) | 2 | **0** | 0 | noise — gone the moment min_LLR ≥ 5 |
| PBMMC_2 (healthy, unstable) | 4.7 | 3.3 | **0** | survives to llr5 (still seed-unstable), only llr10 clears it |
| PBMMC_3 (healthy, LLR 463) | 3 | 2 | **2** | **persists at every threshold** |
| Patient_1 (tumour) | 4 | 4 | 3 | robust |
| Patient_2 (tumour, LLR 558) | 2 | 2 | **0** | real, but **lost at llr10** |

The threshold is **not** a clean universal cut: `min_LLR=5` removes PBM_2's noise while keeping
PBMMC_3 and both tumours, but `min_LLR=10` **over-prunes** — it finally kills PBMMC_2's unstable call
yet also wipes out a *genuine* tumour (Patient_2 → 0 clones). So threshold alone cannot separate
PBMMC_2 (noise) from Patient_2 (real): they cross over between llr5 and llr10.

### PBMMC_3 is a genuine CNV (analysis-03 open question resolved)

PBMMC_3 calls clones at **every** min_LLR including 10, is **perfectly seed-stable (ARI 1.00)** at all,
and carries a median segment **LLR 463** (analysis 03). All three signatures are tumour-like, none are
the low-LLR/seed-unstable signature of a Numbat artefact → PBMMC_3 carries a **real copy-number event**
in a labelled-healthy Caron cord sample (constitutional/mosaic CNV or a sample-quality issue), not a
false positive. It should be excluded from the "healthy over-call" set and looked at on its own.

## Interpretation

Numbat's clone call is **far more reproducible than CopyKAT's** and the instability that exists is
**confined to, and predicted by, the low-LLR calls**: the one seed-unstable sample (PBMMC_2) is the one
flagged weakest by LLR in analysis 03. **Two orthogonal filters fully separate signal from noise:**
- **median segment LLR** (confidence, analysis 03) — removes PBM_2 (LLR 4.6, gone at min_LLR≥5);
- **seed ARI** (reproducibility, this analysis) — removes PBMMC_2 (ARI 0.47), which *no tolerable
  threshold* can drop without also losing a real tumour.

**Operational recommendation:** keep the caller at `max_entropy=0.8` for sensitivity but set the
**operating point at `min_LLR=5`** (not 3: drops PBM_2 noise; not 10: loses Patient_2), and add a
**seed-stability gate** — run ≥3 seeds and report mean clone ARI, treating ARI < ~0.9 as "low-confidence
/ review" regardless of LLR. PBMMC_3-type calls (high LLR, ARI 1.0, persist at high threshold) are
genuine; PBMMC_2-type (modest LLR, ARI < 0.5) are noise. This replaces the single-run, single-threshold
clone membership in the current output contract with a confidence- and reproducibility-aware call.

## Limitations / caveats

- **3 seeds** detects gross instability (per the power calc, R=3 → ~80 % power at flip-rate f≥0.5;
  [[2026-06-23_numbat-sweep]]); subtle instability (f≈0.2) would need ~8 seeds. PBMMC_2's ARI 0.47 is
  gross and safely detected; samples at ARI 1.0 across 3 seeds *could* still have rare-cell instability
  a larger R would surface. Escalating seeds only on borderline samples is the cheap next step.
- **Patient_1 outputs came from a timed-out task** (`run_numbat` wrote clone_post + segs before the 6 h
  kill; `_sweep_status` backfilled). Files verified complete (15,098 cells; 57/57/57 segs across seeds,
  ARI 1.0) but the run was not re-done clean — re-confirm if Patient_1 becomes load-bearing.
- **Active set fixed from the production single run.** A sample that was silent in production but would
  call on some seed is not in the active grid; the 6 silent controls were only swept at min_LLR=3.
- **min_LLR is the only threshold swept** (max_entropy fixed at 0.8); the entropy axis is untested.
- Aneuploid fraction uses `p_cnv>0.5` vs the most-populous (reference) clone, as in analysis 03.

## Links

- Produced by: [[2026-06-23_numbat-sweep]] (worklog: design, power calc, OOM/timeout debugging, job IDs)
- Resolves open items in: [[03_numbat_specificity]] (PBMMC_3 real-CNV flag; the "no sweep yet" caveat)
- Contrasts: [[02_copykat_robustness]] (CopyKAT: 18–98 % cells flip, no usable axis; Numbat: instability
  confined to low-LLR, removable by LLR + ARI)
- Feeds: [[2026-06-09_patient-cohort-runs]] (min_LLR=5 + seed-ARI gate on patient clones; Patient_1 now
  has clone calls)
- Code/figures: `bin/numbat_sweep.R`, `bin/numbat_sweep_manifest.py`, `bin/numbat_reproducibility.py`,
  `jobs/numbat_sweep*.sh`; `results_controls/numbat_robustness/_analysis/`;
  `docs/lab_book/assets/04_numbat_reproducibility/`

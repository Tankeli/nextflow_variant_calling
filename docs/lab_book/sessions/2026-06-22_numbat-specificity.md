---
date: 2026-06-22
project: DDE_33
type: session
status: done
tags: [numbat, cnv, specificity, false-positive, healthy-baseline, controls, llr, analysis]
related: ["[[03_numbat_specificity]]", "[[02_copykat_robustness]]", "[[2026-06-09_patient-cohort-runs]]", "[[2026-06-08_controls-full-cohort]]"]
---

# 2026-06-22 — Numbat specificity on healthy controls (CopyKAT-robustness counterpart)

> **Objective:** turn the 2026-06-12 meeting observation *"Numbat — a lot of the time it's not finding
> anything in healthy"* into a quantified specificity analysis, the counterpart to the CopyKAT
> robustness work ([[02_copykat_robustness]]). Off the existing single Numbat runs — no Viking re-run.

## Context

CopyKAT was shown to over-call aneuploidy on every healthy control ([[02_copykat_robustness]]). The
open question was whether Numbat (the primary clonal axis) is more specific. The Numbat joint runs for
the 9 controls and one patient already existed on disk from [[2026-06-08_controls-full-cohort]] /
[[2026-06-09_patient-cohort-runs]]; this session only parses + characterises them.

## Work done

### 1. Data inventory — what Numbat actually emitted

Walked `results_{controls,patients,patients_dde22}/numbat_joint/*/numbat_out/`:

- **6/9 healthy controls have a truncated `numbat_out/`** (only `bulk_subtrees_1`, `hc.rds`,
  `gexp_roll_wide`, `log.txt`) — no `segs_consensus_*` / `clone_post_*` / `tree_final_*`. `log.txt`:
  *"No CNV remains after filtering by LLR in pseudobulks"*, all arms diploid. Correct negatives.
- **3/9 healthy (PBM_2, PBMMC_2, PBMMC_3)** completed to `clone_post`/`segs_consensus` — false positives.
- **Tumour comparator:** only **Patient_2** has a finished `numbat_out`. **Patient_1 is pileup-only**
  (`numbat_out/` absent); **all 6 DDE_22 patients with a `numbat_joint/<p>/` dir are pileup-only**
  (`numbat_out` empty). So `NUMBAT_RUN` only ever completed for Patient_2 — positive class n=1.

### 2. `bin/numbat_specificity.py`

New analysis script (`aml_scrna` conda), house-style match to `copykat_stability.py`. Auto-discovers
`numbat_joint/*/numbat_out`, classifies SILENT vs CALLED, and for CALLED reads the latest
`segs_consensus_*.tsv` (per-CNV `cnv_state_post`/`LLR`/`seg_length`) + `clone_post_*.tsv`
(`clone_opt`/`p_cnv`) → n_clones, aneuploid fraction, per-CNV LLR, CNV span, and an erythroid fraction
of the aneuploid cells from `reference_mapping/<s>/<s>_celltypes.csv` (the CopyKAT confound test).

```bash
source /opt/apps/eb/software/Miniconda3/23.5.2-0/etc/profile.d/conda.sh; conda activate aml_scrna
python3 bin/numbat_specificity.py --out-dir results_controls/numbat_specificity \
    --results healthy=results_controls --results tumour=results_patients
```

Outputs: `results_controls/numbat_specificity/numbat_specificity_summary.csv` + 3 figures
(`_status`, `_llr`, `_erythroid`), copied to `docs/lab_book/assets/03_numbat_specificity/`.

## Results / outcome

- **Specificity:** Numbat silent on **6/9** healthy controls (CopyKAT over-called all 9).
- **The 3 false-positives grade cleanly by LLR:** PBM_2 median LLR 4.6, PBMMC_2 8.1 (8 clones over
  1093 Mb = low-confidence over-segmentation) vs PBMMC_3 462.7 and Patient_2 (tumour) 558.4 — real
  events sit 50–100× above the noise. `min_LLR=3` is too permissive; post-filter on **median segment
  LLR** rather than tightening the caller (a hard LLR=20 would clip one of Patient_2's real segments).
- **Not the CopyKAT confound:** aneuploid cells are not erythroid-enriched (frac 0.01–0.14; top types
  CD4 T / CLP / PreB) — Numbat's failure mode is statistical over-segmentation, not lineage expression.
- **Flag:** PBMMC_3 (median LLR 463 on focal bamp/del) is tumour-grade in a "healthy" cord sample —
  needs segment-level follow-up; may be a genuine non-malignant CNV rather than a Numbat over-call.

Distilled into [[03_numbat_specificity]].

## Decisions & rationale

- **Decision:** analyse the existing single runs now rather than wait for a seed/threshold sweep.
  **Why:** the specificity question (silent vs called on healthy) is fully answerable off what's on
  disk and was the meeting ask; the CopyKAT-style reproducibility sweep is Viking-gated and a separate
  track. (rejected: block on a sweep — would have delayed a clean result with no new data.)
- **Decision:** keep Patient_2 as the lone positive comparator, caveated, rather than omit a positive
  class. **Why:** anchoring the tumour LLR scale (even at n=1) is what makes the 50–100× separation
  legible; flagged for re-confirmation as `NUMBAT_RUN` completes.

## Issues / blockers

- **Positive class n=1.** `NUMBAT_RUN` only finished for Patient_2; Patient_1 + DDE_22 are pileup-only.
  Re-run/finish those to confirm the tumour LLR scale ([[2026-06-09_patient-cohort-runs]]).
- **PBMMC_3 unexplained** (tumour-grade LLR on a healthy sample) — open item.

## Next steps

- [ ] Inspect PBMMC_3 `segs_consensus_2.tsv` + `bulk_clones` plots — genuine CNV vs artefact?
- [ ] Re-confirm tumour LLR scale once Patient_1 / DDE_22 `NUMBAT_RUN` complete (raise positive n).
- [ ] (Viking) Numbat seed/`min_LLR` sweep = the reproducibility counterpart to the CopyKAT sweep
      (`bin/numbat_specificity.py` already ingests any results dir, so the sweep just adds runs).
- [ ] Consider wiring an LLR-confidence post-filter / two-tier clone band into the Numbat output contract.

## Environment

- git: DDE_33 `main`@2521d03 + uncommitted (`bin/numbat_specificity.py`, lab book)
- conda `aml_scrna` (pandas/matplotlib); inputs from `numbat.sif` `run_numbat()` runs
  (`numbat_max_entropy=0.8`, `numbat_min_llr=3`, hg38) · account biol-stem-2022

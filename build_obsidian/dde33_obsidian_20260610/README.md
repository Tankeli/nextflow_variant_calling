# DDE_33 lab book

nf-core-style Nextflow pipeline formalising DDE_32's variant-calling scripts (Cell Ranger +
Numbat/CopyKAT/souporcell + scanpy annotation). Standard:
[`_workspace_admin/LAB_BOOK_STANDARD.md`](../../../_workspace_admin/LAB_BOOK_STANDARD.md).

> Entries reconstructed 2026-06-06 from `README.md`, `CLAUDE.md` and `scratchpad.md` (the build
> plan). `scratchpad.md` remains the live build tracker; the lab book records *what happened*.

> **Obsidian export** (generated 2026-06-10). Open this folder as an Obsidian vault. Lab-book
> notes are at the root (`analyses/`, `sessions/`, `assets/`); all pipeline output figures are
> embedded in browsable galleries under **[Figures index](Figures/INDEX.md)** (289 PNGs).

## Analyses

| # | Title | Status | Updated |
|---|---|---|---|
| [02](analyses/02_copykat_robustness.md) | CopyKAT robustness / reliability on healthy controls | active | 2026-06-10 |
| [01](analyses/01_pipeline_build_and_validation.md) | Pipeline build & stub validation (Phases 1–6) | active | 2026-06-03 |

## Sessions

| Date | Topic | Status |
|---|---|---|
| [2026-06-10](sessions/2026-06-10_copykat-robustness-sweep.md) | CopyKAT robustness sweep + downstream analysis (9 controls) | done |
| [2026-06-09](sessions/2026-06-09_patient-cohort-runs.md) | Patient cohorts — paired Dx/Rel AML (DDE_32 + DDE_22) | in-progress |
| [2026-06-08](sessions/2026-06-08_controls-full-cohort.md) | Healthy-control cohort expanded to 9 (HD_BM + PBM) | done |
| [2026-06-04](sessions/2026-06-04_first-real-run-caron-controls.md) | First real run — Caron 2020 healthy PBMMC controls | done |

## Headline outputs

- Stub-validated DAG: 4 CELLRANGER + 4 SCANPY_QC + 4 REFERENCE_MAPPING + 4 COPYKAT +
  2 NUMBAT_PILEUP/RUN + 2 SOUPORCELL_PREP + 4 SOUPORCELL (exit 0).
- Healthy controls (9/9, full DAG): `results_controls/` — CopyKAT calls 22–82% "aneuploid" on
  true-normal cells (the baseline behind analysis 02).
- CopyKAT robustness sweep (9/9, 180 runs): `results_controls/copykat_robustness/` — call not
  reproducible (18–98% of cells flip across 20 param×seed combos; seed-only switch up to 0.83),
  drivers confounded with lineage (8/9 sig.), summary fig
  `docs/lab_book/assets/02_copykat_robustness/robustness_summary.png`.
- Patient cohorts (in progress): `results_patients/` (DDE_32 prototype, 2 patients) +
  `results_patients_dde22/` (DDE_22 Vivobank, 10 patients) — first paired Dx/Rel malignant runs.
- Build tracker: [`../../scratchpad.md`](../../scratchpad.md) · architecture: [`../../CLAUDE.md`](../../CLAUDE.md)

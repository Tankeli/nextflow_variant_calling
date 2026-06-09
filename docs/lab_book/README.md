# DDE_33 lab book

nf-core-style Nextflow pipeline formalising DDE_32's variant-calling scripts (Cell Ranger +
Numbat/CopyKAT/souporcell + scanpy annotation). Standard:
[`_workspace_admin/LAB_BOOK_STANDARD.md`](../../../_workspace_admin/LAB_BOOK_STANDARD.md).

> Entries reconstructed 2026-06-06 from `README.md`, `CLAUDE.md` and `scratchpad.md` (the build
> plan). `scratchpad.md` remains the live build tracker; the lab book records *what happened*.

## Analyses

| # | Title | Status | Updated |
|---|---|---|---|
| [01](analyses/01_pipeline_build_and_validation.md) | Pipeline build & stub validation (Phases 1–6) | active | 2026-06-03 |

## Sessions

| Date | Topic | Status |
|---|---|---|
| [2026-06-04](sessions/2026-06-04_first-real-run-caron-controls.md) | First real run — Caron 2020 healthy PBMMC controls | in-progress |

## Headline outputs

- Stub-validated DAG: 4 CELLRANGER + 4 SCANPY_QC + 4 REFERENCE_MAPPING + 4 COPYKAT +
  2 NUMBAT_PILEUP/RUN + 2 SOUPORCELL_PREP + 4 SOUPORCELL (exit 0).
- First live results (PBMMC batch1): `results_controls/` (CellRanger, QC, refmap, CopyKAT,
  Numbat, souporcell — full DAG end-to-end).
- Build tracker: [`../../scratchpad.md`](../../scratchpad.md) · architecture: [`../../CLAUDE.md`](../../CLAUDE.md)

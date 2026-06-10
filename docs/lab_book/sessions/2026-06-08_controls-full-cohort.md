---
date: 2026-06-08
project: DDE_33
type: session
status: done
tags: [controls, caron, healthy-baseline, viking-profile, cellranger, copykat, numbat, souporcell, full-dag]
related: ["[[2026-06-04_first-real-run-caron-controls]]", "[[01_pipeline_build_and_validation]]", "[[02_copykat_robustness]]"]
---

# 2026-06-08 — Healthy-control cohort expanded to 9 samples (batch2: HD_BM + PBM)

> **Objective:** grow the healthy-control set from the 3 Caron PBMMC ([[2026-06-04_first-real-run-caron-controls]])
> to **9 GEX-only controls** spanning three sources (PBMMC, healthy BM, PBMC), so CopyKAT/Numbat
> behaviour on true-normal data has a real baseline cohort — the inputs to the CopyKAT robustness
> follow-up ([[02_copykat_robustness]]). Every stage must run end-to-end on the `viking` profile.

## Context

Batch1 (PBMMC_1, PBMMC_2) completed in [[2026-06-04_first-real-run-caron-controls]]; PBMMC_3 followed
once its download finished. CopyKAT over-called aneuploidy on PBMMC_2 (2,828 aneuploid / 663 diploid)
despite healthy input — to tell "tool over-call" from "real signal" apart we need a bigger normal
panel from independent datasets, not just one cohort. Batch2 adds 6 more controls from two further
public 10x sources.

## Data / provenance

Each control = its own "patient", `timepoint=Dx`, GEX-only (10x 3'). Raw FASTQ → `data/controls/`;
results → `results_controls/` (push to Longship for retention — scratch is purged after 90 d idle).

| Sample | Source (SRA) | Dataset | Notes |
|---|---|---|---|
| PBMMC_1 | SRR9264351 | Caron 2020 (GSE132509) | batch1 — paediatric BM mononuclear |
| PBMMC_2 | SRR9264353 | Caron 2020 | batch1 |
| PBMMC_3 | SRR9264354 | Caron 2020 | added after batch1 |
| HD_BM_1 | SRR12185508 | healthy adult BM (10x) | original 10x BAM → FASTQ (`0064.bam`) |
| HD_BM_2 | SRR12185509 | healthy BM | `3958.bam` |
| HD_BM_3 | SRR12185510 | healthy BM | `5286.bam` |
| HD_BM_4 | SRR12185511 | healthy BM | `7903.bam` |
| PBM_1 | SRR12338699–706 (8 runs) | healthy PBMC (10x) | multi-run, concatenated; R1 28 bp / R2 91 bp |
| PBM_2 | SRR12338707–714 (8 runs) | healthy PBMC | multi-run, concatenated |

## Work done

### 1. Fetch batch2 FASTQ (login/transfer node — outbound internet)
- **What:** HD_BM_* pulled as *original 10x BAM* from `sra-pub-src-1` then re-emitted to cellranger
  FASTQ; PBM_* via `prefetch` + `fasterq-dump --split-files --include-technical` over 8 runs each,
  concatenated per sample. Read-length classifier confirmed (I1≤12, R1≤30, else R2).
- **Command:**
  ```bash
  bash jobs/fetch_controls_batch2.sh   # bin/fetch_sra_10x.sh data/controls HD_BM_1:SRR12185508 ... PBM_1:SRR12338699,...,706
  ```
- **Log:** `logs/fetch_controls_batch2_20260608_115841.log` — all 6 samples `done: .../R{1,2}_001.fastq.gz`.

### 2. Live run — batch2 (HD_BM_1-4, PBM_1-2)
- **What:** full caller + annotation DAG, `viking` profile, `-resume` (batch1/PBMMC_3 cellranger
  work reused). Samplesheet `assets/controls_samplesheet_batch2.csv`; combined 9-sample sheet is
  `assets/controls_samplesheet_all9.csv`.
- **Command:**
  ```bash
  sbatch jobs/run_controls.sh assets/controls_samplesheet_batch2.csv   # orchestrator off the login node
  # → nextflow run . -profile viking -params-file params-controls.yaml \
  #     --input assets/controls_samplesheet_batch2.csv -work-dir work -resume
  ```
- **Job:** orchestrator Slurm **34532455** · 8G · `2026-06-08T22:33 → 2026-06-09T07:03` ·
  **8h 30m**, **Succeeded 67**, exit 0 · log `logs/controls_orchestrator_34532455.log`. (Per-stage
  CellRanger/CopyKAT/Numbat/souporcell jobs are dispatched by Nextflow as their own SLURM jobs.)
- **Result:** ✅ full DAG end-to-end for all 6 batch2 samples. With batch1 + PBMMC_3 the cohort is
  **9/9 complete** under `results_controls/` (cellranger, qc, reference_mapping, copykat,
  numbat_joint, souporcell all present for every sample).

## Results / outcome

**CopyKAT aneuploid/diploid split across the full 9-control cohort** (`results_controls/copykat/<s>/<s>_copykat_prediction.txt`):

| Sample | total | aneuploid | diploid | aneu % |
|---|---|---|---|---|
| PBMMC_1 | 1001 | 223 | 677 | 22% |
| PBMMC_2 | 5035 | 2830 | 661 | 56% |
| PBMMC_3 | 5491 | 1214 | 1829 | 22% |
| HD_BM_1 | 1798 | 1170 | 386 | 65% |
| HD_BM_2 | 1089 | 612 | 289 | 56% |
| HD_BM_3 | 1633 | 1012 | 369 | 62% |
| HD_BM_4 | 1790 | 371 | 1139 | 21% |
| PBM_1 | 6181 | 4422 | 1602 | 72% |
| PBM_2 | 7227 | 1704 | 5180 | 24% |

(Counts exclude "not.defined" cells, so aneu+dip < total per row.) On **healthy** input CopyKAT calls
**21–72% of cells "aneuploid"** with no consistent baseline — the headline that motivates
[[02_copykat_robustness]]. This is tool behaviour, not a pipeline bug: every prediction file was
produced and is well-formed; the variability is in CopyKAT's call, not in the DAG.

## Decisions & rationale

- **Decision:** add HD_BM (adult BM) + PBM (PBMC) from independent datasets rather than only more
  PBMMC. **Why:** a robustness baseline needs cross-dataset normals so the over-call pattern can't be
  blamed on one cohort's chemistry/donor. (rejected: scale up Caron PBMMC only.)
- **Decision:** keep CloneTracer + Phase-2 integration **off** for controls (`run_clonetracer:false`,
  `run_integration:false` in `params-controls.yaml`). **Why:** no malignant clones to trace on healthy
  data; CloneTracer's default GTF path isn't wired for this run.

## Issues / blockers

- **CARRIED FROM [[2026-06-04_first-real-run-caron-controls]]:** trace/report/timeline/dag +
  co2footprint still disabled in the `viking` profile (`numbat.sif` lacks `procps`/`ps`). Add
  `procps` to the image to re-enable. Not blocking.
- No new failures in batch2.

## Next steps

- [x] CopyKAT robustness follow-up over the 9-control baseline → [[02_copykat_robustness]].
- [ ] Push `results_controls/` to Longship for retention (`../snapshot_dde_to_longship.sh`).
- [ ] Add `procps` to `numbat.sif` → re-enable reporting plugins on the `viking` profile.

## Environment

- git: DDE_33 `main`@7b57c25
- profile: `viking` — CellRanger module 9.0.0, SAMtools module, conda `snv` (CopyKAT) + `aml_scrna`
  (scanpy 1.11.5), `numbat.sif`, `souporcell_release.sif`
- params: `params-controls.yaml` (refs under `/mnt/scratch/users/hbp534/references/refdata-gex-GRCh38-2024-A`;
  atlas = DDE_32 paediatric BM h5ad) · account biol-stem-2022 · outdir `results_controls/`, work-dir `work/`

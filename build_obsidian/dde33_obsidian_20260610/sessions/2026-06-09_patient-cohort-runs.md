---
date: 2026-06-09
project: DDE_33
type: session
status: in-progress
tags: [patients, aml, dx-rel, paired, cite-seq, numbat, souporcell, copykat, integration, lsc, viking-profile]
related: ["[[2026-06-08_controls-full-cohort]]", "[[01_pipeline_build_and_validation]]", "[[02_copykat_robustness]]", "[[DDE_32]]", "[[DDE_22]]"]
---

# 2026-06-09 — First malignant-data runs: paired Dx/Rel AML patient cohorts

> **Objective:** run the pipeline on **real paired diagnosis–relapse AML** (CITE-seq, GEX+ADT) for
> the first time — the use case it exists for. Two concurrent runs: the DDE_32 4-sample prototype
> (2 patients) and the larger DDE_22 Vivobank cohort (16 samples, 10 patients, mixed paired/standalone).
> Success = full caller set + annotation + Phase-2 integration (per-patient master table + Dx→Rel
> clonal-tracing Sankeys) on genuine malignant clones.

## Context

All prior live runs were healthy controls ([[2026-06-08_controls-full-cohort]]) — biologically
uninformative for clones, only a DAG exercise. This is the first time joint-per-patient Numbat +
souporcell-paired calling runs on data with real clonal structure, and the first exercise of the new
integration/LSC stages. Patient→sample/timepoint grouping comes entirely from the samplesheet
`patient`/`timepoint` columns (no hardcoded maps).

The two runs use **isolated launch dirs + separate work dirs** (`work_patients` / `work_patients_dde22`)
so they run concurrently without sharing `work/`/`.nextflow` caches.

## Data / provenance

**Cohort A — DDE_32 prototype** (`assets/patients_samplesheet.csv`), CITE-seq GEX+ADT, FASTQ under
`/mnt/scratch/users/hbp534/DDE_32_paediatric_snv_analysis/data/{GEX,CSP}`:

| Patient | Dx | Rel |
|---|---|---|
| Patient_1 | Sample_2395 | Sample_3001 |
| Patient_2 | Sample_2977 | Sample_0109 |

**Cohort B — DDE_22 Vivobank** (`assets/patients_dde22_samplesheet.csv`, 16 libraries → 10 patients;
mapping from `assets/Paed_AML_Vivobank_details.csv`). FASTQ are cellranger-compliant symlinks under
`assets/fastq_dde22/{GEX,CSP}` → `DDE_22_clean_nextflow_basic_pipeline/data`:

| Patient | Samples | Pairing |
|---|---|---|
| AML225 | 0984 + 0984-B (Dx), 4230 (Rel) | paired |
| AML107 | 8178 (Dx), 1894 (Rel) | paired |
| AML066 | 3652 (Dx), 1386 (Rel) | paired (CSV typo 3562→FASTQ 3652) |
| AML079 | 2765 (Dx), 1255 (Rel) | paired |
| AML065 | 5807 + 5807-B (Dx) | standalone |
| AML155 / AML161 / AML124 / AML163 | 2958 / 1199 / 3280 / 1187 (Dx) | standalone |
| AML104 | 8087 (Rel only — Dx 2243 has no FASTQ) | standalone |

> Cohort B deliberately includes single-timepoint patients — per-patient steps + integration figures
> must handle standalone samples, not assume Dx+Rel (see [[cohort-not-all-paired]]).

## Work done

### 1. Cohort A — DDE_32 prototype (4 samples / 2 patients)
- **Command:**
  ```bash
  sbatch jobs/run_patients.sh
  # → nextflow run . -profile viking -params-file params-patients.yaml \
  #     -work-dir work_patients -resume   (isolated launchDir .launch_patients)
  ```
- **Job:** orchestrator Slurm **34532466** · 8G · started `2026-06-08T22:37` · **RUNNING** (~17h47m at
  time of writing) on node053 · log `logs/patients_orchestrator_34532466.log`.
- **Result (progress):** CellRanger 4/4 ✔, SCANPY_QC 4/4 ✔, REFERENCE_MAPPING 4/4 ✔, COPYKAT 4/4 ✔,
  PLOT_COPYKAT 4/4 ✔, LSC_SCORING 4/4 ✔, COHORT_SUMMARY ✔, **NUMBAT_PILEUP 2/2 ✔**. Still running:
  **NUMBAT_RUN 0/2**, **SOUPORCELL_PREP 0/2** → SOUPORCELL (K sweep `2,3,5,12` — the long pole) →
  PHASE0_INTEGRATION + HEADLINE_FIGURES (pending the callers).
  - Numbat joint pileup published per sample, e.g.
    `results_patients/numbat_joint/Patient_1/Patient_1_pileup/Sample_2395_allele_counts.tsv.gz`
    (+ Sample_3001; Patient_2 = Sample_2977 + Sample_0109).
  - LSC scores per sample: `results_patients/lsc_scoring/Sample_2395/Sample_2395_lsc.csv` (×4).

### 2. Cohort B — DDE_22 Vivobank (16 samples / 10 patients)
- **Command:**
  ```bash
  sbatch jobs/run_patients_dde22.sh
  # → nextflow run . -profile viking -params-file params-patients-dde22.yaml \
  #     -work-dir work_patients_dde22 -resume   (launchDir .launch_patients_dde22)
  ```
- **Job:** orchestrator Slurm **34564311** · 8G · started `2026-06-09T14:01` · **RUNNING** (~2h24m) on
  node064 · log `logs/patients_dde22_orchestrator_34564311.log`. (Earlier attempt 34525713.)
- **Result (progress):** early — **CELLRANGER_MULTI 5 of 16** (5 cached from a prior attempt, 11
  counting now as their own SLURM `VARIANTC` jobs), and the 5 ready samples already through COPYKAT /
  QC / REFERENCE_MAPPING / LSC_SCORING. Numbat/souporcell not yet started (await cellranger). Done so
  far: Sample_0984, 1187, 1199, 2958, 8178 (cellranger→copykat→qc→refmap→lsc).

### 3. CopyKAT first-pass calls (malignant data, real signal expected)
`results_patients*/copykat/<sample>/*_copykat_prediction.txt`:

| Sample (patient, tp) | total | aneuploid | diploid |
|---|---|---|---|
| Sample_2395 (P1 Dx) | 5499 | 2706 | 2112 |
| Sample_3001 (P1 Rel) | 9090 | 3346 | 5391 |
| Sample_2977 (P2 Dx) | 1132 | 771 | 252 |
| Sample_0109 (P2 Rel) | 4388 | 2013 | 2138 |
| Sample_0984 (AML225 Dx) | 2232 | 1587 | 456 |
| Sample_8178 (AML107 Dx) | 2002 | 1211 | 709 |
| Sample_1199 (AML161 Dx) | 2387 | 1467 | 774 |
| Sample_2958 (AML155 Dx) | 875 | 659 | 154 |
| Sample_1187 (AML163 Dx) | 609 | 197 | 351 |

Unlike healthy controls, an aneuploid fraction here is *expected* (malignant clones) — but the
control baseline ([[02_copykat_robustness]]) shows the raw split alone can't distinguish true
aneuploidy from CopyKAT's normal over-call. Numbat (primary clonal axis) + souporcell are the
arbiters; both still pending in cohort A.

## Results / outcome

Both patient runs **launched and progressing on real malignant CITE-seq** — neither has reached the
clonal checkpoints yet (Numbat run, souporcell, integration). No interpretation of clones is possible
until souporcell + NUMBAT_RUN complete. This entry is a run-log; clone/Sankey results will be distilled
into a separate analysis once cohort A finishes.

## Decisions & rationale

- **Decision:** run cohorts A and B as two independent orchestrators with separate work dirs. **Why:**
  cohort B is large and slow; isolating it keeps the `-resume` cache and launch state from colliding
  with the smaller prototype run. (rejected: one combined samplesheet — couples their schedules.)
- **Decision:** souporcell K sweep `2,3,5,12` for patients (vs `2,3` for controls). **Why:** real
  clonal substructure may need more clusters; sweeping lets us pick K post-hoc per patient. Cost: K=12
  is the slowest stage and gates integration.
- **Decision:** CloneTracer off for both (`run_clonetracer:false`). **Why:** default GTF not at the
  refs path and the env is still being sorted ([[clonetracer-env-slow]]); enable once pinned.

## Issues / blockers

- **In flight, not failed:** cohort A is gated on NUMBAT_RUN + the K=12 souporcell job; cohort B is
  still in CellRanger. Neither has produced clone calls or Phase-2 Sankeys yet.
- **WATCH — maintenance window YOR796 (2026-06-11 09:00 → 06-15 09:00):** locks ~all nodes for 4 days.
  Anything not finished before 06-11 09:00 must resume after 06-15 via `-resume`. Cohort B (16 samples,
  started 06-09 14:01) is at risk of straddling it.

## Next steps

- [ ] Cohort A: let NUMBAT_RUN + souporcell (K sweep) + PHASE0_INTEGRATION + HEADLINE_FIGURES finish;
      verify `results_patients/` has clone posteriors + per-patient master table + Dx→Rel Sankeys.
- [ ] Distil cohort A clonal results (Numbat clones, souporcell genotype clusters, Sankeys) into a new
      `analyses/` note once complete.
- [ ] Cohort B: monitor across the YOR796 window; `-resume` after 06-15 if it straddles.
- [ ] Push `results_patients*/` to Longship after each cohort completes.
- [ ] Enable CloneTracer once the pinned env/GTF are wired ([[clonetracer-env-slow]]).

## Environment

- git: DDE_33 `main`@7b57c25
- profile: `viking` — CellRanger module 9.0.0, SAMtools module, conda `snv` + `aml_scrna`
  (scanpy 1.11.5), `numbat.sif`, `souporcell_release.sif`
- params: `params-patients.yaml` / `params-patients-dde22.yaml` (CITE-seq feature refs from
  DDE_32 / DDE_22 `config/feature_reference.csv`; atlas = DDE_32 paediatric BM h5ad; frozen seed123
  UMAP via `refmap_umap`) · account biol-stem-2022 · outdirs `results_patients/`, `results_patients_dde22/`

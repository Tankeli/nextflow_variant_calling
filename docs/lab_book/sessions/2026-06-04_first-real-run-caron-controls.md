---
date: 2026-06-04
project: DDE_33
type: session
status: in-progress
tags: [first-run, caron-controls, viking-profile, cellranger, full-dag]
related: ["[[01_pipeline_build_and_validation]]"]
---

# 2026-06-04 — First real run: Caron 2020 healthy PBMMC controls

> **Objective:** exercise the whole pipeline end-to-end on real data for the first time using 3
> healthy paediatric PBMMC controls (GEX-only), on a `viking` profile that uses cluster
> modules/conda/SIFs instead of built containers. Numbat/souporcell are uninformative on healthy
> data but must run the DAG; CopyKAT should call ~all diploid.

> Reconstructed 2026-06-06 from `scratchpad.md` Phase 7. Live job IDs are as recorded there.

## Context

Phases 1–6 were stub-validated only ([[01_pipeline_build_and_validation]]). This is the first run
with real compute, real references, and real environments — the true test of the `viking` profile.
Controls: Caron et al. 2020 (GSE132509 / SRP201012), 10x 3' v2, GEX-only.

## Data / provenance

| Sample | Patient | Modality | Source (SRA) | Notes |
|---|---|---|---|---|
| PBMMC_1 | PBMMC_1 | GEX (10x 3' v2) | SRR9264351 (GSM3872442) | healthy paediatric BM mononuclear |
| PBMMC_2 | PBMMC_2 | GEX | SRR9264353 (GSM3872443) | each control = its own "patient", timepoint=Dx |
| PBMMC_3 | PBMMC_3 | GEX | SRR9264354 (GSM3872444) | still downloading at session end |

Raw FASTQ → `data/controls/`; results → `results_controls/` (push to Longship for retention).

## Work done

### 1. `viking` profile (no container builds)
- **What:** `conf/viking.config` + `-profile viking` — CellRanger module 9.0.0, SAMtools module,
  conda `snv` (CopyKAT), conda `aml_scrna` (scanpy+scrublet), `numbat.sif` (DDE_32) +
  `souporcell_release.sif` (DDE_24) via apptainer.
- **Result:** confirmed `aml_scrna` has scanpy 1.11.5; `snv` has copykat; both R envs 4.5.3.

### 2. GEX-only support + SRA fetch
- **What:** made `fb_reference` optional (`[]` when null; module emits the feature section only
  when given); `bin/fetch_sra_10x.sh` (prefetch + fasterq-dump `--split-files --include-technical`,
  classify reads by length, rename to cellranger convention).
- **Command:**
  ```bash
  sbatch bin/fetch_sra_10x.sh   # writes data/controls/
  ```
- **Job:** Slurm **34329477** · logs `logs/sra_pbmmc_*.log`.
- **Result:** read-length classification correct (I1=8, R1=26 [v2], R2=98). PBMMC_1 + PBMMC_2
  FASTQs complete; PBMMC_3 still downloading.

### 3. Live run — batch1 (PBMMC_1, PBMMC_2)
- **Command:**
  ```bash
  sbatch jobs/run_controls.sh   # orchestrator on a compute node, not the login node
  # → nextflow run . -profile viking -params-file params-controls.yaml \
  #     --input assets/controls_samplesheet_batch1.csv -work-dir work -resume
  ```
- **Job:** orchestrator Slurm **34341569** (node036); CELLRANGER_MULTI **34341670 / 34341672**.
  After the fastqs-path bugfix, resubmit **34347849** → cellranger **34347864 / 34347866**.
- **Result:** ✅ **full DAG end-to-end, exit 0** — 12 succeeded + 6 cached, **3h36m**. Outputs
  under `results_controls/`:
  - CellRanger ✅ `possorted_genome_bam.bam` ~19.6 GB + filtered matrix (both)
  - SCANPY_QC ✅ `*_qc.h5ad` + `_qc_metrics.csv` (aml_scrna conda branch works)
  - REFERENCE_MAPPING ✅ `*_celltypes.csv` + `_mapped.h5ad` (PBMMC_1 ~CLP/T/DC — plausible PBMMC)
  - COPYKAT ✅ `*_prediction.txt` (both)
  - NUMBAT pileup ✅ `*_allele_counts.tsv.gz`; run: PBMMC_1 "No CNV remains after LLR filtering"
    → no clones (correct for healthy); PBMMC_2 called clones ("All done!")
  - SOUPORCELL ✅ `clusters.tsv` + `cluster_genotypes.vcf` per K

## Results / outcome

First real end-to-end run **passes** on the `viking` profile. The biological behaviour on healthy
controls is the expected tool behaviour, not pipeline bugs: CopyKAT over-calls aneuploidy on
**PBMMC_2 (2,828 aneuploid / 663 diploid)** and Numbat called clones on PBMMC_2 — expected when
healthy controls lack a matched-normal baseline + relaxed thresholds (`min_LLR=3`). PBMMC_1
behaves as expected (mostly diploid, no clones). This is exactly the value of a healthy-control run.

## Decisions & rationale

- **Decision:** run controls GEX-only, each as its own "patient". **Why:** Caron data has no
  antibody capture and no Dx/Rel pairs; joint logic degrades to per-sample, which is fine
  mechanically and still exercises the DAG. (rejected: skip Numbat/souporcell on healthy — we want
  to prove they *run*.)
- **Decision:** disable trace/report/timeline/dag + co2footprint in the `viking` profile for now.
  **Why:** the plugin needs `ps`; `numbat.sif` lacks `procps`. (rejected: leave enabled → run fails
  at shutdown.)

## Issues / blockers

- **FIXED:** CellRanger preflight rejected the relative fastqs path ("must be an absolute path:
  fastqs"). Fix: CELLRANGER_MULTI writes an absolute `[libraries]` path via
  `fastqs_dir=$(readlink -f fastqs)`. Resubmitted (34347849) → preflight + chemistry auto-detect
  (v2) + barcode-compat pass.
- **RISK:** SRA prefetch needs outbound internet on the compute node; if blocked, run prefetch on
  a login/transfer node and fasterq-dump on compute. Watch job 34329477 / PBMMC_3.

## Next steps

- [ ] When PBMMC_3 download finishes: `sbatch jobs/run_controls.sh assets/controls_samplesheet.csv`
      (all 3) — `-resume` reuses PBMMC_1/2 cellranger work.
- [ ] Add `procps` to `numbat.sif` (apptainer overlay or rebuild) → re-enable trace/report/
      timeline/dag + co2footprint.
- [ ] For controls-as-baseline: consider a normal reference / stricter thresholds for CopyKAT &
      Numbat (see DDE_32 notebook 08 healthy-ref work).
- [ ] Push `results_controls/` to Longship once PBMMC_3 batch completes.

## Environment

- git: DDE_33 `main`@b771178 (working tree had local mods to configs/modules at session time)
- profile: `viking` — CellRanger module 9.0.0, SAMtools module, conda `snv` + `aml_scrna`
  (scanpy 1.11.5), `numbat.sif`, `souporcell_release.sif`
- account: biol-stem-2022 · outdir `results_controls/`, work-dir `work/`

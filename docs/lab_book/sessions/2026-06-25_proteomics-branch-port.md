---
date: 2026-06-25
project: DDE_33
type: session
status: done
tags: [proteomics, bulk-ms, DESP, deconvolution, R-to-python, limma, port, nextflow, f157, report]
related: ["[[2026-06-24_clonetracer-first-runs]]", "[[proteomics-branch]]", "[[user-prefers-python]]"]
---

# 2026-06-25 — Bulk-proteomics branch: DDE_31 port (R→Python) + first F157 run submitted

> **Objective:** integrate the DDE_31 paediatric bulk-proteomics pipeline (F157, Spectronaut DIA on
> timsTOF) into DDE_33 as an optional, gated branch — **full port, R→Python** except the DESP
> cell-state demix (kept as one R step) — wire DESP's cell-type proportions to this pipeline's scRNA
> reference-mapping output (the plan's Module-C proteogenomic hook), and submit a first real run.

## Context

The scRNA-seq plan's P6 (Module C, proteogenomic anchoring) was zero-code. DDE_31's bulk-MS workflow
(23 R scripts: QC → batch → viz → DE → stage-4 → ML → DESP) is the natural seed. Decisions (user):
**full port**, **DESP kept as one R step**, **proportions derived from our scRNA pipeline**. Input is
a Spectronaut search matrix (NOT FASTQ), so the branch runs parallel to everything and also as a
standalone `-entry PROTEOMICS`. Distinct from the existing CITE-seq ADT `prot_*` modules — this is
`prot_ms_*` (bulk MS).

## Work done

### 1. R→Python port (every stage)

- Shared infra: `bin/prot_ms_utils.py` (config deep-merge + IO + design + imputation + DE filters,
  port of utils.R/imputation.utils.R), `bin/prot_ms_plotting.py` (matplotlib/seaborn port of
  plotting.utils.R), `assets/proteomics_default.yaml` (params only; path-routing dropped — each stage
  takes explicit `--in/--out`).
- Stages: `prot_ms_prep.py` (1a), `prot_ms_batch.py` (1b — **limma `removeBatchEffect` reimplemented**
  + ComBat via inmoose), `prot_ms_de.py` (3a/3b — **limma `lmFit`/`eBayes` reimplemented** incl.
  squeezeVar/fitFDist/trigammaInverse + volcano + per-patient logFC), `prot_ms_viz.py` (2a/2b; UMAP via
  scanpy/DPT — Seurat/slingshot not reproduced), `prot_ms_stage4.py` (4a–4f; pseudotime = DPT+Spearman
  stand-in, tradeSeq omitted), `prot_ms_ml.py` (5a/5b — sklearn DT + RF impurity importance),
  `prot_ms_proportions.py` (scRNA celltypes → cell_type×sample), `prot_ms_desp_viz.py` (6a/6c/6e).
- DESP kept in R: `bin/prot_ms_desp_run.R` wraps `DESP::DESP`, condition-stratified + per-patient.

### 2. Nextflow wiring

9 modules `modules/local/prot_ms_*.nf` + `subworkflows/local/proteomics.nf` (PREP→BATCH spine;
DE/VIZ/STAGE4/ML; PROPORTIONS→DESP→DESP_VIZ; gated `prot_run_*`) + `-entry PROTEOMICS` (main.nf) +
gated `run_proteomics` in `workflows/variantcalling.nf` (feeds `REFERENCE_MAPPING` celltypes to
proportions). Config/schema/containers/envs: `nextflow.config` params (default OFF),
`nextflow_schema.json` `proteomics_options`, `conf/viking.config` env hooks, `conf/modules.config`
publishDir → `results/proteomics/{01_qc..07_desp}`, `containers/{proteomics,desp}/Dockerfile`,
`envs/{proteomics,desp_r}.yml`.

### 3. Environments + DESP verified

```bash
conda env create -f envs/proteomics.yml -n aml_proteomics   # scanpy 1.11.5 + inmoose verified
conda env create -f envs/desp_r.yml     -n desp_r
Rscript -e 'remotes::install_github("davidsjoberg/ggsankey")'   # DESP dep, GitHub-only
Rscript -e 'remotes::install_github("AhmedYoussef95/DESP")'     # the real DESP repo (v1.0)
```

DESP source-confirmed API `DESP(bulk, proportions, lambda, beta, similarities)`: `bulk`=features×
samples, **`proportions`=samples×cell-state**, output=features×cell-state — matches the wrapper's
`t(prop_mat)`. (Earlier placeholder repo path was wrong; fixed.)

### 4. Bugs caught by running on real data (not just stub)

- **Numeric sample IDs** `2977`/`109` read as integers → indexed columns *by position* in R
  (`subscript out of bounds`) → `colClasses="character"` in the wrapper + a `reconcile_sample_ids`
  helper mapping design `0109` ↔ matrix `109` (R's silent int-coercion).
- STAGE4 `RecursionError` on the full-cohort heatmap → `setrecursionlimit` + cap clustermap rows.
- ML `permutation_importance` impractical over ~6500 features → RF impurity importance.
- R colour `grey28` invalid in matplotlib; `ward`+correlation disallowed by scipy → both handled.

### 5. Validation

```bash
nextflow run . -entry PROTEOMICS -profile test -stub -ansi-log false ...   # 8 procs, exit 0
nextflow run . -profile test -stub --run_proteomics ...                    # gated: 9 procs incl.
                                                                            # PROPORTIONS ← REFERENCE_MAPPING
```

Functional run on **real F157 data** (`aml_proteomics` + `desp_r`): PREP **6481 proteins × 12
samples**, BATCH, DE (6481 + per-patient), VIZ, STAGE4, ML (4 s), and the **real DESP = 6481 features
× 54 cell types** over 6 overlap samples + per-patient AML152 — **exactly reproduces DDE_31's
dimensions**.

### 6. First F157 run — submitted, debugged, COMPLETED

No batch correction (single batch, per user) → `batch_correction.method=none` +
`prot_desp_bulk_source=raw`; DESP uses the pre-computed DDE_31 transcriptomics-overlap proportions so
it runs now without waiting on AML-patient reference-mapping.

```bash
sbatch jobs/run_proteomics.sh    # -entry PROTEOMICS -profile viking -params-file params-proteomics.yaml
```

Two **Viking-only** bugs (not caught by local/stub, both staging-related under `bash -ue`):
1. job `35219396` — `PYTHONPATH: unbound variable` (Nextflow runs `bash -ue`; `$PYTHONPATH` unset on
   compute) → `${PYTHONPATH:-}` in all 8 Python modules.
2. job `35273987` — `PROT_MS_VIZ input file name collision: matrix_raw.tsv` (with `bulk_source=raw`
   the `raw` and `corrected` inputs are the same file) → `stageAs` distinct names on VIZ inputs.

**job `35274762` COMPLETED** (exit 0, 5 new tasks + 3 cached, ~no walltime issues). 73 files →
`results_proteomics/proteomics/{01_qc..07_desp}`.

## Results

- Code: `bin/prot_ms_*` (11 files), `modules/local/prot_ms_*.nf` (9), `subworkflows/local/proteomics.nf`,
  entry + gated wiring, config/schema/envs/containers, `params-proteomics.yaml`,
  `assets/proteomics_f157.yaml`, `assets/proteomics_sample_map.tsv`.
- Envs `aml_proteomics` + `desp_r` (DESP 1.0) on the login node.
- **Run 35274762 outputs** (`results_proteomics/proteomics/`): 7603→**6481 proteins × 12 samples**;
  DE **0 proteins at adj.P<0.05** (no-batch baseline; top |logFC| = keratins/LTF = handling
  contamination); RF/DT classifiers (exploratory, n=12); **DESP 6481 × 54 cell states** + Relapse−Dx
  delta + paired AML152. Honest read: no cohort-level Dx/Rel signal without patient-batch modelling;
  DESP is the value-add (cell-state-resolved layer).
- **Obsidian report**: `docs/reports/proteomics_f157/{report.md, figures/}` — self-contained (18
  curated figure copies, relative links), transferable into the vault.

## Notes / open items

- Code ran at `main@2521d03` + uncommitted changes (this whole branch is uncommitted).
- This run uses external (DDE_31) proportions. **Full in-pipeline scRNA-derived proportions** is the
  follow-up: needs AML-patient `REFERENCE_MAPPING` celltypes published, then drop
  `proteomics_proportions` and use `--proteomics_celltypes_glob` + `assets/proteomics_sample_map.tsv`
  (or run inside the main pipeline with `run_reference_mapping=true`).
- Plan P6 *advanced* (fusion-junction peptide DB + bulk↔single-cell consistency) still TODO
  (scratchpad Phase 15b) — needs the fusion branch + raw spectra.
- Watch job 35219396; on completion record DESP/DE/figure outputs + walltime here.

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

`DDE_33_nextflow_variant_calling` is a **formalized nf-core-style Nextflow pipeline** for
single-cell variant calling in paired diagnosis–relapse paediatric AML samples. It replaces a
set of fragmented, ad-hoc SLURM/conda/apptainer wrapper scripts that were developed during the
exploratory analysis in the sibling project `../DDE_32_paediatric_snv_analysis`.

The pipeline takes scRNA-seq (CITE-seq: GEX + Antibody Capture) FASTQ through Cell Ranger and
three variant-calling tools, emitting standardized, resumable per-sample / per-patient variant
checkpoints. The downstream multi-modal integration (Phase-0 per-cell master table, clonal-tracing
Sankey figures) currently stays in DDE_32 and consumes this pipeline's outputs.

> As of the start of this work the directory is empty apart from this file and `scratchpad.md`
> (the build plan). Read `scratchpad.md` for the current implementation phase and task status.

## Scientific context (why these tools)

Goal: identify leukaemic stem cell (LSC) populations and trace clonal evolution across paired
diagnosis (Dx) and relapse (Rel) bone-marrow samples. Clones are defined by **genotype** (CNV +
SNV + mtDNA), and phenotypes (LSC scores, cell types) are overlaid downstream. The prototype
cohort is two patients, each with a Dx + Rel pair:

| Patient | Diagnosis | Relapse |
|---|---|---|
| Patient_1 | Sample_2395 | Sample_3001 |
| Patient_2 | Sample_2977 | Sample_0109 |

Joint-per-patient calling (pooling Dx + Rel cells in one run) is essential for Numbat and
souporcell so clone IDs are comparable across timepoints — per-sample clone IDs are not.

## Pipeline scope (agreed decisions)

- **Boundaries**: FASTQ → variant checkpoints + a parallel cell-annotation branch. Cell Ranger is
  wrapped upstream; the multi-modal integration (Phase-0 master table, Sankeys) stays in DDE_32.
- **Callers included**: Numbat (CNV, primary clonal axis), CopyKAT (expression-based aneuploidy
  gate), souporcell (SNV genotype clusters), and **CloneTracer** (veltenlab) as a *downstream*
  Bayesian clonal-integration branch — it is not a de-novo caller; it consumes per-cell M/N counts
  derived from the callers (Numbat CNVs + souporcell SNVs + a new mtDNA pileup) and emits clone
  trees + per-cell clone posteriors. **mgatk is not a standalone stage** (deselected; silent
  failures on Sample_3001 in DDE_32) but is wired as the *opt-in* mtDNA method for CloneTracer
  (`clonetracer_mtdna_method=mgatk`); the default mtDNA path is cellsnp-lite.
- **Annotation branch (Python/scanpy)**: scanpy/Scrublet QC + reference mapping (`scanpy.tl.ingest`)
  run *parallel* to the callers off the raw cellranger matrices — they do NOT gate caller inputs.
  Ported from the DDE_23 scanpy stack (the user prefers Python over the DDE_32 R/Seurat scripts).
  Atlas is configurable (default DDE_32 paediatric BM `bone_marrow_atlas.h5ad`; Zeng BoneMarrowMap
  or others swap in via `--refmap_atlas`).
- **Structure**: full nf-core layout (`modules/`, `subworkflows/`, `conf/`, `nextflow_schema.json`,
  samplesheet validation).
- **Seqera Platform**: "launchable + config hooks" — `manifest{}`/`tower{}` blocks, launch schema,
  `tower.yml`, retain nf-prov BCO + co2footprint reporting. Workspace/compute-env registration is
  done by the user, not wired here.

## Stages and checkpoints

Each stage is resumable with a stable output contract under `results/`:

1. **CELLRANGER** (`cellranger multi`, GEX + Antibody Capture) → `results/cellranger/<sample>/outs/`
   (`possorted_genome_bam.bam` + `filtered_feature_bc_matrix/`). Vendored module rather than calling
   `nf-core/scrnaseq` as a sub-pipeline.
2. **NUMBAT** → `NUMBAT_PILEUP` (per-sample and joint-per-patient `pileup_and_phase.R`) →
   `*_allele_counts.tsv.gz`; then `NUMBAT_RUN` (`run_numbat()`, relaxed thresholds
   `max_entropy=0.8`, `min_LLR=3`) → `results/numbat_joint/<patient>/numbat_out/`.
3. **COPYKAT** (per-sample) → `results/copykat/<sample>/*_copykat_prediction.txt`.
4. **SOUPORCELL** (optional noNK barcode subset → per-patient CB-retag + merge + sort →
   `souporcell_pipeline.py` over a K sweep) → `results/souporcell/<patient>/k<K>/clusters.tsv`.
5. **CLONETRACER** (downstream, joint per patient) → `MTDNA_PILEUP` (per sample, cellsnp-lite on
   chrM; mgatk opt-in) + `CLONETRACER_BUILD` (synthesise per-cell M/N over CNV/SNV/mtDNA →
   `<patient>.json`) + `CLONETRACER` (`run_clonetracer.py`, pyro) →
   `results/clonetracer/<patient>/<patient>_clone_assignments.csv` + `_out.pickle` + `_tree.pickle`.
   Gated by `run_clonetracer`; uses whatever caller sources are available (mtDNA alone suffices).
6. **ANNOTATION (parallel branch)** — `SCANPY_QC` (per sample) → `results/qc/<sample>/<sample>_qc.h5ad`
   + `_qc_metrics.csv`; then `REFERENCE_MAPPING` (per sample, ingest onto the atlas) →
   `results/reference_mapping/<sample>/<sample>_celltypes.csv` + `_mapped.h5ad`. Gated by
   `run_qc` / `run_reference_mapping` (reference mapping implies QC).
7. **Provenance**: `results/pipeline_info/` — nf-prov BCO/legacy manifests + co2footprint
   trace/report/summary + Nextflow execution report/timeline.

## Containers / environments

- **Numbat**: `docker://pkharchenkolab/numbat-rbase:latest` (runs `/numbat/inst/bin/pileup_and_phase.R`).
- **Souporcell**: `souporcell_release.sif` (runs `souporcell_pipeline.py`). Joint runs require CB
  retagging (`<sample>__<barcode>`), BAM merge + sort + index, combined barcode list; invoked with
  `--no_umi true --skip_remap True --ignore True`.
- **CopyKAT**: no public biocontainer — built locally from `containers/copykat/Dockerfile`
  (R + conda deps per DDE_32 `setup_copykat.sh`, then `remotes::install_github("navinlabcode/copykat")`).
- **scanpy/Scrublet** (annotation branch): `containers/copykat`-style local build from
  `containers/scanpy/Dockerfile` (python3.11 + scanpy + scrublet + leidenalg). `--scanpy_container`.
- Execution: SLURM executor, account `biol-stem-2022`, apptainer enabled (3 h pull timeout),
  `errorStrategy='retry'` / `maxRetries=2`, lenient caching, `afterScript='sleep 30'` to tolerate
  shared-filesystem latency. These are ported from DDE_32 `nextflow.config`.

## Reference data

- Cell Ranger: `references/refdata-gex-GRCh38-2024-A`.
- Numbat 1000G (hg38): SNP VCF `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf`, phasing panel
  `1000G_hg38/`, genetic map `genetic_map_hg38_withX.txt.gz`. One-time download via ported
  `setup_numbat.sh`.
- Souporcell genome: `references/refdata-gex-GRCh38-2024-A/fasta/genome.fa`.
- (Legacy, HoneyBADGER-era) gnomAD common-SNP prep in `download_gnomad_snps.sh` — not in the active
  caller set.

## Source material in DDE_32 / DDE_24 (port from these)

- Numbat: `DDE_32/scripts/run_joint_numbat_pileup.sh`, `run_numbat_pileup.sh`,
  `joint_numbat_analysis.R`, `run_joint_numbat_analysis.sh`, `setup_numbat.sh`.
- CopyKAT: `DDE_32/scripts/copyKAT_profiling.R`, `run_copykat_analysis.sh`, `setup_copykat.sh`.
- Souporcell (in `DDE_24_CITE_seq_reference_map/clean_run_for_grant/scripts/`):
  `08g_souporcell_demux_paired_k2to20_noNK.sh` (paired joint run + retag/merge logic),
  `16b_subset_bam_noNK_batch.sh` (noNK BAM subsetting), `08_souporcell_demux.sh`.
- Cell Ranger driver: `DDE_32/scripts/run_scrnaseq_pipeline.sh`, `config/samplesheet.csv`,
  `config/feature_reference.csv`, `config/params.yaml`.
- Annotation branch (port from DDE_23 scanpy, NOT DDE_32 R): `DDE_23/scripts/stage1_qc.py`
  (scanpy + scrublet), `stage3_reference_mapping_individual.py` (`scanpy.tl.ingest`), `config.py`.
  Atlas h5ad: `DDE_32/references/paediatric_bm_reference/bone_marrow_atlas.h5ad` (default) or
  `DDE_23/zeng_reference_map/BoneMarrowMap_Annotated_Dataset_expandedFeatures.h5ad` (Zeng).

## Samplesheet

Extends the nf-core/scrnaseq format with grouping columns for joint calling:

```
sample,fastq_1,fastq_2,feature_type,expected_cells,patient,timepoint
```

`feature_type` is `gex` or `ab`; `timepoint` is `Dx` or `Rel`. The `patient` column drives
Numbat-joint and souporcell-paired grouping (replaces the hardcoded patient→sample maps in the
original scripts).

## Commands

The pipeline is run with Nextflow on SLURM. Until `main.nf` exists, see `scratchpad.md`. Intended:

```bash
# Main run (FASTQ → variant checkpoints)
nextflow run . -profile apptainer -params-file params.yaml -resume

# Quick DAG/wiring validation without heavy compute (`-stub` is a CLI flag, not a profile)
nextflow run . -profile test -stub -resume

# One-time reference preparation (Numbat 1000G panel etc.) — kept separate from the main run
bash bin/setup_numbat.sh
```

## HPC environment

This runs on a shared **SLURM HPC cluster** (account `biol-stem-2022`). **Do not run compute on the
login node** — submit it. With Nextflow this is largely automatic: the `slurm` executor dispatches
every process as its own SLURM job, so `nextflow run ...` only needs the lightweight orchestrator
process. Keep that orchestrator off the login node too: launch it inside a SLURM job (or `srun`/
interactive allocation), not directly on the login shell. Apply the same rule to any one-off helper
work (reference downloads, container pulls, `samtools`/`bcftools` munging) — wrap it in `sbatch`/
`srun` rather than running it interactively on the login node.

## Data workflow / storage (Viking longship)

Storage tiers on Viking (see https://vikingdocs.york.ac.uk/getting_started/storage_on_viking.html):

| Tier | Path | Notes |
|---|---|---|
| Home | `/users/hbp534` | 100 GB, code/envs only, never deleted. Not for job writes. |
| Scratch (Lustre) | `/mnt/scratch/users/hbp534` | 512 GB, high-performance, compute read/write. **Deleted after 90 days of inactivity.** This project lives here. |
| Longship | `/mnt/longship/users/hbp534` | 2 TB, medium-term warm storage, no deletion policy. **Read-write on login nodes, READ-ONLY on compute nodes.** |

**Recommended flow (push/pull from longship):** stage inputs to longship → jobs read inputs
(longship is readable from compute) and write `work/` + `results/` to scratch → push `results/`
back to longship for retention → back up to Filestore/Vault. Because scratch is purged after 90
days idle, **anything worth keeping must be pushed to longship.**

For this pipeline:
- Run Nextflow with `-work-dir` and `--outdir` on **scratch** (fast Lustre, compute-writable).
- Pull large inputs (FASTQ, Cell Ranger / Numbat / souporcell references) from **longship**; jobs
  may read them directly (read-only on compute is fine) or stage-copy to scratch for hot I/O.
- After a successful run, **push `results/` back to longship.** The existing
  `../snapshot_dde_to_longship.sh` does this: rsync to `/mnt/longship/users/hbp534/dde_snapshots_<stamp>/`,
  zstd-compressing `*.bam/*.cram/*.fastq/*.vcf`, tarring `work/`, and writing `SHA256SUMS.txt`.

**Login-node caveat:** the push to longship must run on a **login / data-transfer node** because
longship is read-only on compute nodes — this is the one sanctioned exception to "keep work off the
login node." Keep it to lightweight `rsync`/`zstd` data staging only; never run analysis there.

**Downloads + environment init go on the login node, not sbatch.** Compute nodes usually cannot
reach the internet, so run SRA `prefetch`/`fasterq-dump`, `apptainer pull`, reference downloads, and
environment setup (conda env creation, `setup_numbat.sh`) interactively on the login/data-transfer
node. Reserve SLURM jobs for actual compute (Cell Ranger, the Nextflow caller processes). This is
the sanctioned login-node exception alongside longship staging.

**Lossless compression — integrate wherever possible.** Genomic data is large and both scratch
(quota) and longship (quota) are finite, so prefer compressed-at-rest formats end-to-end:
- **Alignments**: keep/emit **CRAM** (reference-based, lossless) rather than BAM where a tool
  accepts it; the callers already fall back to `possorted_genome_bam.cram`. Pass the reference
  FASTA so CRAM is decodable.
- **Text/tabular outputs**: write **bgzip** `.gz` (+ `tabix` where applicable) — e.g. Numbat
  `*_allele_counts.tsv.gz`, VCFs as `.vcf.gz`. Avoid publishing uncompressed `.tsv`/`.vcf`.
- **Archival to longship**: **zstd** (`-T0 -8`) for `*.bam/*.cram/*.fastq/*.vcf` and a single
  `work.tar.zst`, as `../snapshot_dde_to_longship.sh` already does.
- In Nextflow, `publishDir` uses `mode: 'link'` (hardlink) here: results are real files that
  survive a `work/` cleanup with no duplication (work/ and results/ share scratch Lustre). Do NOT
  use `'symlink'` — deleting `work/` then destroys the published results. Use `'copy'` only when
  results must live on a different filesystem from `work/`. Compress in the process, not via a copy.
- **Never `rm -rf work` in the project root** while a run's results matter — even with hardlinks the
  `-resume` cache lives in `work/` + `.nextflow/`. Run stub/test validations with an isolated
  `-work-dir work_stub` and a throwaway `--outdir`.

## Conventions / gotchas

- Always pool Dx + Rel into one call for Numbat-joint and souporcell-paired; never compare
  per-sample clone IDs across timepoints.
- Numbat per-sample clone-calling fails at default thresholds (`max_entropy=0.5`, `min_LLR=5`) for
  most of this cohort; the joint run uses relaxed `0.8` / `3`.
- Souporcell paired runs were on noNK-filtered BAMs (~50% of cells, the NK-lineage normals, have no
  call) — expected and fine for malignant clonal tracing.
- Barcode files must be plain text for cellsnp-lite / mgatk-style tools; Cell Ranger ships them
  gzipped (`barcodes.tsv.gz`) so decompress before passing.
- Prefer BAM but fall back to CRAM (`possorted_genome_bam.cram`) — the original scripts all handle
  this; preserve it.

## Lab book

Lab book lives in `docs/lab_book/`: dated worklogs in `sessions/YYYY-MM-DD_slug.md`, distilled
findings in `analyses/NN_slug.md`, index in `docs/lab_book/README.md`. Full standard + templates:
`../_workspace_admin/LAB_BOOK_STANDARD.md`. `scratchpad.md` is the live build tracker; the lab book
records what happened (runs, results, decisions).

When asked to log work, write/update an entry to that standard:
- Use the matching template; keep frontmatter present and `status` accurate.
- No bare claims — every result is a number or an artifact path (job counts, walltime, file paths).
- Commands in fenced ```bash blocks with real flags; Slurm jobs record id · resources · efficiency.
- Tool versions + links; paths repo-relative or absolute; raw/irreplaceable data → Longship path.
- Connect entries with `[[wikilinks]]`; record `branch@commit` (+ container/env) when code ran.
- Session = append a worklog; Analysis = distil + link the sessions that produced it.
- Update `docs/lab_book/README.md` when adding an entry.

# Implementation plan — DDE_33 nf-core variant-calling pipeline

> Working scratchpad. Update task status as you go (`[ ]` → `[~]` in progress → `[x]` done).
> Background and agreed decisions live in `CLAUDE.md`; this file is the build sequence.

## Target layout

```
DDE_33_nextflow_variant_calling/
├── main.nf
├── nextflow.config              # ported DDE_32 config + manifest{} + tower{} hooks
├── nextflow_schema.json         # Seqera launch form
├── params.yaml                  # default params for a real run
├── conf/
│   ├── base.config              # resource labels (process selectors)
│   ├── modules.config           # per-process publishDir / ext.args
│   └── test.config              # tiny inputs + -stub for DAG validation
├── assets/
│   └── schema_input.json        # samplesheet column validation
├── workflows/
│   └── variantcalling.nf        # top-level workflow: wires subworkflows + caller toggles
├── subworkflows/local/
│   ├── input_check.nf           # parse + validate samplesheet, build per-sample/per-patient channels
│   ├── cellranger.nf
│   ├── numbat.nf
│   ├── copykat.nf
│   └── souporcell.nf
├── modules/local/
│   ├── cellranger_multi.nf
│   ├── numbat_pileup.nf         # pileup_and_phase.R (per-sample AND joint-per-patient)
│   ├── numbat_run.nf            # run_numbat()
│   ├── copykat.nf
│   ├── souporcell_prep.nf       # CB retag + merge + sort + index (per patient)
│   └── souporcell.nf            # souporcell_pipeline.py over K sweep
├── bin/                         # executable scripts called by modules
│   ├── pileup_and_phase wrapper
│   ├── run_numbat.R             # from joint_numbat_analysis.R
│   ├── copykat_profiling.R
│   ├── souporcell_retag.sh
│   └── setup_numbat.sh          # one-time reference prep (not part of main run)
├── containers/copykat/Dockerfile
├── README.md
└── CLAUDE.md
```

## Channel design

- `input_check` reads `samplesheet.csv` → emits:
  - `ch_sample` : `[ meta(sample,patient,timepoint), fastq_1, fastq_2, feature_type ]` (per sample, for Cell Ranger)
  - `ch_bam`    : `[ meta, bam, bai, barcodes ]` (per sample, post-Cell Ranger, for per-sample callers)
  - `ch_patient`: `[ meta(patient), [sample metas], [bams], [barcodes] ]` grouped by `patient` (for joint Numbat + paired souporcell)
- Caller subworkflows gated by `params.run_numbat` / `run_copykat` / `run_souporcell`.

## Build phases

### Phase 1 — scaffold + config (review checkpoint A) — DONE
- [x] `nextflow.config`: ported DDE_32 block (SLURM executor, account, apptainer, retry, lenient
      cache, afterScript, resourceLimits, nf-prov, co2footprint) + `manifest{}` + `tower{}` +
      `profiles { apptainer; singularity; docker; test; tower }` + `conf/*.config` includes.
- [x] `conf/base.config` resource labels (process_low/medium/high); `conf/modules.config` default
      symlink publishDir contract.
- [x] `assets/schema_input.json` + `nextflow_schema.json` (launch form, params grouped).
- [x] `subworkflows/local/input_check.nf` (nf-schema `samplesheetToList`) + `params.yaml` defaults.
- [x] `main.nf` includes `workflows/variantcalling.nf` (input parsing + per-sample/per-patient
      channels; no caller logic yet).
- [x] Validated: `nextflow config -profile test` resolves; `nextflow run . -profile test` builds all
      4 sample channels (gex+ab grouped) and correct per-patient grouping; exit 0.

> **Validation note:** `-stub` is a CLI flag, not a profile (use `-profile test -stub`). The
> `nf-co2footprint` plugin throws an NPE at shutdown on **process-free** runs (zero task records);
> harmless and disappears once Phase 2 adds processes. To validate a process-free scaffold cleanly,
> run with `-plugins nf-schema@2.7.2` to skip loading co2footprint/prov. `co2footprint{enabled=false}`
> is NOT honored by the plugin. Placeholder test inputs live in `assets/test/`.

### Phase 2 — Cell Ranger — DONE
- [x] `modules/local/cellranger_multi.nf` — `cellranger multi` (GEX + Antibody Capture); builds the
      multi-config CSV from samplesheet libraries (`expect-cells` from `expected_cells`); normalizes
      output to canonical `<sample>/outs/possorted_genome_bam.bam` + `.bai` +
      `filtered_feature_bc_matrix/`; full `stub:` block. Container via `params.cellranger_container`
      (default `nf-core/cellranger:8.0.0`; licensed, build/pull locally).
- [x] `subworkflows/local/cellranger.nf` — flattens fastqs for staging, derives `fastq_id` per
      library, emits `aln = [meta, bam, bai, filtered_matrix]` + `outs` + `versions`.
- [x] Wired into `workflows/variantcalling.nf`; added per-patient alignment grouping
      (`ch_patient_aln`, ordered Dx→Rel) for Numbat-joint / souporcell-paired.
- [x] publishDir → `results/cellranger/<sample>/outs/` (canonical; matches DDE_32 nf-core layout).
- [x] Validated `nextflow run . -profile test -stub` (exit 0): all 4 samples through CELLRANGER,
      patients grouped Patient_1=[2395 Dx, 3001 Rel], Patient_2=[2977 Dx, 0109 Rel].

> **Notes:** local/login test node has 4 CPUs, so the test profile caps
> `resourceLimits = [cpus:2, memory:6.GB, time:1.h]` (process_high requests 16). co2footprint NPE
> from Phase 1 is gone now that real tasks emit records — run the full plugin set normally.
>
> **Reporting config fixed (was leaking files to the project root):**
> - `nf-co2footprint 1.2.x` changed schema — use nested `report{file}` / `summary{file}` /
>   `trace{file}` blocks, NOT the old flat `reportFile`/`traceFile`/`summaryFile` keys (those are
>   silently ignored → plugin dumps defaults in cwd). This also affects DDE_32 (same stale keys).
> - Config-block paths (`prov`/`co2footprint`/`report`/`timeline`/`trace`/`dag`) resolve
>   `${params.outdir}` at parse time, BEFORE a profile's `outdir` override — so a profile that only
>   sets `params.outdir` still writes reports to the default `results/`. Fine for real runs (outdir
>   comes from `-params-file`/`--outdir`, applied before parsing). The `test` profile therefore
>   DISABLES prov/co2footprint/report/timeline/trace/dag (`prov.enabled=false`, per-file
>   `enabled=false`) — stub reports are meaningless and this keeps test runs leak-free.
> - Execution report/timeline/trace/dag now `enabled = true` in the main config (need it to render).
> - Added `.gitignore` (work/, results*/, .nextflow*, stray report defaults).

### Phase 3 — Numbat (primary axis) (review checkpoint B) — DONE
- [x] `bin/run_numbat.R` ported from `joint_numbat_analysis.R`; cohort comes from CLI args
      (samples/allele/matrix CSVs), thresholds from params (`max_entropy=0.8`, `min_LLR=3`),
      uses built-in `ref_hca`. Executable so Nextflow puts it on PATH via bin/.
- [x] `modules/local/numbat_pileup.nf` — joint multi-sample `pileup_and_phase.R` per patient;
      decompresses barcodes from each matrix dir; SNP VCF/panel/gmap default to in-container paths
      (params). Emits per-sample `*_allele_counts.tsv.gz`. Container `params.numbat_container`.
- [x] `modules/local/numbat_run.nf` — `run_numbat.R` → `<patient>/numbat_out/` (checkpoint
      `clone_post_1.tsv`). Both modules have `stub:` blocks.
- [x] `subworkflows/local/numbat.nf` — pileup → join allele counts with matrices → run; per-patient.
- [x] Wired into `workflows/variantcalling.nf` gated on `params.run_numbat`; publishDir
      `results/numbat_joint/<patient>/{<patient>_pileup, numbat_out}`. Added `numbat_container`
      param + schema; in-container ref paths as defaults.
- [x] **CELLRANGER_MULTI now emits uniquely sample-named handles** (`<id>.bam`/`.bai`/
      `<id>_filtered_feature_bc_matrix`) in `aln` so joint callers stage multiple samples without
      name collisions; publishDir restricted to the canonical `<id>/outs/` only.
- [x] Validated `-profile test -stub` (exit 0): 4 CELLRANGER → 2 NUMBAT_PILEUP (joint) → 2
      NUMBAT_RUN; pileup emits per-sample allele counts (Sample_2395_ + Sample_3001_), run emits
      clone_post_1.tsv.

- [ ] TODO later: `bin/setup_numbat.sh` (one-time reference prep) — deferred to Phase 5 docs.
- [ ] TODO: confirm numbat-rbase container actually ships /data + /Eagle ref paths on the cluster
      (the DDE_32 pileup script assumed it); if not, expose them as stage-mounted reference inputs.

### Phase 4 — CopyKAT + souporcell — DONE
- [x] `bin/copykat.R` (leaner than DDE_32 copyKAT_profiling.R — runs directly on the cellranger
      filtered matrix; the Seurat UMAP/celltype overlays were downstream analysis, dropped here) +
      `containers/copykat/Dockerfile` (micromamba; mirrors setup_copykat.sh, build & push yourself) +
      `modules/local/copykat.nf` + `subworkflows/local/copykat.nf` (per-sample; emits
      `<sample>_copykat_prediction.txt`).
- [x] `bin/souporcell_retag.sh` (CB/CR retag, from DDE_24 08g; CRAM-aware) +
      `modules/local/souporcell_prep.nf` (retag + merge + sort + index + combined barcode list;
      samtools container). **No NK filtering** — full BAMs.
- [x] `modules/local/souporcell.nf` (`souporcell_pipeline.py` per K, `--no_umi --skip_remap
      --ignore`; `samtools faidx` if no .fai) + `subworkflows/local/souporcell.nf` (joint per
      patient × K sweep from `params.souporcell_k`).
- [x] Wired both into `workflows/variantcalling.nf` gated on `run_copykat` / `run_souporcell`;
      container params + schema; publishDir `results/copykat/<sample>/` and
      `results/souporcell/<patient>/k<K>/` (SOUPORCELL_PREP merged BAM not published).
- [x] Validated full DAG `-profile test -stub` (exit 0): 4 CELLRANGER, 4 COPYKAT, 2 NUMBAT_PILEUP +
      2 NUMBAT_RUN, 2 SOUPORCELL_PREP + 4 SOUPORCELL (2 patients × K{2,3}). All checkpoints publish.

> Note: publishDir symlinks the souporcell `k<K>/` *directory*, so `find` without `-L` won't show
> `clusters.tsv` inside — it is present. CopyKAT/souporcell container images must be built/pulled
> before a live run (copykat: build Dockerfile; souporcell: cumulusprod/souporcell:2.5).

### Phase 5 — Seqera hooks + docs — DONE
- [x] `tower.yml` — surfaces execution report / timeline / co2footprint report in the Platform.
- [x] `nextflow.config` tower block: env-based auth (TOWER_ACCESS_TOKEN / TOWER_WORKSPACE_ID /
      TOWER_API_ENDPOINT). **Token is NEVER stored in the repo.**
- [x] `nextflow_schema.json` finalized (all params incl. container images grouped + described).
- [x] `README.md` — pipeline diagram, output table, samplesheet spec, quick start, Viking/storage,
      Seqera (with-tower + tw-agent inside a SLURM alloc; agent connection id
      db7575ab-6cbe-4d38-82fe-7a6d33a47668), caveats.
- [x] `bin/setup_numbat.sh` ported (container pull + optional 1000G ref download).
- [x] Final validation: `nextflow config -profile apptainer` OK; schemas valid JSON;
      `nextflow run . -profile test -stub` exit 0.

> SECURITY: the Tower access token was shared in chat in plaintext — it must be rotated in the
> Seqera Platform. Never commit it; Nextflow reads it from $TOWER_ACCESS_TOKEN.

### Phase 6 — Annotation branch (QC + reference mapping, Python/scanpy) — DONE
Scope add-on (user prefers Python; ported from DDE_23, not DDE_32 R). Runs **parallel** to the
callers off the raw cellranger matrices — does NOT gate caller inputs.
- [x] `bin/qc.py` (scanpy + scrublet, per sample, from DDE_23 stage1_qc.py) →
      `<sample>_qc.h5ad` + `<sample>_qc_metrics.csv` (all cells, pass/fail flag). QC params from
      `params.qc_*`.
- [x] `bin/reference_mapping.py` (`scanpy.tl.ingest`, from DDE_23 stage3) → `<sample>_celltypes.csv`
      + `<sample>_mapped.h5ad` (ref_cell_type, mapping_confidence, poorly_mapped). Atlas-agnostic
      with celltype-key auto-fallback.
- [x] `containers/scanpy/Dockerfile` (python3.11 + scanpy + scrublet + leidenalg).
- [x] `modules/local/scanpy_qc.nf` + `modules/local/reference_mapping.nf` (+stubs);
      `subworkflows/local/annotation.nf` (QC → REFMAP; REFMAP consumes QC h5ad).
- [x] Wired into workflow gated `run_qc` / `run_reference_mapping` (refmap implies QC); params +
      schema (annotation_options group) + publishDir `results/qc/<sample>/`,
      `results/reference_mapping/<sample>/`.
- [x] Atlas: default `../DDE_32.../paediatric_bm_reference/bone_marrow_atlas.h5ad` (exists, h5ad);
      Zeng option `../DDE_23.../zeng_reference_map/...h5ad`; configurable `--refmap_atlas` +
      `--refmap_celltype_key`.
- [x] Validated `-profile test -stub` (exit 0): full DAG now 4 CELLRANGER + 4 SCANPY_QC + 4
      REFERENCE_MAPPING + 4 COPYKAT + 2 NUMBAT_PILEUP/RUN + 2 SOUPORCELL_PREP + 4 SOUPORCELL.

> Build the scanpy-scrublet image before a live run. REFERENCE_MAPPING uses `scanpy.tl.ingest`
> (UMAP-centroid distance → confidence). Confirm the cell-type obs key in bone_marrow_atlas.h5ad
> (default 'cell_type', script auto-falls back to celltype/CellType/annotation).

#### TODO (later) — more rigorous reference mapping
`scanpy.tl.ingest` is the v1 method (quick label transfer + UMAP-centroid confidence). Add better
options as alternatives behind a `--refmap_method` switch:
- **Symphony** — DDE_23's reference is literally `BoneMarrowMap_SymphonyReference.rds`
  (`zeng_reference_map/`), so a Symphony query-mapping path is the natural fit and gives proper
  reference-anchored mapping + per-cell confidence. (R-based; would need a symphony container or
  symphonypy.)
- **scANVI** (scvi-tools) — probabilistic label transfer with calibrated uncertainty; integrates
  well with scanpy/anndata and the Python stack. Heavier (GPU helpful).
Keep `ingest` as the default/fallback; wire method choice + per-method container as params.

### Phase 7 — First real run: Caron 2020 healthy PBMMC controls — IN PROGRESS
Test the pipeline on 3 healthy paediatric controls (Caron et al. 2020, GSE132509 / SRP201012),
10x 3' v2, **GEX-only**. Run everything (Numbat/souporcell uninformative on healthy but exercise the
DAG; CopyKAT should call ~all diploid).

Accessions: PBMMC_1=SRR9264351, PBMMC_2=SRR9264353, PBMMC_3=SRR9264354 (GSM3872442/3/4).

- [x] **`conf/viking.config` + `-profile viking`** — no container builds: CellRanger module 9.0.0,
      SAMtools module (souporcell prep), conda env `snv` (CopyKAT), conda env `aml_scrna`
      (scanpy+scrublet QC/refmap), `numbat.sif` (DDE_32) + `souporcell_release.sif` (DDE_24) via
      apptainer. Confirmed `aml_scrna` has scanpy 1.11.5; `snv` has copykat; both R envs 4.5.3.
- [x] **GEX-only support** — `fb_reference` optional in workflow (`[]` when null) + module emits the
      `[feature]` section only when given. Validated stub with `fb_reference: null`.
- [x] `bin/fetch_sra_10x.sh` — prefetch + fasterq-dump (--split-files --include-technical), classify
      reads by length (I1<=12, R1<=30 [v2 26bp], else R2), rename to cellranger convention.
- [x] `assets/controls_samplesheet.csv` (3 samples, gex, each its own patient, timepoint=Dx) +
      `params-controls.yaml` (absolute ref paths; atlas = DDE_32 paediatric BM h5ad; souporcell_k=2,3).
- [x] Submitted SRA download SLURM job **34329477** (logs/sra_pbmmc_*.log; writes data/controls/).
- [ ] **NEXT (after download):** `nextflow run . -profile viking -params-file params-controls.yaml -resume`
      inside a SLURM alloc. First confirm fastqs landed + read-length classification was right
      (v2: R1=26bp, R2=98bp).

> RISK: SRA prefetch needs outbound internet on the compute node — if Viking compute nodes block it,
> run prefetch on a login/transfer node and fasterq-dump on compute. Watch job 34329477.
> Numbat/souporcell are single-sample per control (each control = its own "patient"); joint logic
> degrades to per-sample, which is fine mechanically.

## Remaining before a live run (not blockers for stub)
- Build/push the CopyKAT image (containers/copykat/Dockerfile) and pull cellranger + numbat +
  souporcell images into apptainer.
- Verify the numbat-rbase container ships /data + /Eagle reference paths on the cluster; if not,
  expose them as stage-mounted reference inputs.
- Provide real references in params.yaml (cellranger_index, fb_reference, souporcell_fasta).
- Confirm Cell Ranger multi-config generation matches the live CITE-seq setup (fastq_id derivation).
- First real run: start with one patient + `--run_numbat` only, then enable the rest.

## Open questions to resolve while building
- Cell Ranger module: confirm `cellranger multi` config CSV generation from the samplesheet
  (GEX + Antibody Capture + fb_reference) matches the CITE-seq setup.
- Numbat per-sample mode: include as an optional path or joint-only? (default joint).
- Souporcell K: single K or sweep range as a param? (DDE_24 swept k2–20; default to a list param).

## Done-but-verify
- mgatk excluded by decision — do not add.
- NK filtering excluded for now — souporcell runs on full BAMs; do not port the DDE_24 noNK step.
- Lossless compression integrated throughout: emit CRAM over BAM where tools allow (pass ref FASTA),
  bgzip+tabix text/VCF outputs, zstd for longship archival. No uncompressed genomic files published.
- Patient→sample mapping comes from the samplesheet `patient`/`timepoint` columns, never hardcoded.
```

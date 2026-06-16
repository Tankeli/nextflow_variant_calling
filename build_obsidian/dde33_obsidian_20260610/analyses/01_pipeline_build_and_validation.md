---
title: Pipeline build & stub validation (Phases 1–6)
project: DDE_33
type: analysis
status: active
updated: 2026-06-03
tags: [nextflow, nf-core, cellranger, numbat, copykat, souporcell, scanpy, pipeline]
related: ["[[2026-06-04_first-real-run-caron-controls]]", "[[DDE_32]]", "[[DDE_24]]"]
---

# 01 — Pipeline build & stub validation (Phases 1–6)

> Reconstructed 2026-06-06 from `scratchpad.md` (the build plan, all phases marked DONE) and
> `CLAUDE.md`. This is the standardised record of *what was built and how it was validated*;
> `scratchpad.md` stays the live task tracker.

## Question / goal

Replace DDE_32's fragmented SLURM/conda/apptainer wrapper scripts (and DDE_24's souporcell) with
a single resumable, nf-core-style Nextflow pipeline: CITE-seq FASTQ → standardised per-sample /
per-patient variant checkpoints + a parallel cell-annotation branch. Validate the full DAG with
`-stub` (no heavy compute) before any live run.

## Data & provenance

| Input | Source | Notes |
|---|---|---|
| Caller/annotation source scripts | [[DDE_32]] (`scripts/`), [[DDE_24]] (`clean_run_for_grant/scripts/08g_*`) | ported, not rewritten |
| Test inputs | `assets/test/` (placeholder FASTQ) | `-profile test -stub` only |
| Cell Ranger ref | `references/refdata-gex-GRCh38-2024-A` | via `/references` symlink |
| Numbat 1000G hg38 | SNP VCF + phasing panel + genetic map | one-time `bin/setup_numbat.sh` |

## Method

nf-core layout (`modules/`, `subworkflows/`, `conf/`, `nextflow_schema.json`, samplesheet
validation). Samplesheet extends nf-core/scrnaseq with `patient` + `timepoint` columns driving
joint grouping. Callers toggled by `--run_{numbat,copykat,souporcell}`;
annotation by `--run_{qc,reference_mapping}`. Build proceeded in reviewed phases, each validated
with `nextflow run . -profile test -stub`.

```bash
# wiring/DAG validation without heavy compute (-stub is a CLI flag, not a profile)
nextflow run . -profile test -stub
```

## Results

### Phases (all DONE, stub-validated exit 0)
| Phase | Delivered | Stub check |
|---|---|---|
| 1 — scaffold + config | `nextflow.config` (SLURM/apptainer/retry/nf-prov/co2footprint), `conf/{base,modules,test}`, schemas, `input_check` | 4 sample channels built, gex+ab grouped, per-patient grouping correct |
| 2 — Cell Ranger | `cellranger_multi.nf` (`cellranger multi`, GEX+AB, multi-config from samplesheet) | 4 samples → CELLRANGER; Patient_1=[2395,3001], Patient_2=[2977,0109] |
| 3 — Numbat (primary axis) | `numbat_pileup.nf` (joint per patient) + `numbat_run.nf` (`max_entropy=0.8`, `min_LLR=3`) | 4 CELLRANGER → 2 NUMBAT_PILEUP → 2 NUMBAT_RUN; per-sample allele counts |
| 4 — CopyKAT + souporcell | `copykat.nf` (per sample), `souporcell_prep.nf` (CB retag+merge, **no NK filter**) + `souporcell.nf` (K sweep) | 4 COPYKAT, 2 SOUPORCELL_PREP + 4 SOUPORCELL (2 patients × K{2,3}) |
| 5 — Seqera hooks + docs | `tower.yml`, env-based `tower{}` auth, finalized schema, README, `setup_numbat.sh` | `nextflow config -profile apptainer` OK |
| 6 — annotation branch | `qc.py` (scanpy+scrublet), `reference_mapping.py` (`scanpy.tl.ingest`), containers, gated `run_qc`/`run_reference_mapping` | full DAG: +4 SCANPY_QC +4 REFERENCE_MAPPING |

**Full stub DAG (exit 0):** 4 CELLRANGER + 4 SCANPY_QC + 4 REFERENCE_MAPPING + 4 COPYKAT +
2 NUMBAT_PILEUP + 2 NUMBAT_RUN + 2 SOUPORCELL_PREP + 4 SOUPORCELL.

### Key build decisions / fixes
- **mgatk excluded** by decision (also had silent Sample_3001 failures in DDE_32).
- **Souporcell on full BAMs** (no NK filtering ported from DDE_24) by decision.
- **Numbat joint-per-patient only** for stable cross-timepoint clone IDs; relaxed thresholds.
- **nf-co2footprint 1.2.x schema change** — nested `report{file}`/`summary{file}`/`trace{file}`
  blocks, not the old flat keys (silently ignored → files leak to cwd). **Same stale keys affect
  DDE_32.** Test profile disables prov/co2footprint/report to stay leak-free.
- **publishDir `mode:'link'` (hardlink)**, not symlink — results survive a `work/` cleanup
  without duplication. Never `rm -rf work` while results matter (resume cache lives there).

## Interpretation

The pipeline is feature-complete and wiring-correct across the full caller + annotation DAG in
stub. It is the production replacement for DDE_32's scripts; DDE_32 now consumes its checkpoints.
Remaining risk is entirely at the **live-run** boundary (real containers, real references, real
compute) — exercised next in [[2026-06-04_first-real-run-caron-controls]].

## Limitations / caveats

- **Stub proves wiring, not biology** — no real compute ran in Phases 1–6.
- **Cell Ranger + CopyKAT containers are not on public registries** (licensing / none published) —
  must be built/pulled before a live run.
- Numbat-rbase container `/data` + `/Eagle` reference paths **assumed**, not yet confirmed on the
  cluster (TODO: expose as stage-mounted inputs if missing).
- `scanpy.tl.ingest` is the v1 reference-mapping method (UMAP-centroid confidence); Symphony /
  scANVI deferred behind a future `--refmap_method` switch.
- **SECURITY:** a Tower access token was pasted in plaintext during the build — must be rotated;
  Nextflow reads it from `$TOWER_ACCESS_TOKEN`, never committed.

## Links

- Feeds: [[2026-06-04_first-real-run-caron-controls]]
- Depends on: [[DDE_32]] (caller scripts), [[DDE_24]] (souporcell), [[DDE_23]] (scanpy annotation)
- Code: DDE_33 `main`@b771178 · build tracker `scratchpad.md` · architecture `CLAUDE.md`

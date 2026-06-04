# DDE_33 — single-cell variant-calling pipeline

An nf-core-style Nextflow pipeline that takes CITE-seq (GEX + Antibody Capture) FASTQ through
Cell Ranger and three variant callers, producing standardized, resumable per-sample / per-patient
variant checkpoints for paired diagnosis–relapse paediatric AML. It formalizes the fragmented
analysis scripts from `../DDE_32_paediatric_snv_analysis` (and souporcell from
`../DDE_24_CITE_seq_reference_map`). See `CLAUDE.md` for architecture and `scratchpad.md` for build
status.

## Pipeline

```
FASTQ ── CELLRANGER (multi) ──┬─ NUMBAT     pileup_and_phase → run_numbat   (joint per patient) → CNV clones
                              ├─ COPYKAT                                    (per sample)         → aneuploid/diploid
                              ├─ SOUPORCELL prep (retag+merge) → K sweep    (joint per patient)  → SNV clusters
                              └─ ANNOTATION  scanpy QC → reference mapping  (per sample)         → cell-type labels
```

The **annotation branch** (Python/scanpy, ported from DDE_23) runs parallel to the callers off the
raw cellranger matrices — it does not gate caller inputs.

Outputs (`results/`):

| Stage | Path | Key file |
|---|---|---|
| Cell Ranger | `cellranger/<sample>/outs/` | `possorted_genome_bam.bam`, `filtered_feature_bc_matrix/` |
| Numbat | `numbat_joint/<patient>/numbat_out/` | `clone_post_1.tsv`, `segs_consensus_*.tsv` |
| Numbat pileup | `numbat_joint/<patient>/<patient>_pileup/` | `<sample>_allele_counts.tsv.gz` |
| CopyKAT | `copykat/<sample>/` | `<sample>_copykat_prediction.txt` |
| Souporcell | `souporcell/<patient>/k<K>/` | `clusters.tsv` |
| QC (scanpy) | `qc/<sample>/` | `<sample>_qc.h5ad`, `<sample>_qc_metrics.csv` |
| Reference mapping | `reference_mapping/<sample>/` | `<sample>_celltypes.csv`, `<sample>_mapped.h5ad` |
| Provenance | `pipeline_info/` | nf-prov BCO, co2footprint, execution report/timeline/trace/dag |

## Samplesheet

CSV with one row per library (a CITE-seq sample has a `gex` and an `ab` row):

```csv
sample,fastq_1,fastq_2,feature_type,expected_cells,patient,timepoint
Sample_2395,/path/GEX-2395_S1_L001_R1_001.fastq.gz,/path/GEX-2395_S1_L001_R2_001.fastq.gz,gex,5000,Patient_1,Dx
Sample_2395,/path/CSP-2395_S1_L001_R1_001.fastq.gz,/path/CSP-2395_S1_L001_R2_001.fastq.gz,ab,5000,Patient_1,Dx
Sample_3001,...,gex,5000,Patient_1,Rel
...
```

`feature_type` ∈ {`gex`, `ab`}; `timepoint` ∈ {`Dx`, `Rel`}. `patient` groups samples for the joint
Numbat and paired souporcell runs (no hardcoded cohort).

## Quick start

```bash
# 1. One-time: prepare references / containers (inside a SLURM job, not the login node)
sbatch --account=biol-stem-2022 --mem=16G --time=06:00:00 --wrap "bash bin/setup_numbat.sh"
#    build & push the CopyKAT image (no biocontainer exists):
#      docker build -t <registry>/copykat:1.1.0 containers/copykat && push   (or apptainer build)
#    pull souporcell: apptainer pull docker://cumulusprod/souporcell:2.5

# 2. Validate wiring without heavy compute (-stub is a CLI flag, not a profile)
nextflow run . -profile test -stub

# 3. Real run (edit params.yaml first). Keep outdir + work-dir on /mnt/scratch.
nextflow run . -profile apptainer -params-file params.yaml -resume
```

Toggle callers with `--run_numbat`, `--run_copykat`, `--run_souporcell`. Numbat thresholds
(`--numbat_max_entropy`, `--numbat_min_llr`) and the souporcell K sweep (`--souporcell_k 2,3,5,12`)
are parameters. Container images are set via `--{cellranger,numbat,copykat,souporcell,samtools}_container`.

## Viking HPC / storage

- SLURM cluster (account `biol-stem-2022`); the `slurm` executor submits each process as its own job.
  **Keep the Nextflow head process off the login node** — launch it inside a SLURM allocation.
- Run with `--outdir` and `-work-dir` on **/mnt/scratch** (fast Lustre, compute-writable; purged
  after 90 days idle). Pull large inputs/references from **longship** (`/mnt/longship/...`,
  read-only on compute, read-write on login). After a run, push `results/` back to longship for
  retention — `../snapshot_dde_to_longship.sh` (rsync + zstd). Longship is the only sanctioned
  login-node write. See `CLAUDE.md` for the full storage policy and the lossless-compression rules.

## Seqera Platform (Tower)

The pipeline is launchable from the Platform (`nextflow_schema.json` renders the launch form) and
surfaces reports via `tower.yml`. Auth is **environment-based — never commit your token**:

```bash
export TOWER_ACCESS_TOKEN=...        # from Seqera Platform → User → Access tokens (rotate if leaked)
export TOWER_WORKSPACE_ID=...        # your workspace ID
# optional, defaults to https://api.cloud.seqera.io
# export TOWER_API_ENDPOINT=...

# A) Monitor a CLI run:
nextflow run . -profile apptainer -params-file params.yaml -with-tower -resume

# B) Tower Agent — let the Platform launch runs onto Viking. Run the agent inside a SLURM
#    allocation (so the Nextflow head lands on a compute node, not the login node):
srun --account=biol-stem-2022 --mem=8G --time=24:00:00 --pty bash
export TOWER_ACCESS_TOKEN=...
/mnt/scratch/users/hbp534/tw-agent db7575ab-6cbe-4d38-82fe-7a6d33a47668 --work-dir=/mnt/scratch/users/hbp534/DDE_33_nextflow_variant_calling/work
```

`db7575ab-6cbe-4d38-82fe-7a6d33a47668` is this environment's agent connection ID; set the matching
"Tower Agent" compute environment + workspace in the Platform.

Handy SLURM alias:

```bash
alias sq='squeue --format="%.18i %.9P %.90j %.8u %.8T %.10M %.9l %.6D %.6C %R" --me'
```

## Notes / caveats

- **Cell Ranger and CopyKAT containers are not on public registries** (licensing / none published) —
  build/pull them yourself before a live run.
- Numbat per-sample clone calling fails at default thresholds for this cohort; the joint run uses
  relaxed `max_entropy=0.8`, `min_LLR=3`.
- Souporcell runs on full BAMs (no NK filtering in this pipeline).
- mgatk (mtDNA) is intentionally excluded.

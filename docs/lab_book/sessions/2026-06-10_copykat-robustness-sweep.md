---
date: 2026-06-10
project: DDE_33
type: session
status: done
tags: [copykat, robustness, seed-stability, controls, sweep, viking-profile, drivers, crossref]
related: ["[[02_copykat_robustness]]", "[[2026-06-08_controls-full-cohort]]", "[[2026-06-04_first-real-run-caron-controls]]"]
---

# 2026-06-10 — CopyKAT robustness sweep + downstream analysis (9 healthy controls)

> **Objective:** run the CopyKAT robustness track (built per plan
> `.claude/plans/ticklish-herding-wilkinson.md`) end-to-end on the 9-control cohort and produce the
> stability / driver / cross-reference / cell-type results that answer "how trustworthy is a CopyKAT
> aneuploid call". Distilled findings → [[02_copykat_robustness]].

## Context

The track is a **hybrid**: a gated Nextflow parameter×seed sweep (`COPYKAT_SWEEP`) + four standalone
Python analyses (`copykat_{stability,drivers,crossref,celltype_matrix}.py`). Implemented and
stub-validated earlier; this session is the first live run on real data — the 9 GEX-only controls
from [[2026-06-08_controls-full-cohort]] (`results_controls/`).

## Work done

### 1. Sweep orchestrator — 180 CopyKAT runs
`jobs/run_controls_robustness.sh` (new) on the `viking` profile: `--run_copykat_robustness` with the
other callers disabled so `-resume` only adds the new tasks; combined 9-sample sheet so all controls
run in one go. Defaults = `KS.cut {0.05,0.1,0.15,0.2} × seeds {1–5}` = 20 combos/sample × 9 = **180**.

```bash
sbatch jobs/run_controls_robustness.sh assets/controls_samplesheet_all9.csv
# nextflow run . -profile viking -params-file params-controls.yaml \
#   --input assets/controls_samplesheet_all9.csv --run_copykat_robustness \
#   --run_numbat false --run_souporcell false --run_clonetracer false \
#   --run_qc false --run_reference_mapping false --run_integration false \
#   -plugins nf-schema@2.7.2,nf-prov@1.7.0 -work-dir work -resume
```
Result: 180/180 combos published to `results_controls/copykat_robustness/<sample>/sweep/<combo>/`.

### 2. Downstream analysis
`jobs/run_copykat_robustness.sh` (`aml_scrna` conda) over the sweep + production CopyKAT + the atlas.
Outputs → `results_controls/copykat_robustness/_analysis/` (per-sample CSVs + PNGs; atlas anchor/marker
gene sets cached once to JSON). Ran twice: first pass (pre-sweep) gave drivers/crossref/celltype but
**skipped stability** (no sweep dir yet); re-run after the sweep added stability/boundary/ARI.

### Slurm accounting
- `34565386` dde33_ck_sweep · 2c/8G orchestrator · **9h37m** · COMPLETED (180 `COPYKAT_SWEEP` tasks,
  each process_medium 8c/64G/12h; CellRanger reused from `work/` cache).
- `34565164` dde33_ck_robust · FAILED 9s, **exit 127** — `copykat_*.py: command not found`.
- `34565441` dde33_ck_robust · 4c/32G · **38m** · COMPLETED — drivers/crossref/celltype (pre-sweep).
- `34605842` dde33_ck_robust · 4c/32G · **4m46s** · COMPLETED — stability re-run (crossref hit the
  atlas cache, hence fast).
Logs: `logs/copykat_sweep_orchestrator_34565386.log`, `logs/copykat_robustness_*.log`.

## Findings (summary — full distillation in [[02_copykat_robustness]])

- **Not reproducible:** 18–98% of cells flip class across the 20 combos; seed-only switch-rate up to
  **0.83** (PBM_1); mean ARI 0.21–0.84. Boundary curves non-monotonic in `KS.cut`.
- **Accidental replicate:** the sweep's `-resume` re-drew the default production calls; PBMMC_2 swung
  **81% → 22%** aneuploid on identical params (see caveats).
- **Over-call:** 22–82% "aneuploid" on true-normal cells (truth ≈ 0).
- **Lineage confound:** drivers enriched for cell-type markers in 8/9 (p down to 6.1e-28) and anchor
  genes; PBMMC_2 overlap = haemoglobin/erythroid/ribosomal (`HBA1/2, HBB, GYPA/B, KLF1, RPS*`).
- **Not LSC-like:** drivers vs pLSC6/LSC17 overlap = 0 → safe to keep CopyKAT as a gate, not a caller.
- Summary figure: `docs/lab_book/assets/02_copykat_robustness/robustness_summary.png`.

## Decisions & rationale

- **Other callers off + `-work-dir work -resume`** for the sweep: reuse the cached CellRanger (re-running
  it for 9 samples would cost hours), add only `COPYKAT_SWEEP`. Rejected a fresh work-dir (would re-run
  CellRanger) and a standalone (non-Nextflow) sweep (loses resume/provenance).
- **Orchestrator walltime cut 48h→40h**: SLURM wouldn't start a 48h job because it overran the
  maintenance reservation **YOR796** (2026-06-11 09:00 → 06-15). 40h clears it; anything unfinished
  resumes after maintenance via `-resume`. The sweep finished in 9h37m, so this was moot in the end.
- **PATH fix**: `bin/` is only on PATH inside Nextflow tasks, so the standalone driver now does
  `export PATH=$PROJECT/bin:$PATH` (the exit-127 cause).

## Caveats / things that bit

- **The sweep overwrote the production CopyKAT calls.** With the combined 9-sample samplesheet the
  production `COPYKAT` task missed `-resume` cache and re-ran for all 9, rewriting
  `results_controls/copykat/<s>/*_copykat_prediction.txt` (mtimes 2026-06-09 19:56–23:56) with fresh
  stochastic draws. The 06-08 cohort numbers (21–72%, [[2026-06-08_controls-full-cohort]]) no longer
  match the files; current numbers are the re-draws. Not data loss (regenerable, and CopyKAT is
  non-deterministic so neither draw is "correct") — but to avoid silent re-draws on patient runs,
  fix `set.seed()` in `bin/copykat.R` or pin the COPYKAT cache.

## Next steps

- [ ] Run the **`norm.cell.names` arm** (`copykat_robustness_use_norm_ref=[true,false]`) to test whether
      a supplied diploid baseline removes the over-call.
- [ ] Add a fixed `set.seed()` to the production `bin/copykat.R` so default calls are reproducible.
- [ ] Apply the same track to the **patient cohort** once their CellRanger completes
      ([[2026-06-09_patient-cohort-runs]]) — sample-agnostic; point the driver at `results_patients*`.
- [ ] Commit the robustness track (code only) — currently uncommitted working tree.

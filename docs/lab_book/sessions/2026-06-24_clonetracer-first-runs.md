---
date: 2026-06-24
project: DDE_33
type: session
status: done
tags: [clonetracer, clonal-inference, pyro, mtdna, cellsnp, figures, patients, debugging, resume]
related: ["[[2026-06-09_patient-cohort-runs]]", "[[2026-06-22_numbat-specificity]]", "[[clonetracer-mutation-cap]]", "[[clonetracer-env-slow]]"]
---

# 2026-06-24 — CloneTracer: first real runs + downstream figures (both patients)

> **Objective:** get CloneTracer (veltenlab Bayesian clonal inference) actually producing output and
> generate the downstream-analysis figures from its `clonal_inference` repo. CloneTracer had never
> completed a model run in DDE_33 — only the per-patient input JSON had ever been built.

## Context

The CloneTracer branch was code-complete and wired (`CLONETRACER_WF`: `MTDNA_PILEUP` →
`CLONETRACER_BUILD` → `CLONETRACER`) but disabled on every prior run (`run_clonetracer:false`): the
JSONs existed for the controls (`results_controls/clonetracer/HD_BM_{3,4}/*.json`) but **no
`*_out.pickle`** anywhere. Memory [[clonetracer-env-slow]] blamed a slow conda env; that turned out to
be a misdiagnosis.

## Work done

### 1. Diagnosed why the model never finished

Ran the pyro model on the existing HD_BM_3 JSON. Symptom on every attempt (100-mut GPU, 8-mut CPU):
**timeout with ~2–3 s CPU used** — i.e. blocked, not computing. Ruled out, by controlled tests:

- **Env** — `clonetracer_gpu` (pyro 1.8.4 / torch 1.13.1+cu117) imports cleanly with
  `PYTHONNOUSERSITE=1` (numpy 2.x was a `~/.local` leak) and runs fine **on CPU**. Not the blocker.
- **GPU** — every a40 was allocated (3/3 all nodes, ~16 h queue); moved to the `nodes` CPU partition.
- **OpenMP threads / stdout redirect** — tested and excluded as the *cause* (a no-redirect copy
  `bin/run_clonetracer_noredir.py` behaved identically).

**Real cause:** `infer_hierarchy` (helper_functions.py) adds mutations one at a time and re-fits SVI
to **every** candidate tree; the candidate count explodes super-exponentially. Measured live (`bin/ct_diag.py`,
num_iter=60): 2 muts→4 trees, 3→12, 4→**57** (7.8 min), 5→**178** (~30+ min). CloneTracer is built for a
handful of curated drivers (~3–6); the build step emitted up to 100 (`clonetracer_max_snvs=50` +
`clonetracer_mtdna_max_sites=50`). Also self-inflicted: `-t<100` makes `init=num_iter-100` negative and
crashes `print_elbo`. → [[clonetracer-mutation-cap]].

First successful pickle: 3-mut subset of HD_BM_3 (`-t` default 300), **6 min wall / 18m54s CPU**
(job 35140461) → `*_out.pickle` + `_tree.pickle` + assignments.

### 2. Downstream figures ported to Python

`bin/clonetracer_figures.py` (runs in `aml_scrna`) — port of the R vignette `funct_clonal_analysis.R`:
clonal-hierarchy **trees** (networkx), **ELBO**/model-evidence, single-cell **VAF + clone-posterior
heatmap**. UMAP overlay stays `bin/plot_clonetracer_umap.py`. Validated on HD_BM_3 (all 4 figure types).

### 3. Pipeline hardening

- **Total-mutation cap** (priority CNV>SNV>mtDNA) in `bin/clonetracer_build_json.py` via
  `--max-total-muts` + new param `clonetracer_max_total_muts=6`; per-source caps 50→6 (config + schema).
- `conf/viking.config` `CLONETRACER`: both paths use the pinned `clonetracer_gpu` env (CPU default),
  with `OMP/MKL/OPENBLAS_NUM_THREADS=${task.cpus}`.
- **Wired trees/ELBO/heatmap into the pipeline** — threaded `CLONETRACER_WF.out.trees` through
  `VISUALIZATION` into an extended `PLOT_CLONETRACER`. Validated with `-profile test -stub` (both
  patients produced all 5 figures, exit 0).

### 4. Ran both patients

Full-pipeline `-resume` was **unusable**: Cell Ranger task hash changed since the cached June run
(workdir `22/3bde50…`→`e9/89ce6f…`, revision `f8057f0cd1`→`3765530376`) → would re-run
cellranger/numbat/souporcell from scratch (days). Caught + cancelled twice. Pivoted to a **standalone
run off existing `results_patients/` outputs** (mirrors `CLONETRACER_WF`), capped at **4** mutations
(Patient_2 has CNVs → 500 iters):

```bash
# Stage 1 — mtDNA pileup (cellsnp-lite via numbat.sif), array over the 4 patient BAMs
sbatch jobs/patients_mtdna_pileup.sh                       # job 35182324[0-3]
# Stage 2 — build JSON -> model -> figures, both patients (dependency on Stage 1)
sbatch --dependency=afterok:35182324 jobs/patients_clonetracer_run.sh   # job 35182415
```

- **mtDNA pileup** 35182324 — 4/4 COMPLETED, 00:32–01:15 each; `results_patients/clonetracer/mtdna/<sample>_mtdna/`.
- **Build+model+figures** 35182415 — both patients, exit 0 (~17:32 / 18:39).

## Results

`results_patients/clonetracer/<patient>/` — `<patient>.json`, `_out.pickle`, `_tree.pickle`,
`_clone_assignments.csv`, `figures/` (trees, elbo, heatmap, clonetracer_umap, composition · png+pdf).

- **Patient_2** (Sample_2977 Dx + Sample_0109 Rel) — CNV(chr21 del) + 3 souporcell SNVs (chr9/chr4/chr8),
  **5887 cells**. **Branching** hierarchy: Healthy → chr4+chr8-SNV clone → {chr21-del subclone,
  chr9-SNV subclone} (tree 8, lowest ELBO). 4 clones; Dx→Rel dynamics, e.g. clone 1 **222→2470** cells,
  clone 3 671→1530.
- **Patient_1** (Sample_2395 Dx + Sample_3001 Rel) — **mtDNA-only** (no `numbat_out`/souporcell exist
  for it, per [[2026-06-22_numbat-specificity]]): 4 mtDNA variants, **15058 cells**, 4 clones, Dx/Rel shifts.

## Notes / open items

- Code ran at `main@2521d03` **+ uncommitted changes** (clonetracer figures, caps, env, figure wiring).
  Standalone outputs are not in the Nextflow `-resume` cache; a future clean pipeline run would
  regenerate from scratch.
- Caron healthy-control CloneTracer figures (`run_controls_clonetracer.sh`) not regenerated this
  session — the controls are degenerate (no real clones); the patient runs are the deliverable.
- `clonetracer_max_total_muts` default is 6 but Patient_2 (CNV → 500 iters) was run at 4 for runtime;
  consider lowering the default for CNV-bearing samples.

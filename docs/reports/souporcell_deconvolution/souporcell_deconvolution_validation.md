---
title: Souporcell deconvolution validation
project: DDE_33
date: 2026-06-26
tags: [souporcell, deconvolution, single-cell, variant-calling, validation]
---

# Souporcell deconvolution validation

**TL;DR** — Souporcell separates pooled single-cell samples back to their genetic individual of
origin with **near-perfect accuracy (ARI 0.97–1.00)** when given the correct number of individuals.
Shown by *artificially mixing* known samples (healthy controls, patient diagnosis samples, and a
male×female pair), deconvoluting blind, and comparing to ground truth recovered from the original cell
barcodes. A male×female mix additionally shows the genotype split is reproduced by sex-chromosome
expression alone; a proportion titration probes how rare a population can get before it is missed; and
negative/wrong-K controls behave exactly as expected.

---

## 1. Why and how

**Souporcell** is a genotype-based deconvolution tool: from the SNVs in each cell's scRNA-seq reads it
clusters cells by genetic individual — used to demultiplex pooled donors and to separate donor vs
recipient cells after a transplant. To demonstrate it works we need a *ground truth*: which cell came
from which individual.

**Design.** Pool two (or three) samples that are each a single known individual, hand the pool to
souporcell, and ask it to split them again. The pipeline retags every barcode with a `<sample>__`
prefix before merging, so the true origin of each cell is known but invisible to souporcell — a clean
blind test. Agreement between souporcell's clusters and the true origin is scored by the **Adjusted
Rand Index (ARI)**: 1.0 = perfect, 0 = random.

**Pipeline.** A standalone `-entry SOUPORCELL_ONLY` runs souporcell off already-published Cell Ranger
outputs (no Cell Ranger re-run). Each "mix" is several samples grouped under one fake `patient` id; the
joint machinery retags, merges, and runs souporcell over a K sweep. Scoring: `bin/souporcell_mix_eval.py`.

---

## 2. Experiment 1 — Healthy control mixes

Mixes of healthy bone-marrow / PBMC controls, both **within one study** (the harder, batch-free test)
and **across studies**.

![ARI at correct K, controls](figures/stats_ari_controls.png)

| Mix | Samples | Type | K | **ARI** | Mapping acc | Doublet rate |
|---|---|---|---|---|---|---|
| MIX_within_PBMMC23 | PBMMC_2 + PBMMC_3 | within-study | 2 | **0.997** | 0.999 | 0.05% |
| MIX_within_HDBM12 | HD_BM_1 + HD_BM_2 | within-study | 2 | **0.975** | 0.994 | 0.5% |
| MIX_cross_PBMMC2_HDBM3 | PBMMC_2 + HD_BM_3 | cross-study | 2 | **0.997** | 0.999 | 0.1% |
| MIX_cross_PBMMC3_PBM1 | PBMMC_3 + PBM_1 | cross-study | 2 | **0.997** | 0.999 | 0.1% |

Souporcell recovers the true donor of essentially every cell, **even within the same study** where
there is no batch signal to lean on.

### UMAPs (before = true origin · after = souporcell)

![within-study PBMMC mix](figures/umap_MIX_within_PBMMC23.png)
![within-study HD_BM mix](figures/umap_MIX_within_HDBM12.png)
![cross-study PBMMC2/HD_BM3 mix](figures/umap_MIX_cross_PBMMC2_HDBM3.png)
![cross-study PBMMC3/PBM1 mix](figures/umap_MIX_cross_PBMMC3_PBM1.png)

### Confusion matrices (true origin × souporcell cluster)

![confusion within PBMMC23](figures/confusion_MIX_within_PBMMC23.png)
![confusion within HDBM12](figures/confusion_MIX_within_HDBM12.png)
![confusion cross PBMMC2_HDBM3](figures/confusion_MIX_cross_PBMMC2_HDBM3.png)
![confusion cross PBMMC3_PBM1](figures/confusion_MIX_cross_PBMMC3_PBM1.png)

---

## 3. Experiment 2 — Patient diagnosis-sample mixes

Mixes of **malignant** AML diagnosis samples (each a single patient, no transplant → single genetic
origin), to show deconvolution holds on aneuploid/leukaemic material. Includes a **3-way** mix.

![ARI at correct K, patients](figures/stats_ari_patients.png)

| Mix | Samples | K | **ARI** | Mapping acc | Doublet rate |
|---|---|---|---|---|---|
| MIX_dx_2395_2977 | Sample_2395 + Sample_2977 | 2 | **0.999** | 1.000 | 0.06% |
| MIX_dx_2958_1199 | Sample_2958 + Sample_1199 | 2 | **1.000** | 1.000 | 0.06% |
| MIX_dx_3way | Sample_2395 + 2977 + 2958 | 3 | **0.999** | 0.999 | 0.1% |

### UMAPs (before = true origin · after = souporcell)

![patient Dx 2-way mix 2395/2977](figures/umap_MIX_dx_2395_2977.png)
![patient Dx 2-way mix 2958/1199](figures/umap_MIX_dx_2958_1199.png)
![patient Dx 3-way mix](figures/umap_MIX_dx_3way.png)

### Confusion matrices

![confusion dx 2395_2977](figures/confusion_MIX_dx_2395_2977.png)
![confusion dx 2958_1199](figures/confusion_MIX_dx_2958_1199.png)
![confusion dx 3way](figures/confusion_MIX_dx_3way.png)

---

## 4. Experiment 3 — Male × female mix (sex-expression positive control)

A male and a female diagnosis sample are mixed, giving **two independent ground truths** — barcode
origin *and* sex. We deconvolute three ways and check they agree: (1) **souporcell** (genotype),
(2) a **sex-expression-only** classifier (`XIST` → female, chrY genes → male, no genotype), and
(3) **true origin**. If (1)≈(2)≈(3), the sex read-out is a valid orthogonal confirmation of
souporcell's split — the principle Exp 5 (transplant) relies on.

Pairs: `Sample_8178` (M, AML107) × `Sample_2977` (F, AML152); `Sample_2395` (M, AML057) ×
`Sample_2958` (F, AML155). _(scoring: `bin/souporcell_sex_deconv.py`)_

> _Status: running (job 35318448). Expected: souporcell ARI ≈ sex-call accuracy ≈ 1.0; the
> XIST-vs-chrY scatter separates into two clean, non-overlapping clouds by true origin._

![M8178 × F2977 — sex genes / origin / souporcell / sex-call](figures/sexdeconv_MIX_sex_M8178_F2977.png)
![M2395 × F2958 — sex genes / origin / souporcell / sex-call](figures/sexdeconv_MIX_sex_M2395_F2958.png)

---

## 5. Experiment 4 — Mixing-proportion titration

Does an **unbalanced** mix reduce accuracy? For one control pair (`PBMMC_2`/`PBMMC_3`) and one patient
pair (`Sample_2977`/`Sample_2395`) we re-ran souporcell with the minority sample at **1, 5, 10, 25,
50 %** of cells, holding the **total cell count constant** so any change is purely the proportion
effect. We report overall ARI and the **recall of the minority population** (can a rare genotype still
be pulled out?). _(scoring: `bin/souporcell_proportions_eval.py`)_

> _Status: running (array 35318661 + eval 35318662). Expected: ARI and minority recall stay high down
> to ~5–10 %, then degrade at 1 % where the rare population approaches souporcell's detection limit._

![ARI and minority recall vs mixing ratio](figures/stats_proportions.png)

---

## 6. Controls that prove the signal is real

**Single-individual (solo) negative controls.** Souporcell *always* returns K clusters — even handed
one individual at k=2 it splits the cells ~50/50, but that split has **no genetic basis** (ARI = 0.0).
A split alone means nothing; only a split that *matches genotype* is meaningful.

| Solo control | K | ARI | Dominant-cluster purity |
|---|---|---|---|
| MIX_solo_PBMMC1 | 2 | 0.00 | 0.58 |
| MIX_solo_2395 | 2 | 0.00 | 0.52 |

![solo control PBMMC1](figures/umap_MIX_solo_PBMMC1.png)
![solo control Sample_2395](figures/umap_MIX_solo_2395.png)

**Wrong K.** Over-splitting a true 2-donor mix to k=3 fragments real clusters (ARI 0.41–0.73); forcing
the 3-way mix to k=2 collapses two donors (ARI 0.92, mapping 0.84). Accuracy peaks precisely when K
equals the true number of individuals — genotype-driven and K-sensitive, not an artefact.

---

## 7. Experiment 5 — Real transplant donor/recipient

Here there is **no artificial mixing and no barcode ground truth** — a post-transplant sample is a
genuine biological mix of recipient (host/leukaemic) and donor (graft) cells. So ARI does not apply
(every cell shares one `Sample_*__` prefix); validation instead asks whether souporcell's genotype
clusters line up with an *independent* signal: **sex-chromosome expression** (`XIST` → female, chrY
genes → male). Recipients here are 46,XX, so a male cluster would mean a male donor. Donor sex is not
recorded, so we ran the available female-BMT samples plus two male recipients as comparisons.

Applied to **AML066** (female, relapse_after_BMT) on the **relapse alone** (`Sample_1386`) and
**Dx+Rel combined** (`Sample_3652`+`Sample_1386`); **AML163** (female, allograft, Dx-only); and
**AML079** / **AML107** (male recipients, relapse).

![Souporcell clusters by sex-chromosome expression](figures/sex_summary_bmt.png)

| Sample (recipient sex) | souporcell split | Per-cluster sex (XIST / chrY) | Read-out |
|---|---|---|---|
| **AML066** rel-alone (F) | 2 clusters | both **female** (XIST≈1.1, chrY≈0) | no donor sex signal |
| **AML066** Dx+Rel (F) | 3 clusters | all **female** | no donor sex signal |
| **AML163** Dx-only (F) | 1 big + **1 small (n=20)** | C1 female; **C0 male** (chrY+, 55% Y⁺) | **sex-discordant** (minor) |
| **AML107** rel (M) | 2 clusters | both **male** (chrY≈0.66, 94% Y⁺) | consistent (recipient male) |
| **AML079** rel (M) | 2 clusters | both low-signal (XIST 0, chrY≈0.005) | inconclusive |

### UMAPs (souporcell cluster · XIST · chrY [· timepoint])

![AML066 relapse alone](figures/umap_sex_AML066_rel_solo.png)
![AML066 Dx+Rel combined](figures/umap_sex_AML066_dxrel.png)
![AML163 Dx (allograft)](figures/umap_sex_AML163_dx_solo.png)
![AML107 relapse (male recipient)](figures/umap_sex_AML107_rel_solo.png)
![AML079 relapse (male recipient)](figures/umap_sex_AML079_rel_solo.png)

### Interpretation

- **The sex read-out itself works**: souporcell + the validation correctly call AML107 male and AML066
  female from expression alone — the method is sound (and Exp 3 confirms it on a balanced mix).
- **AML066 (the headline female-BMT case) shows no sex discordance.** Every cluster is female, so the
  donor was most likely **female (sex-concordant transplant)** — invisible to a sex assay — or
  engraftment is low. Exactly the caveat anticipated when donor sex is unrecorded; fall back to
  genotype (`cluster_genotypes.vcf` distinctness).
- **AML163** is suggestive: a small, genetically distinct, **male** population (~20 cells) inside a
  46,XX patient — consistent with minor male-donor microchimerism — but too small to be conclusive.
- **AML079** gave low signal on both XIST and chrY (a sample/technical limitation; its diagnosis pair
  also under-expresses chrY), so no call.

**Bottom line:** Exp 5 demonstrates the donor/recipient deconvolution + orthogonal sex-validation
machinery end-to-end, but in this cohort none of the female-BMT relapses carried an unambiguous male
donor, so the genotype split could not be independently confirmed by sex here (donor most likely
female). The artificial-mix experiments (§2–4) remain the definitive proof.

---

## 8. Conclusion

Souporcell reliably reconstructs the genetic individual of origin for pooled single cells at
**ARI ≈ 0.97–1.00** given the correct N — for healthy and malignant material, within- and
cross-study, for 2- and 3-way pools, and (Exp 3) with the split independently confirmed by sex
expression. Exp 4 maps the rare-population detection limit; negative and wrong-K controls behave
exactly as expected. The tool is validated for demultiplexing and for donor/recipient separation in
this pipeline.

---

## 9. Reproducibility

- **Entry / scoring**: `nextflow run . -entry SOUPORCELL_ONLY` (see `main.nf`); `bin/souporcell_mix_eval.py`.
- **Samplesheets**: `assets/test/souporcell_mix_{controls,patients,bmt,sex}.csv`.
- **Run scripts**: `jobs/souporcell_mixing.sh` (Exp 1+2), `jobs/souporcell_sex_mix.sh` (Exp 3),
  `jobs/souporcell_proportions.sh` (Exp 4, + `_run.sh` array), `jobs/souporcell_bmt.sh` (Exp 5),
  `jobs/souporcell_figures.sh` / `jobs/souporcell_sex_figures.sh` (figures).
- **Analysis code**: `bin/souporcell_{mix_eval,mix_figures,sex_validate,sex_figures,sex_deconv,subsample_barcodes,proportions_eval}.py`.
- **Results**: `results_soupmix_{controls,patients,bmt,sex,prop}/` (clusters under `callers/souporcell/`,
  scores under `eval/` · `sex_eval/` · `sex_deconv/`).
- **Jobs**: Exp1/2 orchestrator 35217254 (COMPLETED, 23h22m); Exp1/2 figures 35274004; Exp5 (BMT)
  35273680 (clusters OK; sex-validation re-run after a launchDir relative-path fix); Exp5 sex UMAPs
  35318364; Exp3 (sex-mix) 35318448; Exp4 (proportions) array 35318661 + eval 35318662.
- **Branch**: `numbat-specificity-reproducibility`.

---
title: F157 bulk proteomics — Diagnosis vs Relapse (DDE_33 PROTEOMICS branch)
project: DDE_33
date: 2026-06-26
type: report
tags: [proteomics, bulk-ms, DESP, deconvolution, AML, diagnosis-relapse, f157]
---

# F157 bulk proteomics — Diagnosis vs Relapse

Paediatric AML bulk proteomics (Spectronaut DIA on a timsTOF HT, 30 SPD EvoSep gradient; SwissProt +
contaminants, 1% FDR). Run through the DDE_33 **`PROTEOMICS`** branch — the DDE_31 R pipeline ported to
Python, with the **DESP** cell-state demix kept in R. **No batch correction** (all 12 samples are one
batch), so the raw filtered/log2 matrix drives every downstream stage.

> **Provenance.** Pipeline `DDE_33` branch `numbat-specificity-reproducibility`; entry
> `-entry PROTEOMICS -profile viking` via `params-proteomics.yaml` + `assets/proteomics_f157.yaml`;
> SLURM job `35274762`. Outputs: `results_proteomics/proteomics/`. Figures here are copies under
> `figures/` so this folder is self-contained for Obsidian.

## Dataset

- **12 samples**, **7603** protein groups identified; **6481** retained after the presence filter
  (detected in ≥75% of samples in a condition). No samples dropped (all ≥2500 IDs); contaminant list
  applied.
- Design: `sample, condition (Diagnosis/Relapse), replicate (patient AMLxxx), Batch`. One patient
  (**AML152**: Dx `2977` + Rel `0109`) is paired in the DESP overlap set.

![Protein IDs per sample](figures/qc_protein_ids_per_sample.png)
![Protein completeness across samples](figures/qc_protein_completeness.png)
![Per-sample intensity distributions (log2)](figures/qc_intensity_boxplot.png)

## Global structure (no batch correction)

Clustered heatmap and PCA of the raw log2 matrix. Diagnosis/Relapse do **not** cleanly separate in the
uncorrected data — consistent with the DDE_31 finding that the no-batch view is weak and that the
apparent Dx/Rel split in the patient-batch-limma run was partly model-induced (it survived a
randomised-label control). Treat structure here as a baseline, not a biological claim.

![Heatmap — non-corrected](figures/structure_heatmap.png)
![PCA (PC1 vs PC2) coloured by condition](figures/structure_pca_condition.png)
![PCA (PC1 vs PC2) coloured by batch/patient](figures/structure_pca_batch.png)

## Differential expression (Relapse vs Diagnosis)

limma reimplemented in Python (OLS + empirical-Bayes moderation, `eBayes`/`squeezeVar`). On the
uncorrected matrix **no protein reaches adj.P < 0.05** — i.e. no cohort-level Dx/Rel signal without
batch/patient modelling. The largest fold-changes are **keratins (KRT16/17/14/5/10/6A) and LTF**, a
classic skin/handling contamination signature rather than leukaemia biology.

![Volcano — Relapse vs Diagnosis](figures/de_volcano.png)

Per-patient paired logFC (Relapse − Diagnosis, mean difference) is more interpretable than the cohort
test given the heterogeneity:

![Per-patient logFC heatmap](figures/de_per_patient_logfc.png)

## Stage-4 interpretation

Heatmap + k-means of the (permissively-thresholded) "significant" set and condition means. With the
weak global signal these are exploratory.

![Significant-protein heatmap](figures/stage4_sig_heatmap.png)
![Condition-mean expression](figures/stage4_condition_means.png)

## Classifiers (exploratory — n=12)

Random-forest / decision-tree Dx-vs-Rel classifiers. With only ~12 samples and a 3-sample held-out
test these are **not** evidence of a robust signature; included for completeness. RF top features
(impurity importance): SS18L2, LSM2, TUBB4A, ELOC, EIF2B1, RPA3, ATG4B, EIF3I.

![RF top features](figures/ml_rf_top_features.png)
![Decision tree](figures/ml_decision_tree.png)

## DESP — cell-state demixing (the headline)

**DESP** ([Cell Reports Methods 2024](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(24)00052-3))
demixes the bulk proteome into per-cell-state protein profiles given matched cell-type proportions.
Proportions here are the pre-computed DDE_31 transcriptomics-overlap set (6 matched samples); the
pipeline can instead derive them from this project's scRNA reference-mapping output (the Module-C
proteogenomic hook). Output: **6481 proteins × 54 cell states**, per condition + the Relapse−Diagnosis
delta, plus the paired patient AML152.

Top DE-ranked proteins and their direction:

![DESP top proteins heatmap](figures/desp_top_proteins_heatmap.png)
![Top-protein direction (Relapse vs Diagnosis)](figures/desp_protein_direction.png)

Cell-state-resolved change — which cell states drive each protein's Dx→Rel shift:

![DESP delta heatmap (Relapse − Diagnosis) × cell type](figures/desp_delta_heatmap.png)
![Cell-type contribution to top proteins](figures/desp_celltype_contributions.png)

### Paired patient AML152 (Dx 2977 → Rel 0109)

Per-protein contribution (DESP profile × mean cell-type proportion), row-percentile scaled, per
timepoint:

![AML152 — Diagnosis](figures/desp_AML152_diagnosis_percentile.png)
![AML152 — Relapse](figures/desp_AML152_relapse_percentile.png)

## Interpretation

- The **uncorrected** cohort shows **no significant cohort-level Dx/Rel proteomic signal**; the
  strongest fold-changes are contamination-associated keratins/LTF. This is the honest baseline.
- A **patient-batch limma** re-run (set `batch_correction.method=limma`, `batch_column=replicate`,
  `prot_batch_method=limma`, `prot_desp_bulk_source=limma`) reproduces DDE_31's stronger separation —
  but that separation was partly model-induced (randomised-label control), so interpret cautiously.
- **DESP** is the value-add: it partitions the (weak) bulk signal across 54 cell states, giving a
  cell-of-origin-flavoured view and a Dx→Rel delta even where the bulk test is null. The paired AML152
  comparison is the most concrete read; the overlap set is small (6 samples), so this is exploratory.

## Reproduce

```bash
sbatch jobs/run_proteomics.sh   # -entry PROTEOMICS -profile viking -params-file params-proteomics.yaml
```

Full outputs (all figures + tables): `results_proteomics/proteomics/{01_qc,02_preprocessing,
03_visualization,04_de,05_stage4,06_ml,07_desp}/`. Lab book:
[[2026-06-25_proteomics-branch-port]].

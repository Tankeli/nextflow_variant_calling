# Implementation Plan — DNA-free clonal tracing of AML from diagnosis to relapse

**Goal:** identify leukemic stem cells (LSCs) and trace clones from diagnosis to relapse in pediatric AML using single-cell CITE-seq (RNA + ADT) and matched bulk proteomics, **without matched DNA-seq**, on a whitelisted server under an approved UKRI ethics protocol.

This document is the engineering and analysis plan for the repository. It assumes the science context established previously: pediatric AML is fusion-driven and copy-number-quiet, adult LSC signatures transfer imperfectly, and the relapse samples are the strongest available substitute for a functional LSC assay.

---

## 0. Scope, assumptions, and data governance

**In scope:** ingestion → QC → malignant-vs-normal classification → clonal structure → LSC identification → diagnosis→relapse clonal tracing, plus the four novel modules from the design discussion (multi-evidence classifier, surface-proteome anomaly detection, proteogenomic anchoring, relapse-supervised fate labeling).

**Assumptions to confirm before coding (these change tool choices):**
- CITE-seq chemistry (10x 3′ vs 5′; 5′ gives better variant/fusion coverage and is strongly preferred for this project).
- Whether any patients have multiome/ATAC or long-read data (changes mito and fusion feasibility dramatically).
- Bulk proteomics acquisition mode (DDA vs DIA) and depth — determines whether proteogenomic fusion-peptide detection (Module C) is viable at all.
- Per-patient known driver lesions/fusions from the diagnostic workup (these become your supervised malignancy anchors).
- Availability of a pediatric healthy bone-marrow reference (most public atlases are adult).

**Data governance (non-negotiable, given human patient data):**
- **No patient data in Git, ever** — not raw, not processed, not derived count matrices that could be re-identifying. Code in Git; data on the whitelisted server only.
- Enforce with a strict `.gitignore` (`data/`, `*.h5`, `*.h5ad`, `*.bam`, `*.fastq*`, `*.mtx`, proteomics raw) and a pre-commit hook that blocks large/binary files (`pre-commit` + `check-added-large-files`).
- Track data with **DVC** or **git-annex** pointing at server storage, so the *pipeline* is versioned without the *data* leaving the enclave.
- Keep a `PROVENANCE.md` linking each dataset to the ethics reference, consent scope, and the approved use. Plan eventual controlled-access deposition (e.g., EGA, which fits UK/UKRI governance) rather than open repositories.
- CI runs only on **synthetic fixtures**, never on patient data.

---

## 1. Repository and infrastructure

### Repository layout

```
aml-clonal-tracing/
├── README.md
├── LICENSE                      # e.g. MIT for code; data separate
├── CITATION.cff
├── PROVENANCE.md                # ethics ref, consent scope, data lineage
├── environment.yml              # conda/mamba (Python side)
├── renv.lock                    # R dependencies (numbat, Seurat, etc.)
├── .gitignore                   # blocks all data/binaries
├── .pre-commit-config.yaml
├── config/
│   ├── config.yaml              # master config (paths, params, patients)
│   └── samples.tsv              # sample sheet: patient, timepoint, subtype, modality
├── workflow/
│   ├── Snakefile
│   └── rules/                   # one rule file per phase
│       ├── 01_ingest_qc.smk
│       ├── 02_reference_map.smk
│       ├── 03_barcodes.smk
│       ├── 04_classifier.smk
│       ├── 05_surface_anomaly.smk
│       ├── 06_proteogenomics.smk
│       ├── 07_lsc.smk
│       ├── 08_relapse_ot.smk
│       └── 09_validation.smk
├── src/aml_ct/                  # installable package (pip install -e .)
│   ├── io/                      # loaders, sample-sheet parsing
│   ├── qc/
│   ├── reference/               # projection, fusion calling wrappers
│   ├── barcodes/                # mito, expressed variants
│   ├── classifier/              # Module A: multi-evidence malignancy model
│   ├── surface/                 # Module B: ADT anomaly detection
│   ├── proteomics/              # Module C: bulk MS + proteogenomic anchor
│   ├── lsc/                     # stemness scoring, trajectory
│   ├── fate/                    # Module D: relapse-supervised OT
│   └── viz/
├── notebooks/                   # exploratory only; outputs stripped before commit
├── tests/                       # pytest, synthetic fixtures
├── containers/                  # Apptainer/Singularity defs
└── docs/                        # mkdocs-material
```

### Environment and reproducibility

- **Package/env:** `mamba` for the Python stack; `renv` for R (numbat, Seurat, BoneMarrowMap, scDblFinder live in R). Pin everything; commit `environment.yml` and `renv.lock`.
- **Containers:** build **Apptainer/Singularity** images (HPC-friendly, rootless) per major stage so the pipeline is portable across the server's nodes. Keep one image for the Python/scvi-tools/moscot GPU stack and one for the R stack.
- **Workflow manager:** **Snakemake** (Python-native, integrates with the package and config). If you prefer community-maintained sub-pipelines, wrap `nf-core/scrnaseq` (CITE-seq) and `nf-core/quantms` (proteomics) as upstream steps. One rule file per phase, all parameters from `config/config.yaml`.
- **Determinism:** set and log a global seed; record tool versions and the git SHA into every output's metadata; never hardcode paths (all via config).
- **Compute:** GPU node for `totalVI`, `scArches`, and `moscot` (JAX/OTT is GPU-accelerated); CPU for the rest. Snakemake resource directives per rule.
- **CI (GitHub Actions):** lint (`ruff`, `black`), type-check, `pytest` on synthetic fixtures, build docs. No data.

---

## 2. Phase roadmap

| Phase | Output | Novel module |
|---|---|---|
| P1 Ingest & QC | Clean CITE-seq objects + proteomics matrix | — |
| P2 Reference & malignant backbone | Per-cell normal-projection + fusion/CNV evidence | — |
| P3 DNA-free barcodes | Mito + expressed-variant evidence per cell | — |
| P4 Multi-evidence classifier | Per-cell malignancy posterior + clone labels | **Module A** |
| P5 Surface anomaly | ADT out-of-distribution / LAIP-analogue score | **Module B** |
| P6 Proteogenomic anchor | Sample-level malignant protein anchor + consistency | **Module C** |
| P7 LSC identification | LSC-like compartment within malignant cells | — |
| P8 Relapse fate labeling | Diagnosis cells scored by relapse-seeding potential | **Module D (capstone)** |
| P9 Validation | Cross-patient, cross-modal, LAIP/MRD bridge | — |
| P10 Packaging | Installable package, figures, docs | — |

The phases are layered: P2–P6 each produce an independent evidence channel; P4 fuses the malignancy channels; P7 finds stemness within malignant cells; P8 converts "malignant + stem-like" into "relapse-fated = LSC."

---

## 3. Phase detail

### P1 — Ingestion and QC

**CITE-seq.** Run Cell Ranger `multi` (RNA + Feature Barcoding) if starting from FASTQ; otherwise ingest filtered matrices. Then:
- Ambient RNA removal with **CellBender** (`remove-background`).
- Doublets with **scDblFinder** (per sample, before integration).
- RNA QC: gene/UMI/mito-% thresholds set per sample from the distribution, not fixed cutoffs.
- ADT: normalize with **DSB** (uses empty droplets + isotype controls to remove ambient + technical noise — important because ambient ADT badly distorts surface phenotypes); keep isotype controls as features for the anomaly model later.
- Store as `AnnData` (`.h5ad`) with RNA in `X`, ADT in `obsm`/a paired layer (`muon`/`MuData`).

**Bulk proteomics.** From raw: **DIA-NN** (DIA) or **FragPipe**/**MaxQuant** (DDA). Then log-transform, normalize (median/quantile), assess and handle missingness (MNAR-aware), batch-correct (track batch in the sample sheet). Keep the search settings — Module C will re-search with a custom database.

**Harmonization.** A single `samples.tsv` keyed by `patient × timepoint × modality × subtype` is the backbone every downstream rule reads.

**Exit criteria:** QC'd per-sample CITE-seq objects, a proteomics abundance matrix, and a validated sample sheet.

### P2 — Normal reference and malignant-vs-normal backbone

Three independent evidence streams, plus integration:

1. **Reference projection (RNA + protein).** Build the joint CITE-seq embedding with **totalVI** (denoises ADT, gives a joint latent), then map onto a normal hematopoietic reference. Use **BoneMarrowMap** (a hematopoietic atlas built to nominate normal-marrow equivalents of AML LSC states) and the **Triana et al. 132-antibody CITE-seq atlas** for the protein side; project with **Symphony**/**scArches**. Cells deviating from the normal differentiation manifold are malignant candidates. **Build/adapt a pediatric reference** where possible — flag age-mismatch as a known limitation otherwise.
2. **Fusion detection (primary malignant tag in pediatric AML).** Short-read scRNA fusion calling is weakly tooled; use a two-pronged approach: pseudobulk per cluster → **Arriba**/**STAR-Fusion** to nominate fusions, then interrogate each cell for junction-spanning reads at the known breakpoint to assign fusion⁺ status per cell. If long-read exists, use **CTAT-LR-fusion**/**scNanoGPS** — far better. A fusion⁺ cell is a near-unambiguous malignancy call requiring no DNA.
3. **Allele/CNV inference.** Run **numbat** (needs `cellsnp-lite` pileup + population phasing with Eagle/1000G) for patients with copy-number lesions. Expect informative results only where cytogenetics are abnormal (monosomy 7, trisomy 8, complex karyotype); treat as null elsewhere.

**Exit criteria:** per-cell projection-deviation score, per-cell fusion status, per-cell CNV-clone assignment where available.

### P3 — DNA-free clonal barcodes

- **Mitochondrial — feasibility first.** Before relying on it, *quantify* informative mtDNA variant recovery from your CITE-seq libraries with **maegatk/mgatk**; call informative variants with **MQuad**; cluster with **vireoSNP**. Standard 3′ CITE-seq often yields too few variants — if so, treat mito as weak/soft evidence and consider adding **MAESTER** enrichment or an mtscATAC arm on the relapse pairs (highest-value samples).
- **Expressed nuclear variants.** Genotype known driver SNVs at the per-cell level with **VarTrix**/**cb_sniffer** at sites nominated from the diagnostic workup. De novo calling from sparse data is unreliable — restrict to known/recurrent sites and highly expressed genes.

**Exit criteria:** per-cell mito variant matrix (with a documented confidence assessment) and per-cell genotype calls at known sites.

### P4 — Module A: multi-evidence malignancy/clone classifier

**Gap addressed:** in CN-quiet, SNV-poor pediatric AML, every DNA-free signal is individually weak; CloneTracer fuses nuclear+mito but assumes usable variant signal.

**Design:** a probabilistic model producing a per-cell posterior over {normal, pre-leukemic, leukemic} and a clone label, fusing independent weak evidence channels:
- fusion⁺ junction reads (strong when present),
- expressed-variant genotype at known sites,
- mito heteroplasmy (soft, given the ReDeeM low-heteroplasmy caveat — propagate uncertainty, don't hard-call),
- RNA reference-projection deviation (P2),
- ADT surface-anomaly score (Module B),
- a normal-clonality null term (deviation from the clonal-architecture statistics that the hematopoiesis lineage-tracing literature established for healthy HSCs).

Implement as a graphical/Bayesian model (e.g., `numpyro`) or a calibrated gradient-boosted ensemble over the channels, anchored by the unambiguous fusion⁺ cells as semi-supervision. Output per-cell posteriors with credible intervals; the explicit "uncertain boundary" class is a feature, not a bug.

**Exit criteria:** per-cell malignancy posterior + clone assignment, with channel-wise attribution for interpretability.

### P5 — Module B: surface-proteome anomaly detection

**Gap addressed:** existing work nominates aberrant markers post hoc; none trains the normal surface proteome as a null and flags malignancy as out-of-distribution on protein.

**Design:** train a one-class/OOD model (Mahalanobis distance on the `totalVI` ADT latent, or a normal-only VAE's reconstruction error, or `scArches` reference-query distance) on normal-marrow CITE-seq. Score query cells; HSPC-like cells that are anomalous on the protein manifold carry an aberrant combinatorial surface state — the single-cell analogue of the leukemia-associated immunophenotype (LAIP) that clinical MRD flow already uses. This both feeds Module A and provides an external validation bridge (LAIP aberrancy is functionally tied to leukemia-initiating cells). Works even when genetics are silent.

**Exit criteria:** per-cell surface-anomaly score and the aberrant marker combinations driving it.

### P6 — Module C: proteogenomic anchoring (bulk MS)

**Gap addressed:** the matched bulk proteome is usually relegated to validation; it can provide an orthogonal, protein-level malignancy anchor.

**Design (feasibility-gated — confirm acquisition depth first):**
- Build a **custom search database** augmented with fusion-junction peptides (from P2 fusions) and known variant peptides; re-search the bulk MS (FragPipe open/custom search) for fusion/aberrant peptides as a sample-level malignant anchor.
- Derive a protein signature of the malignant clone and enforce a **consistency constraint**: the inferred malignant-cell fraction × their pseudo-protein profile should reconstruct the bulk proteome (NNLS/`BayesPrism`-style, with explicit caveats that protein deconvolution is harder than RNA).
- Optionally exploit **RNA–protein discordance** in CITE-seq (stem-like transcriptome vs committed surface phenotype) as a candidate active-LSC feature.

**Exit criteria:** sample-level malignant protein anchor; a bulk↔single-cell consistency metric.

### P7 — LSC identification within malignant cells

Restricting to cells called malignant by Module A:
- Score stemness with **pLSC6** (pediatric) and LSC17 (adult, for comparison) via **AUCell**/**UCell** (rank-based, robust to dropout).
- Order the hierarchy with **CellRank 2** (which consumes moscot kernels) and **CytoTRACE2**; identify the HSC-like/progenitor-like apex (van Galen-style states).
- Cross-check surface LSC phenotype from ADT, **without** hard-gating CD34⁺CD38⁻ (CD34-negative LSCs occur in pediatric subtypes). Stratify by molecular subtype.

**Exit criteria:** an LSC-like compartment per patient/subtype with multi-modal support.

### P8 — Module D (capstone): relapse-supervised fate labeling

**Gap addressed:** paired studies *observe* stem-cell enrichment at relapse but don't use the pair as a supervised label. The relapse sample is effectively a free prospective fate readout.

**Design with moscot:**
- For each patient with a pair, set up a **`TemporalProblem`** in **moscot** mapping diagnosis → relapse on the joint CITE-seq embedding, using **unbalanced** OT (tune the marginal-relaxation `tau` parameters) because therapy is a mass-destroying bottleneck.
- Compute ancestor probabilities: which diagnosis cells are the likely ancestors of the relapse-persistent leukemic population. Those high-ancestry, malignant, stem-like diagnosis cells are the **LSC candidates with quasi-functional support**.
- Because relapse cells have *evolved*, target the shared persistent **core** state (the convergent-evolution signal), not an identity match; encode this in the OT cost (e.g., emphasize a conserved subspace).
- **Cross-patient generalization:** learn a "relapse-competent LSC" signature pooled across pairs (subtype-stratified), then apply it to **diagnosis-only** patients who lack a relapse sample. Validate by held-out pairs.

**Exit criteria:** per-diagnosis-cell relapse-seeding score; a transferable LSC signature applied to diagnosis-only patients. **This is the project end goal.**

### P9 — Validation framework

- **Cross-patient hold-out** for Module D and the LSC signature (train on some pairs, predict held-out pairs).
- **Cross-modal consistency:** RNA-based, ADT-based, and proteomic malignancy calls should agree on high-confidence cells; quantify and report disagreement on the boundary population rather than hiding it.
- **External bridge:** compare the surface-anomaly (Module B) calls against clinical LAIP/MRD flow markers where available.
- **Negative controls:** healthy-donor marrow should yield ~no malignant calls and normal clonal architecture.
- **Ablations:** show each evidence channel's marginal contribution to Module A.
- **Statistics:** mixed models with patient as a random effect (cells are not independent); multiple-testing control throughout.

### P10 — Outputs and packaging

Installable `aml_ct` package; per-patient HTML reports; publication figures (UMAPs of malignancy posterior, clone trees, OT ancestor flows, LSC signatures); `mkdocs` docs; `CITATION.cff`. Plan controlled-access data deposition (EGA) at publication.

---

## 4. Software stack

| Layer | Tools |
|---|---|
| Counting | Cell Ranger `multi`; (nf-core/scrnaseq) |
| Ambient / doublets | CellBender; scDblFinder |
| ADT normalization | DSB; isotype controls retained |
| Joint embedding | totalVI / scvi-tools; muon/MuData |
| Reference mapping | BoneMarrowMap; Triana CITE-seq atlas; Symphony; scArches |
| Fusion | Arriba / STAR-Fusion (pseudobulk); CTAT-LR-fusion (long-read); per-cell junction counting |
| CNV / allele | numbat; cellsnp-lite; Eagle phasing |
| Mito barcodes | maegatk/mgatk; MQuad; vireoSNP; (MAESTER if enriched) |
| Expressed variants | VarTrix; cb_sniffer |
| Classifier (Mod A) | numpyro / scikit-learn / xgboost |
| Surface anomaly (Mod B) | totalVI latent + Mahalanobis / one-class / VAE |
| Proteomics | DIA-NN / FragPipe / MaxQuant; limma/MSstats; BayesPrism (caveated) |
| LSC / trajectory | AUCell/UCell; CellRank 2; CytoTRACE2 |
| Relapse OT (Mod D) | moscot (TemporalProblem, unbalanced); CellRank 2 kernels |
| Infra | mamba; renv; Snakemake; Apptainer; DVC; pytest; GitHub Actions; mkdocs |

---

## 5. Suggested sequencing of work

1. **Foundations:** repo scaffold, env/containers, data governance, sample sheet, CI on fixtures.
2. **P1–P2:** ingestion, QC, reference projection, fusion detection — gets you a first-pass malignant-vs-normal call quickly.
3. **Feasibility checks:** quantify mito recovery (P3) and proteomics depth (P6 gate) early, since both can kill a downstream module; reprioritize before investing.
4. **P4 + P5:** the malignancy classifier and surface-anomaly module (most of the methodological value, runs on existing data).
5. **P7:** LSC identification.
6. **P8:** the moscot relapse capstone on the pairs, then transfer to diagnosis-only patients.
7. **P9–P10:** validation hardening, packaging, docs.

---

## 6. Risk register

| Risk | Likelihood | Mitigation |
|---|---|---|
| Mito signal too sparse in CITE-seq | High | Quantify first; demote to soft evidence; add MAESTER/mtscATAC on relapse pairs only |
| CN-quiet karyotypes defeat numbat | High (pediatric) | Lean on fusion + surface anomaly; treat numbat as a per-patient bonus |
| Proteogenomic fusion peptides undetectable at acquisition depth | Medium | Gate Module C on a feasibility test; fall back to consistency/validation use only |
| No pediatric healthy reference | Medium | Build/adapt one; report age-mismatch as a limitation; subtype-stratify |
| OT ancestor inference confounded by relapse evolution | Medium | Target conserved core subspace; unbalanced OT; held-out validation |
| Boundary cells unclassifiable without DNA | Certain | Explicit "uncertain" class; report rather than force-call |
| Patient data leakage to Git | Low/Catastrophic | `.gitignore` + pre-commit large-file block; DVC; CI on synthetic only |

---

## Appendix — starter config and workflow stubs

**`config/config.yaml`**
```yaml
project: aml-clonal-tracing
seed: 1312
paths:
  data_root: /srv/whitelisted/aml_ct/data      # server only, never in git
  outputs: /srv/whitelisted/aml_ct/outputs
  reference_bmm: /srv/whitelisted/refs/bonemarrowmap
  phasing_panel: /srv/whitelisted/refs/1000G_eagle
samples_sheet: config/samples.tsv
qc:
  min_genes: 200
  max_mito_pct: 15        # set per-sample from distributions in practice
adt:
  normalization: dsb
  isotype_controls: [IgG1, IgG2a, IgG2b]
classifier:
  classes: [normal, preleukemic, leukemic]
  anchor_on_fusion: true
ot:
  unbalanced_tau_a: 0.9
  unbalanced_tau_b: 0.7
  conserved_subspace: true
```

**`config/samples.tsv`**
```
patient	timepoint	subtype	modality	fusion
P01	diagnosis	KMT2A-r	citeseq	KMT2A-MLLT3
P01	relapse	KMT2A-r	citeseq	KMT2A-MLLT3
P02	diagnosis	CBFA2T3-GLIS2	citeseq	CBFA2T3-GLIS2
```

**`workflow/Snakefile` (skeleton)**
```python
configfile: "config/config.yaml"
import pandas as pd
samples = pd.read_table(config["samples_sheet"])

include: "rules/01_ingest_qc.smk"
include: "rules/02_reference_map.smk"
include: "rules/03_barcodes.smk"
include: "rules/04_classifier.smk"
include: "rules/05_surface_anomaly.smk"
include: "rules/06_proteogenomics.smk"
include: "rules/07_lsc.smk"
include: "rules/08_relapse_ot.smk"
include: "rules/09_validation.smk"

rule all:
    input:
        expand("{out}/fate/{patient}_relapse_seeding.h5ad",
               out=config["paths"]["outputs"],
               patient=samples.query("timepoint=='diagnosis'").patient.unique())
```

**`.gitignore` (critical lines)**
```
data/
*.h5
*.h5ad
*.bam
*.fastq*
*.mtx
*.raw            # proteomics
outputs/
.dvc/cache/
```

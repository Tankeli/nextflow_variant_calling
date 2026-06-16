//
// Top-level workflow. Wires the samplesheet through Cell Ranger and the variant callers.
// Phase 1: input parsing + channel construction.   Phase 2: Cell Ranger.
// Caller subworkflows (Numbat / CopyKAT / souporcell) land in Phases 3-4.
//

include { INPUT_CHECK   } from '../subworkflows/local/input_check'
include { CELLRANGER    } from '../subworkflows/local/cellranger'
include { NUMBAT        } from '../subworkflows/local/numbat'
include { COPYKAT_WF    } from '../subworkflows/local/copykat'
include { COPYKAT_ROBUSTNESS_WF } from '../subworkflows/local/copykat_robustness'
include { SOUPORCELL_WF } from '../subworkflows/local/souporcell'
include { CLONETRACER_WF } from '../subworkflows/local/clonetracer'
include { ANNOTATION    } from '../subworkflows/local/annotation'
include { RNA_CORE        } from '../subworkflows/local/rna_core'
include { RNA_INTEGRATION } from '../subworkflows/local/rna_integration'
include { RNA_ADVANCED    } from '../subworkflows/local/rna_advanced'
include { PROTEIN         } from '../subworkflows/local/protein'
include { PLOT_COPYKAT    } from '../modules/local/plot_copykat'
include { PLOT_SOUPORCELL } from '../modules/local/plot_souporcell'
include { PLOT_CLONETRACER } from '../modules/local/plot_clonetracer'
include { COHORT_SUMMARY  } from '../modules/local/cohort_summary'
include { INTEGRATION     } from '../subworkflows/local/integration'

workflow VARIANTCALLING {

    ch_versions = Channel.empty()

    INPUT_CHECK( file(params.input, checkIfExists: true) )
    ch_samples = INPUT_CHECK.out.samples

    //
    // Cell Ranger multi (FASTQ -> per-sample BAM + filtered matrix)
    //
    // fb_reference is optional (GEX-only samples have no Antibody Capture)
    def fb_ref = params.fb_reference ? file(params.fb_reference, checkIfExists: true) : []
    CELLRANGER(
        ch_samples,
        file(params.cellranger_index, checkIfExists: true),
        fb_ref
    )
    ch_aln      = CELLRANGER.out.aln   // [ meta, bam, bai, filtered_feature_bc_matrix ]
    ch_versions = ch_versions.mix(CELLRANGER.out.versions)

    // Per-patient grouping of alignments (drives Numbat-joint + souporcell-paired).
    // Sort Dx before Rel so sample order is deterministic across callers.
    ch_patient_aln = ch_aln
        .map { meta, bam, bai, mtx -> tuple( meta.patient, [ meta: meta, bam: bam, bai: bai, mtx: mtx ] ) }
        .groupTuple()
        .map { patient, members -> tuple( patient, members.sort { it.meta.timepoint == 'Dx' ? 0 : 1 } ) }

    //
    // Numbat (CNV / clones), joint per patient
    //
    ch_numbat = Channel.empty()
    if (params.run_numbat) {
        NUMBAT( ch_patient_aln )
        ch_numbat   = NUMBAT.out.numbat
        ch_versions = ch_versions.mix(NUMBAT.out.versions)
    }

    //
    // CopyKAT (aneuploid/diploid gate), per sample
    //
    ch_copykat = Channel.empty()
    if (params.run_copykat) {
        COPYKAT_WF( ch_aln )
        ch_copykat  = COPYKAT_WF.out.prediction
        ch_versions = ch_versions.mix(COPYKAT_WF.out.versions)
    }

    //
    // Souporcell (SNV clusters), joint per patient over a K sweep
    //
    ch_souporcell = Channel.empty()
    if (params.run_souporcell) {
        SOUPORCELL_WF( ch_patient_aln, file(params.souporcell_fasta, checkIfExists: true) )
        ch_souporcell = SOUPORCELL_WF.out.clusters
        ch_versions   = ch_versions.mix(SOUPORCELL_WF.out.versions)
    }

    //
    // CloneTracer (downstream clonal integration), joint per patient. Synthesises per-cell M/N
    // over CNV (Numbat) + nuclear-SNV (souporcell) + a new per-sample mtDNA pileup, builds the
    // per-patient JSON and runs the Bayesian model -> clone trees + per-cell clone posteriors.
    //
    ch_clonetracer = Channel.empty()
    if (params.run_clonetracer) {
        // GTF is optional: it only powers the CNV pseudo-counts. If absent, CloneTracer simply
        // builds its JSON from the SNV + mtDNA axes (build script logs the skipped source).
        def ct_gtf = []
        if (params.clonetracer_gtf) {
            def gtf_f = file(params.clonetracer_gtf)
            if (gtf_f.exists()) { ct_gtf = gtf_f }
            else { log.warn "clonetracer_gtf not found (${params.clonetracer_gtf}); CloneTracer CNV axis disabled" }
        }
        CLONETRACER_WF( ch_aln, ch_patient_aln, ch_numbat, ch_souporcell, ct_gtf )
        ch_clonetracer = CLONETRACER_WF.out.assignments
        ch_versions    = ch_versions.mix(CLONETRACER_WF.out.versions)
    }

    //
    // Annotation branch (scanpy QC -> reference mapping), parallel to the callers
    //
    ch_qc = Channel.empty()
    ch_celltypes = Channel.empty()
    ch_mapped = Channel.empty()
    if (params.run_qc || params.run_reference_mapping) {
        def atlas = params.run_reference_mapping ? file(params.refmap_atlas, checkIfExists: true) : []
        def refumap = (params.run_reference_mapping && params.refmap_umap) ? file(params.refmap_umap, checkIfExists: true) : []
        ANNOTATION( ch_aln, params.run_reference_mapping, atlas, refumap )
        ch_qc        = ANNOTATION.out.qc
        ch_celltypes = ANNOTATION.out.celltypes
        ch_mapped    = ANNOTATION.out.mapped
        ch_versions  = ch_versions.mix(ANNOTATION.out.versions)
    }

    //
    // RNA downstream best-practices stack (ported from DDE_27), parallel to the callers and to
    // the lightweight ANNOTATION branch. Runs off the Cell Ranger filtered matrices:
    //   RNA_CORE        QC -> normalize -> feature-select -> dimred -> cluster -> annotate (per sample)
    //   RNA_INTEGRATION cohort integration (scVI / scANVI / BBKNN + scib)            [needs run_rna_core]
    //   RNA_ADVANCED    pseudotime / velocity / DE / composition (gated individually)[needs run_rna_core]
    //   PROTEIN         surface-protein / ADT branch (CITE-seq only)
    // All default OFF; enable per run. None of these gate the variant callers.
    //
    ch_rna_mtx = ch_aln.map { meta, bam, bai, mtx -> tuple( meta, mtx ) }

    ch_rna_annotated = Channel.empty()
    if (params.run_rna_core) {
        RNA_CORE( ch_rna_mtx )
        ch_rna_annotated = RNA_CORE.out.annotated
        ch_versions      = ch_versions.mix(RNA_CORE.out.versions)
    }

    ch_rna_integrated = Channel.empty()
    if (params.run_rna_integration) {
        if (!params.run_rna_core) {
            error "run_rna_integration requires run_rna_core (it consumes the annotated objects)"
        }
        RNA_INTEGRATION( ch_rna_annotated )
        ch_rna_integrated = RNA_INTEGRATION.out.integrated
        ch_versions       = ch_versions.mix(RNA_INTEGRATION.out.versions)
    }

    if (params.run_pseudotime || params.run_velocity || params.run_de || params.run_composition) {
        if (!params.run_rna_core) {
            error "RNA advanced stages (pseudotime/velocity/de/composition) require run_rna_core"
        }
        if (params.run_composition && !params.run_rna_integration) {
            error "run_composition requires run_rna_integration (it operates on the integrated object)"
        }
        // Velocity looms (optional): ${params.velocity_loom_dir}/<sample>.loom; samples without
        // one are dropped by the join in RNA_ADVANCED.
        ch_loom = Channel.empty()
        if (params.run_velocity) {
            if (!params.velocity_loom_dir) {
                error "run_velocity requires velocity_loom_dir (this pipeline does not generate velocyto looms)"
            }
            ch_loom = ch_rna_mtx
                .map { meta, mtx -> tuple( meta.id, file("${params.velocity_loom_dir}/${meta.id}.loom") ) }
                .filter { id, loom -> loom.exists() }
        }
        RNA_ADVANCED(
            ch_rna_annotated,
            ch_loom,
            ch_rna_integrated,
            params.run_pseudotime,
            params.run_velocity,
            params.run_de,
            params.run_composition
        )
        ch_versions = ch_versions.mix(RNA_ADVANCED.out.versions)
    }

    if (params.run_protein) {
        PROTEIN( ch_rna_mtx )
        ch_versions = ch_versions.mix(PROTEIN.out.versions)
    }

    //
    // CopyKAT robustness sweep (separate analysis track; standalone Python runs downstream over the
    // published combos via jobs/run_copykat_robustness.sh). Needs the cellranger matrices + the
    // reference-mapped cell types (for the optional known-normal baseline arm).
    //
    if (params.run_copykat_robustness) {
        if (!params.run_copykat) {
            error "run_copykat_robustness requires run_copykat"
        }
        if (params.copykat_robustness_use_norm_ref.any { it as boolean } && !params.run_reference_mapping) {
            error "copykat_robustness_use_norm_ref=true requires run_reference_mapping (for the normal baseline)"
        }
        ch_ck_mtx = ch_aln.map { meta, bam, bai, mtx -> tuple( meta, mtx ) }
        COPYKAT_ROBUSTNESS_WF( ch_ck_mtx, ch_celltypes )
        ch_versions = ch_versions.mix(COPYKAT_ROBUSTNESS_WF.out.versions)
    }

    //
    // Visualisation — diagnostic figures over the published checkpoints. Each reads an
    // existing output in the relevant container; the reference-space overlays need the
    // mapped h5ads, so they are gated on run_reference_mapping.
    //
    if (params.run_reference_mapping && params.run_copykat) {
        // Per-sample CopyKAT call overlaid on the reference-map UMAP (join by sample id).
        ch_ck_plot = ch_mapped
            .map { meta, h5ad -> tuple( meta.id, meta, h5ad ) }
            .join( ch_copykat.map { meta, pred -> tuple( meta.id, pred ) } )
            .map { id, meta, h5ad, pred -> tuple( meta, h5ad, pred ) }
        PLOT_COPYKAT( ch_ck_plot )
        ch_versions = ch_versions.mix(PLOT_COPYKAT.out.versions)
    }

    if (params.run_reference_mapping && params.run_souporcell) {
        // Per-patient souporcell clones on a joint reference-space UMAP (Dx/Rel or standalone).
        ch_mapped_by_patient = ch_mapped
            .map { meta, h5ad -> tuple( meta.patient, h5ad ) }
            .groupTuple()
        ch_soup_plot = SOUPORCELL_WF.out.clusters
            .filter { meta, k, dir -> k == params.souporcell_plot_k }
            .map { meta, k, dir -> tuple( meta.id, meta, k, dir ) }
            .join( ch_mapped_by_patient )
            .map { pid, meta, k, dir, h5ads -> tuple( meta, k, dir, h5ads ) }
        PLOT_SOUPORCELL( ch_soup_plot )
        ch_versions = ch_versions.mix(PLOT_SOUPORCELL.out.versions)
    }

    if (params.run_reference_mapping && params.run_clonetracer) {
        // Per-patient CloneTracer clones + posterior confidence on the joint reference-space UMAP.
        ch_ct_mapped_by_patient = ch_mapped
            .map { meta, h5ad -> tuple( meta.patient, h5ad ) }
            .groupTuple()
        ch_ct_plot = ch_clonetracer
            .map { meta, csv -> tuple( meta.id, meta, csv ) }
            .join( ch_ct_mapped_by_patient )
            .map { pid, meta, csv, h5ads -> tuple( meta, csv, h5ads ) }
        PLOT_CLONETRACER( ch_ct_plot )
        ch_versions = ch_versions.mix(PLOT_CLONETRACER.out.versions)
    }

    if (params.run_qc || params.run_reference_mapping) {
        // Cohort-level QC summary across all samples.
        COHORT_SUMMARY( ch_qc.collect() )
        ch_versions = ch_versions.mix(COHORT_SUMMARY.out.versions)
    }

    //
    // Integration (Phase 2): per-patient master table + headline clonal-tracing Sankeys.
    // Needs the joint clone axes (Numbat + souporcell) plus the reference-mapped phenotype layers.
    //
    ch_cells = Channel.empty()
    if (params.run_integration) {
        if (!(params.run_numbat && params.run_souporcell && params.run_reference_mapping && params.run_copykat)) {
            error "run_integration requires run_numbat, run_souporcell, run_copykat and run_reference_mapping"
        }
        INTEGRATION(
            ch_mapped,
            ch_celltypes,
            ch_copykat,
            ch_numbat,
            ch_souporcell,
            params.souporcell_plot_k
        )
        ch_cells    = INTEGRATION.out.cells
        ch_versions = ch_versions.mix(INTEGRATION.out.versions)
    }

    emit:
    aln         = ch_aln
    cells       = ch_cells
    patient_aln = ch_patient_aln
    numbat      = ch_numbat
    copykat     = ch_copykat
    souporcell  = ch_souporcell
    clonetracer = ch_clonetracer
    qc          = ch_qc
    celltypes   = ch_celltypes
    versions    = ch_versions
}

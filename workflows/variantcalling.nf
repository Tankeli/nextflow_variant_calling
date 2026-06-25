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
include { RNA_DOWNSTREAM  } from '../subworkflows/local/rna_downstream'
include { VISUALIZATION   } from '../subworkflows/local/visualization'
include { INTEGRATION     } from '../subworkflows/local/integration'
include { PROTEOMICS      } from '../subworkflows/local/proteomics'

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
    ch_clonetracer       = Channel.empty()
    ch_clonetracer_trees = Channel.empty()
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
        ch_clonetracer       = CLONETRACER_WF.out.assignments
        ch_clonetracer_trees = CLONETRACER_WF.out.trees
        ch_versions          = ch_versions.mix(CLONETRACER_WF.out.versions)
    }

    //
    // Annotation branch (scanpy QC -> reference mapping), parallel to the callers
    //
    ch_qc = Channel.empty()
    ch_celltypes = Channel.empty()
    ch_mapped = Channel.empty()
    if (params.run_qc || params.run_reference_mapping) {
        if (params.run_reference_mapping && !params.refmap_atlas) {
            error "run_reference_mapping requires refmap_atlas (path to the reference atlas .h5ad); set it in your -params-file"
        }
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
    // the lightweight ANNOTATION branch. Runs off the Cell Ranger filtered matrices; all stages
    // gated by params (default OFF) inside the subworkflow. None of these gate the callers.
    // The same subworkflow backs the standalone DOWNSTREAM entry (main.nf) off published matrices.
    //
    ch_rna_mtx = ch_aln.map { meta, bam, bai, mtx -> tuple( meta, mtx ) }
    RNA_DOWNSTREAM( ch_rna_mtx, CELLRANGER.out.raw )
    ch_versions = ch_versions.mix(RNA_DOWNSTREAM.out.versions)

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
    // Visualisation + cohort reporting — diagnostic figures over the published checkpoints plus the
    // cohort QC summary. Join keys + per-stage gating live in the VISUALIZATION subworkflow.
    //
    VISUALIZATION( ch_mapped, ch_copykat, ch_souporcell, ch_clonetracer, ch_clonetracer_trees, ch_qc )
    ch_versions = ch_versions.mix(VISUALIZATION.out.versions)

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

    //
    // Bulk proteomics branch (DDE_31 port). Starts from a Spectronaut matrix (params), NOT FASTQ, so
    // it runs parallel to everything else. The DESP demix's cell-type proportions are derived from
    // this run's scRNA reference-mapping output (Module-C proteogenomic hook) unless an external
    // proportions TSV is supplied.
    //
    if (params.run_proteomics) {
        if (!params.proteomics_norm || !params.proteomics_design) {
            error "run_proteomics needs --proteomics_norm and --proteomics_design"
        }
        if (params.prot_run_desp && !params.proteomics_proportions && !params.run_reference_mapping) {
            error "prot_run_desp needs --proteomics_proportions or run_reference_mapping (scRNA proportions)"
        }
        def no_file = file("${projectDir}/assets/NO_FILE")
        ch_prot_inputs = Channel.value( tuple(
            params.proteomics_nonnorm ? file(params.proteomics_nonnorm, checkIfExists: true) : file(params.proteomics_norm, checkIfExists: true),
            file(params.proteomics_norm, checkIfExists: true),
            file(params.proteomics_design, checkIfExists: true),
            params.proteomics_contaminants ? file(params.proteomics_contaminants, checkIfExists: true) : no_file
        ) )
        // scRNA cell-type calls -> proportions (only needed for DESP without an external TSV)
        ch_prot_celltypes = (params.prot_run_desp && !params.proteomics_proportions)
            ? ch_celltypes.map { meta, csv -> csv }.collect()
            : Channel.empty()
        PROTEOMICS( ch_prot_inputs, ch_prot_celltypes )
        ch_versions = ch_versions.mix(PROTEOMICS.out.versions)
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

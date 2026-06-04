//
// Top-level workflow. Wires the samplesheet through Cell Ranger and the variant callers.
// Phase 1: input parsing + channel construction.   Phase 2: Cell Ranger.
// Caller subworkflows (Numbat / CopyKAT / souporcell) land in Phases 3-4.
//

include { INPUT_CHECK   } from '../subworkflows/local/input_check'
include { CELLRANGER    } from '../subworkflows/local/cellranger'
include { NUMBAT        } from '../subworkflows/local/numbat'
include { COPYKAT_WF    } from '../subworkflows/local/copykat'
include { SOUPORCELL_WF } from '../subworkflows/local/souporcell'
include { ANNOTATION    } from '../subworkflows/local/annotation'

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
    // Annotation branch (scanpy QC -> reference mapping), parallel to the callers
    //
    ch_qc = Channel.empty()
    ch_celltypes = Channel.empty()
    if (params.run_qc || params.run_reference_mapping) {
        def atlas = params.run_reference_mapping ? file(params.refmap_atlas, checkIfExists: true) : []
        ANNOTATION( ch_aln, params.run_reference_mapping, atlas )
        ch_qc        = ANNOTATION.out.qc
        ch_celltypes = ANNOTATION.out.celltypes
        ch_versions  = ch_versions.mix(ANNOTATION.out.versions)
    }

    emit:
    aln         = ch_aln
    patient_aln = ch_patient_aln
    numbat      = ch_numbat
    copykat     = ch_copykat
    souporcell  = ch_souporcell
    qc          = ch_qc
    celltypes   = ch_celltypes
    versions    = ch_versions
}

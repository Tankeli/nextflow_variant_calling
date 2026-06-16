//
// Protein / ADT branch (MuData). Ported from DDE_27 (notebooks 12-17). Per sample: QC ->
// normalization -> doublet detection -> dimred. Then cohort: batch correction (Harmony) ->
// annotation. Reads the combined Cell Ranger filtered matrix (GEX + Antibody Capture) — enable
// only for CITE-seq samples (GEX-only samples have no 'prot' modality and PROT_QC will error).
//

include { PROT_QC            } from '../../modules/local/prot_qc'
include { PROT_NORMALIZE     } from '../../modules/local/prot_normalize'
include { PROT_DOUBLET       } from '../../modules/local/prot_doublet'
include { PROT_DIMRED        } from '../../modules/local/prot_dimred'
include { PROT_BATCH_CORRECT } from '../../modules/local/prot_batch_correct'
include { PROT_ANNOTATE      } from '../../modules/local/prot_annotate'

workflow PROTEIN {
    take:
    ch_mtx   // [ meta, filtered_feature_bc_matrix ] (Cell Ranger output)

    main:
    ch_versions = Channel.empty()

    PROT_QC( ch_mtx )
    ch_versions = ch_versions.mix(PROT_QC.out.versions)

    PROT_NORMALIZE( PROT_QC.out.h5mu )
    ch_versions = ch_versions.mix(PROT_NORMALIZE.out.versions)

    PROT_DOUBLET( PROT_NORMALIZE.out.h5mu )
    ch_versions = ch_versions.mix(PROT_DOUBLET.out.versions)

    PROT_DIMRED( PROT_DOUBLET.out.h5mu )
    ch_versions = ch_versions.mix(PROT_DIMRED.out.versions)

    // Cohort: collect per-sample dimred objects for batch correction across samples/batches.
    ch_cohort = PROT_DIMRED.out.h5mu
        .map { meta, h5mu -> h5mu }
        .collect()
        .map { h5mus -> tuple( [ id: 'cohort' ], h5mus ) }

    PROT_BATCH_CORRECT( ch_cohort )
    ch_versions = ch_versions.mix(PROT_BATCH_CORRECT.out.versions)

    PROT_ANNOTATE( PROT_BATCH_CORRECT.out.h5mu )
    ch_versions = ch_versions.mix(PROT_ANNOTATE.out.versions)

    emit:
    annotated = PROT_ANNOTATE.out.h5mu
    versions  = ch_versions
}

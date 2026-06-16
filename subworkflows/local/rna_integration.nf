//
// RNA integration (cohort): collect all per-sample annotated objects and integrate across
// samples/batches (scVI / scANVI / BBKNN). Ported from DDE_27 (notebook 07).
//

include { RNA_INTEGRATE } from '../../modules/local/rna_integrate'

workflow RNA_INTEGRATION {
    take:
    ch_annotated   // [ meta, annotated.h5ad ] per sample

    main:
    ch_versions = Channel.empty()

    // Gather every annotated h5ad into one task input (cohort-level).
    ch_cohort = ch_annotated
        .map { meta, h5ad -> h5ad }
        .collect()
        .map { h5ads -> tuple( [ id: 'cohort' ], h5ads ) }

    RNA_INTEGRATE( ch_cohort )
    ch_versions = ch_versions.mix(RNA_INTEGRATE.out.versions)

    emit:
    integrated = RNA_INTEGRATE.out.h5ad
    versions   = ch_versions
}

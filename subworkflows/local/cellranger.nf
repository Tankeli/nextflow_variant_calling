//
// Run Cell Ranger `multi` per sample and emit a normalized alignment channel
// consumed by the variant callers.
//

include { CELLRANGER_MULTI } from '../../modules/local/cellranger_multi'

workflow CELLRANGER {
    take:
    ch_samples       // [ meta, [ {feature_type,fastq_1,fastq_2}, ... ] ]
    cellranger_index // path
    fb_reference     // path

    main:
    ch_versions = Channel.empty()

    // Build per-sample inputs: a flat fastq list for staging + library info for the multi CSV.
    ch_cr_in = ch_samples.map { meta, libraries ->
        def fastqs = libraries.collectMany { [ file(it.fastq_1), file(it.fastq_2) ] }
        def libinfo = libraries.collect { lib ->
            def fastq_id = file(lib.fastq_1).name.replaceAll(/_S\d+_L\d+_R[12]_001\.f(ast)?q\.gz$/, '')
            [ fastq_id: fastq_id, feature_type: lib.feature_type ]
        }
        tuple( meta, libinfo, fastqs )
    }

    CELLRANGER_MULTI( ch_cr_in, cellranger_index, fb_reference )
    ch_versions = ch_versions.mix(CELLRANGER_MULTI.out.versions)

    // [ meta, bam, bai, filtered_feature_bc_matrix ]
    ch_aln = CELLRANGER_MULTI.out.aln

    emit:
    aln      = ch_aln
    raw      = CELLRANGER_MULTI.out.raw   // [ meta, raw_feature_bc_matrix.h5 ] (for SoupX)
    outs     = CELLRANGER_MULTI.out.outs
    versions = ch_versions
}

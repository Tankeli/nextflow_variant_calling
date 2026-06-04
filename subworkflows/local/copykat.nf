//
// CopyKAT (expression-based aneuploid/diploid gate), per sample.
//

include { COPYKAT } from '../../modules/local/copykat'

workflow COPYKAT_WF {
    take:
    ch_aln   // [ meta, bam, bai, matrix_dir ]

    main:
    ch_versions = Channel.empty()

    ch_in = ch_aln.map { meta, bam, bai, mtx -> tuple( meta, mtx ) }

    COPYKAT( ch_in )
    ch_versions = ch_versions.mix(COPYKAT.out.versions)

    emit:
    prediction = COPYKAT.out.prediction
    versions   = ch_versions
}

//
// Annotation branch (parallel to the variant callers): scanpy QC -> reference mapping.
// Runs off the raw cellranger matrices; does NOT filter the callers' inputs.
//

include { SCANPY_QC         } from '../../modules/local/scanpy_qc'
include { REFERENCE_MAPPING } from '../../modules/local/reference_mapping'

workflow ANNOTATION {
    take:
    ch_aln           // [ meta, bam, bai, matrix_dir ]
    run_refmap       // bool
    atlas            // path (only used when run_refmap)
    refmap_umap      // path | [] — frozen reference UMAP for the shared map space

    main:
    ch_versions = Channel.empty()

    SCANPY_QC( ch_aln.map { meta, bam, bai, mtx -> tuple( meta, mtx ) } )
    ch_versions = ch_versions.mix(SCANPY_QC.out.versions)

    ch_celltypes = Channel.empty()
    ch_mapped    = Channel.empty()
    if (run_refmap) {
        REFERENCE_MAPPING( SCANPY_QC.out.h5ad, atlas, refmap_umap )
        ch_celltypes = REFERENCE_MAPPING.out.celltypes
        ch_mapped    = REFERENCE_MAPPING.out.mapped
        ch_versions  = ch_versions.mix(REFERENCE_MAPPING.out.versions)
    }

    emit:
    qc        = SCANPY_QC.out.metrics
    celltypes = ch_celltypes
    mapped    = ch_mapped
    versions  = ch_versions
}

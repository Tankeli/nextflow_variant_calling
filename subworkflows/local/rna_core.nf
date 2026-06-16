//
// RNA core (per sample, chained): QC -> normalization -> feature selection -> dimred ->
// clustering -> annotation. Best-practices scRNA-seq stack ported from DDE_27 (notebooks 01-06);
// each step reads the previous .h5ad. Runs parallel to the variant callers off the Cell Ranger
// filtered matrices — it does NOT gate the callers' inputs.
//

include { RNA_QC             } from '../../modules/local/rna_qc'
include { RNA_NORMALIZE      } from '../../modules/local/rna_normalize'
include { RNA_FEATURE_SELECT } from '../../modules/local/rna_feature_select'
include { RNA_DIMRED         } from '../../modules/local/rna_dimred'
include { RNA_CLUSTER        } from '../../modules/local/rna_cluster'
include { RNA_ANNOTATE       } from '../../modules/local/rna_annotate'

workflow RNA_CORE {
    take:
    ch_mtx       // [ meta, filtered_feature_bc_matrix ] (Cell Ranger output)

    main:
    ch_versions = Channel.empty()

    RNA_QC( ch_mtx )
    ch_versions = ch_versions.mix(RNA_QC.out.versions)

    RNA_NORMALIZE( RNA_QC.out.h5ad )
    ch_versions = ch_versions.mix(RNA_NORMALIZE.out.versions)

    RNA_FEATURE_SELECT( RNA_NORMALIZE.out.h5ad )
    ch_versions = ch_versions.mix(RNA_FEATURE_SELECT.out.versions)

    RNA_DIMRED( RNA_FEATURE_SELECT.out.h5ad )
    ch_versions = ch_versions.mix(RNA_DIMRED.out.versions)

    RNA_CLUSTER( RNA_DIMRED.out.h5ad )
    ch_versions = ch_versions.mix(RNA_CLUSTER.out.versions)

    RNA_ANNOTATE( RNA_CLUSTER.out.h5ad )
    ch_versions = ch_versions.mix(RNA_ANNOTATE.out.versions)

    emit:
    annotated = RNA_ANNOTATE.out.h5ad   // [ meta, rna_06_annotated.h5ad ]
    celltypes = RNA_ANNOTATE.out.celltypes
    qc        = RNA_QC.out.metrics
    versions  = ch_versions
}

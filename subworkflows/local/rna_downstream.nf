//
// RNA downstream best-practices stack (ported from DDE_27), shared by the main VARIANTCALLING
// workflow (off Cell Ranger output, parallel to the callers) and the standalone DOWNSTREAM entry
// (off already-published filtered matrices). All stages gated by params; none gate the callers.
//   RNA_CORE        QC -> normalize -> feature-select -> dimred -> cluster -> annotate (per sample)
//   RNA_INTEGRATION cohort integration (scVI / scANVI / BBKNN + scib)            [needs run_rna_core]
//   RNA_ADVANCED    pseudotime / velocity / DE / composition (gated individually)[needs run_rna_core]
//   PROTEIN         surface-protein / ADT branch (CITE-seq only)
//

include { RNA_CORE        } from './rna_core'
include { RNA_INTEGRATION } from './rna_integration'
include { RNA_ADVANCED    } from './rna_advanced'
include { PROTEIN         } from './protein'

workflow RNA_DOWNSTREAM {
    take:
    ch_mtx   // [ meta, filtered_feature_bc_matrix ]
    ch_raw   // [ meta, raw_feature_bc_matrix.h5 ] — for SoupX; Channel.empty() when unavailable

    main:
    ch_versions   = Channel.empty()
    ch_annotated  = Channel.empty()
    ch_integrated = Channel.empty()

    if (params.run_rna_core) {
        // RNA_QC input: [ meta, filtered_mtx, raw|[] ]. SoupX (run_soupx) needs the per-sample raw
        // matrix; otherwise pass [] so nothing extra is staged.
        ch_qc_in = params.run_soupx
            ? ch_mtx
                .map { meta, mtx -> tuple( meta.id, meta, mtx ) }
                .join( ch_raw.map { meta, raw -> tuple( meta.id, raw ) } )
                .map { id, meta, mtx, raw -> tuple( meta, mtx, raw ) }
            : ch_mtx.map { meta, mtx -> tuple( meta, mtx, [] ) }
        RNA_CORE( ch_qc_in )
        ch_annotated = RNA_CORE.out.annotated
        ch_versions  = ch_versions.mix(RNA_CORE.out.versions)
    }

    if (params.run_rna_integration) {
        if (!params.run_rna_core) {
            error "run_rna_integration requires run_rna_core (it consumes the annotated objects)"
        }
        RNA_INTEGRATION( ch_annotated )
        ch_integrated = RNA_INTEGRATION.out.integrated
        ch_versions   = ch_versions.mix(RNA_INTEGRATION.out.versions)
    }

    if (params.run_pseudotime || params.run_velocity || params.run_de || params.run_composition) {
        if (!params.run_rna_core) {
            error "RNA advanced stages (pseudotime/velocity/de/composition) require run_rna_core"
        }
        if (params.run_composition && !params.run_rna_integration) {
            error "run_composition requires run_rna_integration (it operates on the integrated object)"
        }
        // Velocity looms (optional): ${params.velocity_loom_dir}/<sample>.loom; samples without one
        // are dropped by the join in RNA_ADVANCED.
        ch_loom = Channel.empty()
        if (params.run_velocity) {
            if (!params.velocity_loom_dir) {
                error "run_velocity requires velocity_loom_dir (this pipeline does not generate velocyto looms)"
            }
            ch_loom = ch_mtx
                .map { meta, mtx -> tuple( meta.id, file("${params.velocity_loom_dir}/${meta.id}.loom") ) }
                .filter { id, loom -> loom.exists() }
        }
        RNA_ADVANCED(
            ch_annotated,
            ch_loom,
            ch_integrated,
            params.run_pseudotime,
            params.run_velocity,
            params.run_de,
            params.run_composition
        )
        ch_versions = ch_versions.mix(RNA_ADVANCED.out.versions)
    }

    if (params.run_protein) {
        PROTEIN( ch_mtx )
        ch_versions = ch_versions.mix(PROTEIN.out.versions)
    }

    emit:
    annotated  = ch_annotated
    integrated = ch_integrated
    versions   = ch_versions
}

//
// Visualisation + cohort reporting over the published checkpoints. Diagnostic figures that read an
// existing caller output and overlay it on the shared reference-map UMAP (so they need the mapped
// h5ads — gated on run_reference_mapping), plus the cohort-level QC summary. Pulled out of the top
// workflow so the per-sample/per-patient join keys live in one testable place.
//

include { PLOT_COPYKAT     } from '../../modules/local/plot_copykat'
include { PLOT_SOUPORCELL  } from '../../modules/local/plot_souporcell'
include { PLOT_CLONETRACER } from '../../modules/local/plot_clonetracer'
include { COHORT_SUMMARY   } from '../../modules/local/cohort_summary'

workflow VISUALIZATION {
    take:
    ch_mapped       // [ meta, mapped.h5ad ]            reference-mapped per sample (may be empty)
    ch_copykat      // [ meta, prediction ]             per sample
    ch_souporcell   // [ meta, k, dir ]                 per patient, per K
    ch_clonetracer  // [ meta, assignments.csv ]        per patient
    ch_ct_trees     // [ meta, out.pickle, tree.pickle ] per patient (CloneTracer model pickles)
    ch_qc           // [ meta, qc_metrics.csv ]         per sample

    main:
    ch_versions = Channel.empty()

    // Reference-mapped h5ads grouped by patient — reused by the souporcell + clonetracer overlays,
    // which are per-patient (joint Dx+Rel) and so need every sample's frame for that patient.
    ch_mapped_by_patient = ch_mapped
        .map { meta, h5ad -> tuple( meta.patient, h5ad ) }
        .groupTuple()

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
        ch_soup_plot = ch_souporcell
            .filter { meta, k, dir -> k == params.souporcell_plot_k }
            .map { meta, k, dir -> tuple( meta.id, meta, k, dir ) }
            .join( ch_mapped_by_patient )
            .map { pid, meta, k, dir, h5ads -> tuple( meta, k, dir, h5ads ) }
        PLOT_SOUPORCELL( ch_soup_plot )
        ch_versions = ch_versions.mix(PLOT_SOUPORCELL.out.versions)
    }

    if (params.run_reference_mapping && params.run_clonetracer) {
        // Per-patient CloneTracer: trees/ELBO/heatmap from the model pickle + clone overlay on the
        // joint reference-space UMAP. Join assignments + out.pickle + mapped h5ads by patient id.
        ch_ct_plot = ch_clonetracer
            .map { meta, csv -> tuple( meta.id, meta, csv ) }
            .join( ch_ct_trees.map { meta, outp, treep -> tuple( meta.id, outp ) } )
            .join( ch_mapped_by_patient )
            .map { pid, meta, csv, outp, h5ads -> tuple( meta, csv, outp, h5ads ) }
        PLOT_CLONETRACER( ch_ct_plot )
        ch_versions = ch_versions.mix(PLOT_CLONETRACER.out.versions)
    }

    if (params.run_qc || params.run_reference_mapping) {
        // Cohort-level QC summary across all samples.
        COHORT_SUMMARY( ch_qc.collect() )
        ch_versions = ch_versions.mix(COHORT_SUMMARY.out.versions)
    }

    emit:
    versions = ch_versions
}

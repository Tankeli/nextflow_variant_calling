//
// RNA advanced (gated, fan-out off the annotated objects + integrated cohort object). Ported
// from DDE_27 (notebooks 08-11):
//   - pseudotime    (per sample, notebook 08)
//   - RNA velocity  (per sample, gated on a velocyto loom supplied via params.velocity_loom_dir)
//   - DE            (per patient: Dx vs Rel, notebook 10)
//   - composition   (cohort, on the integrated object, notebook 11)
//
// Velocity note: this pipeline does not produce velocyto looms. When run_velocity is enabled,
// looms are looked up at ${params.velocity_loom_dir}/<sample>.loom; samples without one are
// skipped (mirrors DDE_27's "only samples that supplied a loom").
//

include { RNA_PSEUDOTIME  } from '../../modules/local/rna_pseudotime'
include { RNA_VELOCITY    } from '../../modules/local/rna_velocity'
include { RNA_DE          } from '../../modules/local/rna_de'
include { RNA_COMPOSITION } from '../../modules/local/rna_composition'

workflow RNA_ADVANCED {
    take:
    ch_annotated     // [ meta, annotated.h5ad ] per sample
    ch_loom          // [ id, loom ] per sample with a loom | empty
    ch_integrated    // [ meta, integrated.h5ad ] (cohort) | empty
    run_pseudotime
    run_velocity
    run_de
    run_composition

    main:
    ch_versions = Channel.empty()

    if (run_pseudotime) {
        RNA_PSEUDOTIME( ch_annotated )
        ch_versions = ch_versions.mix(RNA_PSEUDOTIME.out.versions)
    }

    if (run_velocity) {
        // Join annotated objects to their loom by sample id; only samples that supplied a loom.
        ch_vel = ch_annotated
            .map { meta, h5ad -> tuple( meta.id, meta, h5ad ) }
            .join( ch_loom )
            .map { id, meta, h5ad, loom -> tuple( meta, h5ad, loom ) }
        RNA_VELOCITY( ch_vel )
        ch_versions = ch_versions.mix(RNA_VELOCITY.out.versions)
    }

    if (run_de) {
        // Group annotated objects per patient (Dx + Rel) for paired DE.
        ch_patient = ch_annotated
            .map { meta, h5ad -> tuple( meta.patient, meta, h5ad ) }
            .groupTuple()
            .map { patient, metas, h5ads ->
                tuple( [ id: patient, patient: patient ], h5ads ) }
        RNA_DE( ch_patient )
        ch_versions = ch_versions.mix(RNA_DE.out.versions)
    }

    if (run_composition) {
        RNA_COMPOSITION( ch_integrated )
        ch_versions = ch_versions.mix(RNA_COMPOSITION.out.versions)
    }

    emit:
    versions = ch_versions
}

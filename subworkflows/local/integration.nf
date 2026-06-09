//
// Integration branch (Phase 2): LSC scoring -> Phase-0 master table -> headline Sankeys.
// Per patient, joins the joint Numbat + souporcell clones with CopyKAT / cell type / LSC.
//

include { LSC_SCORING        } from '../../modules/local/lsc_scoring'
include { PHASE0_INTEGRATION } from '../../modules/local/phase0_integration'
include { HEADLINE_FIGURES   } from '../../modules/local/headline_figures'

workflow INTEGRATION {
    take:
    ch_mapped      // [ meta, mapped_h5ad ]      per sample
    ch_celltypes   // [ meta, celltypes.csv ]    per sample
    ch_copykat     // [ meta, prediction.txt ]   per sample
    ch_numbat      // [ meta{id:patient}, numbat_out ] per patient
    ch_souporcell  // [ meta{id:patient}, k, clusters_dir ] per patient/K
    plot_k         // int — souporcell K to use

    main:
    ch_versions = Channel.empty()

    // Per-sample LSC scores.
    LSC_SCORING( ch_mapped )
    ch_versions = ch_versions.mix(LSC_SCORING.out.versions)

    // Per-patient meta (sample + timepoint lists, Dx first) derived from the mapped samples.
    ch_pmeta = ch_mapped
        .map { meta, h5ad -> tuple( meta.patient, [ meta.id, meta.timepoint ] ) }
        .groupTuple()
        .map { patient, pairs ->
            def ordered = pairs.sort { it[1] == 'Dx' ? 0 : 1 }
            tuple( patient, [ id: patient,
                              samples:    ordered.collect { it[0] },
                              timepoints: ordered.collect { it[1] } ] )
        }

    // Group per-sample phenotype layers by patient.
    ch_ct  = ch_celltypes.map      { m, f -> tuple( m.patient, f ) }.groupTuple()
    ch_ck  = ch_copykat.map        { m, f -> tuple( m.patient, f ) }.groupTuple()
    ch_lsc = LSC_SCORING.out.lsc.map { m, f -> tuple( m.patient, f ) }.groupTuple()

    // Per-patient clone axes.
    ch_nb = ch_numbat.map { m, d -> tuple( m.id, d ) }
    ch_sp = ch_souporcell.filter { m, k, d -> k == plot_k }.map { m, k, d -> tuple( m.id, d ) }

    // Join everything on patient (inner — a patient needs all modalities to integrate).
    ch_p0 = ch_pmeta
        .join( ch_nb ).join( ch_sp ).join( ch_ct ).join( ch_ck ).join( ch_lsc )
        .map { patient, pmeta, nb, sp, ct, ck, lsc -> tuple( pmeta, nb, sp, ct, ck, lsc ) }

    PHASE0_INTEGRATION( ch_p0 )
    ch_versions = ch_versions.mix(PHASE0_INTEGRATION.out.versions)

    // Headline figures need the master table + the Numbat output (CNV labels).
    ch_headline = PHASE0_INTEGRATION.out.cells
        .map { meta, cells -> tuple( meta.id, meta, cells ) }
        .join( ch_nb )
        .map { id, meta, cells, nb -> tuple( meta, cells, nb ) }

    HEADLINE_FIGURES( ch_headline )
    ch_versions = ch_versions.mix(HEADLINE_FIGURES.out.versions)

    emit:
    cells    = PHASE0_INTEGRATION.out.cells
    versions = ch_versions
}

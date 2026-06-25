//
// PROTEOMICS — bulk mass-spec branch (ported from DDE_31, R -> Python; DESP kept as one R step).
// Starts from a Spectronaut search matrix (NOT FASTQ), so it runs as an optional gated branch and a
// standalone `-entry PROTEOMICS`. Cell-type proportions for the DESP demix are derived from THIS
// pipeline's scRNA reference-mapping output (Module-C proteogenomic hook) or an external TSV
// (--proteomics_proportions).
//
//   PREP (1a) -> BATCH (1b) -> { DE (3), VIZ (2), STAGE4 (4), ML (5) }
//   PROPORTIONS (scRNA cell-type fractions) -> DESP (6 run, R) -> DESP_VIZ (6 viz)
//
// PREP+BATCH are the spine (always run when run_proteomics); other stages gated by prot_run_*.
// All single-item handles are converted to value channels (.first()) so they fan out to multiple
// downstream processes.
//

include { PROT_MS_PREP        } from '../../modules/local/prot_ms_prep'
include { PROT_MS_BATCH       } from '../../modules/local/prot_ms_batch'
include { PROT_MS_DE          } from '../../modules/local/prot_ms_de'
include { PROT_MS_VIZ         } from '../../modules/local/prot_ms_viz'
include { PROT_MS_STAGE4      } from '../../modules/local/prot_ms_stage4'
include { PROT_MS_ML          } from '../../modules/local/prot_ms_ml'
include { PROT_MS_PROPORTIONS } from '../../modules/local/prot_ms_proportions'
include { PROT_MS_DESP        } from '../../modules/local/prot_ms_desp'
include { PROT_MS_DESP_VIZ    } from '../../modules/local/prot_ms_desp_viz'

workflow PROTEOMICS {
    take:
    ch_inputs      // value: tuple( nonnorm, norm, design, contaminants )
    ch_celltypes   // channel emitting ONE collected list of <sample>_celltypes.csv, or Channel.empty()

    main:
    ch_versions = Channel.empty()
    def no_file = file("${projectDir}/assets/NO_FILE")

    // --- spine: QC/prep -> batch correction ---
    PROT_MS_PREP( ch_inputs )
    ch_versions = ch_versions.mix(PROT_MS_PREP.out.versions)
    ch_matrix = PROT_MS_PREP.out.matrix.first()
    ch_design = PROT_MS_PREP.out.design.first()

    PROT_MS_BATCH( ch_matrix.combine(ch_design) )
    ch_versions    = ch_versions.mix(PROT_MS_BATCH.out.versions)
    ch_raw         = PROT_MS_BATCH.out.raw.first()
    ch_design_corr = PROT_MS_BATCH.out.design.first()
    ch_analysis    = ( (params.prot_desp_bulk_source == 'raw')    ? PROT_MS_BATCH.out.raw
                     : (params.prot_desp_bulk_source == 'combat') ? PROT_MS_BATCH.out.combat
                     :                                              PROT_MS_BATCH.out.limma ).first()
    ch_norm = ch_inputs.map { nonnorm, norm, design, cont -> norm }.first()

    // --- differential expression (on the filtered_log2 matrix, as in 3a) ---
    ch_de = Channel.empty()
    if (params.prot_run_de) {
        PROT_MS_DE( ch_matrix.combine(ch_design) )
        ch_de       = PROT_MS_DE.out.de.first()
        ch_versions = ch_versions.mix(PROT_MS_DE.out.versions)
    }

    // --- visualisation ---
    if (params.prot_run_viz) {
        PROT_MS_VIZ( ch_raw.combine(ch_analysis).combine(ch_design_corr) )
        ch_versions = ch_versions.mix(PROT_MS_VIZ.out.versions)
    }

    // --- stage-4 interpretation (needs DE) ---
    if (params.prot_run_stage4) {
        if (!params.prot_run_de) { error "prot_run_stage4 requires prot_run_de" }
        PROT_MS_STAGE4( ch_analysis.combine(ch_design_corr).combine(ch_de).combine(ch_norm) )
        ch_versions = ch_versions.mix(PROT_MS_STAGE4.out.versions)
    }

    // --- ML classifiers ---
    if (params.prot_run_ml) {
        PROT_MS_ML( ch_analysis.combine(ch_design_corr) )
        ch_versions = ch_versions.mix(PROT_MS_ML.out.versions)
    }

    // --- DESP cell-state demixing (needs proportions + DE) ---
    if (params.prot_run_desp) {
        if (!params.prot_run_de) { error "prot_run_desp requires prot_run_de (DE ranks the figures)" }

        def ch_props
        if (params.proteomics_proportions) {
            ch_props = Channel.value(file(params.proteomics_proportions, checkIfExists: true))
        } else {
            def smap = params.proteomics_sample_map ? file(params.proteomics_sample_map) : no_file
            PROT_MS_PROPORTIONS( ch_celltypes, Channel.value(smap) )
            ch_props    = PROT_MS_PROPORTIONS.out.proportions.first()
            ch_versions = ch_versions.mix(PROT_MS_PROPORTIONS.out.versions)
        }

        PROT_MS_DESP( ch_analysis.combine(ch_props).combine(ch_design_corr) )
        ch_versions = ch_versions.mix(PROT_MS_DESP.out.versions)

        PROT_MS_DESP_VIZ(
            ch_analysis.combine(ch_de).combine(ch_props).combine(ch_design_corr)
                .combine(PROT_MS_DESP.out.desp_dir.first())
        )
        ch_versions = ch_versions.mix(PROT_MS_DESP_VIZ.out.versions)
    }

    emit:
    matrix   = ch_matrix
    de       = ch_de
    versions = ch_versions
}

//
// CopyKAT robustness sweep (separate analysis track, NOT a caller checkpoint).
// Fans out CopyKAT over a parameter x seed grid per sample so the downstream Python
// (stability / drivers / cross-ref / celltype matrix) can quantify where the aneuploid/diploid
// boundary lies and how stable the per-cell call is. Gated by params.run_copykat_robustness.
//

include { COPYKAT_NORM_BARCODES } from '../../modules/local/copykat_norm_barcodes'
include { COPYKAT_SWEEP         } from '../../modules/local/copykat_sweep'

workflow COPYKAT_ROBUSTNESS_WF {
    take:
    ch_mtx        // [ meta, matrix_dir ]    per sample
    ch_celltypes  // [ meta, celltypes.csv ] per sample (empty when reference mapping is off)

    main:
    ch_versions = Channel.empty()

    // Cross-product of the sweep dimensions -> one combo map per CopyKAT run.
    def combos = []
    for (ks in params.copykat_robustness_ks_cut)
      for (win in params.copykat_robustness_win_size)
        for (ng in params.copykat_robustness_ngene_chr)
          for (dist in params.copykat_robustness_distance)
            for (nr in params.copykat_robustness_use_norm_ref)
              for (sd in params.copykat_robustness_seeds) {
                def normFlag = (nr as boolean) ? 1 : 0
                combos << [ ks_cut: ks, win_size: win, ngene_chr: ng, distance: dist,
                            norm: normFlag, seed: sd,
                            id: "ks${ks}_win${win}_ng${ng}_${dist}_norm${normFlag}_seed${sd}" ]
              }
    log.info "COPYKAT_ROBUSTNESS: ${combos.size()} CopyKAT runs per sample"

    def needNorm    = params.copykat_robustness_use_norm_ref.any { it as boolean }
    def placeholder = file("${projectDir}/assets/no_norm_barcodes.txt", checkIfExists: true)

    // [ id, meta(+combo), matrix_dir ] — one row per (sample x combo)
    ch_combo = ch_mtx
        .combine(Channel.fromList(combos))
        .map { meta, mtx, combo -> tuple(meta.id, meta + [combo: combo], mtx) }

    if (needNorm) {
        // Per-sample known-normal baseline from the reference-mapped cell types.
        COPYKAT_NORM_BARCODES(ch_celltypes, params.copykat_norm_celltypes)
        ch_versions = ch_versions.mix(COPYKAT_NORM_BARCODES.out.versions)
        ch_norm = COPYKAT_NORM_BARCODES.out.barcodes.map { meta, f -> tuple(meta.id, f) }

        // combine(by:0) fans the per-sample norm file across every combo of that sample.
        ch_sweep_in = ch_combo
            .combine(ch_norm, by: 0)
            .map { id, meta, mtx, nf -> tuple(meta, mtx, nf) }
    } else {
        ch_sweep_in = ch_combo.map { id, meta, mtx -> tuple(meta, mtx, placeholder) }
    }

    COPYKAT_SWEEP(ch_sweep_in)
    ch_versions = ch_versions.mix(COPYKAT_SWEEP.out.versions)

    emit:
    sweep    = COPYKAT_SWEEP.out.prediction   // [ meta(+combo), prediction.txt, CNA_results.txt ]
    versions = ch_versions
}

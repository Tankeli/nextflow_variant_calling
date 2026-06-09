//
// CloneTracer (downstream clonal integration), joint per patient.
// Synthesises per-cell M/N over CNV (Numbat) + nuclear-SNV (souporcell) + mtDNA (new pileup),
// builds the per-patient JSON, then runs the pyro inference model.
//
//   per-sample MTDNA_PILEUP ─┐
//   Numbat numbat_out ───────┼─> CLONETRACER_BUILD (<patient>.json) ─> CLONETRACER
//   souporcell k<K> ─────────┘
//
// Numbat / souporcell inputs are optional (remainder joins) so a patient is still processed
// from whatever mutation sources are available (mtDNA alone is enough).

include { MTDNA_PILEUP      } from '../../modules/local/mtdna_pileup'
include { CLONETRACER_BUILD } from '../../modules/local/clonetracer_build'
include { CLONETRACER       } from '../../modules/local/clonetracer'

workflow CLONETRACER_WF {
    take:
    ch_aln          // [ meta(sample), bam, bai, matrix_dir ]   (per sample)
    ch_patient_aln  // [ patient, [ {meta,bam,bai,mtx} sorted Dx->Rel ] ]
    ch_numbat       // [ meta(patient), numbat_out ]   (may be empty)
    ch_souporcell   // [ meta(patient), k, k<K> ]       (may be empty)
    gtf             // path or []

    main:
    ch_versions = Channel.empty()

    // mtDNA pileup per sample, then grouped by patient.
    MTDNA_PILEUP( ch_aln )
    ch_versions = ch_versions.mix(MTDNA_PILEUP.out.versions)

    ch_mt_by_patient = MTDNA_PILEUP.out.mtdna
        .map { meta, dir -> tuple( meta.patient, [ id: meta.id, dir: dir ] ) }
        .groupTuple()

    // Per-patient base: ordered meta + matrices + mtDNA dirs (mtDNA ordered to match samples).
    ch_base = ch_patient_aln
        .join( ch_mt_by_patient )
        .map { patient, members, mtlist ->
            def meta = [ id        : patient,
                         samples   : members.collect { it.meta.id },
                         timepoints: members.collect { it.meta.timepoint } ]
            def mtx   = members.collect { it.mtx }
            def mtmap = mtlist.collectEntries { [ (it.id): it.dir ] }
            def mt    = members.collect { mtmap[it.meta.id] }
            tuple( patient, meta, mtx, mt )
        }

    // Optional caller inputs keyed by patient.
    ch_numbat_k = ch_numbat.map { meta, dir -> tuple( meta.id, dir ) }
    ch_soup_k   = ch_souporcell
        .filter { meta, k, dir -> k.toString() == params.clonetracer_k.toString() }
        .map { meta, k, dir -> tuple( meta.id, dir ) }

    ch_build_in = ch_base
        .join( ch_numbat_k, remainder: true )
        .join( ch_soup_k,   remainder: true )
        .filter { it[1] != null }   // drop any right-only remainder artefacts (no base meta)
        .map { patient, meta, mtx, mt, numbat, soup ->
            tuple( meta, mtx, mt, numbat ?: [], soup ?: [], gtf )
        }

    CLONETRACER_BUILD( ch_build_in )
    ch_versions = ch_versions.mix(CLONETRACER_BUILD.out.versions)

    CLONETRACER( CLONETRACER_BUILD.out.json )
    ch_versions = ch_versions.mix(CLONETRACER.out.versions)

    emit:
    assignments = CLONETRACER.out.assignments
    trees       = CLONETRACER.out.trees
    json        = CLONETRACER_BUILD.out.json
    versions    = ch_versions
}

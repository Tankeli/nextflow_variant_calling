//
// Numbat (allele-based CNV / clones), joint per patient.
// pileup_and_phase (SNP pileup + phasing) -> run_numbat (clone calling).
//

include { NUMBAT_PILEUP } from '../../modules/local/numbat_pileup'
include { NUMBAT_RUN    } from '../../modules/local/numbat_run'

workflow NUMBAT {
    take:
    ch_patient_aln   // [ patient, [ {meta,bam,bai,mtx} sorted Dx->Rel ] ]

    main:
    ch_versions = Channel.empty()

    // Pileup input: patient meta carries the ordered sample list; bams/bais/matrices as lists.
    ch_pileup_in = ch_patient_aln.map { patient, members ->
        def meta = [ id: patient, samples: members.collect { it.meta.id } ]
        tuple( meta,
               members.collect { it.bam },
               members.collect { it.bai },
               members.collect { it.mtx } )
    }

    NUMBAT_PILEUP( ch_pileup_in )
    ch_versions = ch_versions.mix(NUMBAT_PILEUP.out.versions)

    // run_numbat input: join the patient's allele counts with its 10X matrices.
    ch_mtx = ch_patient_aln.map { patient, members -> tuple( patient, members.collect { it.mtx } ) }

    ch_run_in = NUMBAT_PILEUP.out.allele
        .map { meta, allele -> tuple( meta.id, meta, allele ) }
        .join( ch_mtx )
        .map { id, meta, allele, mtxs -> tuple( meta, allele, mtxs ) }

    NUMBAT_RUN( ch_run_in )
    ch_versions = ch_versions.mix(NUMBAT_RUN.out.versions)

    emit:
    allele   = NUMBAT_PILEUP.out.allele
    numbat   = NUMBAT_RUN.out.numbat
    versions = ch_versions
}

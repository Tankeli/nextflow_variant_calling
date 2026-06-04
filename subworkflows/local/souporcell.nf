//
// Souporcell (SNV genotype clusters), joint per patient over a K sweep.
//

include { SOUPORCELL_PREP } from '../../modules/local/souporcell_prep'
include { SOUPORCELL      } from '../../modules/local/souporcell'

workflow SOUPORCELL_WF {
    take:
    ch_patient_aln   // [ patient, [ {meta,bam,bai,mtx} sorted Dx->Rel ] ]
    fasta            // path

    main:
    ch_versions = Channel.empty()

    ch_prep_in = ch_patient_aln.map { patient, members ->
        def meta = [ id: patient, samples: members.collect { it.meta.id } ]
        tuple( meta,
               members.collect { it.bam },
               members.collect { it.bai },
               members.collect { it.mtx } )
    }

    SOUPORCELL_PREP( ch_prep_in )
    ch_versions = ch_versions.mix(SOUPORCELL_PREP.out.versions)

    // Fan out the merged BAM over the K sweep.
    ch_k = Channel.fromList( params.souporcell_k.toString().split(',').collect { it.trim() as int } )
    ch_soup_in = SOUPORCELL_PREP.out.merged.combine( ch_k )

    SOUPORCELL( ch_soup_in, fasta )
    ch_versions = ch_versions.mix(SOUPORCELL.out.versions)

    emit:
    clusters = SOUPORCELL.out.clusters
    versions = ch_versions
}

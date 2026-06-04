//
// Parse + validate the samplesheet and build per-sample / per-patient channels.
//

include { samplesheetToList } from 'plugin/nf-schema'

workflow INPUT_CHECK {
    take:
    samplesheet   // path to samplesheet CSV

    main:
    // Each row -> [ meta(id,patient,timepoint), fastq_1, fastq_2, feature_type, expected_cells ]
    ch_rows = Channel.fromList(
        samplesheetToList(samplesheet, "${projectDir}/assets/schema_input.json")
    )

    // Group the per-library rows (gex + ab) into one entry per sample.
    ch_samples = ch_rows
        .map { meta, fastq_1, fastq_2, feature_type, expected_cells ->
            def key = meta + [ expected_cells: expected_cells ]
            tuple( key, [ feature_type: feature_type, fastq_1: fastq_1, fastq_2: fastq_2 ] )
        }
        .groupTuple()
        .map { meta, libraries -> tuple( meta, libraries ) }

    emit:
    samples = ch_samples   // [ meta(id,patient,timepoint,expected_cells), [ {feature_type,fastq_1,fastq_2}, ... ] ]
}

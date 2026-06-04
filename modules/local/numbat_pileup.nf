// Numbat step 1: joint per-patient SNP pileup + phasing (cellsnp-lite + Eagle2).
// Pools all of a patient's samples (Dx + Rel) into one pileup_and_phase.R call so the
// downstream run_numbat clone IDs are comparable across timepoints. SNP VCF / phasing
// panel / genetic map default to paths baked into the numbat container.

process NUMBAT_PILEUP {
    tag "$meta.id"
    label 'process_high'

    container params.numbat_container

    input:
    tuple val(meta), path(bams), path(bais), path(mtx_dirs)

    output:
    tuple val(meta), path("${meta.id}_pileup/*_allele_counts.tsv.gz"), emit: allele
    path "${meta.id}_pileup"                                         , emit: pileup
    path "versions.yml"                                              , emit: versions

    script:
    def samples    = meta.samples
    def sample_csv = samples.join(',')
    def bam_csv    = samples.collect { "${it}.bam" }.join(',')
    def bc_csv     = samples.collect { "${it}.barcodes.tsv" }.join(',')
    def decompress = samples.collect { "zcat ${it}_filtered_feature_bc_matrix/barcodes.tsv.gz > ${it}.barcodes.tsv" }.join('\n    ')
    """
    mkdir -p ${meta.id}_pileup
    ${decompress}

    Rscript /numbat/inst/bin/pileup_and_phase.R \\
        --label ${meta.id} \\
        --samples ${sample_csv} \\
        --bams ${bam_csv} \\
        --barcodes ${bc_csv} \\
        --outdir ${meta.id}_pileup \\
        --gmap ${params.numbat_gmap} \\
        --snpvcf ${params.numbat_snpvcf} \\
        --paneldir ${params.numbat_paneldir} \\
        --ncores ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numbat: \$(Rscript -e 'cat(as.character(packageVersion("numbat")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    def mk = meta.samples.collect { "echo | gzip > ${meta.id}_pileup/${it}_allele_counts.tsv.gz" }.join('\n    ')
    """
    mkdir -p ${meta.id}_pileup
    ${mk}
    echo '"${task.process}": {numbat: stub}' > versions.yml
    """
}

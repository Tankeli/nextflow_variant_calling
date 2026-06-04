// Souporcell prep: joint per-patient BAM assembly.
// Retags CB/CR with "<sample>__" so cells stay unique across Dx + Rel, merges + sorts +
// indexes the per-sample BAMs, and builds the combined (prefixed) barcode list.
// (No NK filtering — runs on full BAMs.)

process SOUPORCELL_PREP {
    tag "$meta.id"
    label 'process_medium'

    container params.samtools_container

    input:
    tuple val(meta), path(bams), path(bais), path(mtx_dirs)

    output:
    tuple val(meta), path("${meta.id}.merged.sorted.bam"),
                     path("${meta.id}.merged.sorted.bam.bai"),
                     path("${meta.id}.barcodes.tsv"), emit: merged
    path "versions.yml"                              , emit: versions

    script:
    def samples    = meta.samples
    def retag      = samples.collect { "souporcell_retag.sh ${it} ${it}.bam ${it}.retag.bam" }.join('\n    ')
    def make_bc    = samples.collect { "zcat ${it}_filtered_feature_bc_matrix/barcodes.tsv.gz | sed 's/^/${it}__/' > ${it}.bc.txt" }.join('\n    ')
    def retag_bams = samples.collect { "${it}.retag.bam" }.join(' ')
    def bc_files   = samples.collect { "${it}.bc.txt" }.join(' ')
    """
    ${retag}
    ${make_bc}

    samtools merge -@ ${task.cpus} -f ${meta.id}.merged.bam ${retag_bams}
    samtools sort  -@ ${task.cpus} -o ${meta.id}.merged.sorted.bam ${meta.id}.merged.bam
    samtools index -@ ${task.cpus} ${meta.id}.merged.sorted.bam
    cat ${bc_files} > ${meta.id}.barcodes.tsv
    rm -f ${meta.id}.merged.bam *.retag.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.sorted.bam ${meta.id}.merged.sorted.bam.bai
    printf '%s\\n' ${meta.samples.collect { "${it}__AAAA-1" }.join(' ')} > ${meta.id}.barcodes.tsv
    echo '"${task.process}": {samtools: stub}' > versions.yml
    """
}

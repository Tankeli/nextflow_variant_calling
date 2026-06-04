// Souporcell: SNV genotype clustering on the joint per-patient merged BAM, one run per K.
// --no_umi/--skip_remap/--ignore mirror the DDE_24 paired-run invocation.
// Checkpoint: <patient>/k<K>/clusters.tsv.

process SOUPORCELL {
    tag "${meta.id}_k${k}"
    label 'process_high'

    container params.souporcell_container

    input:
    tuple val(meta), path(bam), path(bai), path(barcodes), val(k)
    path fasta

    output:
    tuple val(meta), val(k), path("${meta.id}/k${k}"), emit: clusters
    path "versions.yml"                              , emit: versions

    script:
    """
    [ -f ${fasta}.fai ] || samtools faidx ${fasta}
    mkdir -p ${meta.id}/k${k}

    souporcell_pipeline.py \\
        -i ${bam} \\
        -b ${barcodes} \\
        -f ${fasta} \\
        -t ${task.cpus} \\
        -o ${meta.id}/k${k} \\
        -k ${k} \\
        --no_umi true \\
        --skip_remap True \\
        --ignore True

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: \$(souporcell_pipeline.py --version 2>&1 | tail -n1 || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/k${k}
    printf 'barcode\\tstatus\\tassignment\\n' > ${meta.id}/k${k}/clusters.tsv
    echo '"${task.process}": {souporcell: stub}' > versions.yml
    """
}

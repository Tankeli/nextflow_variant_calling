// Per-sample LSC scoring (weighted pLSC6 + LSC17) on the reference-mapped h5ad.
// Feeds the Phase-0 master table and the headline pLSC6 figures.

process LSC_SCORING {
    tag "$meta.id"
    label 'process_low'

    container params.scanpy_container

    input:
    tuple val(meta), path(mapped_h5ad)

    output:
    tuple val(meta), path("${meta.id}_lsc.csv"), emit: lsc
    path "versions.yml"                         , emit: versions

    script:
    """
    lsc_scoring.py ${meta.id} ${mapped_h5ad}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'barcode,sample_id,pLSC6_score,LSC17_score\\nAAAA-1,${meta.id},0.5,0.1\\n' > ${meta.id}_lsc.csv
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

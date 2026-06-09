// Per-sample scanpy QC (CITE-seq): QC metrics + Scrublet doublets + filtering.
// Parallel annotation branch — does NOT gate the variant callers. Checkpoint:
// <sample>_qc.h5ad + <sample>_qc_metrics.csv.

process SCANPY_QC {
    tag "$meta.id"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("${meta.id}_qc.h5ad"), emit: h5ad
    path "${meta.id}_qc_metrics.csv"           , emit: metrics
    path "${meta.id}_qc_panel.{png,pdf}"       , emit: panel
    path "versions.yml"                        , emit: versions

    script:
    """
    qc.py ${matrix_dir} ${meta.id} \\
        ${params.qc_min_genes} ${params.qc_min_umi} ${params.qc_max_mito_pct} ${params.qc_expected_doublet_rate}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'barcode,sample_id,n_genes,pass_qc\\nAAAA-1,${meta.id},1500,True\\n' > ${meta.id}_qc_metrics.csv
    echo stub > ${meta.id}_qc.h5ad
    echo stub > ${meta.id}_qc_panel.png
    echo stub > ${meta.id}_qc_panel.pdf
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

// RNA integration (notebook 07): multi-sample batch integration (scVI / scANVI / BBKNN / Seurat)
// over all per-sample annotated objects. Checkpoint: rna_07_integrated.h5ad.

process RNA_INTEGRATE {
    tag "$meta.id"
    label 'process_high'

    container params.rna_integration_container

    input:
    tuple val(meta), path(h5ads, stageAs: 'input?/*')

    output:
    tuple val(meta), path("rna_07_integrated.h5ad"), emit: h5ad
    path "integration_metrics.csv"                 , emit: metrics, optional: true
    path "versions.yml"                            , emit: versions

    script:
    """
    rna_integrate.py \\
        --inputs . \\
        --methods ${params.integration_methods} \\
        --batch_key ${params.integration_batch_key} \\
        --out rna_07_integrated.h5ad \\
        --metrics integration_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvi-tools: \$(python -c 'import scvi; print(scvi.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_07_integrated.h5ad
    printf 'method,metric,value\\nscvi,kBET,0.5\\n' > integration_metrics.csv
    echo '"${task.process}": {scvi-tools: stub}' > versions.yml
    """
}

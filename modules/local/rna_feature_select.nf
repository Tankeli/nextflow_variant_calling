// RNA feature selection (notebook 03): deviance-based highly variable gene selection (scry).
// Checkpoint: rna_03_feature_selected.h5ad.

process RNA_FEATURE_SELECT {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_preprocessing_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_03_feature_selected.h5ad"), emit: h5ad
    path "versions.yml"                                  , emit: versions

    script:
    """
    rna_feature_select.py \\
        --in ${h5ad} --sample ${meta.id} \\
        --n_top_genes ${params.n_top_genes} \\
        --out rna_03_feature_selected.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_03_feature_selected.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

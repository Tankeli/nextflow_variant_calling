// RNA dimensionality reduction (notebook 04): PCA, t-SNE, UMAP.
// Checkpoint: rna_04_dimred.h5ad.

process RNA_DIMRED {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_clustering_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_04_dimred.h5ad"), emit: h5ad
    path "versions.yml"                        , emit: versions

    script:
    """
    rna_dimred.py \\
        --in ${h5ad} --sample ${meta.id} \\
        --n_pcs ${params.n_pcs} --n_neighbors ${params.n_neighbors} \\
        --out rna_04_dimred.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_04_dimred.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

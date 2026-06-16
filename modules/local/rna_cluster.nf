// RNA clustering (notebook 05): Leiden clustering at multiple resolutions.
// Checkpoint: rna_05_clustered.h5ad.

process RNA_CLUSTER {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_clustering_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_05_clustered.h5ad"), emit: h5ad
    path "versions.yml"                           , emit: versions

    script:
    """
    rna_cluster.py \\
        --in ${h5ad} --sample ${meta.id} \\
        --resolutions ${params.leiden_resolutions} \\
        --out rna_05_clustered.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_05_clustered.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

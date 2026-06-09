// Per-sample CopyKAT UMAP overlay: aneuploid/diploid call + cell type on the
// reference-mapping UMAP. Supplements the native CopyKAT genome heatmap.

process PLOT_COPYKAT {
    tag "$meta.id"
    label 'process_low'

    container params.scanpy_container

    input:
    tuple val(meta), path(mapped_h5ad), path(prediction)

    output:
    path "${meta.id}_copykat_umap.{png,pdf}", emit: umap
    path "versions.yml"                      , emit: versions

    script:
    """
    plot_copykat_umap.py ${meta.id} ${mapped_h5ad} ${prediction}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > ${meta.id}_copykat_umap.png
    echo stub > ${meta.id}_copykat_umap.pdf
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

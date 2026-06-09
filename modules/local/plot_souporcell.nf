// Per-patient souporcell clone UMAP (Dx/Rel overlay) + clone x timepoint composition.
// Recomputes a joint UMAP from the patient's mapped h5ads; clones come from clusters.tsv.

process PLOT_SOUPORCELL {
    tag "${meta.id}_k${k}"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), val(k), path(clusters_dir), path(mapped_h5ads)

    output:
    path "${meta.id}_souporcell_*.{png,pdf}", emit: figures
    path "versions.yml"                      , emit: versions

    script:
    def samples    = meta.samples.join(',')
    def timepoints = meta.timepoints.join(',')
    """
    clusters=\$(find -L ${clusters_dir} -name clusters.tsv | head -n1)
    plot_souporcell_umap.py ${meta.id} \$clusters ${k} \\
        ${samples} ${timepoints} ${mapped_h5ads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > ${meta.id}_souporcell_umap.png
    echo stub > ${meta.id}_souporcell_umap.pdf
    echo stub > ${meta.id}_souporcell_composition.png
    echo stub > ${meta.id}_souporcell_composition.pdf
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

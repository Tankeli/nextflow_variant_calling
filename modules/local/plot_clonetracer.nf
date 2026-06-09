// Per-patient CloneTracer clone UMAP (Dx/Rel overlay + posterior confidence) + clone x timepoint
// composition, on the SHARED reference-map space (X_umap_ref from the mapped h5ads). Clones come
// from CloneTracer's per-cell <patient>_clone_assignments.csv.

process PLOT_CLONETRACER {
    tag "$meta.id"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), path(assignments), path(mapped_h5ads)

    output:
    path "${meta.id}_clonetracer_*.{png,pdf}", emit: figures
    path "versions.yml"                       , emit: versions

    script:
    def samples    = meta.samples.join(',')
    def timepoints = meta.timepoints.join(',')
    """
    plot_clonetracer_umap.py ${meta.id} ${assignments} \\
        ${samples} ${timepoints} ${mapped_h5ads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > ${meta.id}_clonetracer_umap.png
    echo stub > ${meta.id}_clonetracer_umap.pdf
    echo stub > ${meta.id}_clonetracer_composition.png
    echo stub > ${meta.id}_clonetracer_composition.pdf
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

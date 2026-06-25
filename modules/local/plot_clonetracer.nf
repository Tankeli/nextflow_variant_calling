// Per-patient CloneTracer figures. Two complementary sets, both in the aml_scrna env (viking profile):
//   * clonetracer_figures.py (pickle only): clonal-hierarchy trees, ELBO/model-evidence, and the
//     single-cell VAF + clone-posterior heatmap — ports of the veltenlab clonal_inference R vignette.
//   * plot_clonetracer_umap.py (assignments + mapped h5ads): clones + posterior confidence on the
//     SHARED reference-map UMAP (X_umap_ref) and a clone x timepoint composition.

process PLOT_CLONETRACER {
    tag "$meta.id"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), path(assignments), path(out_pickle), path(mapped_h5ads)

    output:
    path "${meta.id}_*.{png,pdf}", emit: figures
    path "versions.yml"          , emit: versions

    script:
    def samples    = meta.samples.join(',')
    def timepoints = meta.timepoints.join(',')
    """
    # Trees / ELBO / VAF+posterior heatmap from the model pickle.
    clonetracer_figures.py ${out_pickle} ${meta.id} .

    # Clone overlay + composition on the shared reference-map UMAP.
    plot_clonetracer_umap.py ${meta.id} ${assignments} \\
        ${samples} ${timepoints} ${mapped_h5ads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    for s in trees elbo heatmap clonetracer_umap clonetracer_composition; do
        echo stub > ${meta.id}_\${s}.png
        echo stub > ${meta.id}_\${s}.pdf
    done
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

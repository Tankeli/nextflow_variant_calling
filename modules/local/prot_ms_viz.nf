// Bulk proteomics 2 — visualisation (ported from DDE_31 2a/2b -> prot_ms_viz.py).
// Clustered heatmaps + PCA (+ loadings) + sample-level UMAP/DPT off the raw & corrected matrices.

process PROT_MS_VIZ {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    // stageAs distinct names: with prot_desp_bulk_source=raw the corrected matrix IS matrix_raw,
    // so raw and corrected would collide on staging.
    tuple path(raw, stageAs: 'matrix_raw.tsv'), path(corrected, stageAs: 'matrix_analysis.tsv'), path(design)

    output:
    path "*.png"        , emit: figures, optional: true
    path "*.csv"        , emit: tables,  optional: true
    path "versions.yml" , emit: versions

    script:
    def cfg = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    def corr = corrected.name != 'NO_FILE' ? "--corrected ${corrected}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_viz.py \\
        --raw ${raw} ${corr} --design ${design} --method ${params.prot_batch_method} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > HM_noncorrected.png
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

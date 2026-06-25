// Bulk proteomics 4 — detailed interpretation (ported from DDE_31 4a-4f -> prot_ms_stage4.py).
// Sig-protein heatmap+clusters, PCA, condition-mean k-means, DPT pseudotime assoc, boxplots, offsets.

process PROT_MS_STAGE4 {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    tuple path(analysis), path(design), path(de), path(norm)

    output:
    path "*.png"        , emit: figures, optional: true
    path "*.csv"        , emit: tables,  optional: true
    path "boxplots/*"   , emit: boxplots, optional: true
    path "versions.yml" , emit: versions

    script:
    def cfg  = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    def norm_arg = norm.name != 'NO_FILE' ? "--norm ${norm}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_stage4.py \\
        --analysis ${analysis} --design ${design} --de ${de} ${norm_arg} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scikit-learn: \$(python -c 'import sklearn; print(sklearn.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > HM_clusters_stub.png
    echo '"${task.process}": {sklearn: stub}' > versions.yml
    """
}

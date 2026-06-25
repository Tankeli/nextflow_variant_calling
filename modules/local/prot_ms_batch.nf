// Bulk proteomics 1b — batch correction (ported from DDE_31 1b -> prot_ms_batch.py).
// limma removeBatchEffect (reimplemented) + ComBat (inmoose) -> matrix_{limma,combat,raw} + corrected design.

process PROT_MS_BATCH {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    tuple path(matrix), path(design)

    output:
    path "matrix_raw.tsv"        , emit: raw
    path "matrix_limma.tsv"      , emit: limma,  optional: true
    path "matrix_combat.tsv"     , emit: combat, optional: true
    path "design_corrected.tsv"  , emit: design
    path "versions.yml"          , emit: versions

    script:
    def cfg = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_batch.py \\
        --matrix ${matrix} --design ${design} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    cp ${matrix} matrix_raw.tsv
    cp ${matrix} matrix_limma.tsv
    cp ${matrix} matrix_combat.tsv
    cp ${design} design_corrected.tsv
    echo '"${task.process}": {pandas: stub}' > versions.yml
    """
}

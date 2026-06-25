// Bulk proteomics 5 — ML classifiers (ported from DDE_31 5a/5b -> prot_ms_ml.py).
// Decision tree + random forest (permutation importance) classifying condition from abundance.

process PROT_MS_ML {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    tuple path(matrix), path(design)

    output:
    path "*.csv"        , emit: tables,  optional: true
    path "*.png"        , emit: figures, optional: true
    path "*.txt"        , emit: rules,   optional: true
    path "versions.yml" , emit: versions

    script:
    def cfg = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_ml.py \\
        --matrix ${matrix} --design ${design} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scikit-learn: \$(python -c 'import sklearn; print(sklearn.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'Protein,Gene,Importance\\nP1,GENE1,0.5\\n' > random_forest_importance.csv
    echo '"${task.process}": {sklearn: stub}' > versions.yml
    """
}

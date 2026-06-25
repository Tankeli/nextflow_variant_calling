// Bulk proteomics 1a — QC & preprocessing (ported from DDE_31 1a.QC_and_prep.R -> prot_ms_prep.py).
// Filters/log2/imputes the Spectronaut matrices -> filtered_log2 matrix + reduced design + QC figs.

process PROT_MS_PREP {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    // stageAs distinct names so nonnorm==norm (fallback) doesn't collide on staging
    tuple path(nonnorm, stageAs: 'nonnorm.tsv'), path(norm, stageAs: 'norm.tsv'), path(design, stageAs: 'design.tsv'), path(contaminants)

    output:
    path "filtered_log2_imputed.tsv", emit: matrix
    path "design_reduced.tsv"       , emit: design
    path "*.png"                    , emit: figures, optional: true
    path "filtering_summary.tsv"    , emit: summary,  optional: true
    path "versions.yml"             , emit: versions

    script:
    def cfg  = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    def cont = contaminants.name != 'NO_FILE' ? "--contaminants ${contaminants}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_prep.py \\
        --nonnorm ${nonnorm} --norm ${norm} --design ${design} ${cont} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'Protein.Ids\\tGenes\\tA\\tB\\nP1\\tGENE1\\t1.0\\t2.0\\n' > filtered_log2_imputed.tsv
    printf 'sample\\tcondition\\treplicate\\tBatch\\nA\\tDiagnosis\\tAML1\\t1\\n' > design_reduced.tsv
    printf 'step\\tn_proteins\\n' > filtering_summary.tsv
    echo '"${task.process}": {pandas: stub}' > versions.yml
    """
}

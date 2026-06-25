// Bulk proteomics 3 — differential expression + volcano (ported from DDE_31 3a/3b -> prot_ms_de.py).
// limma lmFit/eBayes reimplemented (squeezeVar moderation). Uses the filtered_log2 matrix (as in 3a).

process PROT_MS_DE {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    tuple path(matrix), path(design)

    output:
    path "de_results.csv"            , emit: de
    path "de_results_per_patient.csv", emit: de_patient, optional: true
    path "volcano/*"                 , emit: volcano,    optional: true
    path "versions.yml"              , emit: versions

    script:
    def cfg = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_de.py \\
        --matrix ${matrix} --design ${design} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scipy: \$(python -c 'import scipy; print(scipy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'protein,logFC [Relapse vs Diagnosis],P.Value [Relapse vs Diagnosis],adj.P.Val [Relapse vs Diagnosis],Genes\\nP1,1.2,0.01,0.04,GENE1\\n' > de_results.csv
    mkdir -p volcano && echo stub > volcano/Volcano_all_Relapse_vs_Diagnosis.png
    echo '"${task.process}": {scipy: stub}' > versions.yml
    """
}

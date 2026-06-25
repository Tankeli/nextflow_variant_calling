// Bulk proteomics 6 (viz) — DESP visualisations (ported from DDE_31 6a/6c/6e -> prot_ms_desp_viz.py).
// Top-protein heatmap + direction bar, DE-ranked delta heatmap + cell-type contribution, per-patient.

process PROT_MS_DESP_VIZ {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_container

    input:
    tuple path(bulk), path(de), path(proportions), path(design), path(desp_dir)

    output:
    path "*.png"           , emit: figures, optional: true
    path "*.csv"           , emit: tables,  optional: true
    path "per_patient/**"  , emit: per_patient, optional: true
    path "versions.yml"    , emit: versions

    script:
    def cfg = params.proteomics_config ? "--config ${file(params.proteomics_config)}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_desp_viz.py \\
        --bulk ${bulk} --de ${de} --proportions ${proportions} --design ${design} \\
        --desp_dir ${desp_dir} \\
        --default_config ${projectDir}/assets/proteomics_default.yaml ${cfg} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seaborn: \$(python -c 'import seaborn; print(seaborn.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > desp_delta_heatmap.png
    echo '"${task.process}": {seaborn: stub}' > versions.yml
    """
}

// Bulk proteomics 6 (run) — DESP cell-state demixing. The one R step (per design choice): wraps the
// external DESP package. prot_ms_desp_run.R. Outputs land under desp/ for the viz step.

process PROT_MS_DESP {
    tag "proteomics"
    label 'process_medium'

    container params.proteomics_r_container

    input:
    tuple path(bulk), path(proportions), path(design)

    output:
    path "desp"        , emit: desp_dir
    path "versions.yml", emit: versions

    script:
    def cond = params.prot_condition_col ?: 'condition'
    def rep  = params.prot_replicate_col ?: 'replicate'
    def idc  = params.prot_id_col ?: 'Protein.Ids'
    def pp   = params.prot_desp_per_patient ? 'TRUE' : 'FALSE'
    """
    prot_ms_desp_run.R \\
        --bulk ${bulk} --proportions ${proportions} --design ${design} \\
        --id_col ${idc} --condition_col ${cond} --replicate_col ${rep} \\
        --per_patient ${pp} --outdir desp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        DESP: \$(Rscript -e 'cat(as.character(packageVersion("DESP")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p desp/per_patient/AML152
    printf 'feature\\tHSC\\tGMP\\nP1\\t1.0\\t2.0\\n' > desp/desp_diagnosis_cell_state_profiles.tsv
    printf 'feature\\tHSC\\tGMP\\nP1\\t1.5\\t2.5\\n' > desp/desp_relapse_cell_state_profiles.tsv
    printf 'feature\\tHSC\\tGMP\\nP1\\t0.5\\t0.5\\n' > desp/desp_delta_matrix.tsv
    printf 'status=completed\\n' > desp/desp_run_summary.txt
    echo '"${task.process}": {DESP: stub}' > versions.yml
    """
}

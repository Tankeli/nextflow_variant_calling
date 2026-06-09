// Pipeline-level QC summary across all samples. Collects every per-sample
// *_qc_metrics.csv into one panel (cells, gene complexity, %MT, doublet rate).

process COHORT_SUMMARY {
    label 'process_low'

    container params.scanpy_container

    input:
    path qc_metrics

    output:
    path "cohort_summary.{png,pdf}", emit: figures
    path "cohort_summary.csv"      , emit: table
    path "versions.yml"            , emit: versions

    script:
    """
    cohort_summary.py ${qc_metrics}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python -c 'import platform; print(platform.python_version())' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > cohort_summary.png
    echo stub > cohort_summary.pdf
    printf 'sample_id,n_cells_total\\nstub,1\\n' > cohort_summary.csv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

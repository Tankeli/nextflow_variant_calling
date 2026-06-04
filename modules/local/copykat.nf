// CopyKAT: per-sample expression-based aneuploid/diploid malignancy gate.
// Runs directly on the cellranger filtered matrix (Gene Expression). Checkpoint:
// <sample>_copykat_prediction.txt.

process COPYKAT {
    tag "$meta.id"
    label 'process_medium'

    container params.copykat_container

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("${meta.id}_copykat_prediction.txt"), emit: prediction
    path "${meta.id}_copykat_*"                               , emit: results
    path "versions.yml"                                       , emit: versions

    script:
    """
    copykat.R ${matrix_dir} ${meta.id} ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copykat: \$(Rscript -e 'cat(as.character(packageVersion("copykat")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'cell.names\\tcopykat.pred\\n' > ${meta.id}_copykat_prediction.txt
    touch ${meta.id}_copykat_CNA_results.txt ${meta.id}_copykat_clustering_results.rds
    echo '"${task.process}": {copykat: stub}' > versions.yml
    """
}

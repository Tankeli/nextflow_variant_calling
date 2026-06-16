// Protein/ADT batch correction (notebook 16): merge per-sample objects and run Harmony across
// batches (cohort). Checkpoint: prot_05_batch_corrected.h5mu.

process PROT_BATCH_CORRECT {
    tag "$meta.id"
    label 'process_high'

    container params.protein_container

    input:
    tuple val(meta), path(h5mus, stageAs: 'input?/*')

    output:
    tuple val(meta), path("prot_05_batch_corrected.h5mu"), emit: h5mu
    path "versions.yml"                                  , emit: versions

    script:
    """
    prot_batch_correct.py --inputs . --batch_key ${params.integration_batch_key} \\
        --out prot_05_batch_corrected.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        harmonypy: \$(python -c 'import harmonypy; print(harmonypy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_05_batch_corrected.h5mu
    echo '"${task.process}": {harmonypy: stub}' > versions.yml
    """
}

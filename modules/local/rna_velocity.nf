// RNA velocity (notebook 09): scVelo over a velocyto loom + the annotated object (per sample).
// Checkpoint: rna_09_velocity.h5ad.

process RNA_VELOCITY {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_velocity_container

    input:
    tuple val(meta), path(h5ad), path(loom)

    output:
    tuple val(meta), path("rna_09_velocity.h5ad"), emit: h5ad
    path "versions.yml"                          , emit: versions

    script:
    """
    rna_velocity.py --in ${h5ad} --loom ${loom} --sample ${meta.id} --out rna_09_velocity.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scvelo: \$(python -c 'import scvelo; print(scvelo.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_09_velocity.h5ad
    echo '"${task.process}": {scvelo: stub}' > versions.yml
    """
}

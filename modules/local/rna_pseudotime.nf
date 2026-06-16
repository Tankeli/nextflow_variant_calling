// RNA pseudotime (notebook 08): diffusion pseudotime trajectory analysis (per sample).
// Checkpoint: rna_08_pseudotime.h5ad.

process RNA_PSEUDOTIME {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_pseudotime_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_08_pseudotime.h5ad"), emit: h5ad
    path "versions.yml"                            , emit: versions

    script:
    """
    rna_pseudotime.py --in ${h5ad} --sample ${meta.id} --out rna_08_pseudotime.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_08_pseudotime.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

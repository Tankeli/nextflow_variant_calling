// RNA normalization (notebook 02): scran size-factor normalization (+ log1p).
// Checkpoint: rna_02_normalized.h5ad.

process RNA_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_preprocessing_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_02_normalized.h5ad"), emit: h5ad
    path "versions.yml"                            , emit: versions

    script:
    """
    rna_normalize.py --in ${h5ad} --sample ${meta.id} --out rna_02_normalized.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_02_normalized.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

// Protein/ADT dimensionality reduction (notebook 15): PCA + UMAP embeddings.
// Checkpoint: prot_04_dimred.h5mu.

process PROT_DIMRED {
    tag "$meta.id"
    label 'process_medium'

    container params.protein_container

    input:
    tuple val(meta), path(h5mu)

    output:
    tuple val(meta), path("prot_04_dimred.h5mu"), emit: h5mu
    path "versions.yml"                         , emit: versions

    script:
    """
    prot_dimred.py --in ${h5mu} --sample ${meta.id} --out prot_04_dimred.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muon: \$(python -c 'import muon; print(muon.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_04_dimred.h5mu
    echo '"${task.process}": {muon: stub}' > versions.yml
    """
}

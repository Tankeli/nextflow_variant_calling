// Protein/ADT normalization (notebook 13): DSB or CLR normalization of protein counts.
// Checkpoint: prot_02_normalized.h5mu.

process PROT_NORMALIZE {
    tag "$meta.id"
    label 'process_medium'

    container params.protein_container

    input:
    tuple val(meta), path(h5mu)

    output:
    tuple val(meta), path("prot_02_normalized.h5mu"), emit: h5mu
    path "versions.yml"                             , emit: versions

    script:
    """
    prot_normalize.py --in ${h5mu} --sample ${meta.id} \\
        --method ${params.prot_norm_method} --out prot_02_normalized.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muon: \$(python -c 'import muon; print(muon.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_02_normalized.h5mu
    echo '"${task.process}": {muon: stub}' > versions.yml
    """
}

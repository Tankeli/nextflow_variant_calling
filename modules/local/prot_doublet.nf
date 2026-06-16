// Protein/ADT doublet detection (notebook 14): marker-based doublet detection + removal.
// Checkpoint: prot_03_doublet_filtered.h5mu.

process PROT_DOUBLET {
    tag "$meta.id"
    label 'process_medium'

    container params.protein_container

    input:
    tuple val(meta), path(h5mu)

    output:
    tuple val(meta), path("prot_03_doublet_filtered.h5mu"), emit: h5mu
    path "versions.yml"                                   , emit: versions

    script:
    """
    prot_doublet.py --in ${h5mu} --sample ${meta.id} --out prot_03_doublet_filtered.h5mu

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muon: \$(python -c 'import muon; print(muon.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_03_doublet_filtered.h5mu
    echo '"${task.process}": {muon: stub}' > versions.yml
    """
}

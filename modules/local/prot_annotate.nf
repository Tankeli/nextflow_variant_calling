// Protein/ADT annotation (notebook 17): Leiden clustering + cell-type annotation on the
// batch-corrected cohort object. Checkpoint: prot_06_annotated.h5mu.

process PROT_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    container params.protein_container

    input:
    tuple val(meta), path(h5mu)

    output:
    tuple val(meta), path("prot_06_annotated.h5mu"), emit: h5mu
    path "prot_celltypes.csv"                      , emit: celltypes
    path "versions.yml"                            , emit: versions

    script:
    """
    prot_annotate.py --in ${h5mu} \\
        --resolution ${params.prot_leiden_resolution} \\
        --out prot_06_annotated.h5mu --celltypes prot_celltypes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muon: \$(python -c 'import muon; print(muon.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_06_annotated.h5mu
    printf 'barcode,cell_type\\nAAAA-1,CD4 T\\n' > prot_celltypes.csv
    echo '"${task.process}": {muon: stub}' > versions.yml
    """
}

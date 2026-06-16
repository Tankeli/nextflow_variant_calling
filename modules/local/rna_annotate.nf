// RNA annotation (notebook 06): celltypist + marker-based cell-type annotation.
// Checkpoint: rna_06_annotated.h5ad.

process RNA_ANNOTATE {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_annotation_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("rna_06_annotated.h5ad"), emit: h5ad
    path "${meta.id}_celltypes.csv"               , emit: celltypes
    path "versions.yml"                           , emit: versions

    script:
    """
    rna_annotate.py \\
        --in ${h5ad} --sample ${meta.id} \\
        --celltypist_model ${params.celltypist_model} \\
        --resolution ${params.leiden_key_resolution} \\
        --out rna_06_annotated.h5ad \\
        --celltypes ${meta.id}_celltypes.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        celltypist: \$(python -c 'import celltypist; print(celltypist.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_06_annotated.h5ad
    printf 'barcode,cell_type\\nAAAA-1,T cell\\n' > ${meta.id}_celltypes.csv
    echo '"${task.process}": {celltypist: stub}' > versions.yml
    """
}

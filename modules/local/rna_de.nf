// RNA differential expression (notebook 10): MAST / edgeR / pertpy, per patient (Dx vs Rel)
// over the patient's annotated objects. Checkpoint: DE_markers.csv.

process RNA_DE {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_de_container

    input:
    tuple val(meta), path(h5ads, stageAs: 'input?/*')

    output:
    tuple val(meta), path("DE_markers.csv"), emit: markers
    path "versions.yml"                    , emit: versions

    script:
    """
    rna_de.py \\
        --inputs . \\
        --patient ${meta.id} \\
        --method ${params.de_method} \\
        --groupby ${params.de_groupby} \\
        --out DE_markers.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        method: ${params.de_method}
    END_VERSIONS
    """

    stub:
    """
    printf 'cell_type,gene,logFC,pval,padj\\nT cell,CD3D,2.1,0.001,0.01\\n' > DE_markers.csv
    echo '"${task.process}": {method: ${params.de_method}}' > versions.yml
    """
}

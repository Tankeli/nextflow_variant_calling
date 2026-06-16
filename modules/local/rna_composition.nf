// RNA compositional analysis (notebook 11): cell-type composition shifts (scCODA / pertpy)
// over the integrated cohort object. Checkpoint: composition_results.csv.

process RNA_COMPOSITION {
    tag "$meta.id"
    label 'process_medium'

    container params.composition_container

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("composition_results.csv"), emit: results
    path "versions.yml"                             , emit: versions

    script:
    """
    rna_composition.py --in ${h5ad} --out composition_results.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pertpy: \$(python -c 'import pertpy; print(pertpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'cell_type,effect,hdi_low,hdi_high\\nT cell,-0.3,-0.6,0.0\\n' > composition_results.csv
    echo '"${task.process}": {pertpy: stub}' > versions.yml
    """
}

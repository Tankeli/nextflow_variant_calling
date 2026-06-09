// Headline Dx->Rel clonal-tracing figures from the Phase-0 master table:
// Numbat-clone Sankey (CNV-labelled, pLSC6-coloured) + pLSC6-quartile fingerprint Sankey.

process HEADLINE_FIGURES {
    tag "$meta.id"
    label 'process_low'

    container params.scanpy_container

    input:
    tuple val(meta), path(cells_tsv), path(numbat_out)

    output:
    path "${meta.id}_fig*.{pdf,png}", emit: figures
    path "versions.yml"             , emit: versions

    script:
    """
    headline_sankey.py ${meta.id} ${cells_tsv} ${numbat_out}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plotly: \$(python -c 'import plotly; print(plotly.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > ${meta.id}_fig1_numbat_sankey.pdf
    echo stub > ${meta.id}_fig1_numbat_sankey.png
    echo stub > ${meta.id}_fig2_pLSC6_fingerprint.pdf
    echo stub > ${meta.id}_fig2_pLSC6_fingerprint.png
    echo '"${task.process}": {plotly: stub}' > versions.yml
    """
}

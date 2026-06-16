// Protein/ADT QC (notebook 12): load the combined Cell Ranger filtered_feature_bc_matrix/
// directory as MuData (GEX + Antibody Capture; prot_qc.py also accepts a combined .h5), QC
// metrics + outlier filtering on the protein modality. Checkpoint: prot_01_quality_control.h5mu.

process PROT_QC {
    tag "$meta.id"
    label 'process_medium'

    container params.protein_container

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("prot_01_quality_control.h5mu"), emit: h5mu
    path "${meta.id}_prot_qc_metrics.csv"                , emit: metrics
    path "versions.yml"                                  , emit: versions

    script:
    """
    prot_qc.py \\
        --filtered_h5 ${matrix_dir} --sample ${meta.id} \\
        --patient ${meta.patient ?: 'NA'} --timepoint ${meta.timepoint ?: 'NA'} --batch ${meta.batch ?: 'NA'} \\
        --out prot_01_quality_control.h5mu \\
        --metrics ${meta.id}_prot_qc_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        muon: \$(python -c 'import muon; print(muon.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > prot_01_quality_control.h5mu
    printf 'barcode,sample_id,total_protein\\nAAAA-1,${meta.id},1200\\n' > ${meta.id}_prot_qc_metrics.csv
    echo '"${task.process}": {muon: stub}' > versions.yml
    """
}

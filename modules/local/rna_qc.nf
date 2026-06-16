// RNA QC (notebook 01): QC metrics + MAD outlier filtering + scDblFinder doublet scoring +
// gene filtering. Reads the Cell Ranger filtered_feature_bc_matrix/ directory emitted by this
// pipeline's CELLRANGER (rna_qc.py also accepts a combined filtered .h5). SoupX is skipped:
// `cellranger multi` does not emit a per-sample raw matrix here. Checkpoint: rna_01_quality_control.h5ad.

process RNA_QC {
    tag "$meta.id"
    label 'process_medium'

    container params.rna_preprocessing_container

    input:
    tuple val(meta), path(matrix_dir)

    output:
    tuple val(meta), path("rna_01_quality_control.h5ad"), emit: h5ad
    path "${meta.id}_qc_metrics.csv"                     , emit: metrics
    path "versions.yml"                                  , emit: versions

    script:
    def scdbl   = params.run_scdblfinder ? "--scdblfinder" : ""
    """
    rna_qc.py \\
        --filtered_h5 ${matrix_dir} \\
        --sample ${meta.id} \\
        --patient ${meta.patient ?: 'NA'} --timepoint ${meta.timepoint ?: 'NA'} --batch ${meta.batch ?: 'NA'} \\
        --nmads_counts ${params.rna_qc_nmads_counts} \\
        --nmads_mt ${params.rna_qc_nmads_mt} \\
        --max_mito_pct ${params.rna_qc_max_mito_pct} \\
        --min_cells ${params.rna_qc_min_cells} \\
        ${scdbl} \\
        --out rna_01_quality_control.h5ad \\
        --metrics ${meta.id}_qc_metrics.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    echo stub > rna_01_quality_control.h5ad
    printf 'barcode,sample_id,pct_counts_mt,scDblFinder_class\\nAAAA-1,${meta.id},2.1,singlet\\n' > ${meta.id}_qc_metrics.csv
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

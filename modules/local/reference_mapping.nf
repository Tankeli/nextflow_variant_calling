// Per-sample reference mapping via scanpy.tl.ingest onto a configurable atlas
// (default: DDE_32 paediatric BM; Zeng / others swappable). Checkpoint:
// <sample>_celltypes.csv + <sample>_mapped.h5ad.

process REFERENCE_MAPPING {
    tag "$meta.id"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), path(qc_h5ad)
    path atlas

    output:
    tuple val(meta), path("${meta.id}_celltypes.csv"), emit: celltypes
    path "${meta.id}_mapped.h5ad"                    , emit: mapped
    path "versions.yml"                              , emit: versions

    script:
    """
    reference_mapping.py ${qc_h5ad} ${atlas} ${meta.id} \\
        ${params.refmap_celltype_key} ${params.refmap_confidence_threshold} ${params.refmap_n_pcs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'barcode,sample_id,ref_cell_type,mapping_confidence,poorly_mapped\\nAAAA-1,${meta.id},HSC,0.9,False\\n' > ${meta.id}_celltypes.csv
    echo stub > ${meta.id}_mapped.h5ad
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

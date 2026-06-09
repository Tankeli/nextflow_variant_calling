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
    path refmap_umap

    output:
    tuple val(meta), path("${meta.id}_celltypes.csv"), emit: celltypes
    tuple val(meta), path("${meta.id}_mapped.h5ad")  , emit: mapped
    path "${meta.id}_mapping_umap.{png,pdf}"         , emit: umap
    path "versions.yml"                              , emit: versions

    script:
    def umap_arg = refmap_umap ? refmap_umap : 'NONE'
    """
    # v2: reference now log-normalised to match the query before PCA/ingest (raw-count atlas fix)
    reference_mapping.py ${qc_h5ad} ${atlas} ${meta.id} \\
        ${params.refmap_celltype_key} ${params.refmap_confidence_threshold} ${params.refmap_n_pcs} \\
        ${umap_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c 'import scanpy; print(scanpy.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'barcode,sample_id,ref_cell_type,mapping_confidence,poorly_mapped\\nAAAA-1,${meta.id},HSC,0.9,False\\n' > ${meta.id}_celltypes.csv
    echo stub > ${meta.id}_mapped.h5ad
    echo stub > ${meta.id}_mapping_umap.png
    echo stub > ${meta.id}_mapping_umap.pdf
    echo '"${task.process}": {scanpy: stub}' > versions.yml
    """
}

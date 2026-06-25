// Bulk proteomics — cell-type proportions for DESP, derived from THIS pipeline's scRNA
// reference-mapping output (Module-C proteogenomic hook). prot_ms_proportions.py.

process PROT_MS_PROPORTIONS {
    tag "proteomics"
    label 'process_low'

    container params.proteomics_container

    input:
    path celltypes      // collected list of <sample>_celltypes.csv from REFERENCE_MAPPING
    path sample_map     // TSV rna_sample<TAB>prot_sample, or NO_FILE

    output:
    path "celltype_proportions.tsv", emit: proportions
    path "versions.yml"            , emit: versions

    script:
    def map_arg = sample_map.name != 'NO_FILE' ? "--map ${sample_map}" : ""
    def ct_col  = params.prot_celltype_col ? "--celltype_col ${params.prot_celltype_col}" : ""
    """
    export PYTHONPATH=${projectDir}/bin:\${PYTHONPATH:-}
    prot_ms_proportions.py \\
        --celltypes ${celltypes} ${map_arg} ${ct_col} \\
        --out celltype_proportions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'cell_type\\t2977\\t109\\nHSC\\t0.2\\t0.1\\nGMP\\t0.3\\t0.4\\n' > celltype_proportions.tsv
    echo '"${task.process}": {pandas: stub}' > versions.yml
    """
}

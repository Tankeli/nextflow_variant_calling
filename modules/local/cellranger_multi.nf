// Cell Ranger `multi` for CITE-seq (Gene Expression + Antibody Capture).
// Generates the multi-config CSV from the samplesheet libraries, runs cellranger,
// then normalizes the output to the canonical <sample>/outs/ layout the variant
// callers expect: possorted_genome_bam.bam (+ .bai) and filtered_feature_bc_matrix/.

process CELLRANGER_MULTI {
    tag "$meta.id"
    label 'process_high'

    container params.cellranger_container

    input:
    tuple val(meta), val(libraries), path(fastqs, stageAs: 'fastqs/*')
    path cellranger_index
    path fb_reference   // optional: pass [] for GEX-only (no Antibody Capture)

    output:
    // Uniquely sample-named handles so joint callers can stage multiple samples without collision.
    tuple val(meta), path("${meta.id}.bam"),
                     path("${meta.id}.bam.bai"),
                     path("${meta.id}_filtered_feature_bc_matrix"), emit: aln
    path "${meta.id}/outs"   , emit: outs
    path "versions.yml"      , emit: versions

    script:
    def lib_lines = libraries.collect { lib ->
        def ft = lib.feature_type == 'gex' ? 'Gene Expression' : 'Antibody Capture'
        // cellranger requires an ABSOLUTE fastqs path
        "${lib.fastq_id},\${fastqs_dir},${ft}"
    }.join('\n')
    def expect = meta.expected_cells ? "expect-cells,${meta.expected_cells}" : ''
    // GEX-only when no feature/antibody reference is supplied
    def feature_section = fb_reference ? "\n[feature]\nreference,\$(readlink -f ${fb_reference})\n" : ''
    """
    fastqs_dir=\$(readlink -f fastqs)
    cat > multi_config.csv <<EOF
[gene-expression]
reference,\$(readlink -f ${cellranger_index})
create-bam,true
${expect}
${feature_section}
[libraries]
fastq_id,fastqs,feature_types
${lib_lines}
EOF

    cellranger multi \\
        --id=run_${meta.id} \\
        --csv=multi_config.csv \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()}

    # Normalize cellranger multi outputs to the canonical <sample>/outs/ layout (published)
    # plus uniquely sample-named handles for downstream joint callers (emitted in `aln`).
    mkdir -p ${meta.id}/outs
    bam=\$(find run_${meta.id}/outs -name 'sample_alignments.bam' -o -name 'possorted_genome_bam.bam' | head -n1)
    mtx=\$(find run_${meta.id}/outs -type d \\( -name 'sample_filtered_feature_bc_matrix' -o -name 'filtered_feature_bc_matrix' \\) | head -n1)
    abs_bam="\$(readlink -f \$bam)"
    abs_mtx="\$(readlink -f \$mtx)"
    ln -s "\$abs_bam"      ${meta.id}/outs/possorted_genome_bam.bam
    ln -s "\$abs_bam.bai"  ${meta.id}/outs/possorted_genome_bam.bam.bai
    ln -s "\$abs_mtx"      ${meta.id}/outs/filtered_feature_bc_matrix
    ln -s "\$abs_bam"      ${meta.id}.bam
    ln -s "\$abs_bam.bai"  ${meta.id}.bam.bai
    ln -s "\$abs_mtx"      ${meta.id}_filtered_feature_bc_matrix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellranger: \$(cellranger --version 2>&1 | sed 's/^.*cellranger-//; s/[^0-9.].*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/outs/filtered_feature_bc_matrix
    touch ${meta.id}/outs/possorted_genome_bam.bam ${meta.id}/outs/possorted_genome_bam.bam.bai
    echo | gzip > ${meta.id}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
    echo | gzip > ${meta.id}/outs/filtered_feature_bc_matrix/features.tsv.gz
    echo | gzip > ${meta.id}/outs/filtered_feature_bc_matrix/matrix.mtx.gz
    ln -s ${meta.id}/outs/possorted_genome_bam.bam      ${meta.id}.bam
    ln -s ${meta.id}/outs/possorted_genome_bam.bam.bai  ${meta.id}.bam.bai
    ln -s ${meta.id}/outs/filtered_feature_bc_matrix    ${meta.id}_filtered_feature_bc_matrix
    echo '"${task.process}": {cellranger: stub}' > versions.yml
    """
}

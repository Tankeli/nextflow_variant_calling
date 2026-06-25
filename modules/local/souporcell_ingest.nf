// Souporcell ingest: re-create the CELLRANGER `aln` channel contract off ALREADY-PUBLISHED
// Cell Ranger outs, for the standalone SOUPORCELL_ONLY entry (no Cell Ranger re-run).
// Aliases the published possorted_genome_bam.bam (or .cram) + filtered_feature_bc_matrix to the
// per-sample names SOUPORCELL_PREP expects (<id>.bam / <id>.bam.bai / <id>_filtered_feature_bc_matrix)
// so the same sample basename never collides when prep merges several into one mix.

process SOUPORCELL_INGEST {
    tag "$meta.id"
    label 'process_low'

    container params.samtools_container

    input:
    tuple val(meta), path(outs)

    output:
    tuple val(meta), path("${meta.id}.bam"),
                     path("${meta.id}.bam.bai"),
                     path("${meta.id}_filtered_feature_bc_matrix"), emit: aln

    script:
    """
    set -euo pipefail

    # Prefer BAM, fall back to CRAM (mirrors the caller convention in CLAUDE.md).
    bam=\$(find -L ${outs} -maxdepth 1 \\( -name 'possorted_genome_bam.bam' -o -name 'possorted_genome_bam.cram' \\) | head -n1)
    mtx=\$(find -L ${outs} -maxdepth 1 -type d -name 'filtered_feature_bc_matrix' | head -n1)
    [ -n "\$bam" ] || { echo "ERROR: no possorted_genome_bam.{bam,cram} under ${outs}" >&2; exit 1; }
    [ -n "\$mtx" ] || { echo "ERROR: no filtered_feature_bc_matrix under ${outs}" >&2; exit 1; }

    ln -s "\$(readlink -f "\$bam")" ${meta.id}.bam
    if [ -e "\$bam.bai" ]; then
        ln -s "\$(readlink -f "\$bam.bai")" ${meta.id}.bam.bai
    else
        samtools index -@ ${task.cpus} ${meta.id}.bam
    fi
    ln -s "\$(readlink -f "\$mtx")" ${meta.id}_filtered_feature_bc_matrix
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bam.bai
    mkdir -p ${meta.id}_filtered_feature_bc_matrix
    echo | gzip > ${meta.id}_filtered_feature_bc_matrix/barcodes.tsv.gz
    """
}

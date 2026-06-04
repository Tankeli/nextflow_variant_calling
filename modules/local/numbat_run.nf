// Numbat step 2: joint per-patient clone calling (run_numbat).
// Concatenates the per-sample allele counts + 10X matrices into one joint call.
// Relaxed thresholds (max_entropy/min_LLR) come from params; per-sample defaults fail
// on this cohort. Checkpoint output: <patient>/numbat_out/clone_post_1.tsv etc.

process NUMBAT_RUN {
    tag "$meta.id"
    label 'process_high'

    container params.numbat_container

    input:
    tuple val(meta), path(allele_counts), path(mtx_dirs)

    output:
    tuple val(meta), path("${meta.id}/numbat_out"), emit: numbat
    path "versions.yml"                           , emit: versions

    script:
    def samples     = meta.samples
    def sample_csv  = samples.join(',')
    def allele_csv  = samples.collect { "${it}_allele_counts.tsv.gz" }.join(',')
    def mtx_csv     = samples.collect { "${it}_filtered_feature_bc_matrix" }.join(',')
    """
    run_numbat.R \\
        ${meta.id} \\
        ${meta.id}/numbat_out \\
        ${sample_csv} \\
        ${allele_csv} \\
        ${mtx_csv} \\
        ${task.cpus} \\
        ${params.numbat_max_entropy} \\
        ${params.numbat_min_llr} \\
        hg38

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        numbat: \$(Rscript -e 'cat(as.character(packageVersion("numbat")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/numbat_out
    touch ${meta.id}/numbat_out/clone_post_1.tsv \\
          ${meta.id}/numbat_out/segs_consensus_1.tsv \\
          ${meta.id}/numbat_out/joint_post_1.tsv
    echo '"${task.process}": {numbat: stub}' > versions.yml
    """
}

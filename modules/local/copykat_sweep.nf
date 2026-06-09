// CopyKAT robustness sweep: one parameterised + seeded CopyKAT run per (sample x combo).
// Combo = KS.cut x win.size x ngene.chr x distance x normref x seed, carried in meta.combo.
// Outputs are named by the combo id (sam.name=comboid) so a sample's runs never collide.

process COPYKAT_SWEEP {
    tag "${meta.id}:${meta.combo.id}"
    label 'process_medium'

    container params.copykat_container

    input:
    tuple val(meta), path(matrix_dir), path(norm_barcodes)

    output:
    tuple val(meta), path("${meta.combo.id}_copykat_prediction.txt"),
                     path("${meta.combo.id}_copykat_CNA_results.txt"), emit: prediction
    path "${meta.combo.id}_copykat_*"                               , emit: results
    path "versions.yml"                                             , emit: versions

    script:
    def c = meta.combo
    def norm = (c.norm == 1) ? norm_barcodes : 'NONE'   // norm==0 arm ignores the staged file
    """
    copykat_sweep.R ${matrix_dir} ${c.id} ${task.cpus} \\
        ${c.ks_cut} ${c.win_size} ${c.ngene_chr} ${c.distance} ${c.seed} ${norm}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        copykat: \$(Rscript -e 'cat(as.character(packageVersion("copykat")))' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    def c = meta.combo
    """
    printf 'cell.names\\tcopykat.pred\\n' > ${c.id}_copykat_prediction.txt
    printf 'chrom\\tchrompos\\tabspos\\n' > ${c.id}_copykat_CNA_results.txt
    touch ${c.id}_copykat_clustering_results.rds
    echo '"${task.process}": {copykat: stub}' > versions.yml
    """
}

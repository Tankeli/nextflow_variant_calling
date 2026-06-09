// Phase-0 per-patient master table: joins Numbat-joint clones, souporcell clones, CopyKAT,
// cell types and LSC scores into one per-cell table (mgatk axis dropped).

process PHASE0_INTEGRATION {
    tag "$meta.id"
    label 'process_low'

    container params.scanpy_container

    input:
    tuple val(meta), path(numbat_out), path(soup_dir), path(celltype_csvs), path(copykat_txts), path(lsc_csvs)

    output:
    tuple val(meta), path("${meta.id}_cells.tsv"), emit: cells
    path "${meta.id}_clone_QC.{png,pdf}"          , emit: qc
    path "versions.yml"                           , emit: versions

    script:
    def samples    = meta.samples.join(',')
    def timepoints = meta.timepoints.join(',')
    def ct  = (celltype_csvs instanceof List ? celltype_csvs : [celltype_csvs]).join(',')
    def ck  = (copykat_txts  instanceof List ? copykat_txts  : [copykat_txts]).join(',')
    def lsc = (lsc_csvs      instanceof List ? lsc_csvs      : [lsc_csvs]).join(',')
    """
    soup=\$(find -L ${soup_dir} -name clusters.tsv | head -n1)
    phase0_integration.py ${meta.id} ${samples} ${timepoints} ${numbat_out} \\
        "\$soup" ${ct} ${ck} ${lsc}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python -c 'import platform; print(platform.python_version())' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    printf 'barcode\\tsample\\tpatient\\ttimepoint\\tnumbat_clone_joint\\tsouporcell_clone\\tpLSC6\\tmalignant\\n' > ${meta.id}_cells.tsv
    printf 'AAAA-1\\t${meta.samples[0]}\\t${meta.id}\\tDx\\tN1\\tS1\\t0.5\\tTrue\\n' >> ${meta.id}_cells.tsv
    echo stub > ${meta.id}_clone_QC.png
    echo stub > ${meta.id}_clone_QC.pdf
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

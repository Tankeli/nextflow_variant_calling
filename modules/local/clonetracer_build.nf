// Build the per-patient CloneTracer input JSON from this pipeline's caller outputs.
// Synthesises per-cell mutant (M) / reference (N) matrices over CNV (Numbat), nuclear-SNV
// (souporcell) and mtDNA (cellsnp/mgatk) mutations into one <patient>.json. Optional inputs
// (numbat_dir / soup_dir / mt_dirs / gtf) may be empty when their source did not run.

process CLONETRACER_BUILD {
    tag "$meta.id"
    label 'process_medium'

    container params.scanpy_container

    input:
    tuple val(meta), path(mtx_dirs), path(mt_dirs), path(numbat_dir), path(soup_dir), path(gtf)

    output:
    tuple val(meta), path("${meta.id}.json"), emit: json
    path "versions.yml"                      , emit: versions

    script:
    def samples    = meta.samples.join(',')
    def timepoints = meta.timepoints.join(',')
    def mtx_join   = (mtx_dirs instanceof List ? mtx_dirs : [mtx_dirs]).join(',')
    def numbat_arg = numbat_dir ? "--numbat-dir ${numbat_dir}"   : ''
    def soup_arg   = soup_dir   ? "--souporcell-dir ${soup_dir}" : ''
    def gtf_arg    = gtf        ? "--gtf ${gtf}"                  : ''
    def mt_join    = (mt_dirs instanceof List ? mt_dirs : (mt_dirs ? [mt_dirs] : [])).join(',')
    def mt_arg     = mt_join ? "--mtdna-dirs ${mt_join}" : ''
    def pseudobulk = params.clonetracer_pseudobulk ? '--pseudobulk' : ''
    """
    clonetracer_build_json.py \\
        --patient ${meta.id} \\
        --samples ${samples} \\
        --timepoints ${timepoints} \\
        --matrices ${mtx_join} \\
        ${numbat_arg} ${soup_arg} ${mt_arg} ${gtf_arg} ${pseudobulk} \\
        --max-snvs ${params.clonetracer_max_snvs} \\
        --mtdna-min-cells ${params.clonetracer_mtdna_min_cells} \\
        --mtdna-max-sites ${params.clonetracer_mtdna_max_sites} \\
        --output ${meta.id}.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '{"M": [[5]], "N": [[9]], "mut_type": [2], "mut_names": ["mt_3243_A_G"], "r_cnv": [0.0], "cell_barcode": ["${meta.samples[0]}__AAAA-1"]}' > ${meta.id}.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}

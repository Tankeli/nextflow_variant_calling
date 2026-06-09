// mtDNA variant pileup (per sample) — produces per-cell alt/ref counts over chrM variants
// for the CloneTracer mtDNA axis. Method-switchable via params.clonetracer_mtdna_method:
//   'cellsnp' (default) -> cellsnp-lite -p whole-mito pileup (cellsnp-lite ships in numbat.sif)
//   'mgatk'   (opt-in)  -> mgatk, then normalised to cellSNP-style files (mgatk_to_cellsnp.py)
// Emits a normalised <sample>_mtdna/ dir with cellSNP.tag.{AD,DP}.mtx + base.vcf + samples.tsv
// regardless of method, so clonetracer_build_json.py reads one format.

process MTDNA_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    container params.numbat_container

    input:
    tuple val(meta), path(bam), path(bai), path(matrix_dir)

    output:
    tuple val(meta), path("${meta.id}_mtdna"), emit: mtdna
    path "versions.yml"                       , emit: versions

    script:
    def method = params.clonetracer_mtdna_method
    def chrom_override = params.clonetracer_mtdna_chrom ?: ''
    if (method == 'mgatk')
        """
        zcat ${matrix_dir}/barcodes.tsv.gz > barcodes.txt 2>/dev/null || cp ${matrix_dir}/barcodes.tsv barcodes.txt

        mgatk tenx \\
            -i ${bam} \\
            -n ${meta.id} \\
            -o mgatk_out \\
            -c ${task.cpus} \\
            -bt CB -b barcodes.txt

        mgatk_to_cellsnp.py mgatk_out/final ${meta.id}_mtdna

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            mgatk: \$(mgatk --version 2>&1 | tail -n1 || echo NA)
        END_VERSIONS
        """
    else
        """
        zcat ${matrix_dir}/barcodes.tsv.gz > barcodes.txt 2>/dev/null || cp ${matrix_dir}/barcodes.tsv barcodes.txt

        # Mito contig name: explicit override, else auto-detect from the BAM, else fall back to chrM.
        MT='${chrom_override}'
        if [ -z "\$MT" ]; then
            MT=\$(samtools idxstats ${bam} 2>/dev/null | cut -f1 | grep -iE '^(chr)?M(T)?\$' | head -n1)
        fi
        [ -z "\$MT" ] && MT=chrM

        cellsnp-lite \\
            -s ${bam} \\
            -b barcodes.txt \\
            -O ${meta.id}_mtdna \\
            -p ${task.cpus} \\
            --chrom \$MT \\
            --minMAF ${params.clonetracer_mtdna_min_maf} \\
            --minCOUNT ${params.clonetracer_mtdna_min_count} \\
            --cellTAG CB --UMItag UB --genotype

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cellsnp-lite: \$(cellsnp-lite --version 2>&1 | tail -n1 || echo NA)
        END_VERSIONS
        """

    stub:
    """
    mkdir -p ${meta.id}_mtdna
    printf '%%%%MatrixMarket matrix coordinate integer general\\n1 1 1\\n1 1 5\\n' > ${meta.id}_mtdna/cellSNP.tag.AD.mtx
    printf '%%%%MatrixMarket matrix coordinate integer general\\n1 1 1\\n1 1 9\\n' > ${meta.id}_mtdna/cellSNP.tag.DP.mtx
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\nchrM\\t3243\\t.\\tA\\tG\\t.\\t.\\t.\\n' > ${meta.id}_mtdna/cellSNP.base.vcf
    printf '${meta.id}__AAAA-1\\n' > ${meta.id}_mtdna/cellSNP.samples.tsv
    echo '"${task.process}": {cellsnp-lite: stub}' > versions.yml
    """
}

// Build a plain-text known-normal barcode list for the CopyKAT normref sweep arm.
// Selects cells whose reference-mapped ref_cell_type is in params.copykat_norm_celltypes
// (the lineages we treat as confidently diploid, e.g. T / NK). Fed to copykat()'s
// norm.cell.names so the aneuploid/diploid boundary is anchored on a real baseline.

process COPYKAT_NORM_BARCODES {
    tag "$meta.id"
    label 'process_low'

    container params.scanpy_container

    input:
    tuple val(meta), path(celltypes_csv)
    val  norm_celltypes   // list of ref_cell_type labels

    output:
    tuple val(meta), path("${meta.id}_norm_barcodes.txt"), emit: barcodes
    path "versions.yml"                                  , emit: versions

    script:
    def labels = norm_celltypes.join(',')
    """
    python3 - <<'PY'
    import pandas as pd
    labels = {x.strip() for x in "${labels}".split(",") if x.strip()}
    df = pd.read_csv("${celltypes_csv}", index_col=0)
    ct = "ref_cell_type" if "ref_cell_type" in df.columns else df.columns[0]
    keep = df.index[df[ct].astype(str).isin(labels)]
    with open("${meta.id}_norm_barcodes.txt", "w") as fh:
        fh.write("\\n".join(map(str, keep)))
    print(f"${meta.id}: {len(keep)}/{len(df)} cells flagged normal ({sorted(labels)})")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: \$(python3 -c 'import pandas; print(pandas.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_norm_barcodes.txt
    echo '"${task.process}": {pandas: stub}' > versions.yml
    """
}

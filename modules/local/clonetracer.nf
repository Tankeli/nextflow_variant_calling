// CloneTracer (veltenlab) Bayesian clonal inference. Consumes the per-patient JSON built by
// CLONETRACER_BUILD and infers a clonal hierarchy + per-cell clone posterior probabilities.
// run_clonetracer.py / helper_functions.py are vendored in bin/; helper_functions is imported
// as a module, so bin/ is put on PYTHONPATH. --multiple_samples (-s) is set when the patient
// has >1 timepoint (Dx+Rel), matching the joint design. Checkpoint: <patient>_out.pickle.

process CLONETRACER {
    tag "$meta.id"
    label 'process_high'

    container params.clonetracer_container

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("${meta.id}_clone_assignments.csv"), emit: assignments
    tuple val(meta), path("${meta.id}_out.pickle"), path("${meta.id}_tree.pickle"), emit: trees
    path "versions.yml"                                       , emit: versions

    script:
    // -s (multiple_samples) needs bulk_M/bulk_N in the JSON, which CLONETRACER_BUILD only writes
    // under --pseudobulk. Without bulk data the model must run pooled (clones are still joint across
    // Dx+Rel because all cells share one M/N matrix). So gate -s on pseudobulk + >1 timepoint.
    def multi = (params.clonetracer_pseudobulk && meta.timepoints.unique().size() > 1) ? '-s' : ''
    def gpu   = params.clonetracer_gpu ? '-g' : ''
    // run_clonetracer sets init = iters - 100, so iters in 1..99 break the ELBO diagnostic; only
    // pass -t when >=100, otherwise let it auto-pick (300 SNV-only / 500 with CNVs).
    def iters = (params.clonetracer_iterations as int) >= 100 ? "-t ${params.clonetracer_iterations}" : ''
    """
    export PYTHONPATH="${projectDir}/bin:\$PYTHONPATH"

    run_clonetracer.py \\
        -i ${json} \\
        -n ${meta.id} \\
        -o . \\
        ${multi} ${gpu} ${iters}

    clonetracer_assignments.py ${meta.id}_out.pickle ${meta.id}_clone_assignments.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyro: \$(python3 -c 'import pyro; print(pyro.__version__)' 2>/dev/null || echo NA)
        torch: \$(python3 -c 'import torch; print(torch.__version__)' 2>/dev/null || echo NA)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_out.pickle ${meta.id}_tree.pickle
    printf 'barcode,sample,clone,max_prob,tree\\n${meta.samples[0]}__AAAA-1,${meta.samples[0]},0,0.99,1\\n' > ${meta.id}_clone_assignments.csv
    echo '"${task.process}": {pyro: stub}' > versions.yml
    """
}

#!/usr/bin/env nextflow
/*
 * DDE_33 — single-cell variant-calling pipeline
 * Numbat (CNV) + CopyKAT (aneuploidy) + souporcell (SNV) for paired Dx/Rel paediatric AML.
 * See CLAUDE.md (architecture) and scratchpad.md (build plan).
 */

nextflow.enable.dsl = 2

include { VARIANTCALLING } from './workflows/variantcalling'
include { RNA_DOWNSTREAM } from './subworkflows/local/rna_downstream'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    if (params.validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    VARIANTCALLING()
}

//
// Standalone downstream-only entry: run the RNA best-practices stack off ALREADY-PUBLISHED
// Cell Ranger filtered matrices (no Cell Ranger re-run, no variant callers). Input is a small CSV
// (sample,patient,timepoint,matrix) via --downstream_samplesheet; toggle stages with the run_rna_*
// / run_pseudotime / run_de / run_composition / run_protein params. SoupX is unavailable here (the
// raw matrix is not published) — leave run_soupx=false.
//   nextflow run . -entry DOWNSTREAM -profile viking -params-file params-downstream-patients.yaml
//
workflow DOWNSTREAM {
    log.info paramsSummaryLog(workflow)

    ch_mtx = Channel.fromPath(params.downstream_samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample, patient: row.patient, timepoint: row.timepoint ]
            tuple( meta, file(row.matrix, checkIfExists: true) )
        }

    RNA_DOWNSTREAM( ch_mtx, Channel.empty() )
}

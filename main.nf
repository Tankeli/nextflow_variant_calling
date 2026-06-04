#!/usr/bin/env nextflow
/*
 * DDE_33 — single-cell variant-calling pipeline
 * Numbat (CNV) + CopyKAT (aneuploidy) + souporcell (SNV) for paired Dx/Rel paediatric AML.
 * See CLAUDE.md (architecture) and scratchpad.md (build plan).
 */

nextflow.enable.dsl = 2

include { VARIANTCALLING } from './workflows/variantcalling'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    if (params.validate_params) {
        validateParameters()
    }
    log.info paramsSummaryLog(workflow)

    VARIANTCALLING()
}

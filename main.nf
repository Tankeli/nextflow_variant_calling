#!/usr/bin/env nextflow
/*
 * DDE_33 — single-cell variant-calling pipeline
 * Numbat (CNV) + CopyKAT (aneuploidy) + souporcell (SNV) for paired Dx/Rel paediatric AML.
 * See CLAUDE.md (architecture) and scratchpad.md (build plan).
 */

nextflow.enable.dsl = 2

include { VARIANTCALLING } from './workflows/variantcalling'
include { RNA_DOWNSTREAM } from './subworkflows/local/rna_downstream'
include { PROTEOMICS as PROTEOMICS_WF } from './subworkflows/local/proteomics'
include { SOUPORCELL_INGEST } from './modules/local/souporcell_ingest'
include { SOUPORCELL_WF } from './subworkflows/local/souporcell'
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

//
// Standalone bulk-proteomics entry: run the DDE_31 bulk-MS branch off a Spectronaut search matrix
// (NO FASTQ, no callers). Inputs via params: proteomics_norm / proteomics_nonnorm / proteomics_design
// (+ optional proteomics_contaminants). DESP proportions come from --proteomics_proportions, or from
// already-published scRNA celltypes via --proteomics_celltypes_glob. Toggle stages with prot_run_*.
//   nextflow run . -entry PROTEOMICS -profile viking -params-file params-proteomics.yaml
//
workflow PROTEOMICS {
    log.info paramsSummaryLog(workflow)

    if (!params.proteomics_norm || !params.proteomics_design) {
        error "PROTEOMICS entry needs --proteomics_norm and --proteomics_design"
    }
    def no_file = file("${projectDir}/assets/NO_FILE")
    ch_inputs = Channel.value( tuple(
        params.proteomics_nonnorm ? file(params.proteomics_nonnorm, checkIfExists: true) : file(params.proteomics_norm, checkIfExists: true),
        file(params.proteomics_norm, checkIfExists: true),
        file(params.proteomics_design, checkIfExists: true),
        params.proteomics_contaminants ? file(params.proteomics_contaminants, checkIfExists: true) : no_file
    ) )

    ch_celltypes = Channel.empty()
    if (params.prot_run_desp && !params.proteomics_proportions) {
        if (!params.proteomics_celltypes_glob) {
            error "prot_run_desp (standalone) needs --proteomics_proportions or --proteomics_celltypes_glob"
        }
        ch_celltypes = Channel.fromPath(params.proteomics_celltypes_glob, checkIfExists: true).collect()
    }

    PROTEOMICS_WF( ch_inputs, ch_celltypes )
}

//
// Standalone souporcell-only entry: run souporcell deconvolution off ALREADY-PUBLISHED Cell Ranger
// outs (no Cell Ranger re-run, no other callers). Built for the deconvolution-validation experiments
// (artificial multi-individual mixes + real post-BMT donor/recipient samples). Group several samples
// under one `patient` value to make an artificial mix; the souporcell K sweep (--souporcell_k) then
// tries to recover them. Ground truth is preserved via the `<sample>__` barcode prefix that
// SOUPORCELL_PREP adds, so bin/souporcell_mix_eval.py can score the assignment.
// Input CSV (sample,patient,timepoint,outs) via --souporcell_samplesheet, where `outs` points at a
// published results_*/cellranger/<sample>/outs directory.
//   nextflow run . -entry SOUPORCELL_ONLY -profile viking --souporcell_samplesheet assets/test/souporcell_mix_controls.csv
//
workflow SOUPORCELL_ONLY {
    log.info paramsSummaryLog(workflow)

    if (!params.souporcell_samplesheet) {
        error "SOUPORCELL_ONLY entry needs --souporcell_samplesheet (sample,patient,timepoint,outs)"
    }

    ch_in = Channel.fromPath(params.souporcell_samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample, patient: row.patient, timepoint: row.timepoint ]
            // resolve relative `outs` against projectDir, not launchDir (the job uses an isolated launchDir)
            def outs = row.outs.startsWith('/') ? row.outs : "${projectDir}/${row.outs}"
            tuple( meta, file(outs, checkIfExists: true) )
        }

    SOUPORCELL_INGEST( ch_in )

    // Mirror VARIANTCALLING's per-patient grouping (Dx sorted before Rel for determinism).
    ch_patient_aln = SOUPORCELL_INGEST.out.aln
        .map { meta, bam, bai, mtx -> tuple( meta.patient, [ meta: meta, bam: bam, bai: bai, mtx: mtx ] ) }
        .groupTuple()
        .map { patient, members -> tuple( patient, members.sort { it.meta.timepoint == 'Dx' ? 0 : 1 } ) }

    SOUPORCELL_WF( ch_patient_aln, file(params.souporcell_fasta, checkIfExists: true) )
}

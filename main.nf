#!/usr/bin/env nextflow

/*
 * main.nf — transcript_assembly pipeline entry point
 *
 * Usage:
 *   nextflow run main.nf \
 *     --input      samplesheet.csv \
 *     --gtf        Homo_sapiens.GRCh38.95.gtf \
 *     --fasta      GRCh38.fa \
 *     --star_sjouts star_sjouts.txt \
 *     --star_bams   star_bams.txt \
 *     --outdir      results/ \
 *     --espresso_scratch /storage1/scratch/espresso_run1 \
 *     -profile ris
 *
 * See README.md for detailed usage and parameter descriptions.
 */

nextflow.enable.dsl = 2

include { TRANSCRIPT_ASSEMBLY } from './workflows/transcript_assembly'

// ── Parameter validation ──────────────────────────────────────────────────────

def required = ['input', 'gtf', 'fasta', 'star_sjouts', 'star_bams', 'outdir', 'espresso_scratch']
required.each { p ->
    if (!params[p]) {
        log.error "Required parameter '--${p}' is not set."
        System.exit(1)
    }
}

// ── Pipeline summary ──────────────────────────────────────────────────────────

log.info """\
    T R A N S C R I P T   A S S E M B L Y
    ======================================
    samplesheet   : ${params.input}
    gtf           : ${params.gtf}
    fasta         : ${params.fasta}
    star_sjouts   : ${params.star_sjouts}
    star_bams     : ${params.star_bams}
    outdir        : ${params.outdir}
    espresso_scratch: ${params.espresso_scratch}
    """.stripIndent()

// ── Entry workflow ────────────────────────────────────────────────────────────

workflow {

    TRANSCRIPT_ASSEMBLY(
        file(params.input,        checkIfExists: true),
        file(params.gtf,          checkIfExists: true),
        file(params.fasta,        checkIfExists: true),
        file(params.star_sjouts,  checkIfExists: true),
        file(params.star_bams,    checkIfExists: true)
    )
}

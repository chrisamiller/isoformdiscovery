/*
 * subworkflows/local/input_check.nf
 *
 * Validates the samplesheet CSV and emits a channel of (meta, bam) tuples
 * consumed by the main workflow.
 *
 * The meta map carries: [group: 'g1'|'g2', sample: 'AML001']
 *
 * Emit:
 *   reads   — channel of tuple(val(meta), path(bam))
 *   versions
 */

include { SAMPLESHEET_CHECK } from '../../modules/local/SAMPLESHEET_CHECK'

workflow INPUT_CHECK {

    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK(samplesheet)

    // Parse validated CSV into (meta, bam) channel
    ch_reads = SAMPLESHEET_CHECK.out.csv
        .splitCsv(header: true, strip: true)
        .map { row ->
            def meta = [group: row.group, sample: row.sample]
            def bam  = file(row.bam, checkIfExists: true)
            tuple(meta, bam)
        }

    emit:
    reads    = ch_reads
    versions = SAMPLESHEET_CHECK.out.versions
}

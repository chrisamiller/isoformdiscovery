/*
 * SAMPLESHEET_CHECK
 *
 * Validates the input CSV samplesheet using check_samplesheet.py (bin/).
 * Emits a validated CSV and the parsed sample records as a channel of
 * [meta, bam] tuples consumed by the main workflow.
 *
 * Input CSV columns: group, sample, bam
 */

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'
    container 'python:3.11'

    input:
    path samplesheet

    output:
    path 'validated_samplesheet.csv', emit: csv
    path 'versions.yml',              emit: versions

    script:
    """
    check_samplesheet.py \\
        ${samplesheet} \\
        validated_samplesheet.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}

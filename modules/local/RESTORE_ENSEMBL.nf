/*
 * RESTORE_ENSEMBL
 *
 * Restores Ensembl-annotated transcripts (transcript_id starting with 'ENST')
 * that the SQANTI3 ML filter classified as Artifact back into the filtered GTF.
 *
 * Rationale: Ensembl annotations are treated as a ground-truth reference.
 * Even when the ML model flags them as artifact-like (e.g., due to low
 * expression or unusual junction patterns in this cohort), we keep them
 * in the transcript set to avoid systematic loss of known biology.
 *
 * Steps:
 *   1. Extract ENST transcript IDs marked 'Artifact' in ML results
 *   2. Pull their GTF records from the original (unfiltered) merged GTF
 *      using filter_gtf.pl (bin/)
 *   3. Append to the ML-filtered GTF
 *   4. Sort the combined GTF by chromosome then start position
 *
 * Input:
 *   path ml_classification  — *_MLresult_classification.txt from SQANTI3_FILTER
 *   path original_gtf       — merged_N2_R0_updated.gtf (pre-filter, all transcripts)
 *   path filtered_gtf       — *.filtered.gtf from SQANTI3_FILTER
 *
 * Output:
 *   path "*.ensrestore.gtf", emit: gtf
 *   path "versions.yml",     emit: versions
 */

process RESTORE_ENSEMBL {
    label 'process_single'
    container 'chrisamiller/genomic-analysis:0.2'

    input:
    path ml_classification
    path original_gtf
    path filtered_gtf

    output:
    path "*.ensrestore.gtf", emit: gtf
    path "versions.yml",     emit: versions

    script:
    def out_name = filtered_gtf.baseName + '.ensrestore.gtf'
    """
    # Find ENST transcripts flagged as Artifact
    grep '^ENST' ${ml_classification} \\
        | awk '\$NF == "Artifact"' \\
        | cut -f1 \\
        > enst_restore.txt

    # If any exist, extract their GTF records and append to filtered GTF
    cp ${filtered_gtf} restore_combined.gtf
    if [ -s enst_restore.txt ]; then
        filter_gtf.pl enst_restore.txt ${original_gtf} >> restore_combined.gtf
    fi

    # Sort: header lines first, then data sorted by chr + start
    grep '^#' restore_combined.gtf > ${out_name} || true
    grep -v '^#' restore_combined.gtf \\
        | sort -k1,1 -k4,4n \\
        >> ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | grep 'v[0-9]' | sed 's/.*v/v/' | sed 's/ .*//')
    END_VERSIONS
    """
}

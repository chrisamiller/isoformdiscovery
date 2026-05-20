/*
 * FILTER_NORMALIZE
 *
 * Applies a coverage filter to the merged FSM abundance matrix and produces
 * three normalized output tables plus a filtered GTF.
 *
 * Coverage filter: transcripts with total read count < params.filter_min_count
 * are dropped.  The "noism" label in output filenames reflects that the FSM
 * re-quantification round already excludes ISM reads by design.
 *
 * Input:
 *   path abundance  — merged_fsm_N2_R2_abundance.esp (from MERGE_ESPRESSO_FSM)
 *   path gtf        — merged_fsm_N2_R2_updated.gtf   (from MERGE_ESPRESSO_FSM)
 *
 * Output:
 *   path "*.covfilt.noism.esp",     emit: filtered_counts
 *   path "*.covfilt.noism.cpm.esp", emit: cpm
 *   path "*.covfilt.noism.tpm.esp", emit: tpm
 *   path "*.covfilt.noism.gtf",     emit: filtered_gtf
 *   path "versions.yml",            emit: versions
 */

process FILTER_NORMALIZE {
    label 'process_single'
    container 'chrisamiller/genomic-analysis:0.3.1'

    input:
    path abundance
    path gtf

    output:
    path "${abundance.baseName}.covfilt.noism.esp",     emit: filtered_counts
    path "${abundance.baseName}.covfilt.noism.cpm.esp", emit: cpm
    path "${abundance.baseName}.covfilt.noism.tpm.esp", emit: tpm
    path "${gtf.baseName}.covfilt.noism.gtf",           emit: filtered_gtf
    path "versions.yml",                                 emit: versions

    script:
    def abund_prefix = abundance.baseName
    def gtf_prefix   = gtf.baseName
    def min_count    = params.filter_min_count ?: 5
    """
    filter_normalize_abundance.R ${abundance} ${gtf} ${min_count}

    cut -f1 ${abund_prefix}.covfilt.noism.esp | tail -n +2 > kept_transcripts.txt
    filter_gtf.pl kept_transcripts.txt ${gtf} > ${gtf_prefix}.covfilt.noism.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //' | sed 's/ .*//')
        edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}

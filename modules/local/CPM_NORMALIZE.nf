/*
 * CPM_NORMALIZE
 *
 * Converts the merged FSM abundance count matrix to TMM-normalized CPM values
 * using espresso_to_cpm.R (bin/), which applies edgeR's calcNormFactors() +
 * cpm() workflow.
 *
 * Input:
 *   path abundance  — merged_fsm_N2_R2_abundance.esp from MERGE_ESPRESSO_FSM
 *
 * Output:
 *   path "*_abundance.cpm.tsv", emit: cpm
 *   path "versions.yml",        emit: versions
 */

process CPM_NORMALIZE {
    label 'process_single'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'

    input:
    path abundance

    output:
    path "*_abundance.cpm.tsv", emit: cpm
    path "versions.yml",        emit: versions

    script:
    def out_name = abundance.baseName + '.cpm.tsv'
    """
    Rscript -e "if (!requireNamespace('edgeR', quietly=TRUE)) BiocManager::install('edgeR', ask=FALSE, update=FALSE)"
    espresso_to_cpm.R ${abundance} ${out_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //' | sed 's/ .*//')
        edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """
}

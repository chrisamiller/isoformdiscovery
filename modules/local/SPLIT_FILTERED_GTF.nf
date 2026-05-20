/*
 * SPLIT_FILTERED_GTF
 *
 * Re-scatters the genome-wide filtered GTF back to per-chromosome slices for
 * the ESPRESSO_Q round-2 re-quantification step.
 *
 * One invocation per chromosome; the chr value is used as an exact-word grep
 * against the first GTF column.
 *
 * Input:
 *   path filtered_gtf  — merged, Ensembl-restored, filtered GTF
 *   val  chr           — chromosome name to extract (e.g. 'chr1')
 *
 * Output:
 *   tuple val(chr), path("${chr}_filtered.gtf"), emit: chr_gtf
 *   path "versions.yml",                          emit: versions
 */

process SPLIT_FILTERED_GTF {
    tag "$chr"
    label 'process_single'
    container 'chrisamiller/genomic-analysis:0.2'

    input:
    path filtered_gtf
    val  chr

    output:
    tuple val(chr), path("${chr}_filtered.gtf"), emit: chr_gtf
    path  "versions.yml",                         emit: versions

    script:
    """
    # Match lines where the first field (chromosome) equals chr exactly
    grep -P "^${chr}\\t" ${filtered_gtf} > ${chr}_filtered.gtf || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        grep: \$(grep --version | head -n1 | sed 's/grep //')
    END_VERSIONS
    """
}

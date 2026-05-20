/*
 * BAM_MERGE_INDEX
 *
 * Assembles a per-sample genome-wide FSM BAM by merging all per-chromosome
 * FSM BAMs produced by FSM_FILTER and indexes the result.
 *
 * The per-chromosome FSM BAMs reside on params.espresso_scratch; their paths
 * are discovered via glob.  The FSM_FILTER sentinels are declared as input to
 * ensure all per-chr filter jobs have finished before the merge starts.
 *
 * Input:
 *   tuple val(meta), path(fsm_sentinels)  — all FSM_FILTER sentinels for this
 *                                           sample (collected across chromosomes)
 *
 * Output:
 *   tuple val(meta), path("${meta.sample}.bam"),     emit: bam
 *   tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
 *   path "versions.yml",                              emit: versions
 */

process BAM_MERGE_INDEX {
    tag "${meta.sample}"
    label 'process_merge'
    container 'chrisamiller/docker-genomic-analysis:latest'

    input:
    tuple val(meta), path(fsm_sentinels)

    output:
    tuple val(meta), path("${meta.sample}.bam"),     emit: bam
    tuple val(meta), path("${meta.sample}.bam.bai"), emit: bai
    path  "versions.yml",                             emit: versions

    script:
    def sample   = meta.sample
    def scratch  = params.espresso_scratch
    """
    # Collect all per-chromosome FSM BAMs for this sample
    bam_list=\$(ls ${scratch}/*/fsm/bams/${sample}.bam 2>/dev/null | tr '\\n' ' ')

    if [ -z "\$bam_list" ]; then
        echo "ERROR: No FSM BAMs found for sample ${sample}" >&2
        exit 1
    fi

    samtools merge -f ${sample}.bam \$bam_list
    samtools index ${sample}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

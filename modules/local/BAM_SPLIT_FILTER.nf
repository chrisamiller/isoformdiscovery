/*
 * BAM_SPLIT_FILTER
 *
 * Splits a Nanopore cDNA BAM by chromosome and removes secondary/supplementary
 * alignments (samtools view -F 0x900).  Chromosome names are detected from the
 * BAM @SQ headers — both 'chrN' and plain 'N' styles are supported.
 *
 * One output BAM per standard chromosome (1–22, X, Y) is produced; non-standard
 * contigs and random/unplaced sequences are ignored.
 *
 * Input:
 *   tuple val(meta), path(bam)   — meta = [group, sample]
 *
 * Output:
 *   tuple val(meta), path("${meta.sample}.*.bam"), emit: bams
 *   path "versions.yml",            emit: versions
 *
 * Note: output BAMs are named <sample>.<chr>.bam to avoid filename collisions
 * when multiple samples' per-chr BAMs are grouped together downstream.
 */

process BAM_SPLIT_FILTER {
    tag "${meta.sample}"
    label 'process_bam_split'
    container 'chrisamiller/genomic-analysis:0.2'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.*.bam"), emit: bams
    path  "versions.yml",                          emit: versions

    script:
    """
    # Index must exist in the work dir for random-access region extraction
    samtools index ${bam}

    # Extract standard chromosome names from @SQ headers
    # Handles both 'chr1' and '1' naming conventions
    samtools view -H ${bam} \\
        | grep '^@SQ' \\
        | sed 's/.*SN://' \\
        | sed 's/\t.*//' \\
        | grep -E '^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)\$' \\
        > chroms.txt

    # Split and filter per chromosome; prefix with sample name to ensure unique filenames
    while read chr; do
        samtools view -F 0x900 -b -o "${meta.sample}.\${chr}.bam" ${bam} "\${chr}"
    done < chroms.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

/*
 * BUILD_CHR_TSV
 *
 * Assembles a per-chromosome ESPRESSO input TSV from the grouped-by-chromosome
 * channel of (meta, bam) tuples produced by BAM_SPLIT_FILTER.
 *
 * Each TSV row is:   /absolute/path/to/sample_chr.bam <TAB> sample_name
 *
 * Paths are resolved to their true on-disk location via toRealPath() so that
 * Docker-mounted paths inside ESPRESSO containers are correct (not Nextflow
 * work-dir symlinks).
 *
 * This process uses an exec: block (pure Groovy) — no container is launched.
 *
 * Input:
 *   tuple val(chr), val(metas), path(bams)
 *     — chr:   chromosome name (e.g. 'chr1')
 *     — metas: list of meta maps, one per sample ([group, sample])
 *     — bams:  list of per-sample BAM paths for this chromosome
 *
 * Output:
 *   tuple val(chr), path("${chr}.tsv"), emit: tsv
 */

process BUILD_CHR_TSV {
    tag "$chr"

    input:
    tuple val(chr), val(metas), val(bams)

    output:
    tuple val(chr), path("${chr}.tsv"), emit: tsv

    exec:
    def tsv = task.workDir.resolve("${chr}.tsv")
    tsv.withWriter { w ->
        [metas, bams].transpose().each { meta, bam ->
            def realPath = bam.toRealPath().toString()
            w.writeLine("${realPath}\t${meta.sample}")
        }
    }
}

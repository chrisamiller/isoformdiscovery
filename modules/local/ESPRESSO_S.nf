/*
 * ESPRESSO_S
 *
 * ESPRESSO Phase 1: Signal detection (per chromosome).
 *
 * Runs ESPRESSO_S.pl for all samples of one chromosome simultaneously.
 * Output is written directly to params.espresso_scratch/${chr}/ — this is the
 * Option C shared-scratch design (see implementation_plan.md §8 and README).
 *
 * A small sentinel file is emitted as the Nextflow output so downstream
 * processes (ESPRESSO_C, ESPRESSO_Q) can declare a dependency without staging
 * the large ESPRESSO working directory.
 *
 * Resume behaviour: if Nextflow's work cache has the sentinel, the scratch
 * outputs are assumed present.  The scratch directory is wiped at the start of
 * the script to guarantee a clean re-run when the cache is invalidated.
 *
 * Input:
 *   tuple val(chr), path(tsv)  — per-chromosome ESPRESSO input TSV
 *   path  fasta                — reference genome FASTA
 *   path  gtf                  — reference annotation GTF
 *
 * Output:
 *   tuple val(chr), path("done_s_${chr}.txt"), emit: sentinel
 *   path "versions.yml",                        emit: versions
 */

process ESPRESSO_S {
    tag "$chr"
    label 'process_espresso_s'
    container 'sridnona/espresso:v2'

    input:
    tuple val(chr), path(tsv)
    path  fasta
    path  gtf

    output:
    tuple val(chr), path("done_s_${chr}.txt"), emit: sentinel
    path  "versions.yml",                       emit: versions

    script:
    def scratch = "${params.espresso_scratch}/${chr}"
    """
    export PATH="/opt/conda/envs/env/bin:\${PATH}"

    # Remove any stale scratch state so this run starts clean
    rm -rf ${scratch}
    mkdir -p ${scratch}

    perl /bin/espresso/src/ESPRESSO_S.pl \\
        -L ${tsv} \\
        -F ${fasta} \\
        -A ${gtf} \\
        -O ${scratch} \\
        -T ${task.cpus}

    touch done_s_${chr}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        espresso: \$(perl /bin/espresso/src/ESPRESSO_S.pl --version 2>&1 | head -n1 || echo 'v2')
    END_VERSIONS
    """
}

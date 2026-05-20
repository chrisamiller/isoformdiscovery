/*
 * ESPRESSO_C
 *
 * ESPRESSO Phase 2: Error correction (per chromosome × per sample).
 *
 * Runs one ESPRESSO_C.pl job for a single (chromosome, sample_index) pair.
 * Because ESPRESSO_C reads and writes to params.espresso_scratch/${chr}/, all
 * per-sample jobs for the same chromosome share that directory on the cluster
 * filesystem.  Nextflow tracks them via sentinel files only.
 *
 * Parallelism: one Nextflow process per (chr, sample_index) — the maximum
 * parallelism option chosen in the implementation plan (Q2).
 *
 * Input:
 *   tuple val(chr), val(sample_idx), path(done_s)  — sentinel from ESPRESSO_S
 *   path  fasta                                     — reference genome FASTA
 *
 * Output:
 *   tuple val(chr), path("done_c_${chr}_${sample_idx}.txt"), emit: sentinel
 *   path "versions.yml",                                      emit: versions
 */

process ESPRESSO_C {
    tag "${chr} idx=${sample_idx}"
    label 'process_espresso_c'
    container 'sridnona/espresso:v2'

    input:
    tuple val(chr), val(sample_idx), path(done_s)
    path  fasta

    output:
    tuple val(chr), path("done_c_${chr}_${sample_idx}.txt"), emit: sentinel
    path  "versions.yml",                                     emit: versions

    script:
    def scratch = "${params.espresso_scratch}/${chr}"
    """
    export PATH="/opt/conda/envs/env/bin:\${PATH}"

    # Skip correction if ESPRESSO_S produced no reads for this chr/sample.
    if [ ! -s "${scratch}/${sample_idx}/sam.list3" ]; then
        touch done_c_${chr}_${sample_idx}.txt
        printf '"${task.process}":\\n    espresso: v2\\n' > versions.yml
        exit 0
    fi

    perl /bin/espresso/src/ESPRESSO_C.pl \\
        -I ${scratch} \\
        -F ${fasta} \\
        -X ${sample_idx} \\
        -T ${task.cpus}

    touch done_c_${chr}_${sample_idx}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        espresso: \$(perl /bin/espresso/src/ESPRESSO_C.pl --version 2>&1 | head -n1 || echo 'v2')
    END_VERSIONS
    """
}

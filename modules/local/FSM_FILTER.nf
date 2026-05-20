/*
 * FSM_FILTER
 *
 * Filters ESPRESSO intermediate files and per-sample BAMs to retain only
 * Full Splice Match (FSM) reads, using remove_nonfsm_data.pl (bin/).
 *
 * FSM reads are those whose alignment exactly matches all splice junctions of
 * a known or novel transcript.  Restricting to FSM-only reads before the final
 * ESPRESSO_Q round increases count-matrix precision at the cost of some recall.
 *
 * This process operates directly on the shared scratch directory
 * (params.espresso_scratch/${chr}/) — it reads the round-2 compatible_isoform.tsv
 * to identify FSM read IDs, then creates a new fsm/ subdirectory in the
 * chromosome scratch dir with FSM-only versions of all intermediate files and BAMs.
 *
 * Input:
 *   tuple val(chr), path(compatible_isoform)  — *_compatible_isoform.tsv from
 *                                                ESPRESSO_Q_R2 (filtered/ subdir)
 *   path(done_q_r2)                           — ESPRESSO_Q_R2 sentinel (ordering)
 *
 * Output:
 *   tuple val(chr), path("done_fsm_${chr}.txt"), emit: sentinel
 *   path "versions.yml",                          emit: versions
 */

process FSM_FILTER {
    tag "$chr"
    label 'process_fsm_filter'
    container 'chrisamiller/docker-genomic-analysis:latest'

    input:
    tuple val(chr), path(compatible_isoform)
    path  done_q_r2

    output:
    tuple val(chr), path("done_fsm_${chr}.txt"), emit: sentinel
    path  "versions.yml",                         emit: versions

    script:
    def scratch_dir  = "${params.espresso_scratch}/${chr}"
    def tsv_updated  = "${scratch_dir}/${chr}.tsv.updated"
    """
    compat_file="${compatible_isoform}"
    if [[ "${compatible_isoform}" == *.gz ]]; then
        gunzip -c "${compatible_isoform}" > "${chr}_compatible_isoform.tsv"
        compat_file="${chr}_compatible_isoform.tsv"
    fi

    remove_nonfsm_data.pl \\
        \${compat_file} \\
        ${scratch_dir} \\
        ${tsv_updated}

    touch done_fsm_${chr}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | grep 'v[0-9]' | sed 's/.*v/v/' | sed 's/ .*//')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

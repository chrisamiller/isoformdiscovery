/*
 * ESPRESSO_Q
 *
 * ESPRESSO Phase 3: Quantification (per chromosome).
 *
 * Single module reused for all three quantification rounds via include aliases:
 *   ESPRESSO_Q_R0  — round 0 (first-pass against full reference GTF)
 *   ESPRESSO_Q_R2  — round 2 (re-quantification against SQANTI3-filtered GTF)
 *   ESPRESSO_Q_FSM — FSM round (quantification with FSM-filtered intermediate files)
 *
 * Resource labels and ext.args are overridden per alias in nextflow.config:
 *   R0:  label 'process_espresso_q'       (1 CPU, 300 GB, 96 h)
 *   R2:  label 'process_espresso_requant' (8 CPU, 200 GB, 72 h)
 *   FSM: label 'process_espresso_requant' (8 CPU, 200 GB, 72 h)
 *
 * The subdir parameter controls which subdirectory of the chromosome scratch
 * space to write into:
 *   ''         → round 0  (writes directly into scratch/chr/)
 *   'filtered' → round 2  (writes into scratch/chr/filtered/)
 *   'fsm'      → FSM      (writes into scratch/chr/fsm/)
 *
 * Input:
 *   tuple val(chr), path(sentinels)  — all ESPRESSO_C sentinels for this chr
 *                                       (grouped by chr); for rounds 2/FSM these
 *                                       are the round-0 Q sentinel + FSM-filter
 *                                       sentinel respectively
 *   path  gtf                        — annotation GTF (per-chr slice for r2/FSM,
 *                                       full reference for r0)
 *   val   subdir                     — output subdirectory: '', 'filtered', or 'fsm'
 *
 * Output:
 *   tuple val(chr), path("${chr}_abundance.esp"),          emit: abundance
 *   tuple val(chr), path("${chr}_updated.gtf"),            emit: gtf
 *   tuple val(chr), path("${chr}_compatible_isoform.tsv.gz"), emit: compatible_isoform
 *   tuple val(chr), path("${chr}_SJ_simplified.list"),     emit: sj_list   (r0 only)
 *   tuple val(chr), path("done_q_${chr}*.txt"),            emit: sentinel
 *   path "versions.yml",                                    emit: versions
 */

process ESPRESSO_Q {
    tag "${chr}${subdir ? '/' + subdir : ''}"
    label 'process_espresso_q'
    container 'sridnona/espresso:v2'

    input:
    tuple val(chr), path(sentinels)
    path  gtf
    val   subdir

    output:
    tuple val(chr), path("${chr}_abundance.esp"),          emit: abundance
    tuple val(chr), path("${chr}_updated.gtf"),            emit: gtf
    tuple val(chr), path("${chr}_compatible_isoform.tsv.gz"), emit: compatible_isoform
    tuple val(chr), path("${chr}_SJ_simplified.list"),     emit: sj_list,   optional: true
    tuple val(chr), path("done_q_${chr}${subdir ? '_' + subdir : ''}.txt"), emit: sentinel
    path  "versions.yml",                                   emit: versions

    script:
    def scratch_dir = "${params.espresso_scratch}/${chr}"
    def work_dir    = subdir ? "${scratch_dir}/${subdir}" : scratch_dir
    // For FSM round the tsv.updated is inside the fsm/ subdir (written by remove_nonfsm_data.pl)
    def tsv_updated = (subdir == 'fsm')
        ? "${work_dir}/${chr}.tsv.updated"
        : "${scratch_dir}/${chr}.tsv.updated"
    def compat_out  = "${work_dir}/${chr}_compatible_isoform.tsv"
    def out_flag    = subdir ? "-O ${work_dir}" : ""
    def extra_args  = task.ext.args ?: ''
    def sentinel_name = "done_q_${chr}${subdir ? '_' + subdir : ''}.txt"
    """
    export PATH="/opt/conda/envs/env/bin:\${PATH}"

    mkdir -p ${work_dir}

    perl /bin/espresso/src/ESPRESSO_Q.pl \\
        -L ${tsv_updated} \\
        -A ${gtf} \\
        -V ${compat_out} \\
        -T ${task.cpus} \\
        ${out_flag} \\
        ${extra_args}

    # Copy outputs from scratch into the Nextflow work dir for tracking
    cp \$(ls ${work_dir}/${chr}_*abundance*.esp  | head -n1) ${chr}_abundance.esp
    cp \$(ls ${work_dir}/${chr}_*updated*.gtf    | head -n1) ${chr}_updated.gtf
    cp ${work_dir}/${chr}_compatible_isoform.tsv  ${chr}_compatible_isoform.tsv
    gzip "${chr}_compatible_isoform.tsv"

    # SJ_simplified.list only exists after round-0 Q
    [ -f ${work_dir}/SJ_simplified.list ] && cp ${work_dir}/SJ_simplified.list ${chr}_SJ_simplified.list || true

    touch ${sentinel_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        espresso: \$(perl /bin/espresso/src/ESPRESSO_Q.pl --version 2>&1 | head -n1 || echo 'v2')
    END_VERSIONS
    """
}

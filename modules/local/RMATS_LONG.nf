/*
 * RMATS_LONG
 *
 * Differential isoform analysis between group1 and group2 samples using
 * rMATS-long on the FSM-only abundance and GTF files.
 *
 * Group membership is derived from the samplesheet 'group' column: g1 samples
 * form group1, g2 samples form group2.  Sample names are passed as newline-
 * delimited text files (process substitution) to work around argument-length
 * limits with large cohorts.
 *
 * rMATS-long is invoked via conda activation inside the container:
 *   source activate /docker_data/rMATS-long/conda_env && rmats_long.py ...
 *
 * Input:
 *   path abundance     — merged_fsm_N2_R2_abundance.esp
 *   path merged_gtf    — merged_fsm_N2_R2_updated.gtf
 *   path ref_gtf       — reference annotation GTF (Ensembl)
 *   val  g1_samples    — List of group1 sample names
 *   val  g2_samples    — List of group2 sample names
 *
 * Output:
 *   path "rmats_out/", emit: results
 *   path "versions.yml", emit: versions
 */

process RMATS_LONG {
    label 'process_rmats'
    container 'sridnona/rmats_long:v3'

    input:
    path abundance
    path merged_gtf
    path ref_gtf
    val  g1_samples
    val  g2_samples

    output:
    path "rmats_out/",   emit: results
    path "versions.yml", emit: versions

    script:
    // Join sample names comma-separated; read_sample_file() in the R script
    // expects line 1 to be a single comma-delimited string of sample names.
    def g1_str = g1_samples.join(',')
    def g2_str = g2_samples.join(',')
    """
    set +u  # conda activation scripts reference unbound vars (e.g. ADDR2LINE on ARM)
    source activate /docker_data/rMATS-long/conda_env
    set -u

    RMATS_LONG=/docker_data/rMATS-long/scripts/rmats_long.py

    # Write group membership files; each file: one line, comma-separated sample names
    printf '${g1_str}\n' > group1.txt
    printf '${g2_str}\n' > group2.txt

    mkdir -p rmats_out

    python \$RMATS_LONG \\
        --abundance ${abundance} \\
        --updated-gtf ${merged_gtf} \\
        --gencode-gtf ${ref_gtf} \\
        --group-1 group1.txt --group-1-name group1 \\
        --group-2 group2.txt --group-2-name group2 \\
        --out-dir rmats_out \\
        --num-threads ${task.cpus} \\
        --delta-proportion ${params.rmats_delta_proportion} \\
        --adj-pvalue ${params.rmats_adj_pvalue} \\
        --plot-file-type .pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmats_long: \$(python \$RMATS_LONG --version 2>&1 | head -n1 || echo 'v3')
    END_VERSIONS
    """
}

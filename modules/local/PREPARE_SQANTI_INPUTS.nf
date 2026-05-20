/*
 * PREPARE_SQANTI_INPUTS
 *
 * Renames STAR SJ.out.tab files into a numbered sequence (1.SJ.out.tab,
 * 2.SJ.out.tab, …) in a single staging directory for SQANTI3_QC.
 *
 * Background: SQANTI3 -c flag accepts a directory of SJ files.  When the
 * number of samples grows large, passing all paths as a comma-separated
 * argument string exceeds shell limits.  Renaming into a temp directory and
 * passing the directory path sidesteps this issue (matches what run.sh does
 * with the $scratch/tmp/ staging approach).
 *
 * Input:
 *   path sj_list  — text file with one STAR SJ.out.tab path per line
 *
 * Output:
 *   path "sj_dir/", emit: sj_dir
 *   path "versions.yml", emit: versions
 */

process PREPARE_SQANTI_INPUTS {
    label 'process_single'
    container 'chrisamiller/genomic-analysis:0.2'

    input:
    path sj_list

    output:
    path "sj_dir",       emit: sj_dir
    path "versions.yml", emit: versions

    script:
    """
    mkdir sj_dir
    count=1
    while IFS= read -r sj; do
        cp "\${sj}" "sj_dir/\${count}.SJ.out.tab"
        count=\$((count + 1))
    done < ${sj_list}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //' | sed 's/ .*//')
    END_VERSIONS
    """
}

/*
 * SQANTI3_QC
 *
 * Runs SQANTI3 transcript QC on the merged round-0 ESPRESSO GTF.
 *
 * Short-read splice junctions (STAR SJ.out.tab files, numbered in sj_dir/)
 * and short-read BAMs (star_bams_file, newline-delimited list) are used to
 * validate long-read transcript models.
 *
 * SQANTI3 is invoked via conda activation inside the container:
 *   source activate SQANTI3.env && sqanti3_qc.py ...
 *
 * Input:
 *   path merged_gtf      — merged_N2_R0_updated.gtf from MERGE_ESPRESSO_R0
 *   path ref_gtf         — reference annotation GTF (Ensembl)
 *   path fasta           — reference genome FASTA
 *   path sj_dir          — directory of numbered SJ.out.tab files
 *   path star_bams_file  — newline-delimited list of STAR BAM paths
 *
 * Output:
 *   path "*_classification.txt",   emit: classification
 *   path "*_corrected.gtf",        emit: corrected_gtf
 *   path "*_corrected.fasta",      emit: corrected_fasta,     optional: true
 *   path "*_corrected.faa",        emit: corrected_faa,       optional: true
 *   path "*_corrected.genePred",   emit: corrected_genepred,  optional: true
 *   path "*_corrected.gtf.cds.gff",emit: corrected_cds,       optional: true
 *   path "*_junctions.txt",        emit: junctions,           optional: true
 *   path "versions.yml",           emit: versions
 */

process SQANTI3_QC {
    label 'process_sqanti3'
    container 'chrisamiller/sqanti3:v5.2.1'

    input:
    path merged_gtf
    path ref_gtf
    path fasta
    path sj_dir
    path star_bams_file

    output:
    path "*_classification.txt",    emit: classification
    path "*_corrected.gtf",         emit: corrected_gtf
    path "*_corrected.fasta",       emit: corrected_fasta,    optional: true
    path "*_corrected.faa",         emit: corrected_faa,      optional: true
    path "*_corrected.genePred",    emit: corrected_genepred, optional: true
    path "*_corrected.gtf.cds.gff", emit: corrected_cds,      optional: true
    path "*_junctions.txt",         emit: junctions,          optional: true
    path "versions.yml",            emit: versions

    script:
    """
    export PATH="/opt/conda/envs/SQANTI3.env/bin:\${PATH}"

    /app/SQANTI3-5.2.1/sqanti3_qc.py \\
        -t ${task.cpus} \\
        -d sqanti_out \\
        -c ${sj_dir} \\
        --SR_bam ${star_bams_file} \\
        ${merged_gtf} \\
        ${ref_gtf} \\
        ${fasta}

    # Flatten outputs into the work directory
    mv sqanti_out/* ./ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(/app/SQANTI3-5.2.1/sqanti3_qc.py --version 2>&1 | head -n1 || echo 'v5.2.1')
    END_VERSIONS
    """
}

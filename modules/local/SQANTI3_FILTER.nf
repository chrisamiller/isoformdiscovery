/*
 * SQANTI3_FILTER
 *
 * Applies the SQANTI3 machine-learning filter to classify transcripts as
 * Artifact or non-Artifact based on SQANTI3 QC features.
 *
 * The ML filter uses a random-forest model trained on SQANTI3 features.
 * Transcripts classified as non-Artifact form the initial filtered transcript
 * set; Ensembl (ENST-prefixed) transcripts flagged as Artifact are later
 * restored by RESTORE_ENSEMBL.
 *
 * Input:
 *   path classification  — SQANTI3 *_classification.txt output
 *   path corrected_gtf   — SQANTI3 *_corrected.gtf output
 *
 * Output:
 *   path "*_MLresult_classification.txt", emit: ml_classification
 *   path "*.filtered.gtf",                emit: filtered_gtf
 *   path "*_inclusion-list.txt",          emit: inclusion_list,   optional: true
 *   path "*_TN_list.txt",                 emit: tn_list,          optional: true
 *   path "*_TP_list.txt",                 emit: tp_list,          optional: true
 *   path "classifier_variable-importance_table.txt", emit: var_importance, optional: true
 *   path "randomforest.RData",            emit: rf_model,         optional: true
 *   path "testSet_*",                     emit: test_set,         optional: true
 *   path "*_params.txt",                  emit: params_txt,       optional: true
 *   path "versions.yml",                  emit: versions
 */

process SQANTI3_FILTER {
    label 'process_sqanti3'
    container 'chrisamiller/sqanti3:v5.2.1'

    input:
    path classification
    path corrected_gtf

    output:
    path "*_MLresult_classification.txt",         emit: ml_classification
    path "*.filtered.gtf",                        emit: filtered_gtf
    path "*_inclusion-list.txt",                  emit: inclusion_list,  optional: true
    path "*_TN_list.txt",                         emit: tn_list,         optional: true
    path "*_TP_list.txt",                         emit: tp_list,         optional: true
    path "classifier_variable-importance_table.txt", emit: var_importance, optional: true
    path "randomforest.RData",                    emit: rf_model,        optional: true
    path "testSet_*",                             emit: test_set,        optional: true
    path "*_params.txt",                          emit: params_txt,      optional: true
    path "versions.yml",                          emit: versions

    script:
    """
    export PATH="/opt/conda/envs/SQANTI3.env/bin:\${PATH}"

    /app/SQANTI3-5.2.1/sqanti3_filter.py ml \\
        ${classification} \\
        --gtf ${corrected_gtf} \\
        -d filter_ml_default

    # Flatten output directory
    mv filter_ml_default/* ./ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(/app/SQANTI3-5.2.1/sqanti3_filter.py --version 2>&1 | head -n1 || echo 'v5.2.1')
    END_VERSIONS
    """
}

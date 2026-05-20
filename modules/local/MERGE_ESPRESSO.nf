/*
 * MERGE_ESPRESSO
 *
 * Merges per-chromosome ESPRESSO abundance and GTF files into single
 * genome-wide files.  Reused for two merge points:
 *
 *   MERGE_ESPRESSO_R0  — merges round-0 per-chr outputs into merged_N2_R0_*
 *   MERGE_ESPRESSO_FSM — merges FSM per-chr outputs into merged_fsm_N2_R2_*;
 *                        also extracts novel (ESPRESSO:) transcripts
 *
 * Abundance merge: header from the first file, data rows from all files (header
 * skipped for files 2–N).  Rows are sorted by transcript_id.
 *
 * GTF merge: header from the first file; data lines from all files are sorted
 * by chromosome (col 1) then start position (col 4); malformed entries where
 * end < start ($F[3] > $F[4] in 0-based coords) are dropped.
 *
 * The output prefix is controlled by ext.prefix in modules.config:
 *   MERGE_ESPRESSO_R0:  prefix = 'merged_N2_R0'
 *   MERGE_ESPRESSO_FSM: prefix = 'merged_fsm_N2_R2'
 *
 * Input:
 *   path(abundances) — all per-chr *_abundance.esp files (collected)
 *   path(gtfs)       — all per-chr *_updated.gtf files (collected)
 *
 * Output:
 *   path("*_abundance.esp"),        emit: abundance
 *   path("*_updated.gtf"),          emit: gtf
 *   path("*novel_transcripts.gtf"), emit: novel_gtf  (FSM merge only)
 *   path "versions.yml",            emit: versions
 */

process MERGE_ESPRESSO {
    label 'process_merge'
    container 'chrisamiller/genomic-analysis:0.2'

    input:
    path abundances
    path gtfs

    output:
    path "*_abundance.esp",         emit: abundance
    path "*_updated.gtf",           emit: gtf
    path "*novel_transcripts.gtf",  emit: novel_gtf,  optional: true
    path "versions.yml",            emit: versions

    script:
    def prefix = task.ext.prefix ?: 'merged_N2_R0'
    """
    # Capture input file lists before any output files are created, so that the
    # output globs (*_abundance.esp, *_updated.gtf) don't self-referentially pick
    # up the files being written.
    readarray -t _abund_files < <(ls *_abundance.esp | sort)
    readarray -t _gtf_files   < <(ls *_updated.gtf   | sort)

    # ── Merge abundance files ────────────────────────────────────────────────
    head -n1 "\${_abund_files[0]}" > ${prefix}_abundance.esp
    for f in "\${_abund_files[@]}"; do
        tail -n +2 "\$f"
    done | sort >> ${prefix}_abundance.esp

    # ── Merge GTF files ──────────────────────────────────────────────────────
    grep '^#' "\${_gtf_files[0]}" > ${prefix}_updated.gtf || true
    for f in "\${_gtf_files[@]}"; do
        grep -v '^#' "\$f" || true
    done | sort -k1,1 -k4,4n \\
         | perl -nae 'print unless \$F[3] > \$F[4]' \\
         >> ${prefix}_updated.gtf

    # ── Extract novel transcripts (ESPRESSO: prefix) ─────────────────────────
    if grep -q 'transcript_id "ESPRESSO:' ${prefix}_updated.gtf; then
        grep 'transcript_id "ESPRESSO:' ${prefix}_updated.gtf \\
            > ${prefix}_updated.novel_transcripts.gtf
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl -e 'print \$^V')
    END_VERSIONS
    """
}

/*
 * subworkflows/local/espresso_requant.nf
 *
 * Per-chromosome re-quantification after SQANTI3 filtering:
 *   1. SPLIT_FILTERED_GTF: extract chr-specific slice of the filtered GTF
 *   2. ESPRESSO_Q_R2:      re-quantify against filtered GTF (--read_ratio_cutoff 2)
 *   3. FSM_FILTER:         restrict intermediate files and BAMs to FSM reads only
 *   4. ESPRESSO_Q_FSM:     final quantification with FSM-filtered inputs
 *
 * The final FSM abundance/GTF files are collected by the main workflow and
 * merged into the genome-wide merged_fsm_N2_R2_* output files.
 *
 * Take:
 *   chr_sentinel_r0     — tuple(val(chr), path(done_q_r0))  from ESPRESSO_CHROMOSOME
 *   filtered_gtf        — path: genome-wide ensrestore.gtf from RESTORE_ENSEMBL
 *   fasta               — path: reference genome FASTA
 *
 * Emit:
 *   fsm_abundance           — tuple(val(chr), path(*_abundance.esp))
 *   fsm_gtf                 — tuple(val(chr), path(*_updated.gtf))
 *   fsm_compatible_isoform  — tuple(val(chr), path(*_compatible_isoform.tsv))
 *   fsm_sentinel            — tuple(val(chr), path(done_q_fsm_sentinel))
 *   versions
 */

include { SPLIT_FILTERED_GTF             } from '../../modules/local/SPLIT_FILTERED_GTF'
include { ESPRESSO_Q as ESPRESSO_Q_R2    } from '../../modules/local/ESPRESSO_Q'
include { ESPRESSO_Q as ESPRESSO_Q_FSM   } from '../../modules/local/ESPRESSO_Q'
include { FSM_FILTER                     } from '../../modules/local/FSM_FILTER'

workflow ESPRESSO_REQUANT {

    take:
    chr_sentinel_r0  // tuple(chr, done_q_r0_sentinel)
    filtered_gtf     // path (single file, broadcast to all chromosomes)
    fasta

    main:

    // ── Step 1: Split filtered GTF per chromosome ─────────────────────────
    // Broadcast the single filtered_gtf to every chr; extract chr-specific slice
    ch_split_input = chr_sentinel_r0
        .map { chr, _sentinel -> chr }
        .combine(filtered_gtf)  // tuple(chr, filtered_gtf_path)

    SPLIT_FILTERED_GTF(
        ch_split_input.map { chr, gtf -> gtf },   // path input
        ch_split_input.map { chr, gtf -> chr }    // val input
    )

    // ── Step 2: ESPRESSO_Q round 2 ──────────────────────────────────────────
    // Join r0 sentinels with chr-specific filtered GTF slices, keyed on chr
    ch_r2_sentinels = chr_sentinel_r0
        .join(SPLIT_FILTERED_GTF.out.chr_gtf)
        // -> tuple(chr, r0_sentinel, chr_filtered_gtf)

    ch_r2_input = ch_r2_sentinels.map { chr, sentinel, gtf ->
        tuple(chr, [sentinel])  // wrap sentinel in list for groupTuple-compatible input
    }

    ESPRESSO_Q_R2(
        ch_r2_input,
        ch_r2_sentinels.map { chr, _s, gtf -> gtf },
        'filtered'
    )

    // ── Step 3: FSM filtering ─────────────────────────────────────────────
    // Join round-2 compatible_isoform + round-2 sentinel per chr
    ch_fsm_input = ESPRESSO_Q_R2.out.compatible_isoform
        .join(ESPRESSO_Q_R2.out.sentinel)
        // -> tuple(chr, compat_isoform, r2_sentinel)

    FSM_FILTER(
        ch_fsm_input.map { chr, ci, _s -> tuple(chr, ci) },
        ch_fsm_input.map { _chr, _ci, s -> s }
    )

    // ── Step 4: ESPRESSO_Q FSM ───────────────────────────────────────────
    // FSM Q needs the same chr-filtered GTF as round 2 (the ensrestore GTF)
    ch_fsm_q_sentinels = FSM_FILTER.out.sentinel
        .join(SPLIT_FILTERED_GTF.out.chr_gtf)
        // -> tuple(chr, fsm_sentinel, chr_filtered_gtf)

    ch_fsm_q_input = ch_fsm_q_sentinels.map { chr, sentinel, _gtf ->
        tuple(chr, [sentinel])
    }

    ESPRESSO_Q_FSM(
        ch_fsm_q_input,
        ch_fsm_q_sentinels.map { _chr, _s, gtf -> gtf },
        'fsm'
    )

    emit:
    fsm_abundance           = ESPRESSO_Q_FSM.out.abundance
    fsm_gtf                 = ESPRESSO_Q_FSM.out.gtf
    fsm_compatible_isoform  = ESPRESSO_Q_FSM.out.compatible_isoform
    fsm_sentinel            = ESPRESSO_Q_FSM.out.sentinel
    versions = SPLIT_FILTERED_GTF.out.versions
        .mix(ESPRESSO_Q_R2.out.versions)
        .mix(FSM_FILTER.out.versions)
        .mix(ESPRESSO_Q_FSM.out.versions)
}

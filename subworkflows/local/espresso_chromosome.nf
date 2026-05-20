/*
 * subworkflows/local/espresso_chromosome.nf
 *
 * Per-chromosome first-pass ESPRESSO:  S → C (all samples in parallel) → Q round 0
 *
 * This subworkflow encapsulates the scatter-gather for one chromosome:
 *   1. ESPRESSO_S:  signal detection for all samples of this chr simultaneously
 *   2. ESPRESSO_C:  error correction, one job per (chr, sample_index)
 *   3. ESPRESSO_Q_R0: quantification (round 0), produces abundance + GTF +
 *                     compatible_isoform + SJ_simplified.list
 *
 * All ESPRESSO processes use the shared-scratch Option C design: they read/write
 * to params.espresso_scratch/${chr}/ directly and emit only small sentinel files
 * for Nextflow dependency tracking.  See ESPRESSO_S.nf and README for details.
 *
 * Take:
 *   chr_tsv    — channel of tuple(val(chr), path(tsv))
 *                  per-chromosome ESPRESSO input TSV (from BUILD_CHR_TSV)
 *   fasta      — val/path: reference genome FASTA
 *   gtf        — val/path: reference annotation GTF
 *   n_samples  — val: total number of samples (integer)
 *
 * Emit:
 *   chr_sentinel       — tuple(val(chr), path(done_q_sentinel))
 *   abundance          — tuple(val(chr), path(*_abundance.esp))
 *   gtf_r0             — tuple(val(chr), path(*_updated.gtf))
 *   compatible_isoform — tuple(val(chr), path(*_compatible_isoform.tsv))
 *   sj_list            — tuple(val(chr), path(*_SJ_simplified.list))
 *   versions
 */

include { ESPRESSO_S   } from '../../modules/local/ESPRESSO_S'
include { ESPRESSO_C   } from '../../modules/local/ESPRESSO_C'
include { ESPRESSO_Q as ESPRESSO_Q_R0 } from '../../modules/local/ESPRESSO_Q'

workflow ESPRESSO_CHROMOSOME {

    take:
    chr_tsv    // tuple(chr, tsv)
    fasta
    gtf
    n_samples  // integer

    main:

    // ── Phase 1: ESPRESSO_S ──────────────────────────────────────────────────
    ESPRESSO_S(chr_tsv, fasta, gtf)

    // ── Phase 2: ESPRESSO_C (one job per chr × sample_index) ────────────────
    // Expand each chr sentinel to (num_samples) (chr, idx, sentinel) tuples
    ch_c_input = ESPRESSO_S.out.sentinel
        .combine(n_samples)
        .flatMap { chr, sentinel, n ->
            (0..<n).collect { idx -> tuple(chr, idx, sentinel) }
        }

    ESPRESSO_C(ch_c_input, fasta)

    // ── Phase 3: ESPRESSO_Q round 0 ──────────────────────────────────────────
    // Group all per-sample C sentinels back by chr before running Q
    ch_q_input = ESPRESSO_C.out.sentinel
        .groupTuple()  // groups by chr (element 0)

    ESPRESSO_Q_R0(ch_q_input, gtf, '')

    emit:
    chr_sentinel       = ESPRESSO_Q_R0.out.sentinel
    abundance          = ESPRESSO_Q_R0.out.abundance
    gtf_r0             = ESPRESSO_Q_R0.out.gtf
    compatible_isoform = ESPRESSO_Q_R0.out.compatible_isoform
    sj_list            = ESPRESSO_Q_R0.out.sj_list
    versions           = ESPRESSO_S.out.versions
        .mix(ESPRESSO_C.out.versions)
        .mix(ESPRESSO_Q_R0.out.versions)
}

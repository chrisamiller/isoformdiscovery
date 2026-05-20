/*
 * workflows/transcript_assembly.nf
 *
 * Main workflow: Nanopore cDNA transcript discovery, quantification, and
 * differential isoform analysis.
 *
 * Data flow (matches background.md Pipeline Flow Summary):
 *
 *   INPUT_CHECK              → ch_samples: tuple(meta{group,sample}, bam)
 *   BAM_SPLIT_FILTER         → per-chr BAMs per sample
 *   BUILD_CHR_TSV            → per-chr ESPRESSO input TSV
 *   ESPRESSO_CHROMOSOME      → S → C → Q(r0) per chr
 *   MERGE_ESPRESSO_R0        → merged abundance + GTF
 *   PREPARE_SQANTI_INPUTS    → numbered SJ.out.tab directory
 *   SQANTI3_QC               → classification + corrected GTF
 *   SQANTI3_FILTER           → ML-filtered GTF + classification
 *   RESTORE_ENSEMBL          → final filtered GTF (ensrestore)
 *   ESPRESSO_REQUANT         → SPLIT_GTF → Q(r2) → FSM_FILTER → Q(FSM) per chr
 *   MERGE_ESPRESSO_FSM       → merged FSM abundance + GTF + novel transcript GTF
 *   FILTER_NORMALIZE         → coverage-filtered counts + CPM + TPM + filtered GTF
 *   RMATS_LONG               → differential isoform analysis (g1 vs g2)
 *   BAM_MERGE_INDEX          → per-sample indexed FSM BAMs
 */

include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { ESPRESSO_CHROMOSOME } from '../subworkflows/local/espresso_chromosome'
include { ESPRESSO_REQUANT   } from '../subworkflows/local/espresso_requant'

include { BAM_SPLIT_FILTER   } from '../modules/local/BAM_SPLIT_FILTER'
include { BUILD_CHR_TSV      } from '../modules/local/BUILD_CHR_TSV'
include { MERGE_ESPRESSO as MERGE_ESPRESSO_R0  } from '../modules/local/MERGE_ESPRESSO'
include { MERGE_ESPRESSO as MERGE_ESPRESSO_FSM } from '../modules/local/MERGE_ESPRESSO'
include { PREPARE_SQANTI_INPUTS } from '../modules/local/PREPARE_SQANTI_INPUTS'
include { SQANTI3_QC         } from '../modules/local/SQANTI3_QC'
include { SQANTI3_FILTER     } from '../modules/local/SQANTI3_FILTER'
include { RESTORE_ENSEMBL    } from '../modules/local/RESTORE_ENSEMBL'
include { FILTER_NORMALIZE   } from '../modules/local/FILTER_NORMALIZE'
include { RMATS_LONG         } from '../modules/local/RMATS_LONG'
include { BAM_MERGE_INDEX    } from '../modules/local/BAM_MERGE_INDEX'

workflow TRANSCRIPT_ASSEMBLY {

    take:
    ch_input         // path: samplesheet CSV
    ch_gtf           // path: reference annotation GTF
    ch_fasta         // path: reference genome FASTA
    ch_star_sjouts   // path: text file listing STAR SJ.out.tab paths
    ch_star_bams     // path: text file listing STAR BAM paths

    main:

    // ── Samplesheet validation ───────────────────────────────────────────────
    INPUT_CHECK(ch_input)
    ch_samples = INPUT_CHECK.out.reads   // tuple(meta{group,sample}, bam)

    // ── Phase 0: Split BAMs by chromosome ────────────────────────────────────
    BAM_SPLIT_FILTER(ch_samples)

    // Flatten: (meta, [sample.chr1.bam, ...]) → (chr, meta, bam) per bam.
    // chr is extracted from the filename suffix (<sample>.<chr>.bam).
    // Normalize to list: when only one chr BAM exists Nextflow emits a bare Path, not a List.
    ch_chr_bams_flat = BAM_SPLIT_FILTER.out.bams
        .flatMap { meta, bams ->
            def bamList = (bams instanceof List) ? bams : [bams]
            bamList.collect { bam ->
                def chr = bam.name.tokenize('.')[-2]  // <sample>.<chr>.bam → chr
                tuple(chr, meta, bam)
            }
        }

    // ── Build per-chr ESPRESSO input TSVs ────────────────────────────────────
    // Group all per-sample BAMs for each chromosome; chr is already element 0.
    // Result: tuple(chr, [meta...], [bam...])
    ch_chr_grouped = ch_chr_bams_flat
        .groupTuple()  // groups by chr (element 0) → tuple(chr, [metas], [bams])

    BUILD_CHR_TSV(ch_chr_grouped)   // emits tuple(chr, tsv_path)

    // Optional chromosome filter: --chromosomes 'chr1,chr2' restricts which
    // chromosomes flow through ESPRESSO.  null = process all detected chrs.
    ch_chr_tsv = params.chromosomes
        ? BUILD_CHR_TSV.out.tsv.filter { chr, _tsv ->
              params.chromosomes.tokenize(',').contains(chr) }
        : BUILD_CHR_TSV.out.tsv

    // Collect the list of chromosome names for later scatter steps
    ch_chromosomes = ch_chr_tsv.map { chr, _tsv -> chr }

    // Count samples (used for ESPRESSO_C index generation)
    ch_n_samples = ch_samples.count()

    // ── Phases 1–3: ESPRESSO per chromosome ──────────────────────────────────
    ESPRESSO_CHROMOSOME(
        ch_chr_tsv,
        ch_fasta,
        ch_gtf,
        ch_n_samples
    )

    // ── Phase 4a: Merge round-0 outputs ──────────────────────────────────────
    MERGE_ESPRESSO_R0(
        ESPRESSO_CHROMOSOME.out.abundance.map { _chr, f -> f }.collect(),
        ESPRESSO_CHROMOSOME.out.gtf_r0.map   { _chr, f -> f }.collect()
    )

    // ── Phase 4b: SQANTI3 QC and filtering ───────────────────────────────────
    PREPARE_SQANTI_INPUTS(ch_star_sjouts)

    SQANTI3_QC(
        MERGE_ESPRESSO_R0.out.gtf,
        ch_gtf,
        ch_fasta,
        PREPARE_SQANTI_INPUTS.out.sj_dir,
        ch_star_bams
    )

    SQANTI3_FILTER(
        SQANTI3_QC.out.classification,
        SQANTI3_QC.out.corrected_gtf
    )

    // Restore Ensembl transcripts flagged as Artifact by the ML filter
    RESTORE_ENSEMBL(
        SQANTI3_FILTER.out.ml_classification,
        MERGE_ESPRESSO_R0.out.gtf,
        SQANTI3_FILTER.out.filtered_gtf
    )

    ch_final_filtered_gtf = RESTORE_ENSEMBL.out.gtf  // single-value channel

    // ── Phases 5–7: Re-quantification per chromosome ─────────────────────────
    ESPRESSO_REQUANT(
        ESPRESSO_CHROMOSOME.out.chr_sentinel,
        ch_final_filtered_gtf,
        ch_fasta
    )

    // ── Phase 8: Merge FSM outputs ────────────────────────────────────────────
    MERGE_ESPRESSO_FSM(
        ESPRESSO_REQUANT.out.fsm_abundance.map { _chr, f -> f }.collect(),
        ESPRESSO_REQUANT.out.fsm_gtf.map       { _chr, f -> f }.collect()
    )

    // ── Filtering and normalization ───────────────────────────────────────────
    FILTER_NORMALIZE(
        MERGE_ESPRESSO_FSM.out.abundance,
        MERGE_ESPRESSO_FSM.out.gtf
    )

    // ── Phase 8b: rMATS-long differential isoform analysis ───────────────────
    ch_g1_samples = ch_samples
        .filter { meta, _bam -> meta.group == 'g1' }
        .map    { meta, _bam -> meta.sample }
        .collect()

    ch_g2_samples = ch_samples
        .filter { meta, _bam -> meta.group == 'g2' }
        .map    { meta, _bam -> meta.sample }
        .collect()

    RMATS_LONG(
        MERGE_ESPRESSO_FSM.out.abundance,
        MERGE_ESPRESSO_FSM.out.gtf,
        ch_gtf,
        ch_g1_samples,
        ch_g2_samples
    )

    // ── Phase 9: Per-sample FSM BAM assembly ─────────────────────────────────
    // Collect all FSM_FILTER sentinels per sample for merge ordering.
    // The FSM_FILTER sentinels in ESPRESSO_REQUANT are keyed by chr; we group
    // them per sample by joining against the original sample list.
    // Because BAM_MERGE_INDEX reads BAMs directly from scratch via glob, we
    // only need the sentinels — one entry per sample with all-chr sentinels.
    ch_sample_metas = ch_samples
        .map { meta, _bam -> meta.sample }
        .collect()

    // Build (meta{sample}, [all_fsm_sentinels]) for BAM_MERGE_INDEX.
    // combine() flattens a collected list into separate tuple elements, so we
    // re-group with map: [meta, s1, s2, ...] → [meta, [s1, s2, ...]]
    ch_merge_input = ch_samples
        .map { meta, _bam -> meta }
        .combine(
            ESPRESSO_REQUANT.out.fsm_sentinel
                .map { _chr, s -> s }
                .collect()
        )
        .map { items -> [items[0], items[1..-1]] }
        // → tuple(meta, [all_chr_fsm_sentinels])

    BAM_MERGE_INDEX(ch_merge_input)

    emit:
    fsm_abundance    = MERGE_ESPRESSO_FSM.out.abundance
    fsm_gtf          = MERGE_ESPRESSO_FSM.out.gtf
    novel_gtf        = MERGE_ESPRESSO_FSM.out.novel_gtf
    cpm              = FILTER_NORMALIZE.out.cpm
    rmats_results    = RMATS_LONG.out.results
    fsm_bams         = BAM_MERGE_INDEX.out.bam
}

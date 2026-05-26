# transcript_assembly

** WORK IN PROGRESS - check back often for updates **


Nextflow DSL2 pipeline for **Nanopore cDNA transcript discovery, quantification, and differential isoform analysis**.

The pipeline takes per-sample minimap2-aligned BAM files, runs ESPRESSO for transcript detection and quantification, filters the resulting transcript set with SQANTI3, re-quantifies against the filtered transcripts, normalizes to CPM, and reports differential isoform usage with rMATS-long.

---

## Table of Contents

1. [Quick start](#quick-start)
2. [Pipeline overview](#pipeline-overview)
3. [Inputs](#inputs)
   - [Samplesheet](#samplesheet)
   - [STAR auxiliary files](#star-auxiliary-files)
4. [Parameters](#parameters)
5. [Profiles](#profiles)
6. [Resource configuration](#resource-configuration)
7. [The Option C scratch design](#the-option-c-scratch-design)
8. [Outputs](#outputs)
9. [Software requirements](#software-requirements)
10. [Running on WUSTL RIS](#running-on-wustl-ris)
11. [Troubleshooting](#troubleshooting)

---

## Quick start

```bash
nextflow run /path/to/nextflow/ \
  --input      samplesheet.csv \
  --gtf        Homo_sapiens.GRCh38.95.gtf \
  --fasta      GRCh38.fa \
  --star_sjouts star_sjouts.txt \
  --star_bams   star_bams.txt \
  --outdir      results/ \
  --espresso_scratch /scratch/espresso_run1 \
  -profile ris
```

> **Tip:** Pass `-resume` on repeated runs to skip already-completed steps.

---

## Pipeline overview

```
INPUT_CHECK              validate samplesheet; emit (meta, bam) tuples
BAM_SPLIT_FILTER         split BAMs by chromosome; filter secondary/supplementary reads
BUILD_CHR_TSV            build per-chromosome ESPRESSO input TSV with absolute BAM paths
  ┌─────────────────── scatter: one job per chromosome ───────────────────────┐
  │ ESPRESSO_S            signal detection (all samples, one chr)             │
  │ ESPRESSO_C            error correction (one job per chr × sample index)   │
  │ ESPRESSO_Q_R0         round-0 quantification                              │
  └────────────────────────────────────────────────────────────────────────────┘
MERGE_ESPRESSO_R0        gather: merge per-chr abundance tables and GTFs
PREPARE_SQANTI_INPUTS    rename STAR SJ.out.tab files to SQANTI3's expected format
SQANTI3_QC               structural/quality classification of assembled transcripts
SQANTI3_FILTER           ML-based filtering (removes likely artifacts)
RESTORE_ENSEMBL          rescue Ensembl transcripts incorrectly flagged as Artifact
  ┌─────────────────── scatter: one job per chromosome ───────────────────────┐
  │ SPLIT_FILTERED_GTF    extract chr-specific slice of the filtered GTF      │
  │ ESPRESSO_Q_R2         re-quantify against filtered GTF (read_ratio_cutoff)│
  │ FSM_FILTER            restrict reads and BAMs to FSM-only transcripts     │
  │ ESPRESSO_Q_FSM        final quantification with FSM-filtered inputs       │
  └────────────────────────────────────────────────────────────────────────────┘
MERGE_ESPRESSO_FSM       gather: merge FSM abundance tables and GTFs
CPM_NORMALIZE            CPM-normalize the merged FSM abundance matrix
RMATS_LONG               differential isoform analysis (g1 vs g2)
BAM_MERGE_INDEX          assemble and index per-sample FSM BAMs
```

---

## Inputs

### Samplesheet

A CSV file with three required columns:

| Column   | Description                                         |
|----------|-----------------------------------------------------|
| `group`  | Comparison group: must be `g1` or `g2`              |
| `sample` | Unique sample identifier (used in all output names) |
| `bam`    | Absolute path to the minimap2-aligned BAM file      |

**Example** (`assets/samplesheet_example.csv`):

```csv
group,sample,bam
g1,AML_1,/storage1/data/aml_study/AML_1.bam
g1,AML_2,/storage1/data/aml_study/AML_2.bam
g1,AML_3,/storage1/data/aml_study/AML_3.bam
g2,HSC_1,/storage1/data/aml_study/HSC_1.bam
g2,HSC_2,/storage1/data/aml_study/HSC_2.bam
g2,HSC_3,/storage1/data/aml_study/HSC_3.bam
```

Notes:
- BAM files must be indexed (`.bai` or `.csi` alongside each `.bam`).
- Sample names must be unique across the entire samplesheet.
- At least one sample per group is required; rMATS-long is most meaningful with ≥ 3 per group.

### STAR auxiliary files

Two newline-delimited text files pointing to the outputs of a paired short-read STAR alignment run. These files are used by SQANTI3 to cross-validate splice junctions.

**`--star_sjouts`** — one absolute path per line to `SJ.out.tab` files:
```
/storage1/data/star/AML_1_SJ.out.tab
/storage1/data/star/AML_2_SJ.out.tab
...
```

**`--star_bams`** — one absolute path per line to STAR BAM/CRAM files:
```
/storage1/data/star/AML_1.Aligned.sortedByCoord.out.bam
/storage1/data/star/AML_2.Aligned.sortedByCoord.out.bam
...
```

The order of entries in these files does not need to match the samplesheet; SQANTI3 uses them collectively.

---

## Parameters

| Parameter                     | Default          | Description                                                                 |
|-------------------------------|------------------|-----------------------------------------------------------------------------|
| `--input`                     | *required*       | Path to samplesheet CSV                                                     |
| `--gtf`                       | *required*       | Reference annotation GTF (Ensembl format recommended)                       |
| `--fasta`                     | *required*       | Reference genome FASTA                                                      |
| `--star_sjouts`               | *required*       | Text file listing STAR SJ.out.tab paths (one per line)                      |
| `--star_bams`                 | *required*       | Text file listing STAR BAM/CRAM paths (one per line)                        |
| `--outdir`                    | `./results`      | Directory for all published outputs                                         |
| `--espresso_scratch`          | *required*       | Shared cluster filesystem path for ESPRESSO intermediate files              |
| `--chromosomes`               | `null`           | Comma-separated chromosome names to process; `null` = detect from BAM headers |
| `--espresso_read_ratio_cutoff`| `2`              | ESPRESSO `--read_ratio_cutoff` for rounds 2 and FSM                         |
| `--rmats_delta_proportion`    | `0.1`            | rMATS-long `--delta-proportion` threshold                                   |
| `--rmats_adj_pvalue`          | `0.1`            | rMATS-long adjusted p-value cutoff                                          |
| `--publish_dir_mode`          | `copy`           | Nextflow publish mode (`copy`, `symlink`, `link`, `move`)                   |
| `--max_memory`                | `320.GB`         | Hard cap on memory for any single process                                   |
| `--max_cpus`                  | `32`             | Hard cap on CPU count for any single process                                |
| `--max_time`                  | `240.h`          | Hard cap on walltime for any single process                                 |
| `--user_group`                | `null`           | LSF user group (RIS profile only), e.g. `compute-timley`                   |
| `--job_group_name`            | `null`           | LSF job group (RIS profile only), e.g. `/timley/nextflow`                  |
| `--queue_cpu`                 | `general`        | LSF queue name (RIS profile only)                                           |

---

## Profiles

Select a profile with `-profile <name>`.

### `standard` (default)

Runs all processes locally with no resource manager. Useful for small test runs on a workstation. Resource limits default to `max_memory`, `max_cpus`, and `max_time` — lower these if your machine has less RAM (ESPRESSO_Q requires up to 300 GB).

### `ris`

Submits jobs to the WUSTL RIS cluster via LSF with Docker. Requires:

```bash
-profile ris \
--user_group    compute-timley \
--job_group_name /timley/nextflow
```

Each process runs inside the appropriate Docker container (`-a 'docker(...)'` LSF option). See `conf/ris.config` for defaults specific to the RIS environment.

### `test`

Minimal parameter overrides for smoke-testing pipeline structure without real data:

```bash
nextflow run /path/to/nextflow/ -profile test,ris
```

---

## Resource configuration

Process resource labels are defined in `conf/base.config`. Override any label in `nextflow.config` or pass `--max_*` flags to adjust cluster-wide caps.

| Label                   | CPUs | Memory   | Walltime | Used by                          |
|-------------------------|------|----------|----------|----------------------------------|
| `process_espresso_s`    | 8    | 96 GB    | 48 h     | ESPRESSO_S                       |
| `process_espresso_c`    | 8    | 72 GB    | 24 h     | ESPRESSO_C                       |
| `process_espresso_q`    | 1    | 300 GB   | 96 h     | ESPRESSO_Q_R0                    |
| `process_espresso_requant` | 8 | 200 GB   | 72 h     | ESPRESSO_Q_R2, ESPRESSO_Q_FSM    |
| `process_bam_split`     | 4    | 16 GB    | 12 h     | BAM_SPLIT_FILTER                 |
| `process_merge`         | 2    | 32 GB    | 8 h      | MERGE_ESPRESSO_R0/FSM            |
| `process_fsm_filter`    | 4    | 32 GB    | 24 h     | FSM_FILTER                       |
| `process_sqanti3`       | 16   | 64 GB    | 24 h     | SQANTI3_QC, SQANTI3_FILTER       |
| `process_rmats`         | 8    | 64 GB    | 24 h     | RMATS_LONG                       |

> **ESPRESSO_Q memory note:** `ESPRESSO_Q_R0` requires ~300 GB for human whole-genome runs. Set `--max_memory 320.GB` or higher, or ensure your cluster queue has nodes with sufficient RAM.

---

## The Option C scratch design

ESPRESSO produces large intermediate files (tens to hundreds of GB per chromosome). Staging these through Nextflow's work directory would be impractical. Instead, this pipeline uses ESPRESSO's "Option C" design: all ESPRESSO processes read and write directly to a **shared scratch directory** on the cluster filesystem, and Nextflow tracks dependencies via small **sentinel files** (e.g., `done_s_chr1.txt`).

### How it works

1. `BAM_SPLIT_FILTER` splits per-sample BAMs into per-chromosome BAMs and publishes them to `--outdir`.
2. `BUILD_CHR_TSV` writes a per-chromosome TSV listing absolute BAM paths — no file staging needed.
3. `ESPRESSO_S` writes its outputs to `${espresso_scratch}/${chr}/` and emits a sentinel file.
4. `ESPRESSO_C` (one job per chr × sample index) reads from and writes to the same scratch subdir.
5. `ESPRESSO_Q_R0` reads from scratch, copies its output files back to the Nextflow work directory for tracking, then emits abundance/GTF as normal Nextflow outputs.
6. Re-quantification writes to `${espresso_scratch}/${chr}/filtered/` and `${espresso_scratch}/${chr}/fsm/`.

### Implications

- `--espresso_scratch` **must** be on a POSIX filesystem accessible from all compute nodes simultaneously (e.g., a network-mounted cluster scratch volume). It cannot be a local `/tmp`.
- Each `ESPRESSO_S` run **deletes** its scratch subdir before starting — this ensures idempotent re-runs. Do **not** store anything in `${espresso_scratch}/${chr}/` that you cannot regenerate.
- The scratch directory is **not cleaned automatically** after the pipeline completes. After verifying outputs, remove it with:
  ```bash
  rm -rf /path/to/espresso_scratch
  ```
- If the pipeline fails midway, re-run with `-resume`. Nextflow will skip completed processes; any ESPRESSO step that re-runs will wipe its own scratch subdir first.

---

## Outputs

All outputs are published to `--outdir` (default: `results/`). The pipeline also writes timeline, report, trace, and DAG files to `results/pipeline_info/`.

```
results/
├── pipeline_info/
│   ├── execution_timeline.html
│   ├── execution_report.html
│   ├── execution_trace.txt
│   └── pipeline_dag.svg
├── bams/
│   ├── chr/                         # per-chr, per-sample split BAMs
│   └── merged/                      # per-sample merged FSM BAMs + indices
├── espresso/
│   ├── r0/
│   │   ├── merged_N2_R0_abundance.esp      # round-0 merged abundance matrix
│   │   └── merged_N2_R0_updated.gtf        # round-0 merged transcript GTF
│   └── fsm/
│       ├── merged_fsm_N2_R2_abundance.esp  # FSM-filtered abundance matrix
│       ├── merged_fsm_N2_R2_updated.gtf    # FSM-filtered transcript GTF
│       └── merged_fsm_N2_R2_updated.novel_transcripts.gtf  # novel-only GTF
├── sqanti3/
│   ├── qc/                          # SQANTI3_QC outputs (classification, corrected GTF, etc.)
│   └── filter/                      # SQANTI3_FILTER ML outputs + ensrestore GTF
├── cpm/
│   └── *.cpm.tsv                    # CPM-normalized count matrix
└── rmats_long/
    └── *.tsv                        # rMATS-long differential isoform results
```

### Key output files

| File | Description |
|------|-------------|
| `espresso/fsm/merged_fsm_N2_R2_abundance.esp` | Final per-transcript read count matrix (N samples × M transcripts) |
| `espresso/fsm/merged_fsm_N2_R2_updated.gtf` | All retained transcripts (FSM + Ensembl-restored) |
| `espresso/fsm/merged_fsm_N2_R2_updated.novel_transcripts.gtf` | Subset: novel ESPRESSO-detected transcripts only |
| `cpm/*.cpm.tsv` | CPM-normalized version of the FSM abundance matrix |
| `rmats_long/*.tsv` | Differential isoform usage: g1 vs g2 |
| `bams/merged/*.bam` | Per-sample BAMs containing only FSM-assigned reads |

---

## Software requirements

All tools run inside Docker containers. No local installation is needed beyond Nextflow itself.

| Container | Tools |
|-----------|-------|
| `sridnona/espresso:v2` | ESPRESSO_S, ESPRESSO_C, ESPRESSO_Q |
| `chrisamiller/genomic-analysis:0.2` | samtools (BAM split, merge, index), R + edgeR (CPM) |
| `chrisamiller/docker-genomic-analysis:latest` | Perl (FSM_FILTER, RESTORE_ENSEMBL scripts) |
| `sridnona/rmats_long:v3` | rMATS-long (conda env at `/docker_data/rMATS-long/conda_env`) |
| `python:3.11-slim` | samplesheet validation |
| SQANTI3 conda env | SQANTI3 QC and filter (activated inside `chrisamiller/genomic-analysis:0.2`) |

**Nextflow version:** ≥ 23.04.0

---

## Running on Local LSF cluster

### 1. Set up a job group

```bash
bgadd -L 200 /timley/nextflow
```

This caps concurrent LSF jobs to 200 (Nextflow's `queueSize` default).

### 2. Prepare a run script

```bash
#!/bin/bash
#BSUB -G compute-timley
#BSUB -g /timley/nextflow
#BSUB -q general
#BSUB -n 2
#BSUB -M 8GB
#BSUB -W 240:00
#BSUB -o nextflow_%J.log
#BSUB -a 'docker(nextflow/nextflow:23.10.1)'

nextflow run /path/to/nextflow/ \
  --input      /storage1/data/project/samplesheet.csv \
  --gtf        /storage1/ref/Homo_sapiens.GRCh38.95.gtf \
  --fasta      /storage1/ref/GRCh38.fa \
  --star_sjouts /storage1/data/project/star_sjouts.txt \
  --star_bams   /storage1/data/project/star_bams.txt \
  --outdir      /storage1/data/project/results \
  --espresso_scratch /scratch1/timley/espresso_$(date +%Y%m%d) \
  --user_group    compute-timley \
  --job_group_name /timley/nextflow \
  -profile ris \
  -resume
```

### 3. Submit

```bash
bsub < run_pipeline.sh
```

### 4. Monitor

```bash
# Watch LSF job queue
bjobs -g /timley/nextflow

# Watch Nextflow log
tail -f .nextflow.log

# Check process trace
tail -f results/pipeline_info/execution_trace.txt
```

### Notes for running on an LSF cluster:

- The Nextflow head job itself should request at least 8 GB RAM — it holds all channel state in memory.
- Use a dedicated scratch volume (not `/tmp`) for `--espresso_scratch` - it must be persistent across jobs.
- If jobs fail with out-of-memory errors, increase the relevant resource label in `conf/base.config` or override via `--max_memory`.
- SQANTI3 processes activate a conda environment inside the container (`source activate SQANTI3.env`). If the environment path changes in a new image, update `SQANTI3_QC.nf` and `SQANTI3_FILTER.nf`.

---

## Troubleshooting

### Pipeline stalls at ESPRESSO_C

ESPRESSO_C jobs are one per chromosome × sample index, so a run with 24 chromosomes and 6 samples submits 144 jobs. If your job group limit is too low, jobs queue but the pipeline appears stalled. Increase the limit:

```bash
bgmod -L 500 /timley/nextflow
```

### ESPRESSO_Q_R0 fails with out-of-memory

Round-0 quantification of the human genome requires ~300 GB. Confirm that `--max_memory` is at least `320.GB` and that your LSF queue has nodes with sufficient RAM. The `general` queue on RIS supports up to 500 GB nodes.

### "No such variable: chr" in BAM_SPLIT_FILTER

The chromosome names extracted from @SQ headers are matched by the regex `^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y)$`. If your BAM uses non-standard contig names (e.g., `1` vs `chr1`, or `chrM`), check the regex in `modules/local/BAM_SPLIT_FILTER.nf` and adjust as needed.

### SQANTI3 finishes but produces no transcripts

This usually means the merged round-0 GTF has chromosome names that don't match the reference GTF or FASTA. Confirm that ESPRESSO and the reference use consistent naming (both `chr`-prefixed or both not).

### rMATS-long exits with conda environment error

The rMATS-long container activates its environment with `source activate /docker_data/rMATS-long/conda_env`. If the container image changes, verify the conda env path inside the container:

```bash
docker run sridnona/rmats_long:v3 ls /docker_data/rMATS-long/
```

Update `modules/local/RMATS_LONG.nf` if the path has changed.

### Resuming after scratch deletion

If `--espresso_scratch` was deleted and you want to resume from an earlier successful step, you cannot resume through any ESPRESSO step — those processes read from scratch. You must re-run from the beginning (or from the last ESPRESSO merge step if those outputs are still in `--outdir`). Always confirm scratch outputs before deleting them.

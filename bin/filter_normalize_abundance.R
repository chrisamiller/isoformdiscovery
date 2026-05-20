#!/usr/bin/env Rscript
#
# filter_normalize_abundance.R — Filter and normalize ESPRESSO FSM abundance.
#
# Usage:
#   filter_normalize_abundance.R <abundance.esp> <gtf> <min_count>
#
# Arguments:
#   [1] abundance.esp  — ESPRESSO merged FSM abundance file (tab-delimited)
#   [2] gtf            — corresponding merged GTF (for transcript lengths / TPM)
#   [3] min_count      — minimum total read count to retain a transcript
#
# Outputs (written to the current directory using the input basename):
#   <base>.covfilt.noism.esp      filtered raw count matrix
#   <base>.covfilt.noism.cpm.esp  TMM-normalized CPM (edgeR)
#   <base>.covfilt.noism.tpm.esp  TPM (transcript-length normalized)
#
# Coverage filter: transcripts with rowSum < min_count are removed.  The "noism"
# label reflects that the FSM re-quantification round already excludes ISM reads
# by design; this step adds the coverage threshold on top.

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: filter_normalize_abundance.R <abundance.esp> <gtf> <min_count>",
         call. = FALSE)
}

input_file <- args[1]
gtf_file   <- args[2]
min_count  <- as.integer(args[3])
prefix     <- sub("\\.esp$", "", basename(input_file))

a      <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)
a.meta <- a[, 1:3]
a.data <- a[, 4:ncol(a)]

# Coverage filter
keep   <- rowSums(a.data) >= min_count
a.meta <- a.meta[keep, , drop = FALSE]
a.data <- a.data[keep, , drop = FALSE]

# CPM (TMM-normalized via edgeR)
y          <- DGEList(counts = a.data)
y          <- calcNormFactors(y)
counts.cpm <- as.data.frame(cpm(y, normalized.lib.sizes = TRUE))

# TPM: derive transcript lengths from GTF exon records
gtf      <- read.table(gtf_file, sep = "\t", header = FALSE, comment.char = "#",
                        stringsAsFactors = FALSE, quote = "")
exons    <- gtf[gtf[, 3] == "exon", ]
get_tid  <- function(attr) {
    m <- regmatches(attr, regexpr('transcript_id "[^"]+"', attr))
    if (length(m) == 0) return(NA_character_)
    sub('^transcript_id "', "", sub('"$', "", m))
}
exons$tid <- vapply(exons[, 9], get_tid, character(1))
exons$len <- exons[, 5] - exons[, 4] + 1
tx_len    <- tapply(exons$len, exons$tid, sum)

t_ids      <- as.character(a.meta[, 1])
med_len    <- median(tx_len, na.rm = TRUE)
lengths    <- ifelse(t_ids %in% names(tx_len), tx_len[t_ids], med_len)

rpk        <- sweep(a.data, 1, lengths / 1000, FUN = "/")
counts.tpm <- as.data.frame(sweep(rpk, 2, colSums(rpk) / 1e6, FUN = "/"))

write.table(cbind(a.meta, a.data),     paste0(prefix, ".covfilt.noism.esp"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(a.meta, counts.cpm), paste0(prefix, ".covfilt.noism.cpm.esp"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(cbind(a.meta, counts.tpm), paste0(prefix, ".covfilt.noism.tpm.esp"),
            sep = "\t", row.names = FALSE, quote = FALSE)

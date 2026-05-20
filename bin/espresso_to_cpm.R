#!/usr/bin/env Rscript
#
# espresso_to_cpm.R — Convert ESPRESSO abundance counts to TMM-normalized CPM.
#
# Usage:
#   espresso_to_cpm.R <abundance.esp> <output.cpm.tsv>
#
# Arguments:
#   [1] abundance.esp  — ESPRESSO abundance file (tab-delimited); first three columns
#                        are transcript_id, gene_name, support_count_total; remaining
#                        columns are per-sample integer counts.
#   [2] output.cpm.tsv — output path for CPM matrix; same structure as input but
#                        count columns replaced with CPM values.
#
# Normalization: TMM (trimmed mean of M-values) via edgeR calcNormFactors(),
# applied to the full count matrix before CPM calculation.
#
# Called by CPM_NORMALIZE module as the final post-processing step on the
# merged FSM abundance file.

suppressPackageStartupMessages(library(edgeR))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: espresso_to_cpm.R <abundance.esp> <output.cpm.tsv>", call. = FALSE)
}
input_file  <- args[1]
output_file <- args[2]

a      <- read.table(input_file, sep = "\t", header = TRUE, check.names = FALSE)
a.meta <- a[, 1:3]
a.data <- a[, 4:ncol(a)]

y              <- DGEList(counts = a.data)
y              <- calcNormFactors(y)
counts.cpm     <- cpm(y, normalized.lib.sizes = TRUE)

write.table(
    cbind(a.meta, counts.cpm),
    output_file,
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
)

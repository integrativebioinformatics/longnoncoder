#!/usr/bin/env Rscript

# Convert an aligned BAM into a genome-wide coverage track in bigWig format.
#
# Two reasons this step exists rather than reading BAMs directly downstream:
#
#   1. plotgardener's plotSignal does not accept BAM. Its `data` argument takes a
#      bigWig file path, a data frame in BED format, or a GRanges carrying a
#      `score` column -- nothing else.
#   2. Coverage is two to three orders of magnitude smaller than the alignments it
#      derives from. A 30 GB ONT BAM becomes roughly 100-200 MB, which is what
#      makes the tracks publishable and usable away from the cluster.

suppressPackageStartupMessages({
  library(optparse)
  library(Rsamtools)
  library(GenomicAlignments)
  library(rtracklayer)
})

option_list <- list(
  make_option(c("-b", "--bam"), type = "character", default = NULL,
              help = "Indexed BAM file", metavar = "character"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
              help = "Output basename; writes <prefix>.bw", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$bam)) {
  stop("--bam is required.", call. = FALSE)
}

prefix <- if (!is.null(opt$prefix)) opt$prefix else sub("\\.bam$", "", basename(opt$bam))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Nextflow stages the index alongside the BAM. Fail with something readable if that
# wiring ever breaks, rather than letting Rsamtools report it obliquely.
if (!file.exists(paste0(opt$bam, ".bai")) &&
    !file.exists(sub("\\.bam$", ".bai", opt$bam))) {
  stop(sprintf("No index found next to %s. BAM_COVERAGE requires an indexed BAM.",
               opt$bam), call. = FALSE)
}

bam <- Rsamtools::BamFile(opt$bam)

cat("Computing coverage for", basename(opt$bam), "\n")

# coverage() returns an RleList: run-length encoded, so memory scales with the
# number of coverage *changes*, not with genome length or read count. That is why
# a 30 GB BAM does not need anything like 30 GB of RAM here.
cov <- GenomicAlignments::coverage(bam)

out_bw <- file.path(opt$outdir, paste0(prefix, ".bw"))
rtracklayer::export.bw(cov, out_bw)

cat(sprintf("Wrote %s (%d contigs, %.1f MB)\n",
            out_bw, length(cov), file.info(out_bw)$size / 1024^2))

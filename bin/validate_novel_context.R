#!/usr/bin/env Rscript

#' Structural evidence for sense-intronic novel transcripts, from the alignments
#' they were built from.
#'
#' One question: could this model be a fragment of its host's unspliced precursor
#' rather than a transcript in its own right?
#'
#' Only sense-intronic candidates are evaluated, because only for them is that
#' question open. An antisense model cannot come from the host's precursor whatever
#' its reads do, and for a model sharing splice structure with the host, carrying
#' the host's junctions is what its class code already says. Measuring the rest and
#' reporting it invites a reading the measurement does not support.
#'
#' Two measurements:
#'
#'   Overrun    How far do supporting reads extend past the model's ends? A discrete
#'              transcript has reads stopping at its boundaries; a fragment of
#'              something longer has reads running straight through. Reported as the
#'              distribution of distances, separately for the 5' and 3' end.
#'
#'   Junctions  Do supporting reads carry splice junctions belonging to the host
#'              rather than to the candidate? A read that does came from the
#'              molecule that has them.
#'
#' Nothing here is a verdict. Every output is a measurement reported per transcript
#' so that disagreement between the signals stays visible; collapsing them into one
#' score would hide exactly the cases worth looking at.

suppressPackageStartupMessages({
    library(optparse)
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(Rsamtools)
    library(GenomeInfoDb)
})

option_list <- list(
    make_option(c("--metadata"), type = "character", default = NULL,
                help = "Novel transcript metadata CSV", metavar = "character"),
    make_option(c("--annotation"), type = "character", default = NULL,
                help = "Reference annotation GTF", metavar = "character"),
    make_option(c("--bams"), type = "character", default = NULL,
                help = "Comma-separated indexed BAM files", metavar = "character"),
    make_option(c("--junction_tolerance"), type = "integer", default = 10L,
                help = "Bases of slack when matching a read junction to a host intron [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metadata) || is.null(opt$annotation) || is.null(opt$bams)) {
    stop("--metadata, --annotation and --bams are all required.", call. = FALSE)
}

bam_files <- trimws(strsplit(opt$bams, ",", fixed = TRUE)[[1]])
bam_files <- bam_files[nzchar(bam_files)]
if (!length(bam_files)) stop("No BAM files supplied.", call. = FALSE)

TOLERANCE <- opt$junction_tolerance

# ---------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------

cat("Loading novel transcript metadata...\n")
meta <- read.csv(opt$metadata, stringsAsFactors = FALSE)

# An empty candidate set is a legitimate outcome on a small or highly complete
# annotation. Write the empty outputs and stop rather than failing the run.
write_empty <- function() {
    cols <- c("qry_id", "class_code", "strand", "ref_gene_id", "ref_gene_biotype",
              "ref_gene_strand", "same_strand_as_host", "reads_total",
              "reads_same_strand", "median_overrun_5p", "q90_overrun_5p",
              "median_overrun_3p", "q90_overrun_3p",
              "reads_with_host_junction", "reads_spliced_into_host_exon",
              "frac_with_host_junction")
    empty <- as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols))
    write.csv(empty, "novel_context_flags.csv", row.names = FALSE)
    write.csv(data.frame(metric = character(0), value = character(0)),
              "novel_context_summary.csv", row.names = FALSE)
}

if (!nrow(meta)) {
    cat("No novel transcripts to evaluate.\n")
    write_empty()
    quit(save = "no", status = 0)
}

# `cand` is built after the sense-intronic filter below, which needs the host
# strand, which needs the annotation. Nothing between here and there uses it.

# ---------------------------------------------------------------------------
# Reference structure
# ---------------------------------------------------------------------------

cat("Loading reference annotation...\n")
ref <- import(opt$annotation, feature.type = c("exon"))

ref_gene_id <- if (!is.null(ref$gene_id)) as.character(ref$gene_id) else rep(NA_character_, length(ref))

# Host gene extent and strand, keyed on gene id. The metadata already carries
# ref_gene_id from gffcompare's reference match, so the annotation is only
# needed for the geometry.
cat("Deriving host intron structure...\n")
gene_exons <- split(ref, ref_gene_id)
host_exons <- reduce(gene_exons)

# range() yields more than one entry for the rare gene whose exons span two
# seqnames or strands; keeping the first makes the gene-to-range mapping 1:1,
# which the psetdiff below depends on.
host_gene <- unlist(range(gene_exons), use.names = TRUE)
host_gene <- host_gene[!duplicated(names(host_gene))]

# Host introns, per gene: the gaps between that gene's exons after merging
# overlapping exons from all its isoforms. A read junction matching one of these
# is evidence the read came from the host rather than from the candidate.
host_introns <- psetdiff(host_gene, host_exons[names(host_gene)])

#' Harmonise sequence naming between the BAM and the annotation.
#'
#' Ensembl writes "1", GENCODE writes "chr1", and BAM headers follow whichever
#' reference the alignment used. A mismatch produces zero overlaps and no error,
#' so this is checked rather than assumed.
match_seqlevels <- function(gr, target_levels) {
    if (any(seqlevels(gr) %in% target_levels)) return(gr)
    alt <- gr
    if (all(grepl("^chr", target_levels[1]))) {
        seqlevels(alt) <- paste0("chr", seqlevels(alt))
    } else {
        seqlevels(alt) <- sub("^chr", "", seqlevels(alt))
    }
    if (any(seqlevels(alt) %in% target_levels)) alt else gr
}

bam_levels <- names(scanBamHeader(bam_files[1])[[1]]$targets)

# ---------------------------------------------------------------------------
# Strand versus host, then restrict to sense intronic
# ---------------------------------------------------------------------------

cat("Comparing candidate strand with host gene strand...\n")

host_id     <- as.character(meta$ref_gene_id)
host_id[is.na(host_id) | host_id == "-" | !nzchar(host_id)] <- NA_character_

host_strand <- rep(NA_character_, nrow(meta))
known       <- !is.na(host_id) & host_id %in% names(host_gene)
host_strand[known] <- as.character(strand(host_gene[host_id[known]]))

# NA rather than FALSE where there is no host: an intergenic transcript has no
# host to agree or disagree with, and calling that FALSE would put it in the same
# bucket as a genuine antisense call.
same_strand <- ifelse(is.na(host_strand), NA,
                      host_strand == as.character(meta$strand))

# Sense intronic: inside a reference intron (class code i) and on that gene's own
# strand. Everything else is dropped rather than measured, for the reason in the
# header -- the precursor question these tests answer is not open for it.
#
# Derived from class_code and the strand comparison rather than from the
# classification wording, so a change to the label does not silently change which
# transcripts are evaluated.
sense_intronic <- !is.na(meta$class_code) & meta$class_code == "i" &
                  !is.na(same_strand) & same_strand

cat(sprintf("%d of %d novel transcripts are sense intronic; the rest are not evaluated\n",
            sum(sense_intronic), nrow(meta)))

meta        <- meta[sense_intronic, , drop = FALSE]
host_id     <- host_id[sense_intronic]
host_strand <- host_strand[sense_intronic]
same_strand <- same_strand[sense_intronic]

if (!nrow(meta)) {
    cat("No sense-intronic candidates to evaluate.\n")
    write_empty()
    quit(save = "no", status = 0)
}

cand <- GRanges(
    seqnames = meta$seqnames,
    ranges   = IRanges(start = meta$start, end = meta$end),
    strand   = meta$strand
)
names(cand) <- meta$qry_id
cand <- match_seqlevels(cand, bam_levels)

if (!any(seqlevels(cand) %in% bam_levels)) {
    stop("No sequence names shared between the candidates and the BAM header. ",
         "The alignment and the annotation appear to use different references.",
         call. = FALSE)
}

# ---------------------------------------------------------------------------
# Boundary overrun and host junctions
# ---------------------------------------------------------------------------

cat("Querying alignments in", length(cand), "sense-intronic regions across",
    length(bam_files), "BAM file(s)...\n")

n <- length(cand)
reads_total        <- integer(n)
reads_same_strand  <- integer(n)
reads_host_junc    <- integer(n)
reads_into_exon    <- integer(n)

# Overruns are accumulated per read rather than summed, because a median cannot be
# computed from a running total. One entry per (read, candidate) overlap; on the
# sense-intronic subset that is a small fraction of what the full candidate set
# would have produced.
ov_ci <- integer(0)
ov_5p <- integer(0)
ov_3p <- integer(0)

# Host introns and exons restricted to the candidates' hosts, in BAM naming.
host_introns_flat <- unlist(host_introns, use.names = TRUE)
host_exons_flat   <- unlist(host_exons,   use.names = TRUE)
if (length(host_introns_flat)) host_introns_flat <- match_seqlevels(host_introns_flat, bam_levels)
if (length(host_exons_flat))   host_exons_flat   <- match_seqlevels(host_exons_flat,   bam_levels)

for (bam in bam_files) {
    cat("  ", basename(bam), "\n")

    param <- ScanBamParam(which = cand, what = "qname",
                          flag = scanBamFlag(isSecondaryAlignment = FALSE,
                                             isSupplementaryAlignment = FALSE))
    ga <- readGAlignments(bam, param = param, use.names = TRUE)
    if (!length(ga)) next

    # readGAlignments returns the same record once per region it was queried
    # under, so a read spanning two candidates arrives twice. Deduplicate on the
    # read's identity and alignment rather than on the alignment alone: two
    # distinct reads can align identically, and collapsing those would silently
    # undercount support.
    ga <- ga[!duplicated(paste(names(ga), start(ga), cigar(ga), sep = ":"))]
    gr <- granges(ga)

    hits <- findOverlaps(gr, cand, ignore.strand = TRUE)
    if (!length(hits)) next

    ri <- queryHits(hits)
    ci <- subjectHits(hits)

    reads_total <- reads_total + tabulate(ci, nbins = n)

    on_strand <- as.character(strand(gr))[ri] == as.character(strand(cand))[ci]
    reads_same_strand <- reads_same_strand + tabulate(ci[on_strand], nbins = n)

    # How far each read runs past each end of the candidate, in bases.
    #
    # Clamped at zero: a read that stops inside the candidate has not overrun it,
    # and a negative distance would pull a median below zero and read as though the
    # boundary had been respected by more than it was.
    #
    # Split by transcript orientation rather than by genomic left/right, because 5'
    # and 3' overrun mean different things. Reverse transcription falls short at the
    # 5' end and the poly(A) tail anchors the 3', so the two ends have different
    # noise floors and pooling them hides that.
    left_ov  <- pmax(0L, start(cand)[ci] - start(gr)[ri])
    right_ov <- pmax(0L, end(gr)[ri]     - end(cand)[ci])
    on_plus  <- as.character(strand(cand))[ci] != "-"

    ov_ci <- c(ov_ci, ci)
    ov_5p <- c(ov_5p, ifelse(on_plus, left_ov,  right_ov))
    ov_3p <- c(ov_3p, ifelse(on_plus, right_ov, left_ov))

    # Junction tests. junctions() returns the N gaps of each read; a read with no
    # gaps contributes nothing to either test.
    juncs <- junctions(ga)
    jlen  <- lengths(juncs)

    if (any(jlen > 0) && length(host_introns_flat)) {
        jflat <- unlist(juncs, use.names = FALSE)
        jread <- rep(seq_along(juncs), jlen)

        # A read junction that coincides with one of the host's introns came from
        # the host's transcript, not from a candidate that lacks that junction.
        carrying <- integer(0)
        host_hit <- findOverlaps(jflat, host_introns_flat, type = "equal",
                                 maxgap = TOLERANCE, ignore.strand = TRUE)
        if (length(host_hit)) {
            carrying <- unique(jread[queryHits(host_hit)])
            reads_host_junc <- reads_host_junc + tabulate(ci[ri %in% carrying], nbins = n)
        }

        # A junction landing in one of the host's exons means the candidate is
        # spliced into the host -- an unannotated host exon, which is a finding
        # rather than an artifact. Reads already counted as carrying a host
        # intron are excluded so the two outcomes stay distinguishable.
        if (length(host_exons_flat)) {
            into <- findOverlaps(jflat, host_exons_flat, ignore.strand = TRUE)
            if (length(into)) {
                spliced <- setdiff(unique(jread[queryHits(into)]), carrying)
                reads_into_exon <- reads_into_exon + tabulate(ci[ri %in% spliced], nbins = n)
            }
        }
    }
}

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

safe_frac <- function(num, den) ifelse(den > 0, num / den, NA_real_)

#' Summarise one overrun vector per candidate.
#'
#' The median says where the bulk of the reads stop. It is zero whenever most reads
#' respect the boundary, which is the answer for a discrete transcript -- so the 90th
#' percentile is reported beside it. A candidate where a third of the reads run
#' through by a kilobase has a median of zero and a very large q90, and the median
#' alone would hide it.
#'
#' NA where a candidate recovered no reads, which is a different statement from an
#' overrun of zero.
overrun_stats <- function(v) {
    if (!length(ov_ci)) {
        return(list(median = rep(NA_real_, n), q90 = rep(NA_real_, n)))
    }
    f <- factor(ov_ci, levels = seq_len(n))
    list(
        median = as.numeric(tapply(v, f, stats::median)),
        q90    = as.numeric(tapply(v, f, function(x)
                     stats::quantile(x, 0.9, names = FALSE, type = 7)))
    )
}

ov5 <- overrun_stats(ov_5p)
ov3 <- overrun_stats(ov_3p)

flags <- data.frame(
    qry_id                       = meta$qry_id,
    class_code                   = meta$class_code,
    strand                       = as.character(meta$strand),
    ref_gene_id                 = host_id,
    ref_gene_biotype            = meta$ref_gene_biotype,
    ref_gene_strand             = host_strand,
    same_strand_as_host          = same_strand,
    reads_total                  = reads_total,
    reads_same_strand            = reads_same_strand,
    median_overrun_5p            = ov5$median,
    q90_overrun_5p               = ov5$q90,
    median_overrun_3p            = ov3$median,
    q90_overrun_3p               = ov3$q90,
    reads_with_host_junction     = reads_host_junc,
    reads_spliced_into_host_exon = reads_into_exon,
    frac_with_host_junction      = safe_frac(reads_host_junc, reads_total),
    stringsAsFactors             = FALSE
)

write.csv(flags, "novel_context_flags.csv", row.names = FALSE)

# Summary for the report. Every row here is sense intronic by construction, so the
# old same-strand and intronic breakdowns would all equal the total; the counts that
# remain are the ones that still vary.
#
# The overrun quartiles describe the population rather than thresholding it: a run
# where the median 3' overrun is a few bases across the board is a different result
# from one where a quarter of candidates overrun by hundreds.
with_reads <- flags$reads_total > 0

q <- function(v, p) {
    v <- v[with_reads & !is.na(v)]
    if (!length(v)) return(NA_real_)
    round(stats::quantile(v, p, names = FALSE, type = 7), 1)
}

summary_df <- data.frame(
    metric = c(
        "sense_intronic_evaluated",
        "no_reads_recovered",
        "any_read_with_host_junction",
        "any_read_spliced_into_host_exon",
        "median_overrun_5p_q50",
        "median_overrun_5p_q75",
        "median_overrun_3p_q50",
        "median_overrun_3p_q75"
    ),
    value = c(
        nrow(flags),
        sum(!with_reads),
        sum(flags$reads_with_host_junction > 0),
        sum(flags$reads_spliced_into_host_exon > 0),
        q(flags$median_overrun_5p, 0.50),
        q(flags$median_overrun_5p, 0.75),
        q(flags$median_overrun_3p, 0.50),
        q(flags$median_overrun_3p, 0.75)
    ),
    stringsAsFactors = FALSE
)

write.csv(summary_df, "novel_context_summary.csv", row.names = FALSE)

cat("\nDone.\n")
print(summary_df, row.names = FALSE)

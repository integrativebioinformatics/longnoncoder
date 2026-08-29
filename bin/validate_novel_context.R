#!/usr/bin/env Rscript

#' Structural evidence for novel transcripts, from the alignments they were built
#' from.
#'
#' Two questions, both answered from structure alone so that neither needs a
#' threshold tuned to a particular genome:
#'
#'   A  Is the novel transcript on the same strand as the gene whose intron it
#'      lies in? On the opposite strand the host's pre-mRNA cannot explain it, so
#'      the ambiguity does not arise.
#'
#'   D  Do the reads supporting it stop where it stops, and do they carry the
#'      host's splice junctions? A discrete transcript has reads ending at its
#'      boundaries. A fragment of something longer has reads running straight
#'      through, and reads that carry junctions the candidate does not have came
#'      from the molecule that does.
#'
#' Nothing here is a verdict. Every output is a count or a flag, reported per
#' transcript so that disagreement between the signals stays visible; collapsing
#' them into one score would hide exactly the cases worth looking at.

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
    make_option(c("--boundary_margin"), type = "integer", default = 25L,
                help = "Bases a read must extend past a boundary to count as read-through [%default]"),
    make_option(c("--junction_tolerance"), type = "integer", default = 5L,
                help = "Bases of slack when matching a read junction to a host intron [%default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$metadata) || is.null(opt$annotation) || is.null(opt$bams)) {
    stop("--metadata, --annotation and --bams are all required.", call. = FALSE)
}

bam_files <- trimws(strsplit(opt$bams, ",", fixed = TRUE)[[1]])
bam_files <- bam_files[nzchar(bam_files)]
if (!length(bam_files)) stop("No BAM files supplied.", call. = FALSE)

MARGIN    <- opt$boundary_margin
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
              "reads_same_strand", "reads_crossing_boundary",
              "reads_with_host_junction", "reads_spliced_into_host_exon",
              "frac_crossing_boundary", "frac_with_host_junction")
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

cand <- GRanges(
    seqnames = meta$seqnames,
    ranges   = IRanges(start = meta$start, end = meta$end),
    strand   = meta$strand
)
names(cand) <- meta$qry_id

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
cand       <- match_seqlevels(cand, bam_levels)

if (!any(seqlevels(cand) %in% bam_levels)) {
    stop("No sequence names shared between the candidates and the BAM header. ",
         "The alignment and the annotation appear to use different references.",
         call. = FALSE)
}

# ---------------------------------------------------------------------------
# Test A -- strand versus host
# ---------------------------------------------------------------------------

cat("Test A: comparing candidate strand with host gene strand...\n")

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

# ---------------------------------------------------------------------------
# Test D -- boundary read-through and host junctions
# ---------------------------------------------------------------------------

cat("Test D: querying alignments in", length(cand), "candidate regions across",
    length(bam_files), "BAM file(s)...\n")

n <- length(cand)
reads_total        <- integer(n)
reads_same_strand  <- integer(n)
reads_crossing     <- integer(n)
reads_host_junc    <- integer(n)
reads_into_exon    <- integer(n)

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

    # Boundary read-through: the read extends materially past either end of the
    # candidate. MARGIN absorbs the few bases of soft-clip and alignment wobble
    # that would otherwise make almost every read look like read-through.
    crosses <- start(gr)[ri] < (start(cand)[ci] - MARGIN) |
               end(gr)[ri]   > (end(cand)[ci]   + MARGIN)
    reads_crossing <- reads_crossing + tabulate(ci[crosses], nbins = n)

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
    reads_crossing_boundary      = reads_crossing,
    reads_with_host_junction     = reads_host_junc,
    reads_spliced_into_host_exon = reads_into_exon,
    frac_crossing_boundary       = safe_frac(reads_crossing,  reads_total),
    frac_with_host_junction      = safe_frac(reads_host_junc, reads_total),
    stringsAsFactors             = FALSE
)

write.csv(flags, "novel_context_flags.csv", row.names = FALSE)

# Summary for the report. Intronic same-strand candidates are the population the
# whole test exists for, so they are counted separately from the rest.
intronic_same <- flags$class_code == "i" & !is.na(flags$same_strand_as_host) &
                 flags$same_strand_as_host

summary_df <- data.frame(
    metric = c(
        "candidates_evaluated",
        "with_host_gene",
        "same_strand_as_host",
        "opposite_strand_to_host",
        "intronic_same_strand",
        "any_read_crossing_boundary",
        "any_read_with_host_junction",
        "any_read_spliced_into_host_exon",
        "intronic_same_strand_with_host_junction",
        "no_reads_recovered"
    ),
    value = c(
        nrow(flags),
        sum(!is.na(flags$same_strand_as_host)),
        sum(flags$same_strand_as_host %in% TRUE),
        sum(flags$same_strand_as_host %in% FALSE),
        sum(intronic_same),
        sum(flags$reads_crossing_boundary > 0),
        sum(flags$reads_with_host_junction > 0),
        sum(flags$reads_spliced_into_host_exon > 0),
        sum(intronic_same & flags$reads_with_host_junction > 0),
        sum(flags$reads_total == 0)
    ),
    stringsAsFactors = FALSE
)

write.csv(summary_df, "novel_context_summary.csv", row.names = FALSE)

cat("\nDone.\n")
print(summary_df, row.names = FALSE)

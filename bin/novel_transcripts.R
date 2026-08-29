#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(rtracklayer)
    library(GenomicRanges)
    library(optparse)

})

# Shared reference-GTF helpers, staged alongside this script by the module
source("gtf_annotation_utils.R")

# Define command line options
option_list <- list(
    make_option(c("--bambu_gtf"), type="character", default=NULL,
                help="Path to bambu novel transcripts GTF file", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character"),
    make_option(c("--compared_gtf"), type="character", default=NULL,
                help="Path to compared transcriptome annotated GTF file", metavar="character"),
    make_option(c("--tmap_file"), type="character", default=NULL,
                help="Path to tmap results file", metavar="character"),
    make_option(c("--coding_predictions"), type="character", default=NULL,
                help="Path to the coding-potential table (CPC2 or RNAmining)", metavar="character"),
    make_option(c("--tx_counts"), type="character", default=NULL,
                help="Path to transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default=NULL,
                help="Path to gene counts file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$bambu_gtf) || is.null(opt$compared_gtf) || is.null(opt$tmap_file) ||
    is.null(opt$coding_predictions) || is.null(opt$tx_counts) || is.null(opt$gene_counts) ||
    is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

# Reference lookups, used to say what each novel transcript was compared against:
# the biotype and name of the reference gene, and of the reference transcript
# itself. Kept as four named vectors rather than the whole table, which is several
# hundred megabytes for a full GENCODE annotation.
reference <- read_reference_gtf(opt$annotation)
ref_gene_biotype <- reference$gene_biotype
ref_gene_name    <- reference$gene_name
ref_tx_biotype   <- reference$tx_biotype
ref_tx_name      <- reference$tx_name
rm(reference)

# Import GTF files
cat("Loading GTF files...\n")
gtf <- import(opt$bambu_gtf)
tx_gtf <- import(opt$compared_gtf)

tx_table <- as.data.frame(tx_gtf[tx_gtf$type == "transcript"], )

# Load tmap results
cat("Loading tmap results...\n")
tmap <- read_table(opt$tmap_file)

cat("Loading coding-potential predictions...\n")

#' Read a coding-potential table from either predictor.
#'
#' The two formats differ in more than layout. RNAmining writes a preamble and then
#' <id> <coding|non-coding> <score>, where the score is the probability of whichever
#' class it chose. CPC2 writes a "#ID" header with eight tab-separated columns,
#' spells the label "noncoding" without the hyphen, and its coding_probability is
#' the probability of CODING specifically -- so a non-coding row carries a low value
#' where RNAmining carries a high one. Passing either straight through would put two
#' different quantities in the same column.
#'
#' Both are normalised here: prediction in {coding, non-coding}, and coding_prob as
#' P(coding) whichever tool produced it. CPC2 columns are located by header name
#' rather than position, because it inserts ORF_Start when run with --ORF.
read_coding_predictions <- function(path) {
    lines <- readLines(path)
    hdr   <- grep("^#ID\t", lines)

    if (length(hdr)) {
        cols <- strsplit(sub("^#", "", lines[hdr[1]]), "\t", fixed = TRUE)[[1]]
        body <- lines[seq.int(hdr[1] + 1L, length(lines))]
        body <- body[nzchar(body)]
        if (!length(body)) {
            stop("CPC2 table ", path, " has a header but no rows.", call. = FALSE)
        }

        d <- read.table(text = body, sep = "\t", header = FALSE, quote = "",
                        comment.char = "", stringsAsFactors = FALSE)
        if (ncol(d) != length(cols)) {
            stop("CPC2 table ", path, " has ", ncol(d), " columns against a ",
                 length(cols), "-column header.", call. = FALSE)
        }
        names(d) <- cols

        missing <- setdiff(c("ID", "label", "coding_probability"), names(d))
        if (length(missing)) {
            stop("CPC2 table ", path, " is missing column(s): ",
                 paste(missing, collapse = ", "), ". Found: ",
                 paste(names(d), collapse = ", "), call. = FALSE)
        }

        cat(sprintf("Read %d CPC2 predictions\n", nrow(d)))
        return(data.frame(
            transcript_id    = as.character(d$ID),
            prediction       = ifelse(d$label == "coding", "coding", "non-coding"),
            coding_prob      = as.numeric(d$coding_probability),
            coding_predictor = "cpc2",
            stringsAsFactors = FALSE))
    }

    # RNAmining. The preamble line count has changed between versions, and skipping a
    # fixed number silently drops one transcript's call whenever that count is wrong,
    # so the body is identified by matching it instead.
    is_data <- grepl("^[^#[:space:]][^\t]*\t(coding|non-coding)\t", lines)
    if (!any(is_data)) {
        stop("No prediction rows found in ", path, ". Expected either a CPC2 table ",
             "with a #ID header, or RNAmining <id> <coding|non-coding> <score> lines.",
             call. = FALSE)
    }

    d <- read.table(text = lines[is_data], header = FALSE, sep = "\t",
                    stringsAsFactors = FALSE)
    colnames(d) <- c("transcript_id", "prediction", "score")

    cat(sprintf("Read %d RNAmining predictions (%d preamble lines skipped)\n",
                nrow(d), sum(!is_data)))
    data.frame(
        transcript_id    = as.character(d$transcript_id),
        prediction       = as.character(d$prediction),
        coding_prob      = ifelse(d$prediction == "coding", d$score, 1 - d$score),
        coding_predictor = "rnamining",
        stringsAsFactors = FALSE)
}

rnam <- read_coding_predictions(opt$coding_predictions)

# Select relevant information from gtf
tx_table <- dplyr::select(tx_table, seqnames, transcript_id, gene_name, start, end, strand)

# Complete info dataframe
cat("Merging data...\n")
tx_info <- merge(tmap, tx_table, by.x="qry_id", by.y="transcript_id", all.x=TRUE)

# Attach the coding-potential calls
tx_info <- merge(tx_info, rnam, by.x="qry_id", by.y="transcript_id", all.x=TRUE)

# Select relevant info and reorder columns
tx_info <- dplyr::select(tx_info, seqnames, qry_id, ref_id, qry_gene_id, ref_gene_id, gene_name, 
                  class_code, strand, start, end, len, num_exons, prediction, coding_prob,
                  coding_predictor)

# Remove unstranded transcripts
tx_info <- tx_info[tx_info$strand != "*", ]

tx_info$exon <- ifelse(tx_info$num_exons == 1, 'mono-exonic', 'multi-exonic')

#' Reference gene id with gffcompare's placeholders normalised to NA. "-" means
#' no reference match, which is a different statement from a missing value.
normalise_ref_gene <- function(ids) {
    ids <- as.character(ids)
    ids[is.na(ids) | ids == "-" | !nzchar(ids)] <- NA_character_
    ids
}

# What the novel transcript was compared against, recorded on every novel record
# rather than only in the GTF.
#
# "ref" rather than "host": a host is something a transcript sits inside, which is
# true for i, x, m and n but backwards for k and y, where the novel transcript
# CONTAINS the reference. There is no separate host id -- it was a normalised copy
# of ref_gene_id and nothing more -- so the placeholder is cleaned in place instead.
tx_info$ref_gene_id <- normalise_ref_gene(tx_info$ref_gene_id)
tx_info$ref_id      <- normalise_ref_gene(tx_info$ref_id)

tx_info$ref_gene_biotype <- unname(ref_gene_biotype[tx_info$ref_gene_id])
tx_info$ref_gene_name    <- unname(ref_gene_name[tx_info$ref_gene_id])

# Transcript level as well as gene level, because the two disagree in exactly the
# cases worth looking at. A y-class novel transcript containing TARDBP-221 has
# ref_gene_biotype "protein_coding" and ref_transcript_biotype
# "nonsense_mediated_decay", and only the second says what it actually contains.
tx_info$ref_transcript_biotype <- unname(ref_tx_biotype[tx_info$ref_id])
tx_info$ref_transcript_name    <- unname(ref_tx_name[tx_info$ref_id])

# Load counts
cat("Filtering transcripts with zero counts...\n")
tx_counts <- read_table(opt$tx_counts)
gene_counts <- read_table(opt$gene_counts)

# Remove novel transcripts that had 0 transcript counts in all samples
tx_info <- tx_info[tx_info$qry_id %in% tx_counts$TXNAME, ]

# Save the metadata of all novel transcripts
cat("Saving novel transcripts metadata...\n")
write.csv(tx_info, file="novel_transcripts_metadata.csv", row.names = FALSE)

#' Human-readable label for a gffcompare class code. Kept in one place because
#' three outputs need the same mapping and copies of it drifted apart.
CLASS_LABELS <- c(
    u = 'intergenic',
    i = 'intronic',
    x = 'antisense',
    j = 'multiexon SJ match',
    k = 'extends reference',
    o = 'same-strand overlap',
    y = 'reference within intron',
    m = 'total intron retention',
    n = 'partial intron retention'
)

classify_class_code <- function(codes) unname(CLASS_LABELS[as.character(codes)])

#' Build the attribute set written into the novel GTFs.
#'
#' transcript_biotype is carried on the metadata rather than derived here, so the
#' combined GTF keeps each transcript's own category when the three subsets are
#' bound together. At gene level, a novel transcript arising from a known gene
#' keeps that gene's real biotype; only genuinely new loci, which gffcompare
#' reports without a reference gene, are labelled "novel".
novel_attrs <- function(meta) {
    gene_biotype <- meta$ref_gene_biotype
    gene_biotype[is.na(gene_biotype)] <- "novel"

    list(
        transcript_status      = rep("novel", nrow(meta)),
        transcript_biotype     = as.character(meta$transcript_biotype),
        gene_biotype           = gene_biotype,
        class_code             = as.character(meta$class_code),
        classification         = as.character(meta$classification),
        # What it was compared against, carried on the GTF as well as the CSV so the
        # file is self-describing in IGV without the metadata beside it.
        ref_gene_id            = meta$ref_gene_id,
        ref_gene_name          = as.character(meta$ref_gene_name),
        ref_gene_biotype       = as.character(meta$ref_gene_biotype),
        ref_transcript_id      = as.character(meta$ref_id),
        ref_transcript_name    = as.character(meta$ref_transcript_name),
        ref_transcript_biotype = as.character(meta$ref_transcript_biotype),
        gene_name              = as.character(meta$gene_name)
    )
}

#' Subset the Bambu GTF to a set of transcripts, attach the novel attributes and
#' write it out. Returns the exon records, which the exon-length summaries reuse.
write_novel_gtf <- function(meta, path) {
    ids   <- meta$qry_id
    tx    <- subset(gtf, type == "transcript" & transcript_id %in% ids)
    exons <- subset(gtf, type == "exon" & transcript_id %in% ids)

    out <- c(tx, exons)
    out <- out[order(out$transcript_id, out$type == "transcript", decreasing = TRUE)]
    out <- annotate_gtf(out, ids, novel_attrs(meta))
    export(out, path)

    invisible(exons)
}

# Which class codes are eligible to become candidates, and which get their own
# category instead.
#
# k, o and y are in because each can reveal something the others cannot: k means
# the novel model extends past the reference it contains, which is where an
# unannotated exon would show up; y means the model contains a reference inside its
# own intron, which can be a distinct locus rather than an isoform; o is
# same-strand overlap that matches no junction. Excluding them was not a judgement
# that they are uninteresting, only that they had not been looked at.
#
# m and n stay separate on GENCODE's terms: retained_intron is carried alongside
# protein_coding and lncRNA rather than beneath either, and "retained intron"
# describes the structure precisely where "lncRNA" would assert more than the data
# supports. The coding prediction is kept as a column on them, not used as the key.
#
# Note = and c never appear at all: only novel transcripts reach gffcompare, so a
# model matching or contained by a reference was already resolved as annotated.
CANDIDATE_CODES        <- c('u', 'i', 'x', 'j', 'k', 'o', 'y')
INTRON_RETENTION_CODES <- c('m', 'n')

# Select novel lncRNA candidates
cat("Processing lncRNA candidates...\n")
new_lncRNAs <- tx_info[tx_info$class_code %in% CANDIDATE_CODES &
                    tx_info$len >= 200 &
                    tx_info$prediction == 'non-coding', ]
new_lncRNAs$classification     <- classify_class_code(new_lncRNAs$class_code)
new_lncRNAs$transcript_biotype <- "novel_lncRNA"

write.csv(new_lncRNAs, file="novel_lncRNAs_metadata.csv", row.names = FALSE)
new_lncRNAs_exons_gtf <- write_novel_gtf(new_lncRNAs, "novel_lncRNAs.gtf")

# Select novel protein-coding candidates
cat("Processing protein-coding candidates...\n")
new_mRNAs <- tx_info[tx_info$class_code %in% CANDIDATE_CODES &
                         tx_info$len >= 200 &
                         tx_info$prediction == 'coding', ]
new_mRNAs$classification     <- classify_class_code(new_mRNAs$class_code)
new_mRNAs$transcript_biotype <- "novel_protein_coding"

write.csv(new_mRNAs, file="novel_protein-coding_metadata.csv", row.names = FALSE)
new_mRNAs_exons_gtf <- write_novel_gtf(new_mRNAs, "novel_protein-coding.gtf")

# Select intron-retention events. Not split by the coding prediction: it stays as
# a column so nothing is lost, but it is not the routing key here. novel_retained_intron
# mirrors GENCODE, which carries retained_intron alongside protein_coding and
# lncRNA rather than beneath either.
cat("Processing intron retention events...\n")
new_intron_retention <- tx_info[tx_info$class_code %in% INTRON_RETENTION_CODES &
                                    tx_info$len >= 200, ]
new_intron_retention$classification     <- classify_class_code(new_intron_retention$class_code)
new_intron_retention$transcript_biotype <- "novel_retained_intron"

write.csv(new_intron_retention, file="novel_intron_retention_metadata.csv", row.names = FALSE)
write_novel_gtf(new_intron_retention, "novel_intron_retention.gtf")

# Save combined metadata. All three categories: taking m and n out of the lncRNA
# branch must not drop them from the validated GTF and the counts derived from it.
new_mRNA_lncRNA <- rbind(new_mRNAs, new_lncRNAs, new_intron_retention)
write.csv(new_mRNA_lncRNA, "novel_pc_lnc_RNAs_metadata.csv", row.names = FALSE)

# Update novel transcripts GTF with the new classifications
write_novel_gtf(new_mRNA_lncRNA, "novel_transcripts_validated.gtf")

# Get exon lengths
cat("Calculating exon lengths...\n")
new_lncRNA_exon_len <- data.frame(
  transcript_id = new_lncRNAs_exons_gtf$transcript_id,
  exon_number = new_lncRNAs_exons_gtf$exon_number,
  width = width(new_lncRNAs_exons_gtf)
) 

write.csv(new_lncRNA_exon_len, "novel_lncRNA_exon_lengths.csv", row.names = FALSE)

new_mRNA_exon_len <- data.frame(
  transcript_id = new_mRNAs_exons_gtf$transcript_id,
  exon_number = new_mRNAs_exons_gtf$exon_number,
  width = width(new_mRNAs_exons_gtf)
) 

write.csv(new_mRNA_exon_len, "novel_protein-coding_exon_lengths.csv", row.names = FALSE)

# Save counts for novel mRNAs and lncRNAs
tx_ids <- new_mRNA_lncRNA$qry_id
gn_ids <- new_mRNA_lncRNA$qry_gene_id

novel_tx_counts <- subset(tx_counts, TXNAME %in% tx_ids)
write.csv(novel_tx_counts, "bambu_novel_pc_lnc_RNA_tx_counts.csv", row.names = FALSE)

novel_gn_counts <- subset(gene_counts, GENEID %in% gn_ids)
write.csv(novel_gn_counts, "bambu_novel_pc_lnc_RNA_gene_counts.csv", row.names = FALSE)

cat("Analysis completed successfully!\n")
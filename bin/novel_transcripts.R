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
    make_option(c("--rnamining_predictions"), type="character", default=NULL,
                help="Path to rnamining predictions file", metavar="character"),
    make_option(c("--tx_counts"), type="character", default=NULL,
                help="Path to transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default=NULL,
                help="Path to gene counts file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$bambu_gtf) || is.null(opt$compared_gtf) || is.null(opt$tmap_file) ||
    is.null(opt$rnamining_predictions) || is.null(opt$tx_counts) || is.null(opt$gene_counts) ||
    is.null(opt$annotation)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

# Reference gene biotypes, used to give novel transcripts arising from a known
# gene that gene's real biotype rather than a generic label
reference <- read_reference_gtf(opt$annotation)
ref_gene_biotype <- reference$gene_biotype
rm(reference)

# Import GTF files
cat("Loading GTF files...\n")
gtf <- import(opt$bambu_gtf)
tx_gtf <- import(opt$compared_gtf)

tx_table <- as.data.frame(tx_gtf[tx_gtf$type == "transcript"], )

# Load tmap results
cat("Loading tmap results...\n")
tmap <- read_table(opt$tmap_file)

# Load rnamining prediction results
cat("Loading rnamining predictions...\n")
rnamres <- readLines(opt$rnamining_predictions)
rnamres <- rnamres[6:length(rnamres)]

rnam <- read.table(text = rnamres, header = FALSE, sep = "\t")
colnames(rnam) <- c("transcript_id", "prediction", "rnamining_score")

# Select relevant information from gtf
tx_table <- dplyr::select(tx_table, seqnames, transcript_id, gene_name, start, end, strand)

# Complete info dataframe
cat("Merging data...\n")
tx_info <- merge(tmap, tx_table, by.x="qry_id", by.y="transcript_id", all.x=TRUE)

# Add the rnamining results
tx_info <- merge(tx_info, rnam, by.x="qry_id", by.y="transcript_id", all.x=TRUE)

# Select relevant info and reorder columns
tx_info <- dplyr::select(tx_info, seqnames, qry_id, ref_id, qry_gene_id, ref_gene_id, gene_name, 
                  class_code, strand, start, end, len, num_exons, prediction, rnamining_score)

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

# Biotype of the gene a novel transcript arises from, recorded on every novel
# record rather than only in the GTF. For intronic, antisense and retained-intron
# transcripts the host's biotype is the difference between "another isoform of a
# known lncRNA" and "a novel transcript at a protein-coding locus", and hosts are
# frequently neither -- so it cannot be assumed. Genuinely new loci have no host
# and keep NA.
tx_info$host_gene_id      <- normalise_ref_gene(tx_info$ref_gene_id)
tx_info$host_gene_biotype <- unname(ref_gene_biotype[tx_info$host_gene_id])

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
    gene_biotype <- meta$host_gene_biotype
    gene_biotype[is.na(gene_biotype)] <- "novel"

    list(
        transcript_status  = rep("novel", nrow(meta)),
        transcript_biotype = as.character(meta$transcript_biotype),
        gene_biotype       = gene_biotype,
        class_code         = as.character(meta$class_code),
        classification     = as.character(meta$classification),
        ref_gene_id        = meta$host_gene_id,
        gene_name          = as.character(meta$gene_name)
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

# Novel transcripts are sorted by coding potential, but that prediction only
# carries information for some class codes. RNAmining scores length-normalised
# trinucleotide composition; a retained intron is sequence that has never been
# under coding selection, so it pulls the composition toward non-coding in
# proportion to the share of the transcript it occupies. An m or n transcript is
# therefore predicted non-coding largely by construction, and routing it down the
# coding/non-coding branch would file it as a lncRNA candidate on the strength of
# a determination its own structure already made. These get their own category
# instead, with the prediction kept as a column rather than used as the key.
CANDIDATE_CODES        <- c('u', 'i', 'x', 'j')
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
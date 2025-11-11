#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(rtracklayer)
    library(GenomicRanges)
    library(optparse)
})

# Define command line options
option_list <- list(
    make_option(c("--bambu_gtf"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/bambu_novel_transcripts.gtf',
                help="Path to bambu novel transcripts GTF file", metavar="character"),
    make_option(c("--compared_gtf"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/bambu_novel_compared_transcriptome.annotated.gtf',
                help="Path to compared transcriptome annotated GTF file", metavar="character"),
    make_option(c("--tmap_file"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/bambu_novel_compared_transcriptome.bambu_novel_transcripts.gtf.tmap',
                help="Path to tmap results file", metavar="character"),
    make_option(c("--rnamining_predictions"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/predictions.txt',
                help="Path to rnamining predictions file", metavar="character"),
    make_option(c("--tx_counts"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/BambuOutput_counts_transcript.txt',
                help="Path to transcript counts file", metavar="character"),
    make_option(c("--gene_counts"), type="character", default='/mnt/beegfs/scratch/gmdazevedo/rstudio/ratomics/longnoncoder/modules/local/tx_annotation/lncRNA_filter/sandbox-e14748b9-acc0-4427-9298-d7bc49368287/test_data/Bambu_nf/BambuOutput_counts_gene.txt',
                help="Path to gene counts file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if all required arguments are provided
if (is.null(opt$bambu_gtf) || is.null(opt$compared_gtf) || is.null(opt$tmap_file) || 
    is.null(opt$rnamining_predictions) || is.null(opt$tx_counts) || is.null(opt$gene_counts)) {
    print_help(opt_parser)
    stop("All input files must be specified.", call.=FALSE)
}

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

# Save the metadata of all novel transcripts
cat("Saving novel transcripts metadata...\n")
write.csv(tx_info, file="novel_transcripts_metadata.csv", row.names = FALSE)

# Select new lncRNAs
cat("Processing lncRNAs...\n")
new_lncRNAs <- tx_info[tx_info$class_code %in% c('u', 'i', 'x', 'j', 'm', 'n') & 
                    tx_info$len >= 200 &
                    tx_info$prediction == 'non-coding', ]

new_lncRNAs$classification <- ifelse(new_lncRNAs$class_code == 'u', 'intergenic',
                                     ifelse(new_lncRNAs$class_code == 'i', 'intronic',
                                            ifelse(new_lncRNAs$class_code == 'x', 'antisense',
                                                   ifelse(new_lncRNAs$class_code == 'j', 'multiexon SJ match',
                                                          ifelse(new_lncRNAs$class_code == 'm', 'total intron retention',
                                                                 ifelse(new_lncRNAs$class_code == 'n', 'partial intron retention', NA))))))

# Save the metadata of novel lncRNAs
write.csv(new_lncRNAs, file="novel_lncRNAs_metadata.csv", row.names = FALSE)

# Save separate lncRNA gtf file
lnc_ids <- new_lncRNAs$qry_id
new_lncRNAs_tx_gtf <- subset(gtf, type == "transcript" & transcript_id %in% lnc_ids)
new_lncRNAs_exons_gtf <- subset(gtf, type == "exon" & transcript_id %in% lnc_ids)
new_lncRNAs_gtf <- c(new_lncRNAs_tx_gtf, new_lncRNAs_exons_gtf)
new_lncRNAs_gtf <- new_lncRNAs_gtf[order(new_lncRNAs_gtf$transcript_id, new_lncRNAs_gtf$type == "transcript", decreasing = TRUE)]
export(new_lncRNAs_gtf, "novel_lncRNAs.gtf")

# Select new mRNAs
cat("Processing mRNAs...\n")
new_mRNAs <- tx_info[tx_info$class_code %in% c('u', 'i', 'x', 'j', 'm', 'n') & 
                         tx_info$len >= 200 &
                         tx_info$prediction == 'coding', ]

new_mRNAs$classification <- ifelse(new_mRNAs$class_code == 'u', 'intergenic',
                                     ifelse(new_mRNAs$class_code == 'i', 'intronic',
                                            ifelse(new_mRNAs$class_code == 'x', 'antisense',
                                                   ifelse(new_mRNAs$class_code == 'j', 'multiexon SJ match',
                                                          ifelse(new_mRNAs$class_code == 'm', 'total intron retention',
                                                                 ifelse(new_mRNAs$class_code == 'n', 'partial intron retention', NA))))))

# Save the metadata of novel protein-coding RNAs
write.csv(new_mRNAs, file="novel_protein-coding_metadata.csv", row.names = FALSE)

# Save separate mRNA gtf file
mrna_ids <- new_mRNAs$qry_id
new_mRNAs_tx_gtf <- subset(gtf, type == "transcript" & transcript_id %in% mrna_ids)
new_mRNAs_exons_gtf <- subset(gtf, type == "exon" & transcript_id %in% mrna_ids)
new_mRNAs_gtf <- c(new_mRNAs_tx_gtf, new_mRNAs_exons_gtf)
new_mRNAs_gtf <- new_mRNAs_gtf[order(new_mRNAs_gtf$transcript_id, new_mRNAs_gtf$type == "transcript", decreasing = TRUE)]
export(new_mRNAs_gtf, "novel_protein-coding.gtf")

# Get exon lengths
cat("Calculating exon lengths...\n")
new_lncRNA_exon_len <- data.frame(
  transcript_id = new_lncRNAs_exons_gtf$transcript_id,
  exon_number = new_lncRNAs_exons_gtf$exon_number,
  width = width(new_lncRNAs_exons_gtf)
) 

new_mRNA_exon_len <- data.frame(
  transcript_id = new_mRNAs_exons_gtf$transcript_id,
  exon_number = new_mRNAs_exons_gtf$exon_number,
  width = width(new_mRNAs_exons_gtf)
) 

write.csv(new_lncRNA_exon_len, "novel_lncRNA_exon_lengths.csv", row.names = FALSE)
write.csv(new_mRNA_exon_len, "novel_protein-coding_exon_lengths.csv", row.names = FALSE)

# Export transcript and gene counts
cat("Processing counts data...\n")
tx_counts <- read_table(opt$tx_counts)
gene_counts <- read_table(opt$gene_counts)

new_mRNA_lncRNA <- rbind(new_mRNAs, new_lncRNAs)
write.csv(new_mRNA_lncRNA, "novel_pc_lnc_RNAs_metadata.csv", row.names = FALSE)

tx_ids <- new_mRNA_lncRNA$qry_id
gn_ids <- new_mRNA_lncRNA$qry_gene_id

novel_tx_counts <- subset(tx_counts, TXNAME %in% tx_ids)
write.csv(novel_tx_counts, "bambu_novel_pc_lnc_RNA_tx_counts.csv", row.names = FALSE)

novel_gn_counts <- subset(gene_counts, GENEID %in% gn_ids)
write.csv(novel_gn_counts, "bambu_novel_pc_lnc_RNA_gene_counts.csv", row.names = FALSE)

cat("Analysis completed successfully!\n")
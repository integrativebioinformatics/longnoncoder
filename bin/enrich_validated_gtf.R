#!/usr/bin/env Rscript

# Attach biotype and classification attributes to the validated GTFs produced by
# subset_bambu_gtf.sh.
#
# That script filters the Bambu extended annotation down to the transcripts that
# survived count validation, using awk. It is deliberately left untouched: Bambu's
# GTF carries only gene_id, transcript_id and exon_number, so the metadata has to
# come from elsewhere and adding it in awk would be fragile. This script
# post-processes the three outputs instead.
#
# Each validated GTF holds a mixture of known and novel transcripts, so both
# lookups are applied:
#   known -> the reference annotation, via annotated_transcriptome_metadata.csv
#   novel -> the gffcompare/RNAmining results, via novel_pc_lnc_RNAs_metadata.csv

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(rtracklayer)
    library(optparse)
})

# Shared reference-GTF helpers, staged alongside this script by the module
source("gtf_annotation_utils.R")

option_list <- list(
    make_option(c("--annotations_gtf"), type="character", default=NULL,
                help="BambuOutput_annotations_validated.gtf", metavar="character"),
    make_option(c("--fulllength_gtf"), type="character", default=NULL,
                help="BambuOutput_fullLength_validated.gtf", metavar="character"),
    make_option(c("--unique_gtf"), type="character", default=NULL,
                help="BambuOutput_uniquelyMapped_validated.gtf", metavar="character"),
    make_option(c("--known_metadata"), type="character", default=NULL,
                help="annotated_transcriptome_metadata.csv", metavar="character"),
    make_option(c("--novel_metadata"), type="character", default=NULL,
                help="novel_pc_lnc_RNAs_metadata.csv", metavar="character"),
    make_option(c("--annotation"), type="character", default=NULL,
                help="Path to the reference annotation GTF (Ensembl or GENCODE)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

required <- c("annotations_gtf", "fulllength_gtf", "unique_gtf",
              "known_metadata", "novel_metadata", "annotation")
missing <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
if (length(missing) > 0) {
    print_help(opt_parser)
    stop(paste("Missing required arguments:", paste(missing, collapse=", ")), call.=FALSE)
}

# --- Build the attribute lookup ------------------------------------------------

#' read.csv that tolerates the empty placeholder files the upstream steps write
#' when a biotype yielded no transcripts.
read_metadata <- function(path, label) {
    tab <- tryCatch(
        read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) data.frame()
    )
    cat(sprintf("  %s: %d rows\n", label, nrow(tab)))
    tab
}

cat("Reading metadata tables...\n")
known <- read_metadata(opt$known_metadata, "known transcripts")
novel <- read_metadata(opt$novel_metadata, "novel transcripts")

reference <- read_reference_gtf(opt$annotation)
ref_gene_biotype <- reference$gene_biotype
rm(reference)

# The validated GTFs carry Bambu's transcript_id, which for known transcripts is
# whichever id form the reference annotation used. Key the lookup on both forms so
# either matches.
lookup <- list()

if (nrow(known) > 0) {
    known_attrs <- data.frame(
        key                = known$ensembl_transcript_id_version,
        transcript_status  = "known",
        gene_name          = known$external_gene_name,
        gene_biotype       = known$gene_biotype,
        transcript_name    = known$external_transcript_name,
        transcript_biotype = known$transcript_biotype,
        class_code         = NA_character_,
        classification     = NA_character_,
        ref_gene_id        = NA_character_,
        stringsAsFactors   = FALSE
    )
    # duplicate the rows under the unversioned key as well
    known_bare <- known_attrs
    known_bare$key <- known$ensembl_transcript_id
    lookup$known <- rbind(known_attrs, known_bare)
}

if (nrow(novel) > 0) {
    ref_gene <- as.character(novel$ref_gene_id)
    ref_gene[is.na(ref_gene) | ref_gene == "-" | !nzchar(ref_gene)] <- NA_character_

    gene_biotype <- unname(ref_gene_biotype[ref_gene])
    gene_biotype[is.na(gene_biotype)] <- "novel"

    lookup$novel <- data.frame(
        key                = novel$qry_id,
        transcript_status  = "novel",
        gene_name          = as.character(novel$gene_name),
        gene_biotype       = gene_biotype,
        transcript_name    = NA_character_,
        transcript_biotype = ifelse(novel$prediction == "coding",
                                    "novel_protein_coding", "novel_lncRNA"),
        class_code         = as.character(novel$class_code),
        classification     = as.character(novel$classification),
        ref_gene_id        = ref_gene,
        stringsAsFactors   = FALSE
    )
}

if (length(lookup) == 0) {
    stop("Both metadata tables were empty; nothing to attach.", call.=FALSE)
}

attrs <- do.call(rbind, unname(lookup))
attrs <- attrs[!duplicated(attrs$key), ]
cat(sprintf("Attribute lookup built for %d transcripts\n", nrow(attrs)))

attr_cols <- setdiff(names(attrs), "key")

# --- Apply to each validated GTF ----------------------------------------------

enrich <- function(path) {
    if (is.null(path) || !file.exists(path)) {
        cat("Skipping missing file:", path, "\n")
        return(invisible(NULL))
    }

    gr <- import(path)
    if (length(gr) == 0) {
        cat("Empty GTF, re-exporting unchanged:", basename(path), "\n")
        export(gr, basename(path))
        return(invisible(NULL))
    }

    gr <- annotate_gtf(gr, attrs$key, as.list(attrs[attr_cols]))

    matched <- sum(!is.na(gr$transcript_status))
    cat(sprintf("%s: %d of %d records annotated\n",
                basename(path), matched, length(gr)))

    if (matched == 0) {
        warning(sprintf(paste0("No records in %s matched the metadata tables. ",
                               "Check that the transcript identifiers agree."),
                        basename(path)))
    }

    export(gr, basename(path))
    invisible(NULL)
}

cat("Enriching validated GTFs...\n")
enrich(opt$annotations_gtf)
enrich(opt$fulllength_gtf)
enrich(opt$unique_gtf)

cat("Validated GTF enrichment completed successfully!\n")

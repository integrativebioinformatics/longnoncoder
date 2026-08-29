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
        class_code             = NA_character_,
        classification         = NA_character_,
        # A known transcript IS the reference, so it has nothing to have been
        # compared against. NA rather than self-reference, matching class_code.
        ref_gene_id            = NA_character_,
        ref_gene_name          = NA_character_,
        ref_gene_biotype       = NA_character_,
        ref_transcript_id      = NA_character_,
        ref_transcript_name    = NA_character_,
        ref_transcript_biotype = NA_character_,
        stringsAsFactors       = FALSE
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

    # Taken from the metadata, never re-derived from the prediction. The routing
    # upstream uses the prediction AND the reference gene's biotype together: a
    # non-coding call on a model that shares splice structure with a protein-coding
    # gene becomes novel_non_coding, not novel_lncRNA. An ifelse on the prediction
    # alone cannot represent that distinction and silently collapses the two, which
    # is exactly the bug this replaced. The fallback covers a metadata file written
    # before the column existed, and says so rather than pretending otherwise.
    tx_biotype <- if ("transcript_biotype" %in% names(novel)) {
        as.character(novel$transcript_biotype)
    } else {
        warning("novel metadata has no transcript_biotype column; falling back to ",
                "the coding prediction, which cannot distinguish novel_lncRNA from ",
                "novel_non_coding.")
        ifelse(novel$prediction == "coding", "novel_protein_coding", "novel_lncRNA")
    }

    col <- function(name) if (name %in% names(novel)) as.character(novel[[name]]) else NA_character_

    lookup$novel <- data.frame(
        key                    = novel$qry_id,
        transcript_status      = "novel",
        gene_name              = as.character(novel$gene_name),
        gene_biotype           = gene_biotype,
        transcript_name        = NA_character_,
        transcript_biotype     = tx_biotype,
        class_code             = as.character(novel$class_code),
        classification         = as.character(novel$classification),
        ref_gene_id            = ref_gene,
        ref_gene_name          = col("ref_gene_name"),
        ref_gene_biotype       = col("ref_gene_biotype"),
        ref_transcript_id      = col("ref_id"),
        ref_transcript_name    = col("ref_transcript_name"),
        ref_transcript_biotype = col("ref_transcript_biotype"),
        stringsAsFactors       = FALSE
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

#' Synthesise one `gene` feature per gene_id, spanning all of its records.
#'
#' Bambu emits only transcript and exon rows, and subset_bambu_gtf.sh drops any
#' line without a transcript_id, so the validated GTFs reach here with no gene
#' features at all. Everything in this pipeline derives gene extent from the
#' transcript rows and does not need them, but tools fed the published GTF do:
#' IGV cannot collapse a locus without a gene row, and several browsers and
#' converters assume the canonical gene/transcript/exon hierarchy.
add_gene_features <- function(gr) {
    gene_ids <- as.character(mcols(gr)$gene_id)
    ok <- !is.na(gene_ids) & nzchar(gene_ids)
    if (!any(ok)) {
        warning("No gene_id values found; no gene features were added.")
        return(gr)
    }

    sub <- gr[ok]
    gid <- gene_ids[ok]

    # range() over a split keeps seqname and strand with each group, so a gene whose
    # records somehow land on two contigs yields one row per fragment rather than a
    # single span across both. That should not happen, so say so if it does.
    spans <- unlist(range(split(sub, gid)), use.names = TRUE)
    if (anyDuplicated(names(spans))) {
        dups <- unique(names(spans)[duplicated(names(spans))])
        warning(sprintf("%d gene(s) span more than one sequence or strand; one gene row emitted per fragment.",
                        length(dups)))
    }

    # Gene-level values are constant within a gene: annotate_gtf keys on
    # transcript_id, and every transcript of a gene resolves to the same gene
    # metadata. So the first matching record is representative.
    md <- mcols(sub)[match(names(spans), gid), , drop = FALSE]

    # Blank the isoform-level fields. Carrying transcript_id or class_code onto a
    # gene row would assert that the locus has one isoform's identity.
    # `x[] <- NA` rather than `x <- NA`: the former keeps the column's type, the
    # latter would turn a character column logical and then c() below would refuse
    # to bind it against the original.
    for (nm in intersect(c("transcript_id", "transcript_name", "transcript_biotype",
                           "transcript_status", "exon_number", "exon_id",
                           "class_code", "classification", "score", "phase"),
                         colnames(md))) {
        md[[nm]][] <- NA
    }
    md$gene_id <- names(spans)

    # import() returns source and type as factors, neither carrying a "gene" or
    # "pulposeq" level -- assigning to them directly would silently yield NA and
    # export a typeless row. Character on both sides sidesteps the level juggling,
    # and export() coerces either way.
    mcols(gr)$type   <- as.character(mcols(gr)$type)
    mcols(gr)$source <- as.character(mcols(gr)$source)
    md$type   <- "gene"
    md$source <- "pulposeq"

    mcols(spans) <- md
    names(spans) <- NULL

    combined <- c(spans, gr)

    # Canonical GTF ordering: loci by position, and within a locus gene, then
    # transcript, then exons. Without this every gene row would sit in one block at
    # the head of the file -- still valid GTF, but it breaks the contiguity that
    # readers streaming a locus at a time rely on.
    key_gene   <- as.character(mcols(combined)$gene_id)
    gene_start <- start(spans)[match(key_gene, mcols(spans)$gene_id)]
    type_rank  <- match(as.character(mcols(combined)$type),
                        c("gene", "transcript", "exon"))

    combined <- combined[order(as.character(seqnames(combined)),
                               gene_start, key_gene, type_rank, start(combined))]

    cat(sprintf("  added %d gene features\n", length(spans)))
    combined
}

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

    gr <- add_gene_features(gr)

    export(gr, basename(path))
    invisible(NULL)
}

cat("Enriching validated GTFs...\n")
enrich(opt$annotations_gtf)
enrich(opt$fulllength_gtf)
enrich(opt$unique_gtf)

cat("Validated GTF enrichment completed successfully!\n")

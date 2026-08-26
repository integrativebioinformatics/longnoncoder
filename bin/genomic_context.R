#!/usr/bin/env Rscript

# Genomic context figures, drawn with plotgardener. One PNG per gene, stacking:
#
#   title            gene name and the region actually plotted, chromosome included
#   gene structure   the union of the gene's exons, collapsed to one row
#   isoforms         every validated transcript in the window, novel ones highlighted
#   coverage         one track per sample, on a shared y range
#   axis             genomic coordinates
#
# Scope is deliberately narrow: only known genes carrying at least one novel isoform.
# A wholly novel locus has no annotated models to place it against, and that
# comparison is most of what makes the view worth looking at. The bigWigs are
# published, so anything not plotted here can be opened in IGV against the
# validated GTF.

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(plotgardener)
})

option_list <- list(
  make_option(c("-t", "--gtf"), type = "character", default = NULL,
              help = "Validated, enriched transcriptome GTF", metavar = "character"),
  make_option(c("-w", "--bigwigs"), type = "character", default = NULL,
              help = "Comma-separated bigWig files, one per sample", metavar = "character"),
  make_option(c("-s", "--names"), type = "character", default = NULL,
              help = "Comma-separated sample labels, same order as --bigwigs", metavar = "character"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help = "Output directory [default %default]", metavar = "character"),
  make_option(c("-c", "--context_flags"), type = "character", default = NULL,
              help = "novel_context_flags.csv from VALIDATE_NOVEL_CONTEXT", metavar = "character"),
  make_option(c("-a", "--annotation"), type = "character", default = NULL,
              help = "Reference annotation GTF, for host gene intron structure", metavar = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$gtf) || is.null(opt$bigwigs)) {
  stop("--gtf and --bigwigs are both required.", call. = FALSE)
}
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# How many genes to draw, and the isoform count above which a locus stops being
# legible in one panel. Fixed rather than exposed: these are properties of what fits
# on a page, not choices a user needs to make.
N_GENES <- 5L
MAX_TX  <- 12L

# How many flagged intronic candidates to draw. These are the visual control: a
# candidate the flags call host pre-mRNA should look like it when drawn -- coverage
# running flat across the whole intron and continuous with the flanking exons -- and
# if it does not, the flag is what needs revisiting, not the transcript.
N_INTRONIC <- 5L

# No OrgDb is needed, and that is a property of the assembly below rather than an
# oversight. plotgardener consults the OrgDb only to translate one identifier type
# into another: getExons() skips the lookup entirely when gene.id.column equals
# display.column, and only calls AnnotationDbi::select() when they differ. The
# defaults differ -- ENTREZID to SYMBOL -- which is why the documented examples all
# pass org.Hs.eg.db. Setting both to GENEID makes the lookup a no-op, and keeping it
# NULL is what lets this run for any organism. Bambu's novel genes have no entry in
# any OrgDb whatever the species.
ORGDB <- NULL

bw_files <- trimws(strsplit(opt$bigwigs, ",")[[1]])
bw_files <- bw_files[nzchar(bw_files)]
sample_names <- if (!is.null(opt$names)) {
  trimws(strsplit(opt$names, ",")[[1]])
} else {
  sub("\\.bw$", "", basename(bw_files))
}
if (length(sample_names) != length(bw_files)) {
  stop("--names must have the same number of entries as --bigwigs.", call. = FALSE)
}

# --- Read the annotation once -------------------------------------------------

cat("Reading validated annotation:", opt$gtf, "\n")
gtf <- rtracklayer::import(opt$gtf, feature.type = c("transcript", "exon"))
md  <- S4Vectors::mcols(gtf)

col_or_na <- function(name) {
  if (name %in% colnames(md)) as.character(md[[name]]) else NA_character_
}

is_tx  <- as.character(md$type) == "transcript"
tx_gr  <- gtf[is_tx]
tx_md  <- S4Vectors::mcols(tx_gr)

tx_col <- function(name) {
  if (name %in% colnames(tx_md)) as.character(tx_md[[name]]) else NA_character_
}

tx <- data.frame(
  gene_id            = tx_col("gene_id"),
  gene_name          = tx_col("gene_name"),
  transcript_id      = tx_col("transcript_id"),
  transcript_name    = tx_col("transcript_name"),
  transcript_biotype = tx_col("transcript_biotype"),
  transcript_status  = tx_col("transcript_status"),
  chrom              = as.character(GenomeInfoDb::seqnames(tx_gr)),
  start              = BiocGenerics::start(tx_gr),
  end                = BiocGenerics::end(tx_gr),
  stringsAsFactors   = FALSE
)
tx <- tx[!is.na(tx$gene_id) & nzchar(tx$gene_id), , drop = FALSE]

# Label each model by transcript_name where the annotation supplies one, falling
# back to the identifier. Bambu writes no transcript_name at all, so novel isoforms
# always fall back -- which is the desired result, since BambuTx101 is what you would
# search the GTF for. The biotype rides along because it is the thing a reader most
# wants to know about an isoform they have never seen before.
display_name <- ifelse(!is.na(tx$transcript_name) & nzchar(tx$transcript_name),
                       tx$transcript_name, tx$transcript_id)
biotype <- ifelse(is.na(tx$transcript_biotype) | !nzchar(tx$transcript_biotype),
                  "unknown_biotype", tx$transcript_biotype)
tx$label <- paste(display_name, biotype, sep = " | ")

# --- Candidate selection ------------------------------------------------------

# "Known gene" is decided by whether any of its transcripts came from the reference,
# not by an identifier prefix. Prefixes differ between Ensembl, GENCODE and every
# non-model organism, and Bambu assigns novel isoforms at a known locus the reference
# gene_id anyway, so the status column is the only reliable signal.
by_gene <- split(seq_len(nrow(tx)), tx$gene_id)

genes <- do.call(rbind, lapply(names(by_gene), function(gid) {
  idx <- by_gene[[gid]]
  nm  <- tx$gene_name[idx]
  nm  <- nm[!is.na(nm) & nzchar(nm)]
  data.frame(
    gene_id     = gid,
    gene_name   = if (length(nm)) nm[1] else NA_character_,
    chrom       = tx$chrom[idx][1],
    start       = min(tx$start[idx]),
    end         = max(tx$end[idx]),
    n_tx        = length(idx),
    n_known     = sum(tx$transcript_status[idx] == "known", na.rm = TRUE),
    n_novel     = sum(tx$transcript_status[idx] == "novel", na.rm = TRUE),
    n_novel_lnc = sum(tx$transcript_biotype[idx] == "novel_lncRNA", na.rm = TRUE),
    biotypes    = paste(sort(unique(tx$transcript_biotype[idx])), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

cand <- genes[genes$n_known > 0 & genes$n_novel > 0 & genes$n_tx <= MAX_TX, , drop = FALSE]
cat(sprintf("%d genes carry both known and novel isoforms within the %d-isoform limit\n",
            nrow(cand), MAX_TX))

if (nrow(cand) > 0) {
  # Deterministic: novel lncRNA first because that is what the pipeline exists to find,
  # then the richest loci, then gene_id so a re-run selects the same genes.
  cand <- cand[order(-cand$n_novel_lnc, -cand$n_novel, -cand$n_tx, cand$gene_id), , drop = FALSE]
  cand <- head(cand, N_GENES)
  cand$label <- ifelse(is.na(cand$gene_name) | !nzchar(cand$gene_name),
                       cand$gene_id, cand$gene_name)
  cand$pad   <- ceiling((cand$end - cand$start) * 0.05)
  cand$win_s <- pmax(1, cand$start - cand$pad)
  cand$win_e <- cand$end + cand$pad
} else {
  message("No known gene carries a novel isoform.")
}

# --- Flagged intronic candidates ----------------------------------------------

# The second selection, and a different unit: one transcript rather than one gene.
# These are the candidates the structural tests say are hardest -- inside a host
# intron, on the host's strand, and supported by reads that do not stop where the
# transcript does. Drawn on purpose so the report can show what a flagged call
# looks like beside a clean one.
FLAG_COLS <- c("qry_id", "host_gene_id", "host_gene_biotype", "class_code",
               "strand", "reads_total", "reads_crossing_boundary",
               "reads_with_host_junction")

flagged <- data.frame()

if (!is.null(opt$context_flags) && file.exists(opt$context_flags)) {
  flags <- tryCatch(read.csv(opt$context_flags, stringsAsFactors = FALSE),
                    error = function(e) NULL)

  if (!is.null(flags) && nrow(flags) && all(FLAG_COLS %in% names(flags))) {
    sel <- flags$class_code == "i" &
           flags$same_strand_as_host %in% TRUE &
           flags$reads_total > 0
    flagged <- flags[which(sel), , drop = FALSE]

    # Rank by the evidence against, not by expression: the point is to draw the
    # cases where the reads most clearly belong to something longer.
    if (nrow(flagged)) {
      flagged <- flagged[order(-flagged$reads_with_host_junction,
                               -flagged$reads_crossing_boundary,
                               flagged$qry_id), , drop = FALSE]
      flagged <- head(flagged, N_INTRONIC)
    }
  }
}

if (nrow(flagged)) {
  # Window on the host intron the candidate sits in rather than the candidate
  # itself. A window drawn tight around a truncated fragment looks like a discrete
  # transcript no matter what it is; the flanking intron is where the difference
  # shows.
  host_introns_by_gene <- list()
  if (!is.null(opt$annotation) && file.exists(opt$annotation)) {
    ref_ex   <- rtracklayer::import(opt$annotation, feature.type = "exon")
    ref_gid  <- as.character(S4Vectors::mcols(ref_ex)$gene_id)
    wanted   <- unique(flagged$host_gene_id)
    ref_ex   <- ref_ex[ref_gid %in% wanted]
    if (length(ref_ex)) {
      by_g <- split(ref_ex, as.character(S4Vectors::mcols(ref_ex)$gene_id))
      ex_r <- GenomicRanges::reduce(by_g)
      g_r  <- unlist(range(by_g), use.names = TRUE)
      g_r  <- g_r[!duplicated(names(g_r))]
      host_introns_by_gene <- as.list(GenomicRanges::psetdiff(g_r, ex_r[names(g_r)]))
    }
  }

  tx_pos <- tx[match(flagged$qry_id, tx$transcript_id), c("chrom", "start", "end")]
  flagged$chrom <- tx_pos$chrom
  flagged$start <- tx_pos$start
  flagged$end   <- tx_pos$end

  # Exon count is not carried on `tx`, and it is the one number that says whether a
  # candidate is fully unspliced -- worth counting rather than leaving out.
  all_exons <- gtf[as.character(S4Vectors::mcols(gtf)$type) == "exon"]
  n_exons   <- table(as.character(S4Vectors::mcols(all_exons)$transcript_id))
  flagged$num_exons <- as.integer(n_exons[flagged$qry_id])

  # Drop anything the GTF cannot place: the flags come from the metadata, and a
  # transcript missing from the enriched GTF has nothing to draw.
  flagged <- flagged[!is.na(flagged$chrom), , drop = FALSE]
}

if (nrow(flagged)) {
  win <- t(vapply(seq_len(nrow(flagged)), function(i) {
    gid <- flagged$host_gene_id[i]
    ivs <- host_introns_by_gene[[gid]]
    if (!is.null(ivs) && length(ivs)) {
      # The intron containing the candidate; if it straddles more than one, the
      # union of those it touches.
      hit <- which(BiocGenerics::start(ivs) <= flagged$end[i] &
                   BiocGenerics::end(ivs)   >= flagged$start[i])
      if (length(hit)) {
        return(c(min(BiocGenerics::start(ivs)[hit]), max(BiocGenerics::end(ivs)[hit])))
      }
    }
    # No usable intron: fall back to the candidate plus half its length either side,
    # which still shows whether coverage stops at the boundaries.
    pad <- ceiling((flagged$end[i] - flagged$start[i]) * 0.5)
    c(flagged$start[i] - pad, flagged$end[i] + pad)
  }, numeric(2)))

  flagged$win_s <- pmax(1, win[, 1])
  flagged$win_e <- win[, 2]
  flagged$host_gene_name <- genes$gene_name[match(flagged$host_gene_id, genes$gene_id)]
}

cat(sprintf("%d flagged intronic candidates selected for drawing\n", nrow(flagged)))

if (nrow(cand) == 0 && nrow(flagged) == 0) {
  message("Nothing to draw; writing empty candidate tables.")
  write.csv(genes[0, ], file.path(opt$outdir, "genomic_context_candidates.csv"),
            row.names = FALSE)
  write.csv(data.frame(), file.path(opt$outdir, "intronic_context_candidates.csv"),
            row.names = FALSE)
  quit(save = "no", status = 0)
}

# --- Annotation for plotgardener ----------------------------------------------

# plotTranscripts draws from a TxDb and labels each model with its TXNAME, so the
# label has to BE the transcript id as far as the TxDb is concerned. The GTF is
# rewritten accordingly -- but only over the handful of windows being drawn, which
# also keeps makeTxDbFromGFF off the whole-genome annotation. On this data that is a
# few thousand records instead of several million.
# Both selections share one TxDb: makeTxDbFromGFF over the flagged windows as well
# is a few hundred extra records, where building it twice would repeat the whole
# parse.
windows <- c(
  if (nrow(cand)) GenomicRanges::GRanges(
    cand$chrom, IRanges::IRanges(cand$win_s, cand$win_e)) else GenomicRanges::GRanges(),
  if (nrow(flagged)) GenomicRanges::GRanges(
    flagged$chrom, IRanges::IRanges(flagged$win_s, flagged$win_e)) else GenomicRanges::GRanges()
)
sub_gtf <- IRanges::subsetByOverlaps(gtf, windows)

sub_tx_id <- as.character(S4Vectors::mcols(sub_gtf)$transcript_id)
relabel   <- tx$label[match(sub_tx_id, tx$transcript_id)]
# make.unique because transcript_name is not guaranteed unique across an annotation,
# and a TxDb with duplicate transcript names silently loses models.
keep_lab  <- !is.na(relabel)
relabel[keep_lab] <- make.unique(relabel[keep_lab], sep = " #")
S4Vectors::mcols(sub_gtf)$transcript_id <- ifelse(keep_lab, relabel, sub_tx_id)

txdb_src <- file.path(opt$outdir, "genomic_context_regions.gtf")
rtracklayer::export(sub_gtf, txdb_src, format = "gtf")

cat(sprintf("Building TxDb from %d records across %d regions...\n",
            length(sub_gtf), nrow(cand)))
txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
  txdbmaker::makeTxDbFromGFF(txdb_src, format = "gtf")
} else {
  # makeTxDbFromGFF lived in GenomicFeatures before Bioconductor 3.19
  GenomicFeatures::makeTxDbFromGFF(txdb_src, format = "gtf")
}

pulposeq_assembly <- assembly(
  Genome         = "pulposeq_validated",
  TxDb           = txdb,
  OrgDb          = ORGDB,
  gene.id.column = "GENEID",
  display.column = "GENEID"
)

# --- Drawing ------------------------------------------------------------------

exons_gr <- gtf[as.character(S4Vectors::mcols(gtf)$type) == "exon"]
exon_gid <- as.character(S4Vectors::mcols(exons_gr)$gene_id)

#' Draw one region: gene structure, the isoform models over it, and one coverage
#' track per sample.
#'
#' Parameterised rather than branching internally, because the two callers differ
#' only in what they consider the structure row, which transcripts belong in the
#' panel, and what gets highlighted.
#'
#' @param structure_gr GRanges drawn as the single "gene" row above the isoforms
#' @param panel_tx     rows of `tx` in this window; sizes the transcript track
#' @param highlights   labels to colour, as the TxDb now knows them
draw_panel <- function(chrom, win_s, win_e, structure_gr, panel_tx,
                       highlights, title, subtitle, out_png, structure_label = "gene") {

  hl <- if (length(highlights)) {
    data.frame(transcript = highlights, color = "#D95F02", stringsAsFactors = FALSE)
  } else {
    NULL
  }

  # A shared y range across samples, so track heights compare directly instead of
  # each being rescaled to its own maximum.
  peaks <- vapply(bw_files, function(f) {
    sig <- rtracklayer::import.bw(
      f, which = GenomicRanges::GRanges(chrom, IRanges::IRanges(win_s, win_e)))
    if (length(sig) == 0) 0 else max(S4Vectors::mcols(sig)$score, na.rm = TRUE)
  }, numeric(1))
  ymax <- max(c(peaks, 1))

  margin  <- 0.5
  width   <- 8.0
  title_h <- 0.34
  gene_h  <- 0.26
  tx_h    <- max(1.0, nrow(panel_tx) * 0.18)
  sig_h   <- 0.55
  gap     <- 0.12
  label_h <- 0.45
  height  <- margin * 2 + title_h + gene_h + tx_h + gap * 3 +
             length(bw_files) * (sig_h + gap) + label_h

  png(out_png, width = width, height = height, units = "in", res = 300)
  pageCreate(width = width, height = height, default.units = "inches",
             showGuides = FALSE)

  pars <- pgParams(
    chrom = chrom, chromstart = win_s, chromend = win_e,
    assembly = pulposeq_assembly,
    x = margin, width = width - 2 * margin, default.units = "inches"
  )

  # The region is named in full here, chromosome included, so a figure lifted out of
  # the report still says where it came from.
  plotText(
    label = sprintf("%s  |  %s:%s-%s", title, chrom,
                    format(win_s, big.mark = ","), format(win_e, big.mark = ",")),
    x = margin, y = margin, just = c("left", "top"),
    fontsize = 11, fontface = "bold", default.units = "inches"
  )
  plotText(
    label = subtitle,
    x = margin, y = margin + 0.17, just = c("left", "top"),
    fontsize = 8, fontcolor = "grey35", default.units = "inches"
  )

  y <- margin + title_h
  if (length(structure_gr)) {
    plotRanges(
      params = pars, data = structure_gr,
      collapse = TRUE, fill = "#4D4D4D", linecolor = NA,
      y = y, height = gene_h
    )
    plotText(label = structure_label, x = margin - 0.06, y = y + gene_h / 2,
             just = c("right", "center"), fontsize = 7,
             fontcolor = "grey35", default.units = "inches")
  }

  y <- y + gene_h + gap
  plotTranscripts(
    params = pars, y = y, height = tx_h,
    labels = "transcript", fontsize = 6,
    transcriptHighlights = hl
  )

  y <- y + tx_h + gap
  for (i in seq_along(bw_files)) {
    plotSignal(
      params = pars, data = bw_files[i],
      y = y, height = sig_h,
      range = c(0, ymax), scale = TRUE, label = sample_names[i],
      linecolor = "#3B6FB6", fill = "#3B6FB6"
    )
    y <- y + sig_h + gap
  }

  plotGenomeLabel(
    chrom = chrom, chromstart = win_s, chromend = win_e,
    assembly = pulposeq_assembly,
    x = margin, y = y, length = width - 2 * margin,
    default.units = "inches", scale = "bp", fontsize = 8
  )

  pageGuideHide()
  dev.off()
  basename(out_png)
}

# --- Genes carrying novel isoforms --------------------------------------------

if (nrow(cand)) {
  cand$figure <- NA_character_
  for (i in seq_len(nrow(cand))) {
    row <- cand[i, ]
    cat(sprintf("Drawing %s at %s:%d-%d (%d isoforms: %d known, %d novel, %d novel lncRNA)\n",
                row$label, row$chrom, row$win_s, row$win_e,
                row$n_tx, row$n_known, row$n_novel, row$n_novel_lnc))

    gene_tx <- tx[tx$gene_id == row$gene_id, , drop = FALSE]

    # Novel isoforms highlighted against the annotated ones. Without this they are
    # coloured by strand like everything else and the figure stops making its point.
    # Keys must be the rewritten labels, which is what the TxDb now calls them.
    novel_lab <- gene_tx$label[gene_tx$transcript_status == "novel"]

    cand$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      # The union of every exon the gene has, so a reader sees the full exonic
      # footprint before the isoforms start differing from it.
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$gene_id]),
      panel_tx     = gene_tx,
      highlights   = novel_lab[!is.na(novel_lab)],
      title        = row$label,
      subtitle     = sprintf("%d isoforms: %d known, %d novel (%d novel lncRNA candidate%s)",
                             row$n_tx, row$n_known, row$n_novel, row$n_novel_lnc,
                             if (row$n_novel_lnc == 1) "" else "s"),
      out_png      = file.path(opt$outdir, sprintf("genomic_context_%s.png", row$label))
    )
  }

  write.csv(cand[, c("gene_id", "gene_name", "label", "chrom", "start", "end",
                     "win_s", "win_e", "n_tx", "n_known", "n_novel", "n_novel_lnc",
                     "biotypes", "figure")],
            file.path(opt$outdir, "genomic_context_candidates.csv"), row.names = FALSE)
} else {
  write.csv(genes[0, ], file.path(opt$outdir, "genomic_context_candidates.csv"),
            row.names = FALSE)
}
cat("Wrote genomic_context_candidates.csv\n")

# --- Flagged intronic candidates ----------------------------------------------

if (nrow(flagged)) {
  flagged$figure <- NA_character_
  for (i in seq_len(nrow(flagged))) {
    row <- flagged[i, ]
    cat(sprintf("Drawing flagged %s in %s intron at %s:%d-%d\n",
                row$qry_id, ifelse(is.na(row$host_gene_name), row$host_gene_id,
                                   row$host_gene_name),
                row$chrom, row$win_s, row$win_e))

    # Everything in the window, not just the host's isoforms: an intronic candidate
    # carries its own Bambu gene id, so filtering on the host gene would leave the
    # transcript the figure exists to show out of the panel entirely.
    in_window <- tx$chrom == row$chrom & tx$end >= row$win_s & tx$start <= row$win_e
    panel_tx  <- tx[in_window, , drop = FALSE]

    cand_lab <- tx$label[match(row$qry_id, tx$transcript_id)]

    host_lab <- if (is.na(row$host_gene_name) || !nzchar(row$host_gene_name)) {
      row$host_gene_id
    } else {
      row$host_gene_name
    }

    flagged$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$host_gene_id]),
      panel_tx     = panel_tx,
      highlights   = cand_lab[!is.na(cand_lab)],
      title        = sprintf("%s in %s", row$qry_id, host_lab),
      subtitle     = sprintf(
        "intronic, same strand as host (%s) | %d reads: %d cross a boundary, %d carry a host junction",
        ifelse(is.na(row$host_gene_biotype), "unknown biotype", row$host_gene_biotype),
        row$reads_total, row$reads_crossing_boundary, row$reads_with_host_junction),
      out_png      = file.path(opt$outdir,
                               sprintf("intronic_context_%s.png", row$qry_id)),
      structure_label = "host"
    )
  }

  write.csv(flagged[, c("qry_id", "host_gene_id", "host_gene_name", "host_gene_biotype",
                        "chrom", "start", "end", "win_s", "win_e", "strand",
                        "class_code", "num_exons", "reads_total",
                        "reads_crossing_boundary", "reads_with_host_junction",
                        "figure")],
            file.path(opt$outdir, "intronic_context_candidates.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(), file.path(opt$outdir, "intronic_context_candidates.csv"),
            row.names = FALSE)
}
cat("Wrote intronic_context_candidates.csv\n")

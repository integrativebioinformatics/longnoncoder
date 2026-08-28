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

# Flagged candidates are ranked on the fraction of reads carrying a host junction, so
# a candidate needs enough reads for that fraction to mean anything.
MIN_READS_FOR_FRAC <- 20L

# The drawing window for a flagged candidate is the host intron it sits in, padded so
# the flanking exons are on screen. Without them the panel cannot show whether
# coverage stops at the intron boundary or runs straight through it, which is the only
# question it exists to answer.
INTRON_PAD_FRAC <- 0.15
INTRON_PAD_MIN  <- 1000L

# Isoforms are coloured by transcript biotype, and the legend under each panel lists
# the biotypes that panel actually contains. Novel biotypes take the saturated end of
# the palette so a novel model still reads as novel at a glance -- which is what the
# old novel-only highlight did -- while known biotypes stay muted. Anything unlisted
# falls back to grey rather than erroring: an annotation can always carry a biotype
# this list has not seen.
BIOTYPE_COLORS <- c(
  novel_lncRNA                   = "#D95F02",
  novel_protein_coding           = "#7570B3",
  novel_retained_intron          = "#1B9E77",
  protein_coding                 = "#7FBC41",
  lncRNA                         = "#4393C3",
  retained_intron                = "#8DA0CB",
  nonsense_mediated_decay        = "#E7969C",
  protein_coding_CDS_not_defined = "#BDBDBD",
  processed_transcript           = "#D9D9D9",
  processed_pseudogene           = "#C7C7C7",
  unprocessed_pseudogene         = "#C7C7C7",
  unknown_biotype                = "#969696"
)
BIOTYPE_FALLBACK <- "#969696"

biotype_color <- function(bt) {
  out <- unname(BIOTYPE_COLORS[bt])
  out[is.na(out)] <- BIOTYPE_FALLBACK
  out
}

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

# Chromosome lengths for the ideogram, taken from the bigWig header rather than from
# an assembly package: no extra input, no network, and it is guaranteed to agree with
# the coverage being drawn beneath it. If it cannot be read the ideogram is skipped
# and the rest of the panel is unaffected.
chrom_lens <- tryCatch(
  GenomeInfoDb::seqlengths(rtracklayer::BigWigFile(bw_files[1])),
  error = function(e) NULL)
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
# Kept as its own column because the figures colour by it and the legend lists it.
tx$biotype_key <- biotype

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

    # Rank on the FRACTION of supporting reads that carry a host junction, not on the
    # raw count. A count ranks by expression: 225 of 2084 reads (10.8%) outranks 40 of
    # 40 (100%), so the panel fills up with whatever is highly expressed instead of
    # whatever is most clearly host-derived. validate_novel_context.R already writes
    # the fractions; recompute only as a fallback for an older flags table.
    if (nrow(flagged)) {
      # Optional in older flags tables, and the one signal that separates host
      # pre-mRNA from an unannotated host exon -- worth carrying even when absent.
      if (!"reads_spliced_into_host_exon" %in% names(flagged)) {
        flagged$reads_spliced_into_host_exon <- NA_integer_
      }
      frac <- function(num) ifelse(flagged$reads_total > 0, num / flagged$reads_total, 0)
      flagged$frac_host_junction  <- if ("frac_with_host_junction" %in% names(flagged)) {
        flagged$frac_with_host_junction
      } else frac(flagged$reads_with_host_junction)
      flagged$frac_crossing       <- if ("frac_crossing_boundary" %in% names(flagged)) {
        flagged$frac_crossing_boundary
      } else frac(flagged$reads_crossing_boundary)
      flagged$frac_into_host_exon <- frac(flagged$reads_spliced_into_host_exon)

      # A floor on read count, so 3-of-3 does not beat 1700-of-2147 on noise. Skipped
      # if it would empty the panel, which it would on a small test run.
      enough <- flagged$reads_total >= MIN_READS_FOR_FRAC
      if (any(enough)) flagged <- flagged[enough, , drop = FALSE]

      flagged <- flagged[order(-flagged$frac_host_junction,
                               -flagged$frac_crossing,
                               -flagged$reads_total,
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
        i_s <- min(BiocGenerics::start(ivs)[hit])
        i_e <- max(BiocGenerics::end(ivs)[hit])
        # Padded past the intron edges: a window clipped exactly to the intron leaves
        # the flanking exons off screen, so read-through at the boundary is invisible.
        pad <- max(INTRON_PAD_MIN, ceiling((i_e - i_s) * INTRON_PAD_FRAC))
        return(c(i_s - pad, i_e + pad))
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
# Which transcripts each panel is meant to contain. plotTranscripts draws whatever the
# TxDb holds inside the window, so this set -- not the window alone -- is what decides
# a panel's contents. Restricting it stops a neighbouring gene's isoforms from crowding
# out the ones the figure is about: a plain overlap query on one 84 kb window pulled in
# five genes and 31 transcripts to draw a caption that promised eight.
draw_ids <- character(0)
if (nrow(cand)) {
  draw_ids <- c(draw_ids, tx$transcript_id[tx$gene_id %in% cand$gene_id])
}
if (nrow(flagged)) {
  for (i in seq_len(nrow(flagged))) {
    in_win <- tx$chrom == flagged$chrom[i] &
              tx$end   >= flagged$win_s[i] &
              tx$start <= flagged$win_e[i]
    # The host's isoforms for context plus every novel model in the window: an
    # intronic candidate carries its own Bambu gene id, so filtering on the host gene
    # alone would drop the transcript the figure exists to show.
    draw_ids <- c(draw_ids,
                  tx$transcript_id[in_win & (tx$gene_id == flagged$host_gene_id[i] |
                                             tx$transcript_status == "novel")])
  }
}
draw_ids <- unique(draw_ids[!is.na(draw_ids)])

sub_gtf   <- IRanges::subsetByOverlaps(gtf, windows)
row_tx_id <- as.character(S4Vectors::mcols(sub_gtf)$transcript_id)
# Gene-level rows carry no transcript_id and are kept regardless.
sub_gtf   <- sub_gtf[is.na(row_tx_id) | row_tx_id %in% draw_ids]
sub_tx_id <- as.character(S4Vectors::mcols(sub_gtf)$transcript_id)

# make.unique has to run over the transcript -> label mapping, one entry per
# transcript, NOT over GTF rows. sub_gtf carries one row per transcript AND one per
# exon, so uniquifying the row vector renames a transcript's own exons to "... #1",
# "... #2": they stop matching their parent transcript_id, every model loads into the
# TxDb with zero exons, and an 8-exon transcript draws as one featureless bar plus
# eight single-exon decoys. Uniquify the mapping, then apply it to every row.
uniq_tx        <- unique(sub_tx_id[!is.na(sub_tx_id)])
lab            <- tx$label[match(uniq_tx, tx$transcript_id)]
have_lab       <- !is.na(lab)
lab[have_lab]  <- make.unique(lab[have_lab], sep = " #")
lab[!have_lab] <- uniq_tx[!have_lab]
names(lab)     <- uniq_tx
relabel        <- unname(lab[sub_tx_id])
S4Vectors::mcols(sub_gtf)$transcript_id <- ifelse(is.na(relabel), sub_tx_id, relabel)

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

#' How many rows plotTranscripts will need for a set of models.
#'
#' plotgardener packs models into rows by position and writes a label above each, so a
#' row is claimed by whichever is wider, the model or its label. Working the packing
#' out here is what lets the panel height be exact. Sizing on the transcript count
#' alone leaves half the figure blank when models pack together; sizing on the bar
#' height alone clips the overflow, which is how a panel captioned "11 isoforms" came
#' to draw seven and lose three of its four novel candidates.
pack_rows <- function(starts, ends, labels, win_s, win_e, track_w, fontsize) {
  n <- length(starts)
  if (!n) return(1L)
  span    <- max(1, win_e - win_s)
  # Roughly the advance width of one character at this size, in inches, as bp.
  char_in <- fontsize * 0.0075
  lab_bp  <- nchar(labels) * char_in / track_w * span
  mid <- (starts + ends) / 2
  l   <- pmin(starts, mid - lab_bp / 2)
  r   <- pmax(ends,   mid + lab_bp / 2)
  pad <- 0.02 * span                    # plotgardener's default inter-model spacing
  row_end <- numeric(0)
  for (i in order(l)) {
    slot <- which(row_end + pad < l[i])
    k <- if (length(slot)) slot[1] else length(row_end) + 1L
    row_end[k] <- r[i]
  }
  max(1L, length(row_end))
}

#' A chromosome bar with the drawn window marked on it, and the zoom lines down to
#' the tracks.
#'
#' Deliberately not plotgardener's plotIdeogram(): that resolves cytobands through
#' AnnotationHub, which means a network call from a compute node that usually cannot
#' make one, and a failed run at the last step. This carries what the figure needs --
#' where on the chromosome the window sits -- and needs nothing but the chromosome
#' length, which the bigWig header already supplies.
draw_ideogram <- function(chrom, win_s, win_e, chrom_len,
                          x, y, w, h, panel_x0, panel_x1, zoom_y1) {
  plotRect(x = x, y = y, width = w, height = h, just = c("left", "top"),
           default.units = "inches", fill = "#F0F0F0", linecolor = "#BDBDBD")

  f0 <- max(0, min(1, win_s / chrom_len))
  f1 <- max(0, min(1, win_e / chrom_len))
  rx0 <- x + f0 * w
  # Widened to a floor, or a 20 kb window on a 250 Mb chromosome is a hairline.
  rx1 <- max(x + f1 * w, rx0 + 0.03)

  plotRect(x = rx0, y = y, width = rx1 - rx0, height = h, just = c("left", "top"),
           default.units = "inches", fill = "#C6E48B", linecolor = NA)
  plotSegments(x0 = (rx0 + rx1) / 2, y0 = y - 0.05,
               x1 = (rx0 + rx1) / 2, y1 = y + h + 0.05,
               default.units = "inches", linecolor = "#D62728", lwd = 1.2)
  plotText(label = sprintf("Chromosome %s", sub("^chr", "", chrom)),
           x = x + w, y = y + h + 0.07, just = c("right", "top"),
           fontsize = 8, fontcolor = "grey45", default.units = "inches")

  plotSegments(x0 = rx0, y0 = y + h + 0.03, x1 = panel_x0, y1 = zoom_y1,
               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = 0.8)
  plotSegments(x0 = rx1, y0 = y + h + 0.03, x1 = panel_x1, y1 = zoom_y1,
               default.units = "inches", linecolor = "#CCCCCC", lty = 2, lwd = 0.8)
}

#' Draw one region: chromosome position, gene structure, the isoform models over it,
#' and one coverage track per sample.
#'
#' Parameterised rather than branching internally, because the two callers differ
#' only in what they consider the structure row, which transcripts belong in the
#' panel, and what gets highlighted.
#'
#' @param structure_gr GRanges drawn as the single "gene" row above the isoforms
#' @param panel_tx     rows of `tx` this panel will draw; sizes the track and supplies
#'                     the biotype each model is coloured by
draw_panel <- function(chrom, win_s, win_e, structure_gr, panel_tx,
                       title, subtitle, out_png, structure_label = "gene") {

  # Colour every model by its transcript biotype, keyed on the label the TxDb actually
  # ended up using: `lab` is the same transcript -> label mapping the sub-GTF was
  # rewritten with, so a name that genuinely needed a uniquifying suffix still matches.
  panel_lab <- unname(lab[panel_tx$transcript_id])
  keep      <- !is.na(panel_lab)
  hl <- if (any(keep)) {
    data.frame(transcript = panel_lab[keep],
               color      = biotype_color(panel_tx$biotype_key[keep]),
               stringsAsFactors = FALSE)
  } else NULL

  # A shared y range across samples, so track heights compare directly instead of
  # each being rescaled to its own maximum.
  peaks <- vapply(bw_files, function(f) {
    sig <- rtracklayer::import.bw(
      f, which = GenomicRanges::GRanges(chrom, IRanges::IRanges(win_s, win_e)))
    if (length(sig) == 0) 0 else max(S4Vectors::mcols(sig)$score, na.rm = TRUE)
  }, numeric(1))
  ymax <- max(c(peaks, 1))

  margin   <- 0.5
  width    <- 8.0
  track_w  <- width - 2 * margin
  title_h  <- 0.34
  ideo_h   <- 0.17
  ideo_lab <- 0.22
  zoom_h   <- 0.40
  gene_h   <- 0.26
  row_h    <- 0.27          # one packed row: the label plus the model beneath it
  sig_h    <- 0.55
  gap      <- 0.12
  label_h  <- 0.45
  legend_h <- 0.34

  # Height follows the real packing, so the figure is exactly as tall as its contents.
  # One spare row, capped at the worst case of every model on its own row. pack_rows
  # estimates label widths, so it can be off by one either way; erring tall costs a
  # little whitespace, erring short silently drops models off the bottom.
  n_rows <- pack_rows(panel_tx$start, panel_tx$end, panel_tx$label,
                      win_s, win_e, track_w, 6)
  n_rows <- min(max(nrow(panel_tx), 1L), n_rows + 1L)
  tx_h   <- max(0.5, n_rows * row_h)

  chrom_len <- if (!is.null(chrom_lens) && chrom %in% names(chrom_lens)) {
    as.numeric(chrom_lens[[chrom]])
  } else NA_real_
  show_ideo <- is.finite(chrom_len) && chrom_len > 0
  head_h    <- if (show_ideo) ideo_h + ideo_lab + zoom_h else 0

  bts       <- sort(unique(panel_tx$biotype_key))
  show_lgnd <- length(bts) > 0
  lgnd_h    <- if (show_lgnd) legend_h else 0

  height  <- margin * 2 + title_h + head_h + gene_h + tx_h + gap * 3 +
             length(bw_files) * (sig_h + gap) + label_h + lgnd_h

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
  if (show_ideo) {
    draw_ideogram(chrom, win_s, win_e, chrom_len,
                  x = margin, y = y, w = track_w, h = ideo_h,
                  panel_x0 = margin, panel_x1 = margin + track_w,
                  zoom_y1  = y + ideo_h + ideo_lab + zoom_h - 0.04)
    y <- y + ideo_h + ideo_lab + zoom_h
  }

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
    x = margin, y = y, length = track_w,
    default.units = "inches", scale = "bp", fontsize = 8
  )

  # Only the biotypes this panel contains, so the key never explains a colour that is
  # not on screen.
  if (show_lgnd) {
    plotLegend(
      legend = bts, fill = biotype_color(bts), border = FALSE,
      orientation = "h", fontsize = 7,
      x = margin, y = y + label_h, width = track_w, height = legend_h - 0.06,
      just = c("left", "top"), default.units = "inches"
    )
  }

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

    # Everything that will actually be drawn, which after the draw_ids restriction is
    # this gene's isoforms plus any overlapping gene that is itself a candidate.
    # Sizing on the gene's own transcripts under-counts whenever two candidates
    # overlap, and the panel is coloured from these rows as well as sized by them.
    panel_tx <- tx[tx$chrom == row$chrom & tx$end >= row$win_s &
                   tx$start <= row$win_e & tx$transcript_id %in% draw_ids, ,
                   drop = FALSE]

    cand$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      # The union of every exon the gene has, so a reader sees the full exonic
      # footprint before the isoforms start differing from it.
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$gene_id]),
      panel_tx     = panel_tx,
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
    panel_tx  <- tx[in_window & tx$transcript_id %in% draw_ids, , drop = FALSE]

    host_lab <- if (is.na(row$host_gene_name) || !nzchar(row$host_gene_name)) {
      row$host_gene_id
    } else {
      row$host_gene_name
    }

    flagged$figure[i] <- draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$host_gene_id]),
      panel_tx     = panel_tx,
      title        = sprintf("%s in %s", row$qry_id, host_lab),
      # Percentages, because the ranking is on fractions and a bare count invites the
      # same misreading the ranking used to make. The host-exon figure is reported
      # separately: reads splicing into a host exon mean an unannotated exon of the
      # host, which is a finding rather than an artifact.
      subtitle     = sprintf(
        "intronic, same strand as host (%s) | %d reads: %.0f%% cross a boundary, %.0f%% carry a host junction%s",
        ifelse(is.na(row$host_gene_biotype), "unknown biotype", row$host_gene_biotype),
        row$reads_total, 100 * row$frac_crossing, 100 * row$frac_host_junction,
        if (is.na(row$reads_spliced_into_host_exon)) "" else
          sprintf(", %.0f%% splice into a host exon", 100 * row$frac_into_host_exon)),
      out_png      = file.path(opt$outdir,
                               sprintf("intronic_context_%s.png", row$qry_id)),
      structure_label = "host"
    )
  }

  write.csv(flagged[, c("qry_id", "host_gene_id", "host_gene_name", "host_gene_biotype",
                        "chrom", "start", "end", "win_s", "win_e", "strand",
                        "class_code", "num_exons", "reads_total",
                        "reads_crossing_boundary", "reads_with_host_junction",
                        "reads_spliced_into_host_exon", "frac_crossing",
                        "frac_host_junction", "frac_into_host_exon",
                        "figure")],
            file.path(opt$outdir, "intronic_context_candidates.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(), file.path(opt$outdir, "intronic_context_candidates.csv"),
            row.names = FALSE)
}
cat("Wrote intronic_context_candidates.csv\n")

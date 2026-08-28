#!/usr/bin/env Rscript

# Genomic context figures, drawn with plotgardener. Two families of PNG -- one per
# selected gene, one per flagged intronic candidate -- each stacking:
#
#   title       the subject and the region actually plotted, chromosome included
#   chromosome  a bar marking where on the chromosome this window sits
#   structure   the union of the gene's (or the host's) exons, collapsed to one row
#   subject     intronic panels only: the candidate on a row of its own
#   isoforms    every validated transcript in the window, coloured by biotype
#   coverage    one track per sample, on a y range shared across samples
#   axis        genomic coordinates
#   legend      the biotypes this panel contains
#
# Two scope decisions, both deliberate.
#
# Only known genes carrying at least one novel isoform are drawn. A wholly novel locus
# has no annotated model to place it against, and that comparison is most of what makes
# the view worth looking at. No biotype filter applies: a novel lncRNA inside a known
# lncRNA gene competes for a panel on the same terms as one inside a protein-coding
# gene.
#
# The isoform track shows the VALIDATED set -- transcripts supported in these samples --
# and not the reference in full. A reference isoform with no support here is absent by
# design: the panel is a picture of what the data contains, and an unexpressed model
# drawn solid over flat coverage reads as evidence it is not. The bigWigs are published,
# so the complete annotation can be opened against them in IGV.

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
  "biotype not set"              = "#969696"
)
BIOTYPE_FALLBACK <- "#969696"

# Isoform track geometry, passed straight to plotTranscripts.
#
# boxHeight is raised from plotgardener's 2 mm default for contrast: a taller exon box
# against the same thin intron line is the thick-against-thin reading we want, and
# plotTranscripts draws UTR thinner than CDS by itself, so that comes free.
#
# It is NOT the reason labels go missing, which was the earlier guess. plotTranscripts
# centres a label on the transcript's true midpoint and does not draw it when that falls
# outside the window: in the ACOT7 intronic panel exactly two of ten models had an
# in-window midpoint, and exactly those two were labelled. Gene panels never hit this,
# because their window spans min(start) to max(end) over the gene's own transcripts, so
# every midpoint is inside it by construction. Intronic panels are windowed on one host
# intron, so host isoforms running the length of the gene lose their names -- see the
# note on the isoform track in draw_panel().
TX_LABEL_SIZE <- 6
TX_LABEL_H    <- 0.09    # band above each model for its label
TX_BOX_MM     <- 4       # exon box height, against plotgardener's 2 mm default
TX_SPACE_H    <- 0.2     # padding inside a model's own band, as a fraction of the box
TX_STROKE     <- 0.1     # transcript body outline, plotgardener's default

TX_MODEL_H <- (TX_BOX_MM / 25.4) * (1 + TX_SPACE_H)
TX_ROW_H   <- TX_LABEL_H + TX_MODEL_H   # per_model: label band plus the model's own

# Every panel is drawn once per renderer, to a filename carrying the renderer's suffix,
# so the two can be compared on the same data in one run rather than by editing the
# source between runs. The empty suffix is the one the candidate CSVs and the report
# point at; once one renderer wins, drop the other entry and this collapses back to a
# single call.
#
#   packed     one plotTranscripts call. It packs models into rows and labels them, but
#              places a label at the model's midpoint and omits it when that falls
#              outside the window -- so a host isoform spanning a whole gene goes
#              unnamed in a panel windowed on one of its introns.
#   per_model  one call per transcript, filtered to that transcript and placed at a row
#              of our choosing with its own labels off. The label is then ours to draw,
#              clamped inside the track, so an off-window midpoint still gets a name.
#              plotgardener still renders the model, so UTR-vs-CDS thickness is
#              unchanged. Costs one call per model instead of one per panel.
TX_RENDERERS <- c(packed = "", per_model = ".per_model")

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
  class_code         = tx_col("class_code"),
  classification     = tx_col("classification"),
  chrom              = as.character(GenomeInfoDb::seqnames(tx_gr)),
  strand             = as.character(BiocGenerics::strand(tx_gr)),
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
                  "biotype not set", tx$transcript_biotype)
# plotTranscripts labels each model with its TXNAME, so the label has to BE the
# transcript_id as far as the TxDb is concerned -- one string, used on the figure and
# in the rewritten sub-GTF alike.
#
# Biotype is left out of it: every model is coloured by biotype and the legend names
# it, so repeating it on each label spends a third of the track width restating the
# palette. The GTF keeps transcript_biotype as its own attribute, so nothing is lost.
# The bare class code goes too -- the mapping to its wording is one-to-one, "intronic"
# says everything "i" does, and the letter is still on every row of the metadata CSV.
cls <- ifelse(is.na(tx$classification) | !nzchar(tx$classification), "",
              paste0(" | ", tx$classification))
str_mark <- ifelse(tx$strand %in% c("+", "-"), paste0(" | ", tx$strand), "")

tx$label <- paste0(display_name, cls, str_mark)
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

# plotTranscripts draws from a TxDb and labels each model with its TXNAME, so the label
# has to BE the transcript id as far as the TxDb is concerned. The GTF is rewritten
# accordingly -- but only over the handful of windows being drawn, which also keeps
# makeTxDbFromGFF off the whole-genome annotation. On this data that is a few thousand
# records instead of several million.
#
# Note plotTranscripts also has a transcriptFilter argument that takes transcript names
# to display. Restricting the GTF here rather than filtering per call is still the right
# way round: it is what keeps the TxDb build small, and one build is shared by every
# panel.
# Both selections share one TxDb: makeTxDbFromGFF over the flagged windows as well
# is a few hundred extra records, where building it twice would repeat the whole
# parse.
windows <- c(
  if (nrow(cand)) GenomicRanges::GRanges(
    cand$chrom, IRanges::IRanges(cand$win_s, cand$win_e)) else GenomicRanges::GRanges(),
  if (nrow(flagged)) GenomicRanges::GRanges(
    flagged$chrom, IRanges::IRanges(flagged$win_s, flagged$win_e)) else GenomicRanges::GRanges()
)
# Which transcripts each panel is meant to contain -- this set, not the window alone,
# is what decides a panel's contents. Restricting it stops a neighbouring gene's
# isoforms from crowding
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

exons_gr  <- gtf[as.character(S4Vectors::mcols(gtf)$type) == "exon"]
exon_gid  <- as.character(S4Vectors::mcols(exons_gr)$gene_id)
exon_txid <- as.character(S4Vectors::mcols(exons_gr)$transcript_id)

#' Assign each model to a row, packing models that do not overlap onto the same one.
#'
#' A label sits above every model, so a row is claimed by whichever is wider, the model
#' or its label. The assignment returned here both draws the track and sizes it, which
#' is what makes the panel height exact. Sizing on the transcript count
#' alone leaves half the figure blank when models pack together; sizing on the bar
#' height alone clips the overflow, which is how a panel captioned "11 isoforms" came
#' to draw seven and lose three of its four novel candidates.
pack_rows <- function(starts, ends, labels, win_s, win_e, track_w, fontsize) {
  n <- length(starts)
  if (!n) return(integer(0))
  span    <- max(1, win_e - win_s)
  # Roughly the advance width of one character at this size, in inches, as bp.
  char_in <- fontsize * 0.0075
  lab_bp  <- nchar(labels) * char_in / track_w * span
  # Clamped to the window: a model running off the edge only competes for the space it
  # actually occupies on screen, and its label is centred on the visible part.
  vs  <- pmax(starts, win_s)
  ve  <- pmin(ends,   win_e)
  mid <- (vs + ve) / 2
  l   <- pmin(vs, mid - lab_bp / 2)
  r   <- pmax(ve, mid + lab_bp / 2)
  pad <- 0.02 * span
  row_end <- numeric(0)
  rows    <- integer(n)
  for (i in order(l)) {
    slot <- which(row_end + pad < l[i])
    k <- if (length(slot)) slot[1] else length(row_end) + 1L
    row_end[k] <- r[i]
    rows[i]    <- k
  }
  rows
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
#' @param subject_gr   exons of the one transcript the panel is named after, drawn on
#'                     a dedicated row of its own. plotTranscripts packs by position
#'                     and drops labels that would collide, so the subject can end up
#'                     sharing a row with another model and losing its label to it --
#'                     a panel titled BambuTx1528 showed only BambuTx1193, 8 kb away.
#'                     Its own row is the only way to guarantee the title and the
#'                     highlighted model are the same transcript.
draw_panel <- function(chrom, win_s, win_e, structure_gr, panel_tx, renderer,
                       title, subtitle, out_png, structure_label = "gene",
                       subject_gr = NULL, subject_label = NULL,
                       subject_fill = NULL) {

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
  # per_model owns the row pitch, because it places every model itself. packed leaves the
  # pitch to plotgardener, which is boxHeight plus its spacing fraction.
  row_h    <- if (renderer == "per_model") TX_ROW_H else TX_MODEL_H
  sig_h    <- 0.55
  sig_lab_h <- 0.17        # band above each trace for its scale and sample name
  gap      <- 0.12
  label_h  <- 0.45
  legend_h <- 0.34

  # The same packing that draws the track sizes it, so the height is exact rather than
  # estimated -- no spare row, and nothing to clip.
  tx_rows <- pack_rows(panel_tx$start, panel_tx$end, panel_tx$label,
                       win_s, win_e, track_w, TX_LABEL_SIZE)
  n_rows  <- if (length(tx_rows)) max(tx_rows) else 1L
  # per_model is exact: every model sits at a row we chose, so nothing can overflow.
  # packed needs slack, because pack_rows only estimates what plotgardener will do and
  # being short there means models dropped with a "+" in the corner.
  tx_h    <- if (renderer == "per_model") max(0.4, n_rows * row_h)
             else               max(0.4, (n_rows + 0.5) * row_h)

  chrom_len <- if (!is.null(chrom_lens) && chrom %in% names(chrom_lens)) {
    as.numeric(chrom_lens[[chrom]])
  } else NA_real_
  show_ideo <- is.finite(chrom_len) && chrom_len > 0
  head_h    <- if (show_ideo) ideo_h + ideo_lab + zoom_h else 0

  show_subj <- !is.null(subject_gr) && length(subject_gr) > 0
  subj_h    <- if (show_subj) 0.32 else 0

  # plotLegend lays its entries in a single row and divides the width equally between
  # them, so five entries including protein_coding_CDS_not_defined ran off the page.
  # Chunk into as many rows as the longest label needs, and draw one legend per row.
  bts       <- sort(unique(panel_tx$biotype_key))
  show_lgnd <- length(bts) > 0
  lgnd_rows <- if (show_lgnd) {
    # Calibrated against rendered panels: a 30-character biotype occupies ~1.4 in of
    # text at fontsize 7, so four entries fit across 7 in and five do not.
    per_row <- max(1L, floor(track_w / (max(nchar(bts)) * 0.048 + 0.30)))
    unname(split(bts, ceiling(seq_along(bts) / per_row)))
  } else list()
  lgnd_h    <- if (show_lgnd) length(lgnd_rows) * 0.20 + 0.14 else 0

  height  <- margin * 2 + title_h + head_h + gene_h + subj_h + tx_h + gap * 3 +
             length(bw_files) * (sig_lab_h + sig_h + gap) + label_h + lgnd_h

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

  # The subject on a row of its own, directly under the host structure it is being
  # compared against, and labelled with the id in the title.
  if (show_subj) {
    fill <- if (is.null(subject_fill)) BIOTYPE_FALLBACK else subject_fill
    # Anchored top-down, not bottom-up: a bottom-anchored label grows upward into the
    # gap above and can touch the structure row.
    plotText(label = subject_label, x = margin, y = y, just = c("left", "top"),
             fontsize = 7, fontcolor = fill, default.units = "inches")
    plotRanges(
      params = pars, data = subject_gr,
      collapse = TRUE, fill = fill, linecolor = NA,
      y = y + 0.13, height = 0.16
    )
    y <- y + subj_h
  }

  # --- Isoform track, one of two renderers (see TX_RENDERERS) -----------------
  #
  # Colours come from `lab`, the same transcript -> label mapping the sub-GTF was
  # rewritten with, so a name that genuinely needed a uniquifying suffix still matches.
  # colorbyStrand is off in both: it would otherwise decide the colour of anything the
  # highlight table misses, and strand is already on the label.
  panel_lab <- unname(lab[panel_tx$transcript_id])
  panel_col <- biotype_color(panel_tx$biotype_key)

  # Genomic coordinate to page inches, clamped to the window.
  gx <- function(p) {
    margin + (min(max(p, win_s), win_e) - win_s) / max(1, win_e - win_s) * track_w
  }

  if (renderer == "per_model") {
    for (i in seq_len(nrow(panel_tx))) {
      lb <- panel_lab[i]
      if (is.na(lb)) next
      ry <- y + (tx_rows[i] - 1L) * row_h

      plotTranscripts(
        params = pars, y = ry + TX_LABEL_H, height = TX_MODEL_H,
        # Off, or plotgardener's own label lands on top of the one drawn below it and
        # every name appears twice.
        labels = NULL,
        boxHeight = grid::unit(TX_BOX_MM, "mm"), spaceHeight = TX_SPACE_H,
        stroke = TX_STROKE, limitLabel = FALSE,
        fill = c(BIOTYPE_FALLBACK, BIOTYPE_FALLBACK), colorbyStrand = FALSE,
        transcriptFilter     = lb,
        transcriptHighlights = data.frame(transcript = lb, color = panel_col[i],
                                          stringsAsFactors = FALSE)
      )

      # Centred on the visible extent, then clamped inside the track, so a model whose
      # midpoint is off-window -- exactly when packed drops its label -- still gets one.
      lw <- nchar(lb) * TX_LABEL_SIZE * 0.0075
      lx <- min(max((gx(panel_tx$start[i]) + gx(panel_tx$end[i])) / 2,
                    margin + lw / 2),
                margin + track_w - lw / 2)
      plotText(label = lb, x = lx, y = ry, just = c("center", "top"),
               fontsize = TX_LABEL_SIZE, fontcolor = panel_col[i],
               default.units = "inches")
    }
  } else {
    keep <- !is.na(panel_lab)
    hl <- if (any(keep)) {
      data.frame(transcript = panel_lab[keep], color = panel_col[keep],
                 stringsAsFactors = FALSE)
    } else NULL

    plotTranscripts(
      params = pars, y = y, height = tx_h,
      labels = "transcript", fontsize = TX_LABEL_SIZE,
      # grid::unit rather than a bare unit(): grid is a base package and always
      # resolves, where relying on plotgardener to re-export it is an assumption with
      # no upside.
      boxHeight = grid::unit(TX_BOX_MM, "mm"), spaceHeight = TX_SPACE_H,
      stroke = TX_STROKE,
      fill = c(BIOTYPE_FALLBACK, BIOTYPE_FALLBACK), colorbyStrand = FALSE,
      transcriptHighlights = hl,
      # A "+" in the corner means not everything fitted -- raise tx_h.
      limitLabel = TRUE
    )
  }

  y <- y + tx_h + gap
  for (i in seq_along(bw_files)) {
    # Scale and sample name in a band above the trace rather than inside it.
    # plotSignal's own scale = TRUE / label = draw both over the plotting area, where a
    # peak reaching ymax runs straight through them: "[0 - 121]" came out with the 1
    # struck through by the trace.
    # Named, not just bracketed: a bare "[0 - 130]" beside a transcript track reads as
    # though it might be a length in bp. It is the vertical scale -- reads stacked over
    # a single base -- and saying so costs five characters.
    # Rounded because bigWig scores are numeric: an unrounded ymax prints as 130.457.
    plotText(label = sprintf("depth [0 - %s]", format(round(ymax), scientific = FALSE,
                                                      trim = TRUE)),
             x = margin, y = y, just = c("left", "top"),
             fontsize = 7, fontcolor = "grey35", default.units = "inches")
    plotText(label = sample_names[i], x = margin + track_w, y = y,
             just = c("right", "top"), fontsize = 8, default.units = "inches")
    plotSignal(
      params = pars, data = bw_files[i],
      y = y + sig_lab_h, height = sig_h,
      range = c(0, ymax),
      linecolor = "#3B6FB6", fill = "#3B6FB6"
    )
    y <- y + sig_lab_h + sig_h + gap
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
    per_row <- length(lgnd_rows[[1]])
    ly <- y + label_h
    for (r in lgnd_rows) {
      # Width proportional to the entries on this row, so a short final row keeps the
      # same swatch spacing as a full one instead of stretching across the page.
      plotLegend(
        legend = r, fill = biotype_color(r), border = FALSE,
        orientation = "h", fontsize = 7,
        x = margin, y = ly, width = track_w * length(r) / per_row, height = 0.18,
        just = c("left", "top"), default.units = "inches"
      )
      ly <- ly + 0.20
    }
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

    # Once per renderer. The CSV keeps the unsuffixed file, so the report is unaffected
    # by the comparison being run.
    figs <- vapply(names(TX_RENDERERS), function(rend) draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      # The union of every exon the gene has, so a reader sees the full exonic
      # footprint before the isoforms start differing from it.
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$gene_id]),
      panel_tx     = panel_tx,
      renderer     = rend,
      title        = row$label,
      subtitle     = sprintf("%d isoforms: %d known, %d novel (%d novel lncRNA candidate%s)",
                             row$n_tx, row$n_known, row$n_novel, row$n_novel_lnc,
                             if (row$n_novel_lnc == 1) "" else "s"),
      out_png      = file.path(opt$outdir, sprintf("genomic_context_%s%s.png",
                                                   row$label, TX_RENDERERS[[rend]]))
    ), character(1))
    cand$figure[i] <- figs[["packed"]]
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

    figs <- vapply(names(TX_RENDERERS), function(rend) draw_panel(
      chrom        = row$chrom, win_s = row$win_s, win_e = row$win_e,
      structure_gr = GenomicRanges::reduce(exons_gr[exon_gid == row$host_gene_id]),
      panel_tx     = panel_tx,
      renderer     = rend,
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
                               sprintf("intronic_context_%s%s.png",
                                       row$qry_id, TX_RENDERERS[[rend]])),
      structure_label = "host",
      # Its own row, so the transcript named in the title is the one the eye lands on.
      subject_gr    = exons_gr[exon_txid == row$qry_id],
      subject_label = row$qry_id,
      subject_fill  = biotype_color(
        tx$biotype_key[match(row$qry_id, tx$transcript_id)])
    ), character(1))
    flagged$figure[i] <- figs[["packed"]]
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

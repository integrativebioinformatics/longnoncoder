#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: load-libraries
#| include: false

library(tidyverse)
library(cowplot)
library(scales)
library(RColorBrewer)
library(viridis)
library(showtext)

font_add_google("Source Sans 3", "Source Sans 3")
showtext_auto()

base_theme <- theme(
  strip.text = element_text(face = "bold"),
  strip.background = element_rect(fill = "#f0f0f0", color = "black"),
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
)

theme_set(theme_bw(base_family = "Source Sans 3", base_size = 12) + base_theme)
update_geom_defaults("text", list(family = "Source Sans 3", size = 10 /.pt))
update_geom_defaults("label", list(family = "Source Sans 3", size = 11))
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: get-overall-assembly-stats
#| include: false

# List bambu outputs to account for assembled and quantified genes/transcripts
bambu_outputs <- c(
  genes = params$counts_genes,
  transcripts = params$counts_transcript,
  fullLength = params$fl_counts_transcript,
  uniquelyMapped = params$u_counts_transcript
)

# Check if files exist
files_exist <- sapply(bambu_outputs, file.exists)
if (!all(files_exist)) {
  missing_files <- bambu_outputs[!files_exist]
  stop("Error: The following files were not found: ", paste(missing_files, collapse = ", "))
}

# Store results in a list
results_list <- list()
for (name in names(bambu_outputs)) {
  current_file <- bambu_outputs[[name]]
  tryCatch({
    data <- read.delim(current_file, header = TRUE, sep = "\t")
    id_column <- data[[1]]
    total_count <- nrow(data)
    annotated_count <- sum(startsWith(id_column, "ENS"))
    novel_count <- sum(startsWith(id_column, "Bambu"))
    file_result <- data.frame(Bambu = name, Total = total_count, Annotated = annotated_count, Novel = novel_count)
    results_list[[name]] <- file_result
  }, error = function(e) {
    cat("Warning: Could not process file '", basename(current_file), "' (named '", name, "'). Skipping.\n", sep = "")
    cat("Error message:", e$message, "\n")
  })
}

# Organize results
if (length(results_list) > 0) {
  bambu_df <- do.call(rbind, results_list)
  rownames(bambu_df) <- NULL
  bambu_df <- bambu_df %>%
    mutate(FeatureType = if_else(Bambu == "genes", "Gene", "Transcript")) %>%
    pivot_longer(cols = c(Total, Annotated, Novel), names_to = "Status", values_to = "Count")
  write.csv(bambu_df, "Bambu_assembly_summary.csv", row.names = FALSE)
} else {
  cat("No files were successfully processed from the provided list.\n")
}
#
#
#
#| label: prepare-to-plot-overall-assembly-stats
#| include: false

if (!exists("bambu_df")) {
  stop("Error: The data frame 'bambu_df' was not found.")
}
custom_labels <- c("genes" = "Genes", "transcripts" = "All transcripts", "fullLength" = "Full length", "uniquelyMapped" = "Uniquely mapped")
bambu_df$Status <- factor(bambu_df$Status, levels = c("Total", "Annotated", "Novel"))
bambu_df$Bambu <- factor(bambu_df$Bambu, levels = rev(c("genes", "transcripts", "fullLength", "uniquelyMapped")))
#
#
#
#| label: plot-bambu-gene-transcript-stats
#| fig-cap: "**Figure 1: Summary of Assembled Genes and Transcripts.** This figure displays the total counts of genes (A) and various categories of transcripts (B) identified by Bambu. Counts are broken down into 'Total', 'Annotated' (known), and 'Novel' (newly discovered). Transcript counts are further detailed for all transcripts, only full-length and uniquely mapped transcripts."
#| fig-width: 10
#| fig-height: 5

bambu_gene_data <- bambu_df %>% filter(FeatureType == "Gene")
bambu_gene_plot <- ggplot(bambu_gene_data, aes(x = Status, y = Count, fill = Status)) +
  geom_col(show.legend = FALSE, position = "dodge") +
  geom_text(aes(label = Count), vjust = -0.3) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Genes", x = "\nStatus", y = "No. of genes\n") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

bambu_transcript_data <- bambu_df %>% filter(FeatureType == "Transcript")
bambu_transcript_data$Bambu <- factor(bambu_transcript_data$Bambu, levels = names(custom_labels)[-1])
levels(bambu_transcript_data$Bambu) <- custom_labels[-1]

bambu_transcript_plot <- ggplot(bambu_transcript_data, aes(x = Bambu, y = Count, fill = Status)) +
  scale_fill_brewer(palette = "Set1") +
  geom_col(position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3) +
  labs(title = "Transcripts", x = "\nAssembly", y = "No. of transcripts\n", fill = "Status") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

bambu_plot <- plot_grid(bambu_gene_plot, bambu_transcript_plot, ncol = 2, rel_widths = c(0.7, 1.3), labels = "AUTO", label_y = 0.85)
bambu_plot
#
#
#
#
#
#| label: get-lnc-pc-assembly-stats
#| include: false

# --- 1A. Check if Annotated Transcriptome Metadata exists ---
if (!file.exists(params$transcriptome_meta)) {
  stop("Error: The annotated transcriptome metadata file was not found at: ", params$transcriptome_meta)
}

# --- 1B. Check if Novel lncRNA/PC Metadata exists ---
if (!file.exists(params$novel_transcriptome_meta)) {
  stop("Error: The novel transcript metadata file was not found at: ", params$novel_transcriptome_meta)
}

# --- 2A. Import Annotated Transcriptome Metadata ---
tryCatch({
  ann_txome <- read.csv(params$transcriptome_meta)
}, error = function(e) {
  stop("Critical Error reading annotated transcriptome metadata file: ", e$message)
})

# --- 2B. Import Novel lncRNA/PC Metadata ---
tryCatch({
  novel_tx <- read.csv(params$novel_transcriptome_meta)
}, error = function(e) {
  stop("Critical Error reading novel transcript metadata file: ", e$message)
})

# Standardize biotypes in novel_tx based on prediction
novel_tx$biotype[novel_tx$prediction == "coding"] <- "novel protein_coding"
novel_tx$biotype[novel_tx$prediction == "non-coding"] <- "novel lncRNA"

# Separate the raw transcript IDs first
raw_nov_pc_tx <- unique(novel_tx$qry_id[novel_tx$prediction == "coding"])
raw_nov_lnc_tx <- unique(novel_tx$qry_id[novel_tx$prediction == "non-coding"])

# --- 3. Read Bambu IDs from the output files ---
bambu_ids <- list()
for (name in names(bambu_outputs)) {
  if (name == "genes") next # Skip genes file for this section
  current_file <- bambu_outputs[[name]]
    # Only read the first column (IDs)
    data_ids <- read.delim(current_file, header = TRUE, sep = "\t")[[1]]
    bambu_ids[[name]] <- as.character(data_ids)
}

# Identify if annotated transcript IDs contain versioning suffix
if (any(grepl("\\.", bambu_ids$transcripts))) {
  version_suffix <- TRUE
} else {
  version_suffix <- FALSE
}

# --- 4. Extract ID vectors for each category (Transcripts Only) ---
ids <- list(
  # Protein-coding
  ann_pc_tx = if (version_suffix) {
    unique(ann_txome$ensembl_transcript_id_version[ann_txome$gene_biotype == "protein_coding"])
  } else {
    unique(ann_txome$ensembl_transcript_id[ann_txome$gene_biotype == "protein_coding"])
  },
  
  # lncRNAs
  ann_lnc_tx = if (version_suffix) {
    unique(ann_txome$ensembl_transcript_id_version[ann_txome$gene_biotype == "lncRNA"])
  } else {
    unique(ann_txome$ensembl_transcript_id[ann_txome$gene_biotype == "lncRNA"])
  },

  # Filter to only include true novel protein-coding features (starting with "Bambu")
  nov_pc_tx = raw_nov_pc_tx[startsWith(raw_nov_pc_tx, "Bambu")],
  
  # Filter to only include true novel lncRNA features (starting with "Bambu")
  nov_lnc_tx = raw_nov_lnc_tx[startsWith(raw_nov_lnc_tx, "Bambu")]
)

# Helper function to count overlaps between Bambu outputs and our target IDs
count_overlap <- function(bambu_set, target_ids) {
  if (is.null(bambu_set) || length(bambu_set) == 0) return(0)
  sum(bambu_set %in% target_ids)
}

# --- 5. Quantify across categories ---
results_biotype <- list()

for (bt in c("pc", "lnc")) {
  biotype_label <- ifelse(bt == "pc", "Protein-coding", "lncRNA")
  
  # Transcript Level (Check all three transcript files)
  for (tx_cat in c("transcripts", "fullLength", "uniquelyMapped")) {
    ann_t <- count_overlap(bambu_ids[[tx_cat]], ids[[paste0("ann_", bt, "_tx")]])
    nov_t <- count_overlap(bambu_ids[[tx_cat]], ids[[paste0("nov_", bt, "_tx")]])
    tot_t <- ann_t + nov_t
    
    results_biotype[[paste0(bt, "_", tx_cat)]] <- data.frame(
      FeatureType = "Transcript", Biotype = biotype_label, Bambu = tx_cat,
      Total = tot_t, Annotated = ann_t, Novel = nov_t
    )
  }
}

# --- 6. Organize and Save Results ---
if (length(results_biotype) > 0) {
  biotype_df <- do.call(rbind, results_biotype) %>%
    pivot_longer(cols = c(Total, Annotated, Novel), names_to = "Status", values_to = "Count")

  # Ensure proper factor ordering
  biotype_df$Status <- factor(biotype_df$Status, levels = c("Total", "Annotated", "Novel"))
  biotype_df$Bambu <- factor(biotype_df$Bambu, levels = rev(c("transcripts", "fullLength", "uniquelyMapped")))

  write.csv(biotype_df, "Bambu_lncRNA_PC_summary.csv", row.names = FALSE)
} else {
  cat("No files were successfully processed to build the biotype summary.\n")
}
#
#
#
#| label: plot-bambu-lnc-pc-stats
#| fig-cap: "**Figure 1.1: Summary of Assembled Protein-coding and lncRNA Transcripts.** Total counts of lncRNA and protein-coding transcripts identified by Bambu. Counts are broken down into 'Total', 'Annotated' (known), and 'Novel' (newly discovered). Transcript counts are further detailed for all transcripts, only full-length and uniquely mapped transcripts."
#| fig-width: 10
#| fig-height: 5

# 1. Data Preparation
custom_labels <- c("transcripts" = "Total", "fullLength" = "Full length", "uniquelyMapped" = "Uniquely mapped")

# Transcript Data & Factor Leveling
bambu_transcript_biotype_data <- biotype_df %>% 
  filter(FeatureType == "Transcript")

bambu_transcript_biotype_data$Bambu <- factor(
  bambu_transcript_biotype_data$Bambu, 
  levels = names(custom_labels)
)
levels(bambu_transcript_biotype_data$Bambu) <- custom_labels

# 2. Build Transcript Plot
lnc_pc_bambu_plot <- ggplot(bambu_transcript_biotype_data, aes(x = Bambu, y = Count, fill = Status)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.3) +
  scale_fill_brewer(palette = "Set1") +
  facet_wrap(~ Biotype, scales = "free_y") +
  labs(title = "Transcripts", x = "\nAssembly", y = "No. of transcripts\n", fill = "Status") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(margin = margin(b = 10)),
    strip.text = element_text(face = "bold")
  )

# Final Output
lnc_pc_bambu_plot
#
#
#
#
#
#
#
#
#
#| label: plot-availability
#| include: false

plot_params <- list(
  raw_pca                 = params$raw_pca,
  raw_pca_grouped         = params$raw_pca_grouped,
  raw_heatmap_gene        = params$raw_heatmap_gene,
  raw_heatmap_transcript  = params$raw_heatmap_transcript,
  post_pca                = params$post_pca,
  post_pca_grouped        = params$post_pca_grouped,
  post_heatmap_gene       = params$post_heatmap_gene,
  post_heatmap_transcript = params$post_heatmap_transcript
)

# A plot is only embedded if the file is present and non-empty, so a stub run or
# a partially completed execution does not break the render.
plot_available <- vapply(plot_params, function(p) {
  !is.null(p) && file.exists(p) && file.info(p)$size > 0
}, logical(1))

# include_graphics must be handed a vector: called inside a loop it would not be
# auto-printed and nothing would be embedded.
show_plots <- function(keys) {
  paths <- unlist(plot_params[keys][plot_available[keys]], use.names = FALSE)
  if (length(paths) > 0) knitr::include_graphics(paths) else NULL
}
#
#
#
#
#
#| label: fig-pca
#| echo: false
#| fig-cap: "Sample PCA before (left) and after (right) validation and curation."
#| layout-ncol: 2

show_plots(c("raw_pca", "post_pca"))
#
#
#
#| label: fig-pca-grouped
#| echo: false
#| fig-cap: "Sample PCA coloured by group, before (left) and after (right) validation and curation."
#| layout-ncol: 2

show_plots(c("raw_pca_grouped", "post_pca_grouped"))
#
#
#
#
#
#| label: fig-heatmap-transcript
#| echo: false
#| fig-cap: "Transcript-level expression heatmap, before (left) and after (right) validation and curation."
#| layout-ncol: 2

show_plots(c("raw_heatmap_transcript", "post_heatmap_transcript"))
#
#
#
#| label: fig-heatmap-gene
#| echo: false
#| fig-cap: "Gene-level expression heatmap, before (left) and after (right) validation and curation."
#| layout-ncol: 2

show_plots(c("raw_heatmap_gene", "post_heatmap_gene"))
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: data-prep-annotated
#| include: false

# Load and prepare data
summary_data <- ann_txome %>% 
  group_by(gene_biotype, transcript_biotype) %>% 
  summarise(count = n_distinct(ensembl_transcript_id_version), .groups = 'drop')

# Categorize biotypes
protein_coding <- summary_data %>% filter(gene_biotype == "protein_coding")
lncRNAs <- summary_data %>% filter(gene_biotype == "lncRNA")
pseudogene_RNAs <- summary_data %>% filter(str_detect(gene_biotype, "pseudogene")) %>% mutate(category = "Pseudogene")
small_RNAs <- summary_data %>% filter(gene_biotype %in% c("snoRNA", "snRNA", "scaRNA", "sRNA", "rRNA", "misc_RNA", "miRNA"))
Tcell_RNAs <- summary_data %>% filter(gene_biotype %in% c("TR_V_gene", "TR_J_gene", "TR_D_gene", "TR_C_gene")) %>% mutate(category = "TCR elements")
IG_RNAs <- summary_data %>% filter(gene_biotype %in% c("IG_V_gene", "IG_J_gene", "IG_D_gene", "IG_C_gene")) %>% mutate(category = "Immunoglobulin elements")

# Gather all biotypes that have already been classified above
classified_biotypes <- unique(c(
  protein_coding$gene_biotype,
  lncRNAs$gene_biotype,
  pseudogene_RNAs$gene_biotype,
  small_RNAs$gene_biotype,
  Tcell_RNAs$gene_biotype,
  IG_RNAs$gene_biotype
))

# Dynamically filter for ANY biotype that isn't in the classified list
other_RNAs_explicit <- summary_data %>% 
  filter(!gene_biotype %in% classified_biotypes) %>% 
  mutate(category = "all other biotypes")

# Group all other RNAs for the aggregated plot
all_other_RNAs <- bind_rows(pseudogene_RNAs, Tcell_RNAs, IG_RNAs, other_RNAs_explicit) %>% 
  group_by(category) %>% 
  summarise(count = sum(count), .groups = 'drop')

all_other_RNAs$category <- factor(all_other_RNAs$category, levels = c("Pseudogene", "TCR elements", "Immunoglobulin elements", "all other biotypes"))

# Plotting function
create_plot <- function(data, x, fill, title) {
  ggplot(data, aes(x = {{x}}, y = count, fill = {{fill}})) +
    geom_bar(stat = "identity", position = "dodge2") +
    labs(x = "\nGene biotypes", y = "Count (log10 scale)\n", fill = "Transcript biotypes", title = title) +
    theme(
      panel.background = element_rect(fill = "white"), 
      axis.text.x = element_text(hjust = 0.5), 
      legend.position = "right",
      legend.text = element_text(size = 12), 
      legend.title = element_text(face = "bold", size = 12), 
      plot.title = element_text(hjust = 0.5, face = "bold"), 
      axis.title.y = element_text(face = "bold"), 
      axis.title.x = element_text(face = "bold"), 
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = comma) +
    geom_text(aes(label = count), position = position_dodge2(width = 0.9), vjust = -0.5, hjust = 0.5)
}

# Generate plots
p1 <- create_plot(protein_coding, gene_biotype, transcript_biotype, "Protein-coding RNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 <- create_plot(lncRNAs, gene_biotype, transcript_biotype, "LncRNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p3 <- create_plot(small_RNAs, gene_biotype, transcript_biotype, "Small RNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p4 <- create_plot(all_other_RNAs, category,  category, "All Other RNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p5 <- create_plot(pseudogene_RNAs, gene_biotype, transcript_biotype, "Pseudogene RNA diversity") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p6 <- create_plot(other_RNAs_explicit, gene_biotype, transcript_biotype, "Other RNAs (detailed)") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p7 <- create_plot(Tcell_RNAs, gene_biotype, transcript_biotype, "TCR element RNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p8 <- create_plot(IG_RNAs, gene_biotype, transcript_biotype, "Immunoglobulin element RNAs") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
#
#
#
#
#
#
#
#| label: plot-panel-1
#| fig-cap: "**Figure 2: Transcript Counts for Major RNA Categories.** This panel shows the distribution of transcript biotypes within four main gene biotype categories: (A) Protein-coding RNAs, (B) Long non-coding RNAs (lncRNAs), (C) Small RNAs, and (D) a consolidated group of all other RNA types."
#| fig-width: 16
#| fig-height: 10

panel1 <- plot_grid(p1, p2, p3, p4, ncol = 2, rel_widths = c(1.3, 1), labels = "AUTO", label_y = 0.85)
panel1
#
#
#
#
#
#
#
#| label: plot-panel-2
#| fig-cap: "**Figure 3: Detailed Breakdown of Other RNA Categories.** This figure provides a closer look at specific RNA families: (A) diversity within Pseudogene RNAs, (B) T-cell receptor (TCR) elements, (C) Immunoglobulin (IG) elements, and (D) other variety of RNA types."
#| fig-width: 16
#| fig-height: 10

top_row_panel2 <- plot_grid(p5, ncol = 1, labels = "A", label_y = 0.85)
bottom_row_panel2 <- plot_grid(p7, p8, p6, ncol = 3, labels = c("B", "C", "D"), label_y = 0.85)
panel2 <- plot_grid(top_row_panel2, bottom_row_panel2, nrow = 2, rel_heights = c(1, 1))
panel2
#
#
#
#
#
#
#
#
#
#| label: prep-chr-data
#| include: false

chr_gb <- ann_txome %>%
  mutate(gene_biotype = if_else(str_detect(gene_biotype, "pseudogene"), "pseudogene", gene_biotype)) %>%
  mutate(gene_biotype = if_else(gene_biotype %in% c("snoRNA", "snRNA", "scaRNA", "sRNA", "rRNA", "misc_RNA", "miRNA"), "small_RNA", gene_biotype)) %>%
  mutate(gene_biotype = if_else(gene_biotype != "pseudogene" & gene_biotype != "small_RNA" & gene_biotype != "protein_coding" & gene_biotype != "lncRNA", "other", gene_biotype)) %>%
  group_by(chromosome_name, gene_biotype) %>%
  summarise(count = n(), .groups = "drop")
chr_gb$gene_biotype <- factor(chr_gb$gene_biotype, levels = c("protein_coding", "lncRNA", "small_RNA", "pseudogene", "other"))
chromosome_names <- unique(chr_gb$chromosome_name)
numeric_chromosomes <- grep("^[0-9]+$", chromosome_names, value = TRUE)
non_numeric_chromosomes <- setdiff(chromosome_names, numeric_chromosomes)
numeric_chromosomes_sorted <- as.character(sort(as.numeric(numeric_chromosomes)))
chromosome_order <- c(numeric_chromosomes_sorted, sort(non_numeric_chromosomes))
chr_gb$chromosome_name <- factor(chr_gb$chromosome_name, levels = chromosome_order)
#
#
#
#| label: plot-gb-across-chr
#| fig-cap: "**Figure 4: Chromosomal Distribution of Gene Biotypes.** Each panel represents a chromosome, showing the count of transcript isoforms for each biotype found there. This provides a genome-wide view of how different functional classes of genes are organized and transcribed across chromosomes."
#| fig-height: 30
#| fig-width: 35

p9 <- ggplot(chr_gb, aes(x = gene_biotype, y = count)) +
  geom_bar(stat = "identity", aes(fill = gene_biotype), alpha = 1, width = 0.8) +
  facet_wrap(~ chromosome_name, scales = "free_x", ncol = 5) +
  scale_fill_viridis_d(option = "viridis", guide = "none") +
  scale_y_log10(labels = scales::comma_format(), expand = expansion(mult = c(0, 0.2))) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey90"), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 5, b = 10)),
    axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", color = "black", margin = margin(t = 8, b = 8, l = 5, r = 5)),
    strip.background = element_rect(fill = "grey85", color = "grey70", linewidth = 0.5),
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 40, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"), panel.spacing = unit(1.2, "cm")
  ) +
  labs(title = "Gene Biotype Diversity Across Human Chromosomes", subtitle = "Distribution of different gene types per chromosome (log10 scale)", x = "Gene Biotype", y = "Count (log10 scale)") +
  geom_text(aes(label = count), vjust = -0.5, color = "black")
p9
#
#
#
#
#
#
#
#| label: prep-tx-length-data
#| include: false

selected_tx <- ann_txome %>% filter(gene_biotype %in% c("protein_coding", "lncRNA", "miRNA", "snoRNA", "snRNA"))
selected_tx$gene_biotype <- factor(selected_tx$gene_biotype, levels = c("protein_coding", "lncRNA", "miRNA", "snoRNA", "snRNA"))
median_stx <- selected_tx %>% group_by(gene_biotype) %>% summarise(median_length = median(transcript_length), .groups = 'drop')
median_pc <- selected_tx %>% filter(gene_biotype == "protein_coding") %>% group_by(transcript_biotype) %>% summarise(median_length_txb = median(transcript_length), .groups = 'drop')
median_lnc <- selected_tx %>% filter(gene_biotype == "lncRNA") %>% group_by(transcript_biotype) %>% summarise(median_length_txb = median(transcript_length), .groups = 'drop')

p10 <- ggplot(selected_tx, aes(x = gene_biotype, y = transcript_length, fill = factor(gene_biotype))) +
  geom_boxplot() + labs(title = "Transcript length across gene biotypes", x = "Gene biotype", y = "Transcript length (log10 scale)", fill = "Gene biotype\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_stx$gene_biotype, " (", round(median_stx$median_length, 2), ")"))

p11 <- selected_tx %>% filter(gene_biotype == "protein_coding") %>%
  ggplot(aes(x = gene_biotype, y = transcript_length, fill = transcript_biotype)) +
  geom_boxplot() + labs(title = "Transcript length across protein-coding biotypes", x = "Gene biotype", y = "Transcript length (log10 scale)", fill = "Transcript biotypes\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_pc$transcript_biotype, " (", round(median_pc$median_length_txb, 2), ")"))

p12 <- selected_tx %>% filter(gene_biotype == "lncRNA") %>%
  ggplot(aes(x = gene_biotype, y = transcript_length, fill = transcript_biotype)) +
  geom_boxplot() + labs(title = "Transcript length across lncRNA biotypes", x = "Gene biotype", y = "Transcript length (log10 scale)", fill = "Transcript biotypes\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_lnc$transcript_biotype, " (", round(median_lnc$median_length_txb, 2), ")"))
#
#
#
#| label: plot-panel3
#| fig-cap: "**Figure 5: Transcript Length Distributions.** This panel of boxplots compares transcript lengths (on a log10 scale) across different biotypes. (A) Compares major gene biotypes. (B) Details the lengths of different transcript types within the protein-coding category. (C) Does the same for lncRNAs. Median lengths for each category are shown in the legend."
#| fig-width: 8
#| fig-height: 15

panel3 <- plot_grid(p10, p11, p12, nrow = 3, labels="AUTO", label_y = 1, rel_widths = c(1.5,2,1.5))
panel3
#
#
#
#
#
#
#
#| label: prep-exon-data
#| include: false

ann_pc_txome <- read.csv(params$protein_coding_meta)
ann_lnc_txome <- read.csv(params$lncrna_meta)
ann_pc_lnc_txome <- rbind(ann_pc_txome, ann_lnc_txome)

median_exnum <- ann_pc_lnc_txome %>% group_by(gene_biotype) %>% summarise(median_exnum = median(num_exons), .groups = 'drop')
p13 <- ggplot(ann_pc_lnc_txome, aes(x = gene_biotype, y = num_exons, fill = gene_biotype)) +
  geom_boxplot() + labs(title = "Number of exons across protein-coding and lncRNA biotypes", x = "Gene biotype", y = "Exon number (log10 scale)", fill = "Gene biotype\n(median number)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_exnum$gene_biotype, " (", round(median_exnum$median_exnum, 2), ")"))

median_pc_exnum <- ann_pc_txome %>% group_by(transcript_biotype) %>% summarise(median_exnum = median(num_exons), .groups = 'drop')
median_lnc_exnum <- ann_lnc_txome %>% group_by(transcript_biotype) %>% summarise(median_exnum = median(num_exons), .groups = 'drop')
p14 <- ggplot(ann_pc_txome, aes(x = gene_biotype, y = num_exons, fill = transcript_biotype)) + geom_boxplot() + theme(plot.margin = margin(1, 1, 1, 1, "cm")) + labs(title = "Number of exons across protein-coding transcript biotypes", x = "Gene biotype", y = "Exon number (log10 scale)", fill = "Transcript biotype\n(total median number)") + scale_y_log10() + scale_fill_discrete(labels = paste0(median_pc_exnum$transcript_biotype, " (", round(median_pc_exnum$median_exnum, 2), ")"))
p15 <- ggplot(ann_lnc_txome, aes(x = gene_biotype, y = num_exons, fill = transcript_biotype)) + geom_boxplot() + theme(plot.margin = margin(1, 1, 1, 1, "cm")) + labs(title = "Number of exons across lncRNA transcript biotypes", x = "Gene biotype", y = "Exon number (log10 scale)", fill = "Transcript biotype\n(total median number)") + scale_y_log10() + scale_fill_discrete(labels = paste0(median_lnc_exnum$transcript_biotype, " (", round(median_lnc_exnum$median_exnum, 2), ")"))
#
#
#
#
#
#| label: plot-panel4
#| fig-cap: "**Figure 6: Exon Number Distributions.** This panel examines the number of exons per transcript. (A) A direct comparison between protein-coding and lncRNA gene biotypes. (B) A detailed view of exon numbers for different protein-coding transcript types. (C) A similar detailed view for lncRNA transcript types."
#| fig-width: 8
#| fig-height: 15

panel4 <- plot_grid(
  p13, p14, p15, 
  nrow = 3, 
  labels = "AUTO", 
  label_y = 1, 
  rel_heights = c(1, 1.4, 1), 
  align = "v",                
  axis = "lr"                 
)

panel4
#
#
#
#
#
#| label: plot-mono-multi-exon-analysis
#| include: false

ann_pc_lnc_txome$class_exon <- ifelse(ann_pc_lnc_txome$num_exons == 1, "mono-exonic", "multi-exonic")
count_momu <- ann_pc_lnc_txome %>% group_by(gene_biotype, transcript_biotype, class_exon) %>% summarise(count = n_distinct(ensembl_transcript_id_version), .groups = 'drop')

p16 <- ggplot(count_momu, aes(x = transcript_biotype, y = count, fill = transcript_biotype)) +
  geom_bar(stat = "identity") + facet_grid(gene_biotype ~ class_exon, scales = "free_x", space = "free_x") +
  labs(title = "Distribution of mono and multi-exonic transcripts", x = "\nTranscript biotype", y = "Count (log10 scale)", fill = "Transcript biotype") +
  theme(,
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_log10() + geom_text(aes(label = count), vjust = -0.5)

p17 <- ggplot(ann_pc_lnc_txome, aes(x = gene_biotype, y = transcript_length, fill = class_exon)) +
  geom_boxplot() + labs(title = "Mono and multi-exonic transcript length", x = "Gene biotype", y = "Transcript length (log10 scale)", fill = "Exon class") +
  scale_y_log10()
#
#
#
#| fig-cap:
#|   - "**Figure 7: Counts of Mono- vs. Multi-exonic Transcripts.** This plot shows the number of mono- and multi-exonic transcripts for each transcript biotype, separated by protein-coding and lncRNA gene categories."
#| fig-width: 10
#| fig-height: 10 

p16
#
#
#
#| fig-cap:
#|   - "**Figure 8: Length of Mono- vs. Multi-exonic Transcripts.** This boxplot compares the final transcript length of mono- and multi-exonic transcripts within the protein-coding and lncRNA categories."

p17
#
#
#
#
#
#
#
#| label: plot-density
#| fig-cap: "**Figure 9: Relationship Between Exon Count and Transcript Length.** This plot shows the density of transcripts based on their exon count (x-axis) and total length (y-axis). Higher density (brighter color) indicates a common combination of exon number and length. The analysis is shown separately for protein-coding genes and lncRNAs."
#| fig-width: 8
#| fig-height: 6

p18 <- ggplot(ann_pc_lnc_txome, aes(x = num_exons, y = transcript_length)) +
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.8) +
  scale_x_log10() + scale_y_log10() + scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~gene_biotype, scales = "free_x") +
  labs(x = "Number of exons (log10 scale)", y = "Transcript length (log10 scale)", title = "Relationship between number of exons and transcript length")
p18
#
#
#
#
#
ann_pc_exon <- read.csv(params$protein_coding_exonlength)
ann_pc_exon$gene_biotype <- "protein_coding"

ann_lnc_exon <- read.csv(params$lncrna_exonlength)
ann_lnc_exon$gene_biotype <- "lncRNA"

ann_pc_lnc_exon <- rbind(ann_pc_exon, ann_lnc_exon)

ann_pc_lnc_exon <- merge(
  ann_pc_lnc_exon,
  ann_pc_lnc_txome[, c("ensembl_transcript_id_version", "class_exon")],
  by = "ensembl_transcript_id_version",
  all.x = TRUE
)
#
#
#
#| label: plot-indiv-exlen
#| fig-cap: "**Figure 10: Length of individual exons.** This boxplot compares the individual exon length of protein-coding and lncRNA transcripts, separeted by mono and multi-exonic structures."
#| fig-width: 8
#| fig-height: 6

p19 <- ggplot(ann_pc_lnc_exon, aes(x = gene_biotype, y = width,
                                   fill = class_exon)) +
  geom_boxplot() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(title = "Mono and multi-exon individual length across protein-coding\nand lncRNA gene biotypes",
       x = "Gene biotype",
       y = "Exon length (log10 scale)",
       fill = "Exon class") +  # Set legend title
  scale_y_log10() # Add log10 scale to y-axis

p19
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: load-data-novel
#| include: false

class_counts <- novel_tx %>% group_by(biotype, classification) %>% summarise(count = n(), .groups = 'drop')
#
#
#
#| label: plot-novel-biotypes-and-classes
#| fig-cap: "**Figure 11: GFFCompare Classification of Novel Transcripts.** This bar chart shows the counts of novel protein-coding and lncRNA transcripts according to their GFFCompare classification code."
#| fig-width: 10
#| fig-height: 7

p20 <- ggplot(class_counts, aes(x = classification, y = count, fill = classification)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count, y = count), vjust = -0.5, position = position_dodge(width = 0.9), angle = 0) +
  facet_grid(~ biotype, scales = "free_x", space = "free_x") +
  labs(title = "GFFCompare classification stats overview of novel protein-coding and lncRNAs", x = "\n\nClassification", y = "Count (log10)", fill = "Classification") +
  theme(,
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust =1),
        legend.position = "right") +
scale_y_log10()
p20
#
#
#
#
#
#
#
#
#
#| label: prep-novel-chr-data
#| include: false

chr_tb_novel <- novel_tx %>% group_by(seqnames, biotype) %>% summarise(count = n(), .groups = 'drop')
chromosome_names_novel <- unique(chr_tb_novel$seqnames)
numeric_chromosomes_novel <- grep("^[0-9]+$", chromosome_names_novel, value = TRUE)
non_numeric_chromosomes_novel <- setdiff(chromosome_names_novel, numeric_chromosomes_novel)
numeric_chromosomes_sorted_novel <- as.character(sort(as.numeric(numeric_chromosomes_novel)))
chromosome_order_novel <- c(numeric_chromosomes_sorted_novel, sort(non_numeric_chromosomes_novel))
chr_tb_novel$seqnames <- factor(chr_tb_novel$seqnames, levels = chromosome_order_novel)
#
#
#
#| label: plot-novel-across-chr
#| fig-cap: "**Figure 12: Chromosomal Distribution of Novel Transcripts.** Each panel represents a chromosome and shows the number of novel protein-coding and novel lncRNA transcripts discovered there."
#| fig-width: 25
#| fig-height: 30

p21 <- ggplot(chr_tb_novel, aes(x = biotype, y = count)) +
  geom_bar(stat = "identity", aes(fill = biotype), alpha = 1, width = 0.8) +
  facet_wrap(~ seqnames, scales = "free_x", nrow = 5) +
  scale_fill_viridis_d(option = "viridis", guide = "none") +
  scale_y_log10(labels = scales::comma_format(), expand = expansion(mult = c(0, 0.2))) +
  theme(
    panel.background = element_rect(fill = "white", color = "grey90"), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, margin = margin(t = 5, b = 10)),
    axis.title.x = element_text(face = "bold"), axis.title.y = element_text(face = "bold"),
    strip.text = element_text(face = "bold", color = "black", margin = margin(t = 8, b = 8, l = 5, r = 5)),
    strip.background = element_rect(fill = "grey85", color = "grey70", linewidth = 0.5),
    plot.title = element_text(size = 50, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 40, hjust = 0.5),
    plot.margin = margin(1, 1, 1, 1, "cm"), panel.spacing = unit(1.0, "cm"), axis.ticks = element_blank()
  ) +
  labs(title = "Novel RNAs Diversity Across Chromosomes", subtitle = "Distribution of novel RNA biotypes per chromosome (log10 scale)", x = "Biotype", y = "Count (log10 scale)") +
  geom_text(aes(label = count), vjust = -0.5, color = "black")
p21
#
#
#
#
#
#
#
#| label: prep-novel-length-exon-data
#| include: false

median_txlen_novel <- novel_tx %>% group_by(biotype) %>% summarise(median_length = median(len), .groups = 'drop')
p22 <- ggplot(novel_tx, aes(x = biotype, y = len, fill = factor(biotype))) +
  geom_boxplot() + labs(title = "Transcript length across novel transcripts", x = "Biotype", y = "Transcript length (log10 scale)", fill = "Biotype\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_txlen_novel$biotype, " (", round(median_txlen_novel$median_length, 2), ")"))

median_pc_novel <- novel_tx %>% filter(biotype == "novel protein_coding") %>% group_by(classification) %>% summarise(median_length_novel = median(len), .groups = 'drop')
p23 <- novel_tx %>% filter(biotype == "novel protein_coding") %>%
  ggplot(aes(x = biotype, y = len, fill = classification)) +
  geom_boxplot() + labs(title = "Transcript length across novel protein-coding sub-biotypes", x = "Biotype", y = "Transcript length (log10 scale)", fill = "Transcript sub-biotypes\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_pc_novel$classification, " (", round(median_pc_novel$median_length_novel, 2), ")"))

median_lnc_novel <- novel_tx %>% filter(biotype == "novel lncRNA") %>% group_by(classification) %>% summarise(median_length_novel = median(len), .groups = 'drop')
p24 <- novel_tx %>% filter(biotype == "novel lncRNA") %>%
  ggplot(aes(x = biotype, y = len, fill = classification)) +
  geom_boxplot() + labs(title = "Transcript length across novel lncRNA sub-biotypes", x = "Biotype", y = "Transcript length (log10 scale)", fill = "Transcript sub-biotypes\n(median length)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_lnc_novel$classification, " (", round(median_lnc_novel$median_length_novel, 2), ")"))
#
#
#
#| label: plot-panel5
#| fig-cap: "**Figure 13: Length Distribution of Novel Transcripts.** (A) Compares the overall length of novel protein-coding vs. novel lncRNA transcripts. (B) and (C) provide a more detailed view, showing length distributions for each GFFCompare classification within the novel protein-coding and lncRNA categories, respectively."
#| fig-width: 10
#| fig-height: 10

panel5 <- plot_grid(p22, p23, p24, nrow = 3, labels="AUTO", label_y = 1)
panel5
#
#
#
#
#
#
#
#
#
#| label: prep-novel-exon-structure
#| include: false

# Classify transcripts
novel_tx$class_exon <- ifelse(novel_tx$num_exons == 1, "mono-exonic", "multi-exonic")

# Count mono- and multi-exonic transcripts
count_momu_novel <- novel_tx %>%
  group_by(biotype, classification, class_exon) %>%
  summarise(count = n_distinct(qry_id), .groups = 'drop')

p25 <- ggplot(count_momu_novel, aes(x = classification, y = count, fill = classification)) +
  geom_bar(stat = "identity") +
  facet_grid(biotype ~ class_exon, scales = "free_x", space = "free_x") +
  labs(title = "Distribution of mono and multi-exonic novel transcripts", x = "Classification",
       y = "Count (log10 scale)", fill = "Classification") +
  theme( 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_y_log10() + geom_text(aes(label = count), vjust = -0.5)

# plot exon numbers

median_exnum_novel <- novel_tx %>% group_by(biotype) %>% summarise(median_exnum = median(num_exons), .groups = 'drop')
p26 <- ggplot(novel_tx, aes(x = biotype, y = num_exons, fill = biotype)) +
  geom_boxplot() + labs(title = "Number of exons across novel transcripts", x = "Biotype", y = "Exon number (log10 scale)", fill = "Biotype\n(median number)") +
  scale_y_log10() + scale_fill_discrete(labels = paste0(median_exnum_novel$biotype, " (", round(median_exnum_novel$median_exnum, 2), ")"))

median_exnum_class_novel <- novel_tx %>% group_by(biotype, classification) %>% summarise(median_exnum = median(num_exons), .groups = 'drop')
p27 <- ggplot(novel_tx, aes(x = classification, y = num_exons, fill = paste(biotype, classification))) +
  geom_boxplot() + 
  facet_wrap(~biotype, nrow = 2, scales = "free_x") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(
    title = "Number of exons across novel transcripts within subtypes", 
    x = "Classification", 
    y = "Exon number (log10 scale)", 
    fill = "Classification\n(total median number)"
  ) +
  scale_y_log10() +
  scale_fill_discrete(labels = paste0(
      median_exnum_class_novel$biotype, " ",
      median_exnum_class_novel$classification, " (",
      round(median_exnum_class_novel$median_exnum, 2), ")"))

# plot exon structure class length

p28 <- ggplot(novel_tx, aes(x = biotype, y = len, fill = class_exon)) +
  geom_boxplot() + labs(title = "Mono and multi-exonic novel transcripts' length", x = "Biotype", y = "Transcript length (log10 scale)", fill = "Exon classification") +
  scale_y_log10()

p29 <- ggplot(novel_tx, aes(x = classification, y = len, fill = classification)) +
  geom_boxplot() + facet_grid(biotype ~ class_exon, scales = "free_x", space = "free_x") +
  labs(title = "Length of mono and multi-exonic novel transcripts", x = "Classification",
       y = "Transcript length (log10 scale)", fill = "Classification") +
  theme( 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "right", 
        plot.margin = margin(1, 1, 1, 1, "cm")) +
  scale_y_log10()
#
#
#
#| fig-cap:
#|   - "**Figure 14: Distribution of mono and multi-exonic structures in Novel Transcripts within its subtypes.** This barplot compares the distribution of mono and multi-exonic structures for novel protein-coding versus novel lncRNA subtypes."
#| fig-width: 10
#| fig-height: 10

p25
#
#
#
#| fig-cap:
#|   - "**Figure 15: Exon Number in Novel Transcripts.** This boxplot compares the number of exons per transcript for novel protein-coding versus novel lncRNA biotypes."


p26
#
#
#
#| fig-cap:
#|   - "**Figure 16: Exon Number in Novel Transcripts within subtypes.** This boxplot compares the number of exons per transcript for novel protein-coding versus novel lncRNA subtypes."
#| fig-width: 10
#| fig-height: 8

p27
#
#
#
#
#
#| fig-cap:
#|   - "**Figure 17: Length of Mono- vs. Multi-exonic Novel Transcripts.** This plot compares the final transcript length of mono- and multi-exonic novel transcripts."

p28
#
#
#
#| fig-cap:
#|   - "**Figure 18: Length of mono and multi-exonic structures in Novel Transcripts within subtypes.** This boxplot compares the length of mono and multi-exonic structures for novel protein-coding versus novel lncRNA subtypes."
#| fig-width: 12
#| fig-height: 8

p29
#
#
#
#
#
#
#
#| label: plot-density-novel
#| fig-cap: "**Figure 19: Exon Count vs. Transcript Length for Novel RNAs.** Similar to Figure 13, this 2D density plot shows the relationship between exon number and transcript length, but this time for the set of newly discovered transcripts."
#| fig-width: 8
#| fig-height: 6

p30 <- ggplot(novel_tx, aes(x = num_exons, y = len)) +
  stat_density2d(aes(fill = after_stat(level)), geom = "polygon", alpha = 0.8) +
  scale_x_log10() + scale_y_log10() + scale_fill_viridis_c(option = "plasma") +
  facet_wrap(~biotype, scales = "free_x") +
  labs(x = "Number of exons (log10 scale)", y = "Transcript length (log10 scale)", title = "Relationship between number of exons and transcript length across novel RNAs")
p30
#
#
#
#
#

exlen_novel_lncRNA <- read.csv(params$novel_lncrna_exonlength)
exlen_novel_pc <- read.csv(params$novel_protein_coding_exonlength)

exlen_novel_lncRNA$biotype <- "novel lncRNA"
exlen_novel_pc$biotype <- "novel protein-coding"

exlen_novel <- rbind(exlen_novel_lncRNA, exlen_novel_pc)

exlen_novel <- merge(
  exlen_novel,
  novel_tx[, c("qry_id", "class_exon")],
  by.x = "transcript_id",
  by.y = "qry_id",
  all.x = TRUE
)
#
#
#
#| label: plot-exlen-novel
#| fig-cap: "**Figure 20: Length of individual exons across novel transcripts.** This boxplot compares the individual exon length of protein-coding and lncRNA novel transcripts, separated by mono and multi-exonic structures."
#| fig-width: 8
#| fig-height: 6

p32 <- ggplot(exlen_novel, aes(x = biotype, y = width,
                                   fill = class_exon)) +
  geom_boxplot() +
  theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
  labs(title = "Mono and multi-exon individual length across novel protein-coding\nand lncRNA biotypes",
       x = "Biotype",
       y = "Exon length (log10 scale)",
       fill = "Exon class") +  # Set legend title
  scale_y_log10() # Add log10 scale to y-axis

p32
#
#
#

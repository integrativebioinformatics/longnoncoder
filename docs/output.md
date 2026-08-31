# Output 

This document describes the file-level reference for outputs produced by each pipeline step, which include:

1.  Read quality control and optional filtering
2.  Mapping/alingment of reads to a reference genome
3.  Transcriptome assembly and quantification
4.  Novel RNA isoform candidates classification and coding-potential prediction
5.  Annotation and summarization of assembled isoforms
6.  Final visual reporting
7.  Run traceability

> *All output files are organized by stage-specific directories.*

## `transcriptome_report/` (final report outputs)

| File | Description |
|----------------------|--------------------------------------------------|
| `Bambu_assembly_summary.csv` | Summary table of transcriptome assembly metrics derived from bambu outputs. |
| `Bambu_lncRNA_PC_summary.csv` | Summary table comparing lncRNA and protein-coding assembly status . |
| `report.html` | Final transcriptome report with integrated summaries and visualizations. |

### Report interface

![pulposeq-report](../figures/pulposeq-report-interface.gif)

## `chopper/` (read filtering and trimming)

| File | Description |
|-------------------|-----------------------------------------------------|
| `filtered_*.fastq.gz` | Filtered/compressed FASTQ reads generated after filtering/trimming. |

## `multiqc/` (aggregated QC)

### `multiqc/`

| File                  | Description                           |
|-----------------------|---------------------------------------|
| `multiqc_report.html` | Main MultiQC interactive HTML report. |

#### Summary statistics

![multiqc_table](../figures/multiqc-nanostat_fastq-qc.png)

#### Phred-score quality plot

![multiqc_plot_qual](../figures/multiqc-nanostat_fastq-qual-dist.png)

### `multiqc/multiqc_data/`

| File | Description |
|-------------------------|------------------------------------------------|
| `multiqc_citations.txt` | Citation information for tools represented in the MultiQC report. |
| `multiqc_data.json` | Structured JSON with parsed metrics used to render MultiQC outputs. |
| `multiqc_general_stats.txt` | Tabular general statistics exported by MultiQC. |
| `multiqc.log` | MultiQC execution log. |
| `multiqc_nanostat.txt` | Parsed NanoStat metrics as imported by MultiQC. |
| `multiqc_software_versions.txt` | Software version summary collected by MultiQC modules. |
| `multiqc_sources.txt` | Source file provenance for metrics included in MultiQC. |
| `nanostat_fastq_stats_table.txt` | NanoStat FASTQ summary table extracted by MultiQC. |
| `nanostat_quality_dist.txt` | Quality distribution table extracted by MultiQC. |

### `multiqc/multiqc_plots/pdf/`

| File | Description |
|----------------------------|--------------------------------------------|
| `general_stats_table.pdf` | PDF export of MultiQC general statistics table. |
| `nanostat_fastq_stats_table.pdf` | PDF export of NanoStat FASTQ statistics table. |
| `nanostat_quality_dist-cnt.pdf` | PDF plot of quality-score counts distribution. |
| `nanostat_quality_dist-pct.pdf` | PDF plot of quality-score percentage distribution. |

### `multiqc/multiqc_plots/png/`

| File | Description |
|----------------------------|--------------------------------------------|
| `general_stats_table.png` | PNG export of MultiQC general statistics table. |
| `nanostat_fastq_stats_table.png` | PNG export of NanoStat FASTQ statistics table. |
| `nanostat_quality_dist-cnt.png` | PNG plot of quality-score counts distribution. |
| `nanostat_quality_dist-pct.png` | PNG plot of quality-score percentage distribution. |

### `multiqc/multiqc_plots/svg/`

| File | Description |
|----------------------------|--------------------------------------------|
| `general_stats_table.svg` | SVG export of MultiQC general statistics table. |
| `nanostat_fastq_stats_table.svg` | SVG export of NanoStat FASTQ statistics table. |
| `nanostat_quality_dist-cnt.svg` | SVG plot of quality-score counts distribution. |
| `nanostat_quality_dist-pct.svg` | SVG plot of quality-score percentage distribution. |

## `nanocomp/` (read quality comparison)

### `nanocomp/raw_reads/`

| File | Description |
|--------------------------------|----------------------------------------|
| `AllNanoComp_lengths_violin.html` | Violin plot report for raw-read length distribution. |
| `AllNanoComp_log_length_violin.html` | Violin plot report for raw-read log-length distribution. |
| `AllNanoComp_N50.html` | N50 summary plot/report for raw reads. |
| `AllNanoComp_number_of_reads.html` | Plot/report of raw read counts per sample. |
| `AllNanoComp_OverlayHistogram.html` | Overlay histogram report for raw-read lengths. |
| `AllNanoComp_OverlayHistogram_Normalized.html` | Normalized overlay histogram report for raw-read lengths. |
| `AllNanoComp_OverlayLogHistogram.html` | Overlay log-histogram report for raw-read lengths. |
| `AllNanoComp_OverlayLogHistogram_Normalized.html` | Normalized overlay log-histogram report for raw-read lengths. |
| `AllNanoComp_quals_violin.html` | Violin plot report for raw-read quality score distributions. |
| `AllNanoComp-report.html` | Main NanoComp HTML summary report for raw reads. |
| `AllNanoComp_total_throughput.html` | Throughput plot/report for raw reads. |
| `AllNanoStats.txt` | Text summary statistics generated by NanoComp for raw reads. |

### `nanocomp/filt_reads/`

| File | Description |
|---------------------------------|---------------------------------------|
| `AllNanoComp_lengths_violin.html` | Violin plot report for read-length distribution. |
| `AllNanoComp_log_length_violin.html` | Violin plot report for log-transformed read lengths. |
| `AllNanoComp_N50.html` | N50 summary plot/report for read lengths. |
| `AllNanoComp_number_of_reads.html` | Plot/report of read counts per sample. |
| `AllNanoComp_OverlayHistogram.html` | Overlay histogram report for read lengths. |
| `AllNanoComp_OverlayHistogram_Normalized.html` | Normalized overlay histogram report for read lengths. |
| `AllNanoComp_OverlayLogHistogram.html` | Overlay log-histogram report for read lengths. |
| `AllNanoComp_OverlayLogHistogram_Normalized.html` | Normalized overlay log-histogram report for read lengths. |
| `AllNanoComp_quals_violin.html` | Violin plot report for base quality score distributions. |
| `AllNanoComp-report.html` | Main NanoComp HTML summary report for filtered reads. |
| `AllNanoComp_total_throughput.html` | Throughput plot/report for cumulative bases per sample. |
| `AllNanoStats.txt` | Text summary statistics generated by NanoComp. |

### `nanocomp/mapping/`

| File | Description |
|-------------------------------|-----------------------------------------|
| `AllNanoComp_lengths_violin.html` | Violin plot report for mapped-read length distribution. |
| `AllNanoComp_log_length_violin.html` | Violin plot report for mapped-read log-length distribution. |
| `AllNanoComp_N50.html` | N50 summary plot/report for mapped reads. |
| `AllNanoComp_number_of_reads.html` | Plot/report of mapped read counts per sample. |
| `AllNanoComp_OverlayHistogram.html` | Overlay histogram report for mapped-read lengths. |
| `AllNanoComp_OverlayHistogram_Identity.html` | Overlay histogram report for mapping identity distribution. |
| `AllNanoComp_OverlayHistogram_Normalized.html` | Normalized overlay histogram report for mapped-read lengths. |
| `AllNanoComp_OverlayHistogram_PhredScore.html` | Overlay histogram report for mapped-read Phred score distribution. |
| `AllNanoComp_OverlayLogHistogram.html` | Overlay log-histogram report for mapped-read lengths. |
| `AllNanoComp_OverlayLogHistogram_Normalized.html` | Normalized overlay log-histogram report for mapped-read lengths. |
| `AllNanoComp_percentIdentity_violin.html` | Violin plot report for mapping percent identity distribution. |
| `AllNanoComp_quals_violin.html` | Violin plot report for mapped-read quality scores. |
| `AllNanoComp-report.html` | Main NanoComp HTML summary report for mapped reads. |
| `AllNanoComp_total_throughput.html` | Throughput plot/report for mapped reads. |
| `AllNanoStats.txt` | Text summary statistics generated by NanoComp for mapped reads. |

### Mapped number of reads

![NanoComp_Num_Reads](../figures/NanoComp-Num-reads.png)

### Mapped read length

![NanoComp_Read_Length](../figures/NanoComp-Read-Length-Log.png)

### Mapped reads' identity to reference

![NanoComp_Ref_Identity](../figures/NanoComp-Ref-Identity.png)

## `minimap2/` (genome mapping)

| File        | Description                                                |
|-------------|------------------------------------------------------------|
| `*.bam`     | Sorted BAM alignment files generated by minimap2/samtools. |
| `*.bam.bai` | BAM index files for the corresponding `*.bam` alignments.  |

## `bambu/` (transcriptome assembly and quantification)

### Raw count-matrix & annotation standard outputs

> *Bambu standard outputs are constructed as an extension of the reference annotation that incorporates the novel RNA isoform candidates. Therefore, the raw files will include reference transcripts/genes that were not assembled. The count-matrix will will contain zero-counts in all samples for these specific extended reference features.*

| File | Description |
|--------------------------------|----------------------------------------|
| `BambuOutput_counts_gene.txt` | Gene-level raw count matrix. |
| `BambuOutput_counts_transcript.txt` | Transcript-level raw count matrix. |
| `BambuOutput_CPM_transcript.txt` | Transcript-level CPM-normalized expression matrix. |
| `BambuOutput_fullLengthCounts_transcript.txt` | Full-length transcript count matrix. |
| `BambuOutput_uniqueCounts_transcript.txt` | Uniquely mapped transcript count matrix. |
| `BambuOutput_extended_annotations.gtf` | Extended reference transcript annotations. |
| `bambu_novel_transcripts.gtf` | Novel isoform candidates identified by Bambu. |

### Standard automatic plots

| File | Description |
|--------------------------------|----------------------------------------|
| `heatmap_gene.png` | Heatmap visualization of gene-level expression patterns. |
| `heatmap_transcript.png` | Heatmap visualization of transcript-level expression patterns. |
| `pca_grouped.png` | PCA plot for samples organized by experimental design groups. |
| `pca.png` | PCA plot for individual samples. |

#### Spearman correlation coefficients

![SpR_Gene_level](../figures/bambu-gene-level-spR.png)

#### PCA analysis 

![PCA](../figures/bambu-pca.png)

### Rdata

> [!NOTE]
> You can access Bambu Novel Discovery Rate (NDR) metric by accessing the Summarized Experiment R object `se_multiSample.rds`. It is also output to the module's execution log `.command.log`, stored at the `work/` directory.

| File | Description |
|--------------------------------|----------------------------------------|
| `seGene_multiSample.rds` | Serialized R object with gene-level summarized experiment data. |
| `se_multiSample.rds` | Serialized R object with transcript-level summarized experiment data. |

### Filtered count-matrices

> *Transcripts/genes from the extended reference annotations that have 0-counts in all samples are removed, keeping only annotated and novel features that were actually assembled and quantified by Bambu from the input dataset.*

| File | Description |
|--------------------------------|----------------------------------------|
| `BambuOutput_counts_transcript_filtered.txt` | Filtered transcript-level raw count matrix. |
| `BambuOutput_counts_gene_filtered.txt` | Filtered gene-level raw count matrix. |
| `BambuOutput_CPM_transcript_filtered.txt` | Filtered transcript-level CPM-normalized expression matrix. |
| `BambuOutput_fullLengthCounts_transcript_filtered.txt` | Filtered full-length transcript count matrix. |
| `BambuOutput_fullLengthCounts_transcript_validated.txt` | Full-length transcript count matrix for validated features. |
| `BambuOutput_uniqueCounts_transcript_filtered.txt` | Filtered uniquely mapped transcript count matrix. |

### Validated transcripts & annotations

> *Filtered count-matrices are submitted to validation of novel transcript isoforms, further keeping only the reference features and selected novel lncRNA and protein-coding RNA isoform candidates that were assembled* by Bambu from the input dataset.

| File | Description |
|--------------------------------|----------------------------------------|
| `BambuOutput_counts_gene_validated.txt` | Gene-level raw count matrix restricted to validated features. |
| `BambuOutput_counts_transcript_validated.txt` | Transcript-level raw count matrix restricted to validated features. |
| `BambuOutput_CPM_transcript_validated.txt` | Transcript-level CPM matrix for validated features. |
| `BambuOutput_uniqueCounts_transcript_validated.txt` | Uniquely mapped transcript count matrix for validated features. |
| `BambuOutput_fullLength_validated.gtf` | GTF of validated full-length transcript isoforms. |
| `BambuOutput_annotations_validated.gtf` | Final validated annotation. |
| `BambuOutput_uniquelyMapped_validated.gtf` | Validated transcript isoforms supported by uniquely mapped reads. |

### Annotated GTF attributes

Bambu writes only `gene_id`, `transcript_id` and `exon_number` into its GTF output. pulposeq adds the attributes below to every `GTF` it generates, so each file is self-describing and can be filtered without cross-referencing the metadata tables.

| Attribute | Description |
|--------------------------------|----------------------------------------|
| `transcript_status` | `known` if the transcript is present in the reference annotation, `novel` if it was assembled by Bambu. |
| `gene_biotype` | For known transcripts, the biotype from the reference annotation. For novel transcripts arising from a known gene, that gene's reference biotype; for novel transcripts at previously unannotated loci, `novel`. |
| `transcript_biotype` | For known transcripts, the biotype from the reference annotation. For novel transcripts, one of `novel_lncRNA`, `novel_protein_coding` or `novel_non_coding` — see [How novel models are routed](#how-novel-models-are-routed). |
| `gene_name` | Gene symbol, where the reference annotation or gffcompare provides one. |
| `transcript_name` | Transcript name from the reference annotation (known transcripts only). |
| `class_code` | gffcompare class code relative to the reference (novel transcripts only). |
| `classification` | Human-readable reading of `class_code`, following gffcompare's own definitions, with `i` qualified `(sense)` or `(antisense)` (novel transcripts only). See [Which class codes become candidates](#which-class-codes-become-candidates). |
| `ref_gene_id` | Reference gene the novel transcript was matched against, where gffcompare found one (novel transcripts only). |

> [!NOTE]
> Novel transcripts frequently arise from genes that are already annotated. In that case `gene_biotype` reports the reference gene's real biotype while `transcript_biotype` records the novel isoform's predicted class, so the two can legitimately differ, e.g. a `novel_lncRNA` transcript within a `protein_coding` gene.

## `gffcompare/` (novel transcript comparison)

> _Check the GffCompare's official [documentation](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml)._

| File | Description |
|-------------------------|------------------------------------------------|
| `compared.annotated.gtf` | GTF output with reference-aware annotation status from gffcompare. |
| `compared.bambu_novel_transcripts.gtf.refmap` | Mapping of reference transcripts to novel isoform candidates assembled. |
| `compared.bambu_novel_transcripts.gtf.tmap` | Mapping/classification table for novel isoform candidatesagainst the reference. |
| `compared.loci` | Loci-level grouping of reference and novel isoform candidates' structures. |
| `compared.stats` | Summary statistics from novel isoform candidates comparison against the reference annotation. |
| `compared.tracking` | Tracking table. |

### GffCompare classification codes

Class_codes description figure below retrieved from the official [documentation](https://ccb.jhu.edu/software/stringtie/gffcompare_codes.png).

![descriptions](../figures/gffcompare_codes.png)

## `gffread/` (sequence extraction)

| File | Description |
|------------------------|------------------------------------------------|
| `bambu_novel_transcripts.fa` | FASTA sequences extracted for bambu novel isoform candidates. |

## `coding_potential/` (coding potential prediction for novel transcripts)

Contents depend on which predictor `coding_potential_pred` selected. Whichever ran,
the call it produced is recorded per transcript in the novel metadata as
`prediction`, with `coding_prob` holding P(coding) and `coding_predictor` naming the
tool — normalised across the two, since they report different quantities natively.

With `cpc2` (default):

| File | Description |
|------------------|------------------------------------------------------|
| `*.txt` | CPC2 table: one row per novel isoform candidate, with transcript and peptide length, Fickett score, isoelectric point, ORF integrity, coding probability and label. |

With `rnamining`:

| File | Description |
|------------------|------------------------------------------------------|
| `codings.txt` | Novel isoform candidates predicted as protein-coding, as FASTA. |
| `noncodings.txt` | Novel isoform candidates predicted as non-coding, as FASTA. |
| `predictions.txt` | Full RNAmining prediction output for all evaluated candidates. |

> RNAmining is currently under review — see [Coding potential](usage.md#coding-potential)
> in the usage documentation before selecting it.

## `annotations_metadata/` (metadata handling)

| File | Description |
|-------------------------|-----------------------------------------------|
| `annotated_lncRNAs_exonlength.csv` | Exon length summary for annotated lncRNA transcripts. |
| `annotated_lncRNAs_metadata.csv` | Metadata table for annotated lncRNA transcripts. |
| `annotated_protein-coding_exonlength.csv` | Exon length summary for annotated protein-coding transcripts. |
| `annotated_protein-coding_metadata.csv` | Metadata table for annotated protein-coding transcripts. |
| `annotated_transcriptome_metadata.csv` | Metadata summary for the entire annotated transcriptome. |
| `bambu_annotated_lncRNAs.gtf` | GTF annotation file containing annotated lncRNA isoforms from bambu outputs. |
| `bambu_annotated_mRNAs.gtf` | GTF annotation file containing annotated mRNA/protein-coding isoforms from bambu outputs. |
| `bambu_annotated_transcriptome_gene_counts.csv` | Gene-level count matrix for annotated transcriptome features. |
| `bambu_annotated_transcriptome.gtf` | GTF file for validated annotated features assembled/quantified by bambu. |
| `bambu_annotated_transcriptome_tx_counts.csv` | Transcript-level count matrix for annotated transcriptome features. |
| `bambu_novel_pc_lnc_RNA_gene_counts.csv` | Gene-level counts for novel isoform candidates classified as protein-coding or lncRNA. |
| `bambu_novel_pc_lnc_RNA_tx_counts.csv` | Transcript-level counts for novel isoform candidates classified as protein-coding or lncRNA. |
| `novel_lncRNA_exon_lengths.csv` | Exon length summary for novel lncRNA isoform candidates. |
| `novel_lncRNAs.gtf` | GTF file containing novel lncRNA isoform candidates. |
| `novel_lncRNAs_metadata.csv` | Metadata table for novel lncRNA isoform candidates. |
| `novel_pc_lnc_RNAs_metadata.csv` | Combined metadata table for novel protein-coding and lncRNA transcript isoform candidates. |
| `novel_protein-coding_exon_lengths.csv` | Exon length summary for novel protein-coding isoform candidates. |
| `novel_protein-coding.gtf` | GTF file containing novel protein-coding RNA isoform candidates. |
| `novel_protein-coding_metadata.csv` | Metadata table for novel protein-coding RNA isoform candidates. |
| `novel_transcripts_metadata.csv` | Metadata table for all novel RNA isoform candidates. |
| `novel_transcripts.gtf` | GTF file containing only the validated set of novel isoform candidates. |
| `novel_non_coding_metadata.csv` | Metadata table for models predicted non-coding that share splice structure with a reference gene which is not an lncRNA. |
| `novel_non_coding.gtf` | GTF file for the same set, with `transcript_biotype "novel_non_coding"`. |
| `novel_context_flags.csv` | Structural evidence per novel transcript: the reference gene, its biotype and strand, `same_strand_as_host`, and the read counts behind the boundary and junction tests. |
| `novel_context_summary.csv` | Aggregate counts from the same tests, as rendered in the report. |

### Which class codes become candidates

Nine gffcompare class codes are eligible. Each is given its wording in the
`classification` column, following gffcompare's own definitions so the letter never
has to be looked up and the phrasing does not drift from what the tool asserts:

| Code | `classification` |
|---|---|
| `u` | unknown or intergenic |
| `i` | fully contained within ref intron |
| `x` | exonic overlap on the opposite strand |
| `j` | multi-exonic matching ref splice junction(s) |
| `k` | contains reference transcript |
| `o` | exonic overlap on the same strand |
| `y` | contains reference within its introns |
| `m` | retained intron (all matched or retained) |
| `n` | retained intron (not all matched or retained) |

`i` is additionally qualified by orientation — `fully contained within ref intron
(sense)` or `(antisense)`. That is the only class code whose strand relative to the
reference varies: `x` is antisense by definition and the rest are same-strand
matches. The qualifier is pulposeq's addition; gffcompare's definition of `i` says
nothing about strand.

> [!NOTE]
> Class codes are reported exactly as gffcompare assigns them — pulposeq never
> recomputes or overrides one. A code is assigned relative to a single matched
> reference transcript, and cases have been observed where it did not match the
> geometry: a model overlapping none of its reference gene's annotated exons, and
> lying wholly within one of its introns, was assigned `x` rather than `i`. Where
> that happens the sense/antisense qualifier under-counts intron-contained models.
> Routing is unaffected, since `i` and `x` are both decided by the coding
> prediction alone.

`=` and `c` never appear. Only novel transcripts reach gffcompare, so a model
matching or contained by a reference was already resolved as annotated. `s`, `e`,
`p` and `r` are not admitted.

### How novel models are routed

Each candidate is assigned one of three `transcript_biotype` values, from its class
code, its coding-potential prediction, and the biotype of the reference gene it was
matched against.

**`u`, `i` and `x` share no splice structure with a reference.** `u` has no
reference at all, `i` lies wholly inside an intron, and `x` overlaps on the opposite
strand. The coding prediction decides alone:

| Prediction | `transcript_biotype` |
|---|---|
| non-coding | `novel_lncRNA` |
| coding | `novel_protein_coding` |

A sense-intronic or antisense lncRNA inside a protein-coding locus is still an
lncRNA, and both Ensembl and GENCODE annotate them that way.

**`j`, `k`, `o`, `y`, `m` and `n` do share structure with a reference transcript.**
A non-coding call carries less weight here: an unproductive isoform of a
protein-coding gene — a retained intron, an NMD target, a truncated model — is
predicted non-coding too, and calling one a novel lncRNA would assert a new
non-coding gene at a locus that already has a coding one.

| Prediction | Reference gene | `transcript_biotype` |
|---|---|---|
| coding | any | `novel_protein_coding` |
| non-coding | lncRNA | `novel_lncRNA` |
| non-coding | anything else | `novel_non_coding` |

`novel_non_coding` states that a model has no coding potential without asserting
what it is, which is as much as the evidence supports for a non-coding isoform of a
coding gene. The prediction and its probability are kept as columns on every record
in all three categories.

### Reference identity on novel records

Every novel record carries what it was compared against: `ref_gene_id`,
`ref_gene_name`, `ref_gene_biotype`, `ref_id`, `ref_transcript_name` and
`ref_transcript_biotype`, on the metadata CSVs and as GTF attributes.

The **transcript**-level biotype is not redundant with the gene-level one, and the
difference is where the interesting cases live. A `y`-class novel transcript
containing `TARDBP-221` has `ref_gene_biotype` `protein_coding` and
`ref_transcript_biotype` `nonsense_mediated_decay` — only the second says what it
actually contains.

These are named `ref_` rather than `host_` because a host is something a transcript
sits *inside*, which holds for `i`, `x`, `m` and `n` but is backwards for `k` and
`y`, where the novel transcript contains the reference. Columns describing what
*reads* do — `same_strand_as_host`, `reads_with_host_junction`,
`reads_spliced_into_host_exon` — keep the word, since the step producing them runs
only on intronic candidates, where it is accurate.

### Reading `novel_context_flags.csv`

| Column | Meaning |
|--------------------------|----------------------------------------------|
| `same_strand_as_host` | `TRUE` where the transcript shares its host gene's strand. Empty for intergenic transcripts, which have no host — this is deliberately distinct from `FALSE`. |
| `reads_crossing_boundary` | Supporting reads extending past either end of the transcript. A discrete transcript has few; a fragment of a longer molecule has many. |
| `reads_with_host_junction` | Supporting reads carrying a splice junction belonging to the host gene. Evidence that the read came from the host, not from the candidate. |
| `reads_spliced_into_host_exon` | Supporting reads spliced from the candidate into one of the host's exons. This is a **finding**, not an artifact: the candidate is most likely an unannotated exon of the host. |
| `frac_*` | The corresponding count over `reads_total`. `NA` where no reads were recovered. |

No transcript is removed on the basis of these columns. They are reported so that
a count is never quoted without the evidence behind it.


## `genomic_context/` (coverage and transcript models at selected loci)

| File | Description |
|--------------------------|----------------------------------------------|
| `genomic_context_<gene>.png` | Gene structure, every isoform at the locus with novel ones highlighted, and one coverage track per sample on a shared scale. Selected from genes carrying both known and novel isoforms. |
| `genomic_context_candidates.csv` | The loci drawn, with their windows and isoform counts. |
| `genomic_context_regions.gtf` | The plotted regions as a small GTF. Built to feed `makeTxDbFromGFF`, but useful on its own — load it in IGV beside the published bigWigs. |
| `intronic_context_<transcript>.png` | The flagged set: intronic candidates on their host's strand whose supporting reads run past the boundaries or carry host junctions. Windowed on the **whole host intron**, because a window drawn tightly around a truncated fragment looks discrete whatever it is. |
| `intronic_context_candidates.csv` | The flagged candidates drawn, with the read counts behind their selection. |

The two sets exist to be read against each other. A candidate that really is host
pre-mRNA shows coverage running flat across the whole intron and continuous with
the flanking exons; a genuine independent transcript shows a discrete block with
quiet intron either side. If a flagged candidate looks discrete, the flag is what
needs revisiting — not the transcript.

## `bam_coverage/` (per-sample coverage tracks)

| File | Description |
|--------------------------|----------------------------------------------|
| `<sample>.bw` | Genome-wide coverage in bigWig format. Two to three orders of magnitude smaller than the alignments they come from, so this is the form of the data that travels off the cluster. Load them in IGV to explore any region beyond the windows the pipeline chose to draw. |

## `pipeline_info/` (workflow run metadata)

| File | Description |
|--------------------------|----------------------------------------------|
| `execution_report_<timestamp>.html` | Nextflow execution report with runtime, resources, and process-level summaries. |
| `execution_timeline_<timestamp>.html` | Nextflow execution timeline visualization. |
| `execution_trace_<timestamp>.txt` | Nextflow trace table with task-level runtime and resource usage. |
| `nf_core_pipeline_software_mqc_versions.yml` | Consolidated software versions used in the run and reported by MultiQC. |
| `params_<timestamp>.json` | Parameters snapshot used for the current pipeline run. |
| `pipeline_dag_<timestamp>.html` | DAG visualization of process dependencies for the current pipeline run. |

### Resource usage

![nf-rep-res](../figures/nextflow-report-resources.gif)

## Nextflow reports

Nextflow provides excellent functionality for generating various [reports](https://www.nextflow.io/docs/latest/reports.html) relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

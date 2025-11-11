process NOVEL_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    tuple val(meta), path(bambu_gtf), path(compared_gtf), path(tmap_file), path(rnamining_predictions), path(tx_counts), path(gene_counts)

    output:
    tuple val(meta), path("novel_transcripts_metadata.csv")          , emit: novel_transcripts_metadata
    tuple val(meta), path("novel_lncRNAs_metadata.csv")              , emit: novel_lncrnas_metadata
    tuple val(meta), path("novel_protein-coding_metadata.csv")       , emit: novel_mrnas_metadata
    tuple val(meta), path("novel_pc_lnc_RNAs_metadata.csv")          , emit: novel_combined_metadata
    tuple val(meta), path("novel_lncRNAs.gtf")                       , emit: novel_lncrnas_gtf
    tuple val(meta), path("novel_protein-coding.gtf")                , emit: novel_mrnas_gtf
    tuple val(meta), path("novel_lncRNA_exon_lengths.csv")           , emit: novel_lncrna_exon_lengths
    tuple val(meta), path("novel_protein-coding_exon_lengths.csv")   , emit: novel_mrna_exon_lengths
    tuple val(meta), path("bambu_novel_pc_lnc_RNA_tx_counts.csv")    , emit: novel_tx_counts
    tuple val(meta), path("bambu_novel_pc_lnc_RNA_gene_counts.csv")  , emit: novel_gene_counts
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/novel_transcripts.R \\
        --bambu_gtf ${bambu_gtf} \\
        --compared_gtf ${compared_gtf} \\
        --tmap_file ${tmap_file} \\
        --rnamining_predictions ${rnamining_predictions} \\
        --tx_counts ${tx_counts} \\
        --gene_counts ${gene_counts} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    touch novel_transcripts_metadata.csv
    touch novel_lncRNAs_metadata.csv
    touch novel_protein-coding_metadata.csv
    touch novel_pc_lnc_RNAs_metadata.csv
    touch novel_lncRNAs.gtf
    touch novel_protein-coding.gtf
    touch novel_lncRNA_exon_lengths.csv
    touch novel_protein-coding_exon_lengths.csv
    touch bambu_novel_pc_lnc_RNA_tx_counts.csv
    touch bambu_novel_pc_lnc_RNA_gene_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
    END_VERSIONS
    """
}

process BAMBU {
    tag "Running Bambu"
    label 'process_bambu'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/longnoncoder:test3':
        'docker.io/itsiaguara/longnoncoder:test3' }"

    input:
    val bam_list
    path reference
    path annotation
    val sample_info


    output:
    path "heatmap_gene.png"                              , emit: h_gene
    path "heatmap_transcript.png"                        , emit: h_transcript
    path "pca_grouped.png"                               , emit: pca_grouped
    path "pca.png"                                       , emit: pca
    path "BambuOutput_counts_transcript.txt"             , emit: tx_counts
    path "BambuOutput_counts_gene.txt"                   , emit: gene_counts
    path "BambuOutput_CPM_transcript.txt"                , emit: CPM
    path "BambuOutput_fullLengthCounts_transcript.txt"   , emit: full_length
    path "BambuOutput_uniqueCounts_transcript.txt"       , emit: unique_counts
    path "bambu_novel_transcripts.gtf"                   , emit: gtf_new_transcripts
    path "BambuOutput_extended_annotations.gtf"          , emit: gtf_all_transcripts
    path "se_multiSample.rds"                            , emit: rds_transcript
    path "seGene_multiSample.rds"                        , emit: rds_gene
    path "*.rds"                                         , emit: rds
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    """
    bambu.R \\
        -g $reference \\
        -a $annotation \\
        -b $bam_list \\
        -n $task.cpus \\
        -s $sample_info \\
        -o . \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-biocparallel: \$(Rscript -e "cat(as.character(packageVersion('BiocParallel')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))") 
    END_VERSIONS
    """

    stub:
    // REMOVED the unused 'def args' line from here
    """
    touch heatmap_gene.png
    touch heatmap_transcript.png
    touch pca_grouped.png
    touch pca.png
    touch BambuOutput_counts_transcript.txt
    touch BambuOutput_counts_gene.txt
    touch BambuOutput_CPM_transcript.txt
    touch BambuOutput_fullLengthCounts_transcript.txt
    touch BambuOutput_uniqueCounts_transcript.txt
    touch bambu_novel_transcripts.gtf
    touch BambuOutput_extended_annotations.gtf
    touch se_multiSample.rds
    touch seGene_multiSample.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-biocparallel: \$(Rscript -e "cat(as.character(packageVersion('BiocParallel')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))") 
    END_VERSIONS
    """
}

process POST_REFINEMENT {
    tag "Regenerating_Validated_Plots"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/longnoncoder:test3':
        'docker.io/itsiaguara/longnoncoder:test3' }"

    input:
    path se_rds
    path se_gene_rds
    path counts_transcript_validated
    path counts_gene_validated
    path r_script

    output:
    path "postrefinement_pca.png"                , emit: pca
    path "postrefinement_pca_grouped.png"        , emit: pca_grouped
    path "postrefinement_heatmap_gene.png"       , emit: h_gene
    path "postrefinement_heatmap_transcript.png" , emit: h_transcript
    path "se_multiSample_validated.rds"          , emit: rds_transcript
    path "seGene_multiSample_validated.rds"      , emit: rds_gene
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript $r_script \\
        --se_rds ${se_rds} \\
        --se_gene_rds ${se_gene_rds} \\
        --counts_transcript ${counts_transcript_validated} \\
        --counts_gene ${counts_gene_validated} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """

    stub:
    """
    touch postrefinement_pca.png
    touch postrefinement_pca_grouped.png
    touch postrefinement_heatmap_gene.png
    touch postrefinement_heatmap_transcript.png
    touch se_multiSample_validated.rds
    touch seGene_multiSample_validated.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        bioconductor-bambu: \$(Rscript -e "cat(as.character(packageVersion('bambu')))")
        bioconductor-summarizedexperiment: \$(Rscript -e "cat(as.character(packageVersion('SummarizedExperiment')))")
    END_VERSIONS
    """
}

process KNOWN_TRANSCRIPTS {
    tag "Processing_Known_Transcripts"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/longnoncoder:test':
        'docker.io/itsiaguara/longnoncoder:test' }"

    input:
    path transcript_counts
    path gene_counts
    path gtf_file

    output:
    path "annotated_transcriptome_metadata.csv"          , emit: transcriptome_metadata
    path "annotated_lncRNAs_metadata.csv"                , emit: lncrna_metadata
    path "annotated_lncRNAs_exonlength.csv"              , emit: lncrna_exonlength
    path "annotated_protein-coding_metadata.csv"         , emit: protein_coding_metadata
    path "annotated_protein-coding_exonlength.csv"       , emit: protein_coding_exonlength
    path "bambu_annotated_transcriptome.gtf"             , emit: annotated_transcriptome_gtf
    path "bambu_annotated_transcriptome_tx_counts.csv"   , emit: annotated_tx_counts
    path "bambu_annotated_lncRNAs.gtf"                   , emit: annotated_lncrna_gtf
    path "bambu_annotated_mRNAs.gtf"                     , emit: annotated_mrna_gtf
    path "bambu_annotated_transcriptome_gene_counts.csv" , emit: annotated_gene_counts
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ? task.ext.args : "--ensembl_organism_dataset '${params.ensembl_organism_dataset}' --ensembl_version ${params.ensembl_version}"
    """
    # Copy the R script to the working directory
    cp ${projectDir}/bin/known_transcripts.R ./

    # Run the R script
    Rscript known_transcripts.R \\
        --transcript_counts ${transcript_counts} \\
        --gene_counts ${gene_counts} \\
        --gtf_file ${gtf_file} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-httr: \$(Rscript -e "cat(as.character(packageVersion('httr')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-biomart: \$(Rscript -e "cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """

    stub:
    def args     = task.ext.args ? task.ext.args : "--ensembl_organism_dataset '${params.ensembl_organism_dataset}' --ensembl_version ${params.ensembl_version}"
    """
    touch annotated_transcriptome_metadata.csv
    touch annotated_lncRNAs_metadata.csv
    touch annotated_lncRNAs_exonlength.csv
    touch annotated_protein-coding_metadata.csv
    touch annotated_protein-coding_exonlength.csv
    touch bambu_annotated_transcriptome.gtf
    touch bambu_annotated_transcriptome_tx_counts.csv
    touch bambu_annotated_lncRNAs.gtf
    touch bambu_annotated_mRNAs.gtf
    touch bambu_annotated_transcriptome_gene_counts.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        r-httr2: \$(Rscript -e "cat(as.character(packageVersion('httr2')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-biomart: \$(Rscript -e "cat(as.character(packageVersion('biomaRt')))")
    END_VERSIONS
    """
}
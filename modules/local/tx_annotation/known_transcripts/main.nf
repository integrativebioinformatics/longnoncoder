process KNOWN_TRANSCRIPTS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-biomart=2.58.0 bioconda::bioconductor-genomicranges=1.54.0 bioconda::bioconductor-rtracklayer=1.62.0 conda-forge::r-dplyr=1.1.4 conda-forge::r-readr=2.1.4"
    
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    //    'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e2/e29cd4dfa406c2f0aaab7cad1662a69f2fe871bc67a85b468a6f45a94d15505f/data' : 
    //    'community.wave.seqera.io/library/bioconductor-biomart_bioconductor-genomicranges_bioconductor-rtracklayer_r-dplyr_r-readr:latest' }"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lfreitasl/bambu:3.8.0':
        'docker.io/lfreitasl/bambu:3.8.0' }"

    input:
    tuple val(meta), path(transcript_counts), path(gene_counts), path(gtf_file)

    output:
    tuple val(meta), path("annotated_transcriptome_metadata.csv")           , emit: transcriptome_metadata
    tuple val(meta), path("annotated_lncRNAs_metadata.csv")                , emit: lncrna_metadata
    tuple val(meta), path("annotated_lncRNAs_exonlength.csv")              , emit: lncrna_exonlength
    tuple val(meta), path("annotated_protein-coding_metadata.csv")         , emit: protein_coding_metadata
    tuple val(meta), path("annotated_protein-coding_exonlength.csv")       , emit: protein_coding_exonlength
    tuple val(meta), path("bambu_annotated_transcriptome.gtf")             , emit: annotated_transcriptome_gtf
    tuple val(meta), path("bambu_annotated_transcriptome_tx_counts.csv")   , emit: annotated_tx_counts
    tuple val(meta), path("bambu_annotated_lncRNAs.gtf")                   , emit: annotated_lncrna_gtf
    tuple val(meta), path("bambu_annotated_mRNAs.gtf")                     , emit: annotated_mrna_gtf
    tuple val(meta), path("bambu_annotated_transcriptome_gene_counts.csv") , emit: annotated_gene_counts
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Copy the R script to the working directory
    cp ${projectDir}/bin/known_transcripts.R ./

    # Run the R script
    Rscript known_transcripts.R \\
        --transcript_counts ${transcript_counts} \\
        --gene_counts ${gene_counts} \\
        --gtf_file ${gtf_file} \\
        --prefix ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-biomart: \$(Rscript -e "cat(as.character(packageVersion('biomaRt')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-biomart: \$(Rscript -e "cat(as.character(packageVersion('biomaRt')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        r-dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
        r-readr: \$(Rscript -e "cat(as.character(packageVersion('readr')))")
    END_VERSIONS
    """
}

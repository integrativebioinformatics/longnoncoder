process SUBSET_BAMBU_COUNTS {
    tag "Subsetting_Counts"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path counts_gene
    path counts_transcript
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_counts_gene_subset.txt"                 , emit: counts_gene_subset
    path "BambuOutput_counts_transcript_subset.txt"           , emit: counts_transcript_subset
    path "BambuOutput_CPM_transcript_subset.txt"              , emit: cpm_transcript_subset
    path "BambuOutput_fullLengthCounts_transcript_subset.txt" , emit: full_length_counts_transcript_subset
    path "BambuOutput_uniqueCounts_transcript_subset.txt"     , emit: unique_counts_transcript_subset
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    subset_bambu_counts.sh \\
        --counts_gene ${counts_gene} \\
        --counts_transcript ${counts_transcript} \\
        --cpm_transcript ${cpm_transcript} \\
        --full_length ${full_length_counts_transcript} \\
        --unique ${unique_counts_transcript} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_gene_subset.txt
    touch BambuOutput_counts_transcript_subset.txt
    touch BambuOutput_CPM_transcript_subset.txt
    touch BambuOutput_fullLengthCounts_transcript_subset.txt
    touch BambuOutput_uniqueCounts_transcript_subset.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mawk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}
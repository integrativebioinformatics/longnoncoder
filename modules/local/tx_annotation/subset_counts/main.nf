process SUBSET_BAMBU_COUNTS {
    tag "Subsetting_Counts"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/longnoncoder:test':
        'docker.io/itsiaguara/longnoncoder:test' }"

    input:
    path counts_gene
    path counts_transcript
    path cpm_transcript
    path full_length_counts_transcript
    path unique_counts_transcript

    output:
    path "BambuOutput_counts_gene_filtered.txt"                 , emit: counts_gene_filtered
    path "BambuOutput_counts_transcript_filtered.txt"           , emit: counts_transcript_filtered
    path "BambuOutput_CPM_transcript_filtered.txt"              , emit: cpm_transcript_filtered
    path "BambuOutput_fullLengthCounts_transcript_filtered.txt" , emit: full_length_counts_transcript_filtered
    path "BambuOutput_uniqueCounts_transcript_filtered.txt"     , emit: unique_counts_transcript_filtered
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
        --unique ${unique_counts_transcript}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch BambuOutput_counts_gene_filtered.txt
    touch BambuOutput_counts_transcript_filtered.txt
    touch BambuOutput_CPM_transcript_filtered.txt
    touch BambuOutput_fullLengthCounts_transcript_filtered.txt
    touch BambuOutput_uniqueCounts_transcript_filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n1 | sed 's/mawk //; s/,.*//')
        bash: \$(bash --version | head -n1 | sed 's/GNU bash, version //; s/ .*//')
    END_VERSIONS
    """
}
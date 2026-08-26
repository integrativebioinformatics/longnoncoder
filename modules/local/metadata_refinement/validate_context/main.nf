process VALIDATE_NOVEL_CONTEXT {
    tag "Validate Novel Transcript Context"
    label 'process_validate_context'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path metadata
    path annotation
    path bams  , stageAs: 'alignments/*'
    path bais  , stageAs: 'alignments/*'
    path r_script

    output:
    path "novel_context_flags.csv"   , emit: flags
    path "novel_context_summary.csv" , emit: summary
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # BAMs and their indexes are staged into one directory so that Rsamtools
    # finds each .bai beside its .bam. Nextflow stages them as symlinks, so this
    # costs no disk regardless of how large the alignments are.
    bam_list=\$(ls alignments/*.bam | paste -sd, -)

    Rscript $r_script \\
        --metadata ${metadata} \\
        --annotation ${annotation} \\
        --bams "\$bam_list" \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-genomicalignments: \$(Rscript -e "cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))")
    END_VERSIONS
    """

    stub:
    """
    printf 'qry_id,class_code,strand,host_gene_id,host_gene_biotype,host_gene_strand,same_strand_as_host,reads_total,reads_same_strand,reads_crossing_boundary,reads_with_host_junction,reads_spliced_into_host_exon,frac_crossing_boundary,frac_with_host_junction\\n' > novel_context_flags.csv
    printf 'metric,value\\ncandidates_evaluated,0\\n' > novel_context_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-genomicranges: \$(Rscript -e "cat(as.character(packageVersion('GenomicRanges')))")
        bioconductor-genomicalignments: \$(Rscript -e "cat(as.character(packageVersion('GenomicAlignments')))")
        bioconductor-rsamtools: \$(Rscript -e "cat(as.character(packageVersion('Rsamtools')))")
    END_VERSIONS
    """
}

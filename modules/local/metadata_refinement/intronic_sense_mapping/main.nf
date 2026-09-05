process INTRONIC_SENSE_MAPPING {
    tag "Intronic Sense Mapping"
    label 'process_intronic_sense_mapping'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path metadata
    path annotation
    // Per-sample transcript counts. Each candidate is measured only in the
    // samples that quantified it, so the reads of a sample where the model was
    // not called cannot be counted against it.
    path counts
    path bams  , stageAs: 'alignments/*'
    path bais  , stageAs: 'alignments/*'
    path r_script

    output:
    path "intronic_sense_flags.csv"   , emit: flags
    path "intronic_sense_summary.csv" , emit: summary
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
        --counts ${counts} \\
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
    printf 'qry_id,class_code,strand,ref_gene_id,ref_gene_biotype,ref_gene_strand,same_strand_as_host,samples_quantified,samples_total,quantified_samples,reads_total,reads_same_strand,median_overrun_5p,q90_overrun_5p,median_overrun_3p,q90_overrun_3p,reads_with_host_junction,reads_spliced_into_host_exon,frac_with_host_junction\\n' > intronic_sense_flags.csv
    printf 'metric,value\\ncandidates_evaluated,0\\n' > intronic_sense_summary.csv

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

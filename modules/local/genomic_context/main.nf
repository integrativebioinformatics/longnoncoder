process GENOMIC_CONTEXT {
    tag "Genomic_Context"
    label 'process_genomic_context'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://itsiaguara/pulposeq:test':
        'docker.io/itsiaguara/pulposeq:test' }"

    input:
    path gtf_file      // validated, enriched transcriptome GTF
    path bigwigs       // every sample's coverage track, collected
    path context_flags // per-transcript strand and read-through flags
    path annotation    // reference annotation, for host gene structure
    path r_script

    output:
    path "genomic_context_*.png"          , emit: figures, optional: true
    path "genomic_context_candidates.csv" , emit: candidates
    // The flagged set: intronic candidates on the host's strand, windowed on the
    // whole host intron so the coverage either side of the candidate is visible.
    path "intronic_context_*.png"          , emit: intronic_figures, optional: true
    path "intronic_context_candidates.csv" , emit: intronic_candidates
    // The plotted regions as a small GTF. Built to feed makeTxDbFromGFF, but useful
    // on its own: load it in IGV beside the published bigWigs to explore the same
    // loci interactively.
    path "genomic_context_regions.gtf"    , emit: regions_gtf, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def bw_list   = bigwigs.collect { bw -> bw.name }.join(',')
    def name_list = bigwigs.collect { bw -> bw.name.replaceAll(/\.bw$/, '') }.join(',')
    """
    Rscript $r_script \\
        --gtf ${gtf_file} \\
        --bigwigs ${bw_list} \\
        --names ${name_list} \\
        --context_flags ${context_flags} \\
        --annotation ${annotation} \\
        --outdir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-plotgardener: \$(Rscript -e "cat(as.character(packageVersion('plotgardener')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-txdbmaker: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('txdbmaker')), error = function(e) 'not installed'))")
    END_VERSIONS
    """

    stub:
    """
    touch genomic_context_STUBGENE.png
    touch genomic_context_regions.gtf
    printf 'gene_id,gene_name,label,chrom,start,end,win_s,win_e,n_tx,n_known,n_novel,n_novel_lnc,biotypes,figure\\n' > genomic_context_candidates.csv
    printf 'ENSG00000000000,STUBGENE,STUBGENE,chr1,1000,9000,950,9050,4,3,1,1,protein_coding;novel_lncRNA,genomic_context_STUBGENE.png\\n' >> genomic_context_candidates.csv

    touch intronic_context_STUBTX.png
    printf 'qry_id,host_gene_id,host_gene_name,host_gene_biotype,chrom,start,end,win_s,win_e,strand,class_code,num_exons,reads_total,reads_crossing_boundary,reads_with_host_junction,figure\\n' > intronic_context_candidates.csv
    printf 'BambuTxSTUB,ENSG00000000000,STUBGENE,protein_coding,chr1,3000,3800,2000,5000,+,i,1,42,31,18,intronic_context_STUBTX.png\\n' >> intronic_context_candidates.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version 2>&1 | sed 's/R version //; s/ (.*//' | head -1)
        r-optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
        bioconductor-plotgardener: \$(Rscript -e "cat(as.character(packageVersion('plotgardener')))")
        bioconductor-rtracklayer: \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))")
        bioconductor-txdbmaker: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('txdbmaker')), error = function(e) 'not installed'))")
    END_VERSIONS
    """
}

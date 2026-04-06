process GFFCOMPARE {
    tag "Comparing_Gtfs"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffcompare:0.12.6--h9f5acd7_0' :
        'biocontainers/gffcompare:0.12.6--h9f5acd7_0' }"

    input:
    path  gtfs
    path  reference_gtf

    output:
    path "*.annotated.gtf", optional: true, emit: annotated_gtf
    path "*.combined.gtf" , optional: true, emit: combined_gtf
    path "*.tmap"         , optional: true, emit: tmap
    path "*.refmap"       , optional: true, emit: refmap
    path "*.loci"         , emit: loci
    path "*.stats"         , emit: stats
    path "*.tracking"     , emit: tracking
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "compared"
    def ref_gtf = reference_gtf ? "-r ${reference_gtf}" : ''
    """
    gffcompare \\
        $args \\
        $ref_gtf \\
        $gtfs \\
        -o $prefix \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(echo \$(gffcompare --version 2>&1) | sed 's/^gffcompare v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "compared"
    """
    touch ${prefix}.annotated.gtf
    touch ${prefix}.combined.gtf
    touch ${prefix}.tmap
    touch ${prefix}.refmap
    touch ${prefix}.loci
    touch ${prefix}.stats
    touch ${prefix}.tracking

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gffcompare: \$(echo \$(gffcompare --version 2>&1) | sed 's/^gffcompare v//')
    END_VERSIONS
    """
}

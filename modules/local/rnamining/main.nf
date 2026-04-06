process RNAMINING {
    tag 'Predicting_Coding_Potential'
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rnamining:1.0.4--pyhdfd78af_0':
        'biocontainers/rnamining:1.0.4--pyhdfd78af_0' }"

    input:
    val fasta

    output:
    path 'codings.txt'        , emit: coding
    path 'noncodings.txt'     , emit: noncoding
    path 'predictions.txt'    , emit: preds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ? task.ext.args : "-organism_name ${params.organism} -prediction_type coding_prediction"
    def prefix   = task.ext.prefix ?: "Coding_Potential"

    """
    rnamining \\
            -f $fasta \\
            $args \\
            -output_folder ./

    """

    stub:
    def args     = task.ext.args ? task.ext.args : "-organism_name ${params.organism} -prediction_type coding_prediction"
    def prefix   = task.ext.prefix ?: "Coding_Potential"
    
    
    """
    touch ${prefix}.txt

    """
}

//
// Subworkflow with utility functions specific to the nf-core pipeline template
//

import org.yaml.snakeyaml.Yaml
import nextflow.extension.FilesEx

/*
========================================================================================
    SUBWORKFLOW DEFINITION
========================================================================================
*/

workflow UTILS_NFCORE_PIPELINE {

    take:
    nextflow_cli_args

    main:
    valid_config = checkConfigProvided()
    checkProfileProvided(nextflow_cli_args)

    emit:
    valid_config
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
//  Warn if a -profile or Nextflow config has not been provided to run the pipeline
//
def checkConfigProvided() {
    valid_config = true
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[$workflow.manifest.name] You are attempting to run the pipeline without any custom configuration!\n\n" +
            "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
            "   (1) Using an existing pipeline profile e.g. `-profile docker` or `-profile singularity`\n" +
            "   (2) Using an existing nf-core/configs for your Institution e.g. `-profile crick` or `-profile uppmax`\n" +
            "   (3) Using your own local custom config e.g. `-c /path/to/your/custom.config`\n\n" +
            "Please refer to the quick start section and usage docs for the pipeline.\n "
        valid_config = false
    }
    return valid_config
}

//
// Exit pipeline if --profile contains spaces
//
def checkProfileProvided(nextflow_cli_args) {
    if (workflow.profile.endsWith(',')) {
        error "The `-profile` option cannot end with a trailing comma, please remove it and re-run the pipeline!\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
    if (nextflow_cli_args[0]) {
        log.warn "The pipeline does not accept positional arguments. The positional argument `${nextflow_cli_args[0]}` has been detected.\n" +
            "HINT: A common mistake is to provide multiple values separated by spaces e.g. `-profile test, docker`.\n"
    }
}

//
// Citation string for pipeline
//
def workflowCitation() {
    return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
        "* The pipeline\n" +
        "  ${workflow.manifest.doi}\n\n" +
        "* The nf-core framework\n" +
        "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
        "* Software dependencies\n" +
        "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}

//
// Generate workflow version string
//
def getWorkflowVersion() {
    String version_string = ""
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }

    if (workflow.commitId) {
        def git_shortsha = workflow.commitId.substring(0, 7)
        version_string += "-g${git_shortsha}"
    }

    return version_string
}

//
// Get software versions for pipeline
//
def processVersionsFromYAML(yaml_file) {
    Yaml yaml = new Yaml()
    versions = yaml.load(yaml_file).collectEntries { k, v -> [ k.tokenize(':')[-1], v ] }
    return yaml.dumpAsMap(versions).trim()
}

//
// Get workflow version for pipeline
//
def workflowVersionToYAML() {
    return """
    Workflow:
        $workflow.manifest.name: ${getWorkflowVersion()}
        Nextflow: $workflow.nextflow.version
    """.stripIndent().trim()
}

//
// Get channel of software versions used in pipeline in YAML format
//
def softwareVersionsToYAML(ch_versions) {
    return ch_versions
                .unique()
                .map { processVersionsFromYAML(it) }
                .unique()
                .mix(Channel.of(workflowVersionToYAML()))
}

//
// Get workflow summary for MultiQC
//
def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    for (group in summary_params.keySet()) {
        def group_params = summary_params.get(group)  // This gets the parameters of that particular group
        if (group_params) {
            summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
            summary_section += "    <dl class=\"dl-horizontal\">\n"
            for (param in group_params.keySet()) {
                summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
            }
            summary_section += "    </dl>\n"
        }
    }

    String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
    yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
    yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
    yaml_file_text        += "plot_type: 'html'\n"
    yaml_file_text        += "data: |\n"
    yaml_file_text        += "${summary_section}"

    return yaml_file_text
}


//
// Return dashed line
//
def dashedLine(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}

//
// ANSII colours used for terminal logging
//
def logColours(monochrome_logs=true) {
    Map colorcodes = [:]

    // Reset / Meta
    colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
    colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
    colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
    colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
    colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
    colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
    colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

    // Regular Colors
    colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
    colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
    colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

    // Bold
    colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
    colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
    colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
    colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
    colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
    colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
    colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
    colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

    // Underline
    colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
    colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
    colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
    colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
    colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
    colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
    colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
    colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

    // High Intensity
    colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
    colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
    colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
    colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
    colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
    colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
    colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
    colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

    // Bold High Intensity
    colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
    colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
    colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
    colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
    colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
    colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
    colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

    return colorcodes
}

//
// Print pipeline summary on completion
//
def completionSummary(monochrome_logs=true) {
    Map colors = logColours(monochrome_logs)
    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
        } else {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.yellow} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
        }
    } else {
        log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
    }
}
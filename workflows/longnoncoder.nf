/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC } from '../modules/nf-core/multiqc/main'
include { QC_FILT } from '../subworkflows/local/qc'
include { ALIGNMENT } from '../subworkflows/local/alignment'
include { TRANSCRIPT_RECONSTRUCTION } from '../subworkflows/local/transcript_reconstruction'
include { CLASSIFICATION_POTENTIAL_CODING } from '../subworkflows/local/classification_codingpotential.nf'
include { paramsSummaryMap } from 'plugin/nf-validation'
include { paramsSummaryMultiqc } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_longnoncoder_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LONGNONCODER {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_gtf_new_transcripts = Channel.empty()

    //
    //Run QC workflow
    //
    if (!params.skip_qc) {
        QC_FILT(
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(QC_FILT.out.multiqc)
        ch_versions = ch_versions.mix(QC_FILT.out.versions)
    }
    //
    // Run alignment workflow
    //
    if (!params.skip_alignment) {
        ALIGNMENT(QC_FILT.out.filt_reads)

        if (!params.skip_alignment_qc) {
            ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.multiqc)
        }
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        TRANSCRIPT_RECONSTRUCTION(
            ALIGNMENT.out.bam,
            params.reference,
            params.annotation,
        )

        TRANSCRIPT_RECONSTRUCTION.out.gtf_new_transcripts.set { ch_gtf_new_transcripts }

        ch_versions = ch_versions.mix(TRANSCRIPT_RECONSTRUCTION.out.versions)
    }

    if (!params.skip_class) {
        CLASSIFICATION_POTENTIAL_CODING(
            ch_gtf_new_transcripts,
            params.annotation,
            params.reference,
        )
        ch_versions = ch_versions.mix(CLASSIFICATION_POTENTIAL_CODING.out.versions)
    }
    
    //
    // Collate and save software versions - FIX: Add sort: false to prevent sort.index issues
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info", 
            name: 'nf_core_pipeline_software_mqc_versions.yml', 
            sort: true,  
            newLine: true
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    
    // FIX: Separate collectFile operations to avoid race conditions
    ch_workflow_summary_file = ch_workflow_summary
        .collectFile(
            name: 'workflow_summary_mqc.yaml',
            sort: true 
        )
    
    ch_methods_description_file = ch_methods_description
        .collectFile(
            name: 'methods_description_mqc.yaml', 
            sort: false  
        )
    
    // Combine all multiqc files
    ch_multiqc_files = ch_multiqc_files
        .mix(ch_workflow_summary_file)
        .mix(ch_collated_versions)
        .mix(ch_methods_description_file)

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions = ch_versions // channel: [ path(versions.yml) ]
}
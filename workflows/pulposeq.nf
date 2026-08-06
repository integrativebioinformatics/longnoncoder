/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { NOVEL_TRANSCRIPTS                 } from '../modules/local/metadata_refinement/novel_transcripts/main'
include { SUBSET_BAMBU_COUNTS               } from '../modules/local/metadata_refinement/subset_counts/main'
include { SUBSET_BAMBU_GTF                  } from '../modules/local/metadata_refinement/subset_bambu_gtf/main'
include { BAMBU_VALIDATE                    } from '../modules/local/metadata_refinement/validate_counts/main'
include { KNOWN_TRANSCRIPTS                 } from '../modules/local/metadata_refinement/known_transcripts/main'
include { POST_REFINEMENT                   } from '../modules/local/post_refinement/main'
include { RENDER_REPORT                     } from '../modules/local/report/main'
include { QC_FILT                           } from '../subworkflows/local/qc'
include { ALIGNMENT                         } from '../subworkflows/local/alignment'
include { TRANSCRIPT_RECONSTRUCTION         } from '../subworkflows/local/transcript_reconstruction'
include { CLASSIFICATION                    } from '../subworkflows/local/classification'
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_pulposeq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PULPOSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions             = channel.empty()
    ch_multiqc_files        = channel.empty()
    ch_gtf_new_transcripts  = channel.empty()
    ch_reads_for_alignment  = ch_samplesheet
    def run_alignment       = !params.skip_alignment
    def run_classification  = run_alignment && !params.skip_class

    //
    // Run QC workflow
    //
    if (!params.skip_qc) {
        QC_FILT (
            ch_samplesheet
        )
        ch_reads_for_alignment = QC_FILT.out.filt_reads
        ch_multiqc_files = ch_multiqc_files.mix(QC_FILT.out.multiqc)
        ch_versions = ch_versions.mix(QC_FILT.out.versions)
    }

    //
    // Run alignment workflow
    //
    if (run_alignment) {
        ALIGNMENT(ch_reads_for_alignment)

        if (!params.skip_alignment_qc) {
            ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.multiqc)
        }
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)

        TRANSCRIPT_RECONSTRUCTION (
            ALIGNMENT.out.bam,
            params.reference,
            params.annotation
        )

        TRANSCRIPT_RECONSTRUCTION.out.gtf_new_transcripts
            .set { ch_gtf_new_transcripts }

        ch_versions = ch_versions.mix(TRANSCRIPT_RECONSTRUCTION.out.versions)
    }

    if (run_classification) {
        CLASSIFICATION (
            ch_gtf_new_transcripts,
            params.annotation,
            params.reference
        )
        ch_versions = ch_versions.mix(CLASSIFICATION.out.versions)

        SUBSET_BAMBU_COUNTS (
            TRANSCRIPT_RECONSTRUCTION.out.gene_counts,
            TRANSCRIPT_RECONSTRUCTION.out.transcript_counts,
            TRANSCRIPT_RECONSTRUCTION.out.CPM,
            TRANSCRIPT_RECONSTRUCTION.out.full_length,
            TRANSCRIPT_RECONSTRUCTION.out.unique_counts
        )

        NOVEL_TRANSCRIPTS (
            ch_gtf_new_transcripts,
            CLASSIFICATION.out.annotated_gtf,
            CLASSIFICATION.out.tmap,
            CLASSIFICATION.out.predictions,
            SUBSET_BAMBU_COUNTS.out.counts_transcript_subset,
            SUBSET_BAMBU_COUNTS.out.counts_gene_subset,
            file("${projectDir}/bin/novel_transcripts.R", checkIfExists: true) // <- ADDED SCRIPT PATH
        )
        ch_versions = ch_versions.mix(NOVEL_TRANSCRIPTS.out.versions)

        BAMBU_VALIDATE (
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            SUBSET_BAMBU_COUNTS.out.counts_gene_subset,
            SUBSET_BAMBU_COUNTS.out.counts_transcript_subset,
            SUBSET_BAMBU_COUNTS.out.cpm_transcript_subset,
            SUBSET_BAMBU_COUNTS.out.full_length_counts_transcript_subset,
            SUBSET_BAMBU_COUNTS.out.unique_counts_transcript_subset
        )

        KNOWN_TRANSCRIPTS (
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.counts_gene_validated,
            TRANSCRIPT_RECONSTRUCTION.out.gtf_all_transcripts,
            file("${projectDir}/bin/known_transcripts.R", checkIfExists: true) // <- ADDED SCRIPT PATH
        )
        ch_versions = ch_versions.mix(KNOWN_TRANSCRIPTS.out.versions)

        SUBSET_BAMBU_GTF (
            TRANSCRIPT_RECONSTRUCTION.out.gtf_all_transcripts,
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.full_length_counts_transcript_validated,
            BAMBU_VALIDATE.out.unique_counts_transcript_validated
        )
        ch_versions = ch_versions.mix(SUBSET_BAMBU_GTF.out.versions)

        //
        // Regenerate the Bambu PCA and heatmaps from the validated transcriptome
        //
        POST_REFINEMENT (
            TRANSCRIPT_RECONSTRUCTION.out.rds_transcript,
            TRANSCRIPT_RECONSTRUCTION.out.rds_gene,
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.counts_gene_validated,
            file("${projectDir}/bin/post_refinement.R", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(POST_REFINEMENT.out.versions)

        RENDER_REPORT (
            BAMBU_VALIDATE.out.counts_gene_validated,
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.full_length_counts_transcript_validated,
            BAMBU_VALIDATE.out.unique_counts_transcript_validated,
            KNOWN_TRANSCRIPTS.out.transcriptome_metadata,
            KNOWN_TRANSCRIPTS.out.protein_coding_metadata,
            KNOWN_TRANSCRIPTS.out.lncrna_metadata,
            KNOWN_TRANSCRIPTS.out.protein_coding_exonlength,
            KNOWN_TRANSCRIPTS.out.lncrna_exonlength,
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            NOVEL_TRANSCRIPTS.out.novel_lncrna_exon_lengths,
            NOVEL_TRANSCRIPTS.out.novel_mrna_exon_lengths,
            TRANSCRIPT_RECONSTRUCTION.out.pca,
            TRANSCRIPT_RECONSTRUCTION.out.pca_grouped,
            TRANSCRIPT_RECONSTRUCTION.out.h_gene,
            TRANSCRIPT_RECONSTRUCTION.out.h_transcript,
            POST_REFINEMENT.out.pca,
            POST_REFINEMENT.out.pca_grouped,
            POST_REFINEMENT.out.h_gene,
            POST_REFINEMENT.out.h_transcript,
            file("${projectDir}/bin/report.qmd", checkIfExists: true) // <- ADDED SCRIPT PATH
        )
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }
//
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'pulposeq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def multiqc_config_file                   = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def multiqc_custom_config_file            = params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : []
    def multiqc_logo_file                     = params.multiqc_logo ? file(params.multiqc_logo, checkIfExists: true) : []
    def summary_params                        = paramsSummaryMap(workflow)
    def ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    def multiqc_configs = multiqc_custom_config_file ? [multiqc_config_file, multiqc_custom_config_file] : [multiqc_config_file]

    MULTIQC (
        ch_multiqc_files.collect().map { files ->
            [ [id: 'multiqc'], files, multiqc_configs, multiqc_logo_file, [], [] ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }.toList()
    versions       = ch_versions

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

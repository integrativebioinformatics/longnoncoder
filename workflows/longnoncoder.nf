/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { NOVEL_TRANSCRIPTS                 } from '../modules/local/tx_annotation/novel_transcripts/main'
include { SUBSET_BAMBU_COUNTS               } from '../modules/local/tx_annotation/subset_counts/main'
include { SUBSET_BAMBU_GTF                  } from '../modules/local/tx_annotation/subset_bambu_gtf/main'
include { BAMBU_VALIDATE                    } from '../modules/local/tx_annotation/validate_counts/main'
include { KNOWN_TRANSCRIPTS                 } from '../modules/local/tx_annotation/known_transcripts/main'
include { RENDER_REPORT                     } from '../modules/local/report/main'
include { QC_FILT                           } from '../subworkflows/local/qc'
include { ALIGNMENT                         } from '../subworkflows/local/alignment'
include { TRANSCRIPT_RECONSTRUCTION         } from '../subworkflows/local/transcript_reconstruction'
include { CLASSIFICATION_POTENTIAL_CODING   } from '../subworkflows/local/classification_codingpotential.nf'
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { paramsSummaryMultiqc              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_longnoncoder_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LONGNONCODER {

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
        CLASSIFICATION_POTENTIAL_CODING (
            ch_gtf_new_transcripts,
            params.annotation,
            params.reference
        )
        ch_versions = ch_versions.mix(CLASSIFICATION_POTENTIAL_CODING.out.versions)

        SUBSET_BAMBU_COUNTS (
            TRANSCRIPT_RECONSTRUCTION.out.gene_counts,
            TRANSCRIPT_RECONSTRUCTION.out.transcript_counts,
            TRANSCRIPT_RECONSTRUCTION.out.CPM,
            TRANSCRIPT_RECONSTRUCTION.out.full_length,
            TRANSCRIPT_RECONSTRUCTION.out.unique_counts
        )

        NOVEL_TRANSCRIPTS (
            ch_gtf_new_transcripts,
            CLASSIFICATION_POTENTIAL_CODING.out.annotated_gtf,
            CLASSIFICATION_POTENTIAL_CODING.out.tmap,
            CLASSIFICATION_POTENTIAL_CODING.out.predictions,
            SUBSET_BAMBU_COUNTS.out.counts_transcript_filtered,
            SUBSET_BAMBU_COUNTS.out.counts_gene_filtered
        )
        ch_versions = ch_versions.mix(NOVEL_TRANSCRIPTS.out.versions)

        BAMBU_VALIDATE (
            NOVEL_TRANSCRIPTS.out.novel_combined_metadata,
            SUBSET_BAMBU_COUNTS.out.counts_gene_filtered,
            SUBSET_BAMBU_COUNTS.out.counts_transcript_filtered,
            SUBSET_BAMBU_COUNTS.out.cpm_transcript_filtered,
            SUBSET_BAMBU_COUNTS.out.full_length_counts_transcript_filtered,
            SUBSET_BAMBU_COUNTS.out.unique_counts_transcript_filtered
        )

        KNOWN_TRANSCRIPTS (
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.counts_gene_validated,
            TRANSCRIPT_RECONSTRUCTION.out.gtf_all_transcripts
        )
        ch_versions = ch_versions.mix(KNOWN_TRANSCRIPTS.out.versions)

        SUBSET_BAMBU_GTF (
            TRANSCRIPT_RECONSTRUCTION.out.gtf_all_transcripts,
            BAMBU_VALIDATE.out.counts_transcript_validated,
            BAMBU_VALIDATE.out.full_length_counts_transcript_validated,
            BAMBU_VALIDATE.out.unique_counts_transcript_validated
        )
        ch_versions = ch_versions.mix(SUBSET_BAMBU_GTF.out.versions)

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
            NOVEL_TRANSCRIPTS.out.novel_mrna_exon_lengths
        )
        ch_versions = ch_versions.mix(RENDER_REPORT.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_config                     = channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config              = params.multiqc_config ? channel.fromPath(params.multiqc_config, checkIfExists: true) : channel.empty()
    def ch_multiqc_logo                       = params.multiqc_logo ? channel.fromPath(params.multiqc_logo, checkIfExists: true) : channel.empty()
    def summary_params = paramsSummaryMap(workflow)
    def ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                          = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                          = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                          = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    def ch_multiqc_all_configs = ch_multiqc_config
        .mix(ch_multiqc_custom_config)
        .collect()

    MULTIQC (
        ch_multiqc_files.collect()
            .combine(ch_multiqc_all_configs.map { configs -> [configs] })
            .combine(ch_multiqc_logo.collect().ifEmpty([[]]))
            .map { files, configs, logo ->
                [ [id: 'multiqc'], files, configs, logo, [], [] ]
            }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { meta, report -> report }.toList()
    versions       = ch_versions
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
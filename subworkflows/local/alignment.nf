//
// MODULE: from the nf-core modules repository
//
include { MINIMAP2_ALIGN } from '../../modules/nf-core/minimap2/align/main'
include { NANOCOMP as NANOCOMP_MAPPING } from '../../modules/nf-core/nanocomp/main'

/*
========================================================================================
    RUN ALIGNMENT WORKFLOW
========================================================================================
*/

workflow ALIGNMENT {
    take:
    reads

    main:
    ch_versions = Channel.empty()
    ch_bam = Channel.empty()
    ch_index = Channel.empty()
    ch_combined_mapping = Channel.empty()
    ch_alignment_qc = Channel.empty()
    ch_reference = Channel.empty()

    // Building metamap for the reference
    Channel.fromPath(params.reference)
        .map { file_path ->
            def basename = file_path.baseName
            // Extract the base name
            [basename, file_path.toString()]
        }
        .collect()
        .set { ch_reference }
    // Alignment with the minimap2 module in case no filtering is applied to read length

    MINIMAP2_ALIGN(
        reads,
        ch_reference,
        params.bam_format,
        params.bam_index_extension,
        params.cigar_paf_format,
        params.cigar_bam,
    )

    MINIMAP2_ALIGN.out.bam.set { ch_bam }

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first().ifEmpty(null))

    if (!params.skip_alignment_qc) {

        // Collect only BAM files for NanoComp analysis
        // The indices will be available through symlinks in the work directory
        ch_bam
            .collect { it[1] }
            .map { filelist -> [[id: "All"], filelist] }
            .set { ch_combined_mapping }

        NANOCOMP_MAPPING(
            ch_combined_mapping
        )

        ch_alignment_qc = ch_alignment_qc.mix(NANOCOMP_MAPPING.out.stats_txt.collect { it[1] }.ifEmpty([]))

        ch_versions = ch_versions.mix(NANOCOMP_MAPPING.out.versions.first().ifEmpty(null))
    }

    emit:
    index = ch_index
    bam = ch_bam
    reference = ch_reference
    multiqc = ch_alignment_qc
    versions = ch_versions
}
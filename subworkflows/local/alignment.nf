//
// MODULES: MINIMAP2_ALIGN and NANOCOMP
//
include { MINIMAP2_ALIGN                  } from '../../modules/nf-core/minimap2/align/main'
include { NANOCOMP as NANOCOMP_MAPPING    } from '../../modules/nf-core/nanocomp/main'

/*
========================================================================================
    RUN ALIGNMENT WORKFLOW
========================================================================================
*/

workflow ALIGNMENT {
    take:
    reads

    main:
    ch_versions         = channel.empty()
    ch_bam              = channel.empty()
    ch_index            = channel.empty()
    ch_combined_mapping = channel.empty()
    ch_alignment_qc     = channel.empty()
    ch_reference        = channel.empty()

    // Building metamap for the reference
    channel
        .fromPath(params.reference)
        .map { file_path ->
            def basename = file_path.baseName
            [basename, file_path.toString()]
        }
        .collect()
        .set { ch_reference }

    // Alignment with the minimap2 module
    MINIMAP2_ALIGN (
        reads,
        ch_reference,
        params.bam_format,
        params.bam_index_extension,
        params.cigar_paf_format,
        params.cigar_bam
    )

    MINIMAP2_ALIGN.out.bam
        .set { ch_bam }

    if (!params.skip_alignment_qc) {

        ch_bam
            .collect { item -> item[1] }
            .map { filelist -> [[id: "All"], filelist] }
            .set { ch_combined_mapping }

        NANOCOMP_MAPPING (
            ch_combined_mapping
        )

        ch_alignment_qc = ch_alignment_qc.mix(NANOCOMP_MAPPING.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

        ch_versions = ch_versions.mix(NANOCOMP_MAPPING.out.versions.ifEmpty(null))
    }

    emit:
    index    = ch_index
    bam      = ch_bam
    reference = ch_reference
    multiqc  = ch_alignment_qc
    versions = ch_versions
}
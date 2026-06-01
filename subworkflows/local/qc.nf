//
// MODULE: Installed directly from nf-core/modules
//

include { NANOCOMP as NANOCOMP_RAW          } from '../../modules/nf-core/nanocomp/main'
include { NANOCOMP as NANOCOMP_FILT         } from '../../modules/nf-core/nanocomp/main'
include { CHOPPER                           } from '../../modules/nf-core/chopper/main'

/*
========================================================================================
    RUN QC_FILT WORKFLOW
========================================================================================
*/


workflow QC_FILT {
    take:
    reads

    main:
    ch_versions          = channel.empty()
    ch_multiqc_raw       = channel.empty()
    ch_multiqc_filt      = channel.empty()
    ch_multiqc_all       = channel.empty()
    ch_combined_raw      = channel.empty()
    ch_combined_filtered = channel.empty()
    ch_reads             = reads

    // Running nanocomp on raw reads

    ch_reads
        .collect { item -> item[1] }
        .map { filelist -> [[id: "All"], filelist] }
        .set { ch_combined_raw }

    NANOCOMP_RAW(ch_combined_raw)

    ch_versions = ch_versions.mix(NANOCOMP_RAW.out.versions)

    // Generating a multiqc file for raw reads report
    ch_multiqc_raw = ch_multiqc_raw.mix(NANOCOMP_RAW.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

    ch_multiqc_all = ch_multiqc_all.mix(ch_multiqc_raw.ifEmpty([]))

    // Putting conditional to whether run filtering on samples
    if (!params.skip_filtering) {

        CHOPPER(ch_reads, [])

        // Removed the broken CHOPPER.out.versions line entirely

        // Running quality check in filtered reads

        // Running quality check in filtered reads
        CHOPPER.out.fastq
            .collect { item -> item[1] }
            .map { filelist -> [[id: "All"], filelist] }
            .set { ch_combined_filtered }

        NANOCOMP_FILT(ch_combined_filtered)

        ch_multiqc_filt = ch_multiqc_filt.mix(NANOCOMP_FILT.out.stats_txt.collect { item -> item[1] }.ifEmpty([]))

        ch_multiqc_all = ch_multiqc_all.mix(ch_multiqc_filt.ifEmpty([]))

        // Putting the output of chopper as new ch_reads
        ch_reads = CHOPPER.out.fastq
    }

    emit:
    filt_reads = ch_reads
    multiqc    = ch_multiqc_all
    versions   = ch_versions
}
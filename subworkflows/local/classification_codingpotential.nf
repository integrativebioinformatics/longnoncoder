//
// MODULE: nf-core modules
//
include { GFFCOMPARE } from '../../modules/nf-core/gffcompare/main'
include { GFFREAD    } from '../../modules/nf-core/gffread/main'

//
// MODULE: Local to the pipeline
//
include { RNAMINING  } from '../../modules/local/rnamining/main'

/*
========================================================================================
    RUN CLASSIFICATION_POTENTIAL_CODING WORKFLOW
========================================================================================
*/

workflow CLASSIFICATION_POTENTIAL_CODING {
    take:
    gtf         // path: GTF from transcript reconstruction
    annotation  // path/string: reference annotation GTF (params.annotation)
    reference   // path/string: reference genome FASTA (params.reference)

    main:
    ch_versions = channel.empty()

    //
    // Prepare inputs for nf-core/gffcompare
    // nf-core module expects: tuple(meta, gtfs), tuple(meta2, fasta, fai), tuple(meta3, reference_gtf)
    //
    def ch_gffcompare_gtf = gtf.map { gtf_file ->
        [ [id: 'gffcompare'], gtf_file ]
    }
    def ch_gffcompare_fasta = [ [id: 'genome'], [], [] ]  // not using fasta input
    def ch_gffcompare_ref   = [ [id: 'reference'], file(annotation) ]

    GFFCOMPARE(
        ch_gffcompare_gtf,
        ch_gffcompare_fasta,
        ch_gffcompare_ref
    )

    // Extract bare paths from tuple outputs for downstream compatibility
    ch_annotated_gtf = GFFCOMPARE.out.annotated_gtf.map { meta, f -> f }
    ch_tmap          = GFFCOMPARE.out.tmap.map { meta, f -> f }

    //
    // Prepare inputs for nf-core/gffread
    // nf-core module expects: tuple(meta, gff), path(fasta)
    //
    def ch_gffread_input = gtf.map { gtf_file ->
        [ [id: 'gffread'], gtf_file ]
    }

    GFFREAD(
        ch_gffread_input,
        reference
    )

    // Extract bare fasta path for downstream compatibility
    ch_gffread_fasta = GFFREAD.out.gffread_fasta.map { meta, f -> f }

    //
    // RNAMINING: coding potential prediction
    //
    RNAMINING(
        ch_gffread_fasta
    )

    ch_predictions = RNAMINING.out.preds

    emit:
    annotated_gtf = ch_annotated_gtf
    tmap          = ch_tmap
    gffread_fasta = ch_gffread_fasta
    predictions   = ch_predictions
    versions      = ch_versions
}
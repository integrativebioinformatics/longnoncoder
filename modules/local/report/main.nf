process RENDER_REPORT {
    tag "Rendering report"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 'oras://community.wave.seqera.io/library/quarto_r-cowplot_r-knitr_r-rcolorbrewer_pruned:9183d70916f57fd4' : 'community.wave.seqera.io/library/quarto_r-cowplot_r-knitr_r-rcolorbrewer_pruned:a02cb51064eb0c90'}"

    input:
        path(notebook)

        // bambu outputs
        path(counts_genes)
        path(counts_transcript)
        path(fl_counts_transcript)
        path(u_counts_transcript)

        // metadata csv tables
        path(transcriptome_meta)
        path(protein_coding_meta)
        path(lncrna_meta)
        path(protein_coding_exonlength)
        path(lncrna_exonlength)
        path(novel_lncrna_meta)
        path(novel_lncrna_exonlength)
        path(novel_protein_coding_exonlength)

        // plots
        path(plot_heatmap_gene)
        path(plot_heatmap_tcpt)
        path(plot_pca_grouped)
        path(plot_pca)

    output:
        path("*.html")                      , emit: report
        path("Bambu_assembly_summary.csv")  , emit: bambu_assembly_summary
        // path("*.png")                       , emit: generated_plots_png
        // path("*.pdf")                       , emit: generated_plots_pdf

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        export XDG_CACHE_HOME=/tmp/quarto_cache_home
        export XDG_DATA_HOME=/tmp/quarto_data_home

        ENV_QUARTO=/opt/conda/etc/conda/activate.d/quarto.sh
        set +u
        if [ -z "\${QUARTO_DENO}" ] && [ -f "\${ENV_QUARTO}" ]; then
            source "\${ENV_QUARTO}"
        fi
        set -u

        cp ${notebook} report.qmd
        quarto render report.qmd \\
            -P counts_genes:${counts_genes} \\
            -P counts_transcript:${counts_transcript} \\
            -P fl_counts_transcript:${fl_counts_transcript} \\
            -P u_counts_transcript:${u_counts_transcript} \\
            -P transcriptome_meta:${transcriptome_meta} \\
            -P protein_coding_meta:${protein_coding_meta} \\
            -P lncrna_meta:${lncrna_meta} \\
            -P protein_coding_exonlength:${protein_coding_exonlength} \\
            -P lncrna_exonlength:${lncrna_exonlength} \\
            -P novel_transcriptome_meta:${novel_transcriptome_meta} \\
            -P novel_lncrna_exonlength:${novel_lncrna_exonlength} \\
            -P novel_protein_coding_exonlength:${novel_protein_coding_exonlength} \\
            -P plot_heatmap_gene:${plot_heatmap_gene} \\
            -P plot_heatmap_tcpt:${plot_heatmap_tcpt} \\
            -P plot_pca_grouped:${plot_pca_grouped} \\
            -P plot_pca:${plot_pca} \\
            --to html
        """
    stub:
        """
        touch report.html
        touch Bambu_assembly_summary.csv
        touch generated_plot.png
        touch generated_plot.pdf
        """

}
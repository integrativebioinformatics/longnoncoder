## integrativebioinformatics/longnoncoder <img src="figures/logo.svg" align=right height="200px"/>

`LongNonCoder` is a Nextflow pipeline designed for isoform-level lncRNA discovery and characterization from long-read RNA-seq data. The workflow encompasses QC, mapping, transcriptome assembly and quantification, followed by a detailed final characterization of the entire transcriptome with particular emphasis on lncRNA structure and isoforms across known annotations and novel candidates.

For more details and further functionality, please refer to the [usage](docs/usage.md) and [output](docs/output.md) documentations.

> [!IMPORTANT]
> LongNonCoder is compatible Ensembl or GENCODE reference genomes and annotations from the following organisms:
> *Homo sapiens, Mus musculus, Danio rerio, Anolis carolinensis*, *Chrysemys picta belli, Eptatetrus burgeri, Gallus gallus, Latimeria chalumnae, Monodelphis domestica, Notechis scutatus, Ornithorhynchus anatinus*, *Petromyzon marinus, Sphenodon punctatus,* and *Xenopus tropicalis.*
>
>  **In the next releases, we plan to update the pipeline workflow to cover more organisms or even more general taxonomic classes.**

**The workflow**

![longnoncoder workflow](figures/LongNonCoder.drawio.png)

We can describe each step of the workflow as follows:

1.  Quality control of reads ([NanoComp](https://github.com/wdecoster/nanocomp "wdecoster/nanocomp"))
2.  Filtering and trimming ([chopper](https://github.com/wdecoster/chopper "wdecoster/chopper"))
3.  Mapping to a genome reference ([minimap2](https://github.com/lh3/minimap2 "lh3/minimap2") and [samtools](https://github.com/samtools/samtools "samtools"))
4.  Quality control of mapped reads ([NanoComp](https://github.com/wdecoster/nanocomp "wdecoster/nanocomp"))
5.  Transcriptome Assembly ([Bambu](https://www.bioconductor.org/packages/release/bioc/html/bambu.html)))
6.  Compare novel transcripts to the annotation reference ([GffCompare](https://github.com/gpertea/gffcompare "gpertea/gffcompare"))
7.  Convert novel transcripts `GTF` file to `FASTA` ([GffRead](https://github.com/gpertea/gffread "gpertea/gffread"))
8.  Predict transcripts as protein-coding or non-coding ([RNAmining](https://gitlab.com/integrativebioinformatics/RNAmining "integrativebioinformatics/RNAmining"))
9.  Gather all data from previous steps and generate informative and re-usable metadata `.csv` and `GTF` files for both novel and annotated transcripts ([tidyverse](https://tidyverse.org/), [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html), [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html), and [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html))
10. Provide a report and data visualization for the full transcriptome, with emphasis on lncRNAs ([Quarto](https://quarto.org/), [tidyverse](https://tidyverse.org/), [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html), [scales](https://cran.r-project.org/web/packages/scales/index.html), etc)
11. Gather all possible QC information from the previous steps ([MultiQC](https://github.com/MultiQC/MultiQC "MultiQC"))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data. The pipeline is compatible with Docker and Singularity/Apptainer.

> [!WARNING]
> You might be able to customize the modules and configuration to run the pipeline with Conda, although this is NOT recommended nor was evaluated during development.

You can run an example test by following the instructions:

Enter the `test_data` folder

``` bash
cd test_data
```

Download and unzip the reference `FASTA` and `GTF` files, and also download the fastq.gz files:

Make the file executable!!

``` bash
chmod +x download-ref-fastq.sh
```

Execute the script

``` bash
./download-ref-fastq.sh
```

Add YOUR full path for the samples in the `samplesheet.csv` ([file](test_data/samplesheet.csv)). For example, your full path for a sample could be:

`home/user/longnoncoder/test_data/thesample.fastq.gz`

Go back to the main directory and execute the test!

``` bash
cd ..
```

``` bash
nextflow run main.nf -profile test,singularity -params-file test_data/testing.yml
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option and input a `yaml` parameters file. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration ***except for parameters***; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

> [!TIP]
> LongNonCoder is configured in the copy mode for publishing the output directories, creating file replicas of your final results originally created at the `work/` directory. This directory not only contains your results, but all other intermediate execution files (e.g. temporary files, `.command.sh` and `.command.log`). After running the pipeline, remember that this setup allows you to safely delete the pipeline's `work/` directory without losing your published results.

> [!NOTE]
> Copying large datasets can significantly take a long time to complete or occupy a large amount of disk space. If you require other setup to optimize resource usage, follow the instructions from the [Nextflow documentation](https://docs.seqera.io/nextflow/reference/process#publishdir) to change the `publishDir` setting in the [`nextflow.config`](nextflow.config) file.

## Citations

If you use LongNonCoder in your research, please consider citing it.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) initative, and reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.

An extensive list of references for the tools incorporated by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

## Acknowledgments

#### Development & Contributions

The `LongNonCoder` pipeline was originally developed by **Bárbara Borges** ([\@borgessbarbara](https://github.com/borgessbarbara)). We extend our sincere thanks to **Lucas Freitas** ([\@lfreitasl](https://github.com/lfreitasl)), **Gleison Azevedo** ([\@gleisonm](https://github.com/gleisonm)) and **João Cavalcante** ([\@jvfe](https://github.com/jvfe)) for their significant contributions and assistance during development.

#### Supervision & Collaborations

This project was developed under the supervision of **Vinicius Maracaja-Coutinho** and **Thaís Gaudencio**, in collaboration with **Rodrigo Dalmolin** and the [Dalmolin Systems Biology Group](https://github.com/dalmolingroup)

#### Computational resources

This project was supported by the High-Performance Computing Center at UFRN (NPAD/UFRN) and the National Laboratory for High Performance Computing (NLHPC) (CCSS210001) at UChile.

#### Funding

CAPES (001), CNPq (MCTI/FNDCT 445067/2024-1), FONDECYT-ANID Postdoctorado (3250452), FONDECYT-ANID (1211731), FONDAP-ANID (15130011 and 1523A0008), and Anillo-ANID (ATE220016).

#### Laboratories & Institutions Involved

Bárbara Borges <sup>1,2,3</sup>, Lucas Freitas <sup>4,5</sup>, Gleison Azevedo <sup>1</sup>, João Cavalcante <sup>1</sup>, Rodrigo Dalmolin <sup>1</sup>, Thaís Gaudencio <sup>1,3</sup>, Vinicius Maracaja-Coutinho <sup>1,2,6,7</sup>

<details>

<summary>Affiliations</summary>

<sup>1</sup> Bioinformatics Multidisciplinary Environment, Instituto Metrópole Digital, Universidade Federal do Rio Grande do Norte, Natal, Brazil

<sup>2</sup> Unidad de Genómica Avanzada, Advanced Center for Chronic Diseases, Facultad de Ciencias Químicas y Farmacéuticas, Universidad de Chile, Santiago, Chile

<sup>3</sup> Artificial Intelligence Applications, Universidade Federal da Paraíba, João Pessoa, Brazil

<sup>4</sup> Faculty of Natural Sciences, Universität Hohenheim, Stuttgart, Germany

<sup>5</sup> Staatliches Museum für Naturkunde Stuttgart, Department of Biodiversity Monitoring, Stuttgart, Germany

<sup>6</sup> Laboratório de Medicina e Saúde Pública de Precisão, Instituto Gonçalo Moniz, Fundação Oswaldo Cruz, Fiocruz, Salvador, Brazil

<sup>7</sup> Instituto Nacional de Ciência e Tecnologia em Saúde Digital (INCT-DigiSaúde)

</details>

<p align="center">

<picture> <source media="(prefers-color-scheme: dark)" srcset="figures/institutional-logos-dark-theme.png" height=150> <img src="figures/institutional-logos-light.png" height="150"/> </picture>

</p>

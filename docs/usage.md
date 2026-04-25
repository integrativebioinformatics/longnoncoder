# Usage

> *Pipeline parameters are described below and should match the current pipeline schema/configuration.*

## Introduction

<!-- TODO integrative/bioinformatics: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the parameter `--input` in the bash command to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

``` bash
--input '[path to samplesheet file]'
```

### Full samplesheet

A final `samplesheet.csv` file consisting of single-end data may look something like the one below. This is for 6 samples, where we have 2 experimental groups and 3 replicates per group.

``` csv
sample,group,fastq
R1_H1975,H1975,home/user/R1_H1975.fastq.gz
R2_H1975,H1975,home/user/R2_H1975.fastq.gz
R3_H1975,H1975,home/user/R3_H1975.fastq.gz
R1_HCC827,HCC827,home/user/R1_HCC827.fastq.gz
R2_HCC827,HCC827,home/user/R2_HCC827.fastq.gz
R3_HCC827,HCC827,home/user/R3_HCC827.fastq.gz
```

| Column   | Description                                                                                     |
|----------|-------------------------------------------------------------------------------------------------|
| `sample` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`).                        |
| `group`  | Experimental group name. For example: `treatment` vs `control` or `cell_line1` vs `cell_line2`  |
| `fastq`  | Full path to FastQ file for ONT or PacBio long-reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                     |

Another [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Pipeline parameters

After seeting up the samplesheet, follow up to set the pipeline parameters.

> [!WARNING]
> Relative paths might cause issues! Always input double-quoted **FULL PATHS** (e.g. `"/full/path/to/file/file.extension"`). 

> [!TIP]
> Stay up to date! Remember to always use the **lastest Ensembl releases** for reference genomes and annotations.

| Parameter        | Description                                                                              |
|------------------|------------------------------------------------------------------------------------------|
| `input`          | Full path to the `samplesheet.csv` file                                                  |
| `outdir`         | Full path to the output directory where you want the results to be saved                 |
| `skip_qc`        | Skip entire QC process                                                                   |
| `skip_filtering` | Skip filtering with `chopper`, just perform `NanoComp` reporting                         |
| `minqual`        | Minimum read average Phred-score quality cut-off                                         |
| `minlen`         | Minimum read length (bp)                                                                 |
| `maxgc`          | Maximum GC content (%)                                                                   |
| `mingc`          | Minimum GC content (%)                                                                   |
| `headcrop`       | Cut `x` number of bases from the beginning of the reads                                  |
| `tailcrop`       | Cut `y` number of bases from the end of the reads                                        |
| `skip_alignment` | Skip mapping/alignment to genome reference (runs MultiQC, then finishes pipeline execution)     |
| `skip_alignment_qc` | Skip mapping/alignment QC (might reduce resource usage and speed up execution when running large datasets)        |
| `reference`      | Full path to a reference genome `FASTA` file from Ensembl                                |
| `annotation`     | Full path to a reference annotation `GTF` file from Ensembl                              |
| `skip_class`     | Skip transcriptome characterization (runs MultiQC, then finishes pipeline execution)                                  |
| `organism`       | Organism scientific name (e.g. `"Homo_sapiens"`)                                         |
| `ensembl_organism_dataset` | Reference BiomaRt dataset (e.g. `"hsapiens_gene_ensembl"`)                    |
| `ensembl_version` | Number of Ensembl release version (e.g. `114`)                                          |

## Running the pipeline {#running-the-pipeline}

The typical command for running the pipeline is as follows:

``` bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --minqual [value] --reference [fasta] --annotation [gtf] --organism [Genus_species] --ensembl_organism_dataset [Gspecies_gene_ensembl] --ensembl_version [release number] -profile [profile: light, medium, large, etc],[executor profile: docker/singularity]
```

Note that the pipeline will create the following files in your working directory:

``` bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

``` bash
nextflow run main.nf -profile docker -params-file params.yaml
```

with `params.yaml` containing:

``` yaml
input: './samplesheet.csv'
outdir: './results/'
<...>
```

> [!TIP] 
> Follow the examples from the [test](../test_data/testing.yml) and the [example run](../examplerun.yml).

### Updating the pipeline

``` bash
git clone https://github.com/integrativebioinformatics/longnoncoder.git
```

When you run the above command, Git automatically clones the pipeline code from GitHub and stores it. When running the pipeline after this, it will always use this version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the commits in the pipeline:

``` bash
git fetch origin main
```
``` bash
git pull origin main
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [integrativebioinformatics/longnoncoder releases page](https://github.com/integrativebioinformatics/longnoncoder/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a *single* hyphen (pipeline parameters use a double-hyphen).


### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is *not* recommended, since it can lead to different results on different machines dependent on the computer enviroment. You can also create your own profile!

-   `test`
    -   A profile with configuration for testing that consumes low resources
-   `light`
    -   A profile for small-scale data, consumes low resources
-   `medium`
    -   A profile for medium-scale data, consumes medium resources
-   `large`
    -   A profile for large-scale data, consumes high resources
-   `docker`
    -   A generic configuration profile to be used with [Docker](https://docker.com/)
-   `singularity`
    -   A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
-   `apptainer`
    -   A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf=core/configs` repository with the addition of your config file, associated documentation file (see examples in [nf-core/configs/docs](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

``` bash
NXF_OPTS='-Xms1g -Xmx16g'
```

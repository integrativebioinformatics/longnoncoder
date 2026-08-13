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

| Column | Description |
|-------------|-----------------------------------------------------------|
| `sample` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `group` | Experimental group name. For example: `treatment` vs `control` or `cell_line1` vs `cell_line2` |
| `fastq` | Full path to FastQ file for ONT or PacBio long-reads. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

Another [example samplesheet] has been provided with the pipeline.

  [example samplesheet]: ../assets/samplesheet.csv

## Pipeline parameters

After seeting up the samplesheet, follow up to set the pipeline parameters.

> [!WARNING]
> Relative paths might cause issues! Always input double-quoted **FULL PATHS** (e.g. `"/full/path/to/file/file.extension"`).

> [!TIP]
> Stay up to date! Remember to always use the **lastest Ensembl releases** for reference genomes and annotations.

| Parameter | Description |
|----------------|--------------------------------------------------------|
| `input` | Full path to the `samplesheet.csv` file |
| `outdir` | Full path to the output directory where you want the results to be saved |
| `skip_qc` | Skip entire QC process |
| `skip_filtering` | Skip filtering with `chopper`, just perform `NanoComp` reporting |
| `minqual` | Minimum read average Phred-score quality cut-off |
| `minlen` | Minimum read length (bp) |
| `maxlen` | Maximum read length (bp), optional — leave unset for no upper limit |
| `maxgc` | Maximum GC content (%) |
| `mingc` | Minimum GC content (%) |
| `headcrop` | Cut `x` number of bases from the beginning of the reads |
| `tailcrop` | Cut `y` number of bases from the end of the reads |
| `skip_alignment` | Skip mapping/alignment to genome reference (runs MultiQC, then finishes pipeline execution) |
| `skip_alignment_qc` | Skip mapping/alignment QC (might reduce resource usage and speed up execution when running large datasets) |
| `reference` | Full path to a reference genome `FASTA` file from Ensembl |
| `annotation` | Full path to a reference annotation `GTF` file from Ensembl |
| `library` | Sequencing library type: `ONT_cDNA`, `ONT_DRS` or `PacBio` |
| `stranded_library` | Whether the library is stranded (`true`/`false`) |
| `skip_class` | Skip transcriptome characterization (runs MultiQC, then finishes pipeline execution) |
| `organism` | Organism scientific name (e.g. `"Homo_sapiens"`) |

### Library type and strandedness

The `library` and `stranded_library` parameters are **required** whenever alignment is not skipped. Together they set the `minimap2` alignment preset and Bambu's `stranded` argument — the two are deliberately driven from a single declaration of library chemistry, so they cannot be configured independently.

| `library` | `stranded_library` | `minimap2` preset |
|--------------|-----------------------|----------------------------|
| `ONT_DRS` | automatically `true` | `-ax splice -uf -k14` |
| `PacBio` | automatically `true` | `-ax splice:hq -uf` |
| `ONT_cDNA` | `true` | `-ax splice -uf` |
| `ONT_cDNA` | `false` | `-ax splice` |

> [!IMPORTANT]
> A single execution must use one library type. Running PacBio and ONT samples together in the same run is not supported, as the difference in error profiles and library chemistry makes joint assembly difficult to interpret.

**PacBio.** The pipeline requires PacBio reads that have already been processed and stranded by PacBio's own independent standard workflows. pulposeq does not perform that preprocessing, and assumes it has already been done.

**ONT_DRS.** Direct RNA sequencing reads the native RNA molecule, so these libraries are considered automatically stranded. Setting `stranded_library: false` for `ONT_DRS` is ignored, and the pipeline emits a warning.

**ONT_cDNA.** This covers both PCR-cDNA and direct-cDNA protocols. Standard ONT cDNA library preparation protocols can sequence either the first or second cDNA strand, leaving read direction mixed or unaligned to the original mRNA's 5'-to-3' orientation. Such libraries must either go through a manual orienting process outside the pipeline, or be treated as unstranded by setting `stranded_library: false`. Tools like [Pychopper](https://github.com/epi2me-labs/pychopper) and [Restrander](https://github.com/mainguyenanhvu/Restrander) are used to detect poly(A) tails and primer signatures to fix this.

### Ensembl vs GENCODE annotations

Transcript metadata is read **directly from the annotation `GTF` you supply** with `annotation`. Nothing is queried over the network, so runs are reproducible, work on offline compute nodes, and are unaffected by Ensembl archiving older releases or retiring service endpoints. It also guarantees the metadata describes exactly the annotation used for assembly, rather than whichever release a remote server happened to serve.

Both Ensembl and GENCODE annotations are accepted, and the differences between them are handled automatically:

| | Ensembl | GENCODE |
|---------------------|-------------------------------------------|--------------------------------------|
| Biotype attributes | `gene_biotype` / `transcript_biotype` | `gene_type` / `transcript_type` |
| Identifiers | unversioned, with separate `gene_version` / `transcript_version` | versioned inline (e.g. `ENSG00000290825.2`) |
| Sequence names | `1` | `chr1` |

Both identifier forms are always produced, so the `ensembl_transcript_id` and `ensembl_transcript_id_version` columns are populated whichever source you use.

Sequence names are **not** rewritten. `chromosome_name` reads `1` with an Ensembl annotation and `chr1` with GENCODE, exactly as the file you supplied names them. Stripping the prefix would have to special-case the GRC accessions (`KI270728.1`, `GL000009.2`) that both sources use for unplaced scaffolds, and it would put the metadata tables out of step with the novel-transcript outputs, which report the sequence names that reach them from the assembly untouched.

> [!IMPORTANT]
> With GENCODE, use the **CHR** or **PRI** annotation build. The **ALL** build additionally contains alternate loci, assembly patches and haplotypes, which duplicate gene and transcript identifiers. pulposeq checks for duplicated identifiers and stops with an error rather than silently double-counting them.

Regarding the pseudoautosomal regions (PAR) of chromosome Y: the gene annotation in these regions is identical between chromosomes X and Y, and until GENCODE release 43 the chrY copies were distinguished by an `_PAR_Y` suffix appended to their identifiers. From GENCODE release 44 (equivalent to Ensembl release 110) onward, the chrY PAR annotation carries its own distinct identifiers, so within the supported release range below this is not a concern. Additionally, the GENCODE GTF contains a number of attributes not present in the Ensembl GTF, but that does not interfere with the pipeline execution nor is required.

Check detailed information at the [GENCODE FAQ][1].

  [1]: https://www.gencodegenes.org/pages/faq.html

> [!WARNING]
> We highly recommend always using the latest stable release available, whether its GENCODE or Ensembl. The pipeline was only tested from Ensembl release 112 onwards. Due to reported differences between older versions, we recommend using, at least, from GENCODE v44 (Ensembl release 110) onwards. Versions older than that might cause issues with the transcriptome characterization steps.

## Running the pipeline {#running-the-pipeline}

The typical command for running the pipeline is as follows:

``` bash
nextflow run main.nf --input ./samplesheet.csv --outdir ./results --minqual [value] --reference [fasta] --annotation [gtf] --organism [Genus_species] --library [ONT_cDNA/ONT_DRS/PacBio] --stranded_library [true/false] -profile [test/medium/large],[container runtime: docker/singularity/apptainer]
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
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications], other infrastructural tweaks (such as output directories), or module arguments (args).

  [tuning process resource specifications]: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources

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
> Follow the examples from the [test] and the [example run].

  [test]: ../test_data/testing.yml
  [example run]: ../examplerun.yml

### Updating the pipeline

``` bash
git clone https://github.com/integrativebioinformatics/pulposeq.git
```

When you run the above command, Git automatically clones the pipeline code from GitHub and stores it. When running the pipeline after this, it will always use this version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the commits in the pipeline:

``` bash
git fetch origin main
```

``` bash
git pull origin main
```

### Reproducibility

IIt is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [integrativebioinformatics/pulposeq releases page] and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

  [integrativebioinformatics/pulposeq releases page]: https://github.com/integrativebioinformatics/pulposeq/releases

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files] to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

  [parameter files]: #running-the-pipeline

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a *single* hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Apptainer, Conda) - see below.

- `test`
  - Calibrated to the bundled chromosome 1 test data. Peaks at roughly 24 GB, so it runs on a workstation.
- `medium`
  - Full reference genome and annotation with modest sample sizes.
- `large`
  - Full reference with large ONT or PacBio samples, tens of GB each. Bambu alone can need several hundred GB; see [Resource requests](#resource-requests).
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)

Please only use [Conda](https://conda.io/docs/) as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, or Apptainer. Note that pulposeq does not support conda for local modules.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.


### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Resource allocations live in `conf/`. One file is always loaded and the profile you choose overrides it:

```
conf/base.config          always loaded, defines every label
    ↓ overridden by
conf/{test,medium,large}.config   only the labels each one defines
```

A label that a profile does not define keeps its `base.config` value. This is worth knowing, because
a process whose label no config matches falls back to the generic defaults **silently** — the usual
symptom is an out-of-memory kill (exit `137`) on a step you thought you had sized. To make unmatched
selectors report themselves, add the `debug` profile:

``` bash
nextflow run main.nf -profile large,apptainer,debug -params-file params.yml
```

#### Retries

`errorStrategy` resubmits on out-of-memory and related exit codes, with `maxRetries = 1` — so at most
one retry. Whether that retry asks for **more** than the first attempt depends on how the profile is
written:

``` groovy
memory = { 16.GB * task.attempt }   // escalates: 16 GB, then 32 GB
memory = 16.GB                      // does not: 16 GB, then 16 GB again
```

Fixed values are appropriate where the input is known and unchanging, such as the `test` profile.
For real data, `task.attempt` scaling lets you allocate close to the measured requirement and still
survive an unusual sample.

#### `resourceLimits`

`resourceLimits` caps every request **after** the `task.attempt` multiplier, and it clamps silently.
Two consequences:

- It must be at or above the largest request in the file, or that label is quietly cut back.
- If you rely on retry escalation, it must be at or above **twice** the largest first attempt, or the
  retry is clamped to the value that already failed.

#### Tuning for your own system

Prefer a site config passed with `-c` over editing the repository, so your settings survive
`git pull`:

``` groovy
// ~/.nextflow/config or site.config, used with -c
process {
    executor       = 'slurm'
    queue          = 'main'
    resourceLimits = [ cpus: 256, memory: 700.GB, time: 24.h ]

    withLabel:process_bambu { memory = 700.GB }   // override a single label
}

apptainer.cacheDir = '/path/to/apptainer_images'
```

### Measured resource requirements

Peak resident memory observed in real runs. Use these to judge whether a dataset fits your hardware.

| Process | test (chr1 subset) | medium (full reference, small samples) | large (8 ONT samples, 25–30 GB each) |
|---|---|---|---|
| `CHOPPER` | 7.5 GB | ~8 GB | **54.7 GB** |
| `MINIMAP2_ALIGN` | 16.3 GB (8 cpu) | 31.1 GB (20 cpu) | 56.1 GB (32 cpu) |
| `NANOCOMP` | 1.8 GB | 1.2 GB | 33.4 GB |
| `NANOCOMP_MAPPING` | 1.5 GB | 3.5 GB | 48.4 GB |
| `BAMBU` | 13.6 GB | 63.6 GB | **639.8 GB** |
| metadata refinement | 0.9 GB | ~4 GB | 4.1 GB |
| `RENDER_REPORT` | 1.0 GB | 1.3 GB | — |

Four rules let you predict your own case rather than guessing:

**Chopper is linear in input size**, at close to **0.38 × bytes read**, verified from 2.5 GB to
143 GB of input. A 30 GB gzipped FASTQ expands to roughly 130 GB, so expect a ~50 GB peak.

**Bambu grows by roughly 54 GB per additional sample** — 532 GB at 6 samples, 640 GB at 8 — because
it holds every sample in a single object.

> [!IMPORTANT]
> This puts a ceiling on samples per run. On a 768 GB node, Bambu becomes unschedulable somewhere
> around **10 samples**. It is a property of the assembly step, not something configuration can
> tune away: beyond that you need a larger-memory node, or to split the run.

**minimap2 and NanoComp scale with thread count**, since their memory is a fixed index plus
per-thread buffers. Raising `cpus` raises memory too — NanoComp went from 23.8 GB on 1 CPU to
33.4 GB on 3.

**The metadata refinement steps scale with the annotation, not the data.** They parse the reference
`GTF`, so a full human annotation costs ~4 GB whether the input is a single small sample or eight
large ones. A chromosome-1 test run understates them by roughly fourfold.

### SLURM and `--mem-per-cpu`

Some clusters schedule by memory per core rather than by total. Nextflow supports this:

``` groovy
executor.perCpuMemAllocation = true   // Nextflow 23.10+, SLURM only
```

With it enabled, jobs are submitted as `--mem-per-cpu <memory / cpus>` instead of `--mem <memory>`.
**The `memory` directive still means total memory** — Nextflow does the division at submission — so
there is no second set of numbers to maintain, and the setting is ignored by the local executor.

What changes is that `cpus` becomes load-bearing. Every memory value you raise also raises the
per-core figure, which must stay under your partition's `MaxMemPerCPU`:

``` bash
scontrol show config | grep -iE "MaxMemPerCPU|DefMemPerCPU"
sinfo -o "%P %c %m" | sort -u
```

Where a step needs more memory than that ratio allows, **raise `cpus` rather than lowering memory**.
Requesting cores a single-threaded tool will not use is the price of getting the memory at all, and
the cores are unusable by other jobs anyway once the node's memory is committed.

This belongs in your site config, not in the pipeline — it describes a SLURM installation rather
than pulposeq.

### Reading the resource reports

Every run writes `pipeline_info/execution_trace_*.txt`, which is more useful than the HTML report for
this purpose: it carries `rchar` and `wchar` alongside memory, so you can relate usage to input size.

A few traps worth knowing:

- **Size from `peak_rss`, never `peak_vmem`.** `RENDER_REPORT` reports over 1 TB of virtual memory
  against about 1 GB resident — address space that R and Quarto reserve without ever touching.
- **`peak_rss` is unreliable for piped commands.** `CHOPPER` runs `zcat | chopper | gzip` and has
  reported `0.00 GB` while genuinely using ~50 GB, because short-lived children escape Nextflow's
  sampler. For those, ask SLURM instead:

  ``` bash
  sacct -j <jobid> --format=JobID,JobName,MaxRSS,ReqMem,State,ExitCode
  ```

- **Nextflow reads high relative to SLURM**, by around 12% on the same task. Its figures are already
  conservative, so there is no need to add a large safety factor on top.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.


### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time. Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory. We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

``` bash
NXF_OPTS='-Xms1g -Xmx16g'
```
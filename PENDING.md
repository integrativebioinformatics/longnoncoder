# pulposeq — pending fixes & pending validations

Carried over from the `dev`-branch review. Tracks what has been applied on the `updating` branch and
what still needs verifying before it merges to `dev`.

**Applied on `updating`:** M4, M5, P1–P6, POST_REFINEMENT, the biomaRt removal and the GTF
enrichment.

## Run status — 2026-08-07

A full run completed on the chr1 test dataset with `-profile test,singularity` and
`library: PacBio` / `stranded_library: true`. **All 24 tasks reached COMPLETED or CACHED; none
failed.** That settles the structural questions:

- the pipeline parses and runs under the strict syntax on Nextflow 26.x (V1)
- `POST_REFINEMENT` and `ENRICH_VALIDATED_GTF` are wired correctly and produce their outputs
  (V11, V12.4)
- `KNOWN_TRANSCRIPTS` completes with no biomaRt query, so the reference-GTF path works end to end
  on an Ensembl annotation (V12.2)
- the `:test3` container serves every local module including `BAMBU_VALIDATE` (V10)
- publishing still works after the `outputDir` removal (V9)

> [!WARNING]
> Completing is not the same as being correct. Everything below that concerns **output content**
> remains unverified — most importantly the metadata parity check (V12.1), the GTF attribute
> contents (V12.3), the group labels after the collectFile change (V5), and the report panels
> (V11.3). A GENCODE run (V12.7) has not been attempted at all.

`conf/test.config` has since been retuned against this run's measured usage; see Part 2.

---

## Part 1 — Validation plan

Validation is done by **running the pipeline end to end** on the test data and reviewing the
artifacts. A single real run exercises most of the changed surface at once. Three cheap gates run
first, because each catches a whole class of failure in seconds rather than after hours of compute.

### Environment

```bash
conda activate nf-core          # provides Nextflow >= 26.04
module load singularity
cd test_data && chmod +x download-ref-fastq.sh && ./download-ref-fastq.sh && cd ..
# then fill the /full/path/to/... placeholders in test_data/samplesheet.csv and test_data/testing.yml
```

> [!IMPORTANT]
> This branch sets `manifest.nextflowVersion = '!>=26.04.0'`. From 26.04 the **strict syntax parser
> is enabled by default**, which is stricter about config files than any version this pipeline has
> run under before. Stage 0 exists specifically to catch that.

---

### Stage 0 — `nextflow lint` (seconds)

```bash
nextflow lint .
```

Statically checks every `.nf` and `.config` file against the strict syntax. This is the highest-value
single command on the branch, because 26.04 enforces rules the pipeline has never been parsed under.

**Known constraint already handled:** the strict parser allows only config assignments, blocks and
includes — top-level `def` declarations and helper functions are rejected. The library/strandedness
logic was originally written as two shared top-level closures in `conf/modules.config`; it is now
inlined inside the `MINIMAP2_ALIGN` and `BAMBU` `ext.args` closures, where local variables are still
permitted. The rule is deliberately stated twice; do not "refactor" it back into a shared helper.

**Most likely remaining offender:** `nextflow.config` sets

```groovy
trace_report_suffix = new java.util.Date().format( 'yyyy-MM-dd_HH-mm-ss')
```

This is nf-core template code written for 25.10 and involves a constructor call in a config
assignment. If lint rejects it, replace with a literal or move the timestamp generation elsewhere.

Anything lint flags in `modules/nf-core/*` is upstream template code, not part of this branch —
note it, don't fix it here.

### Stage 1 — config resolution (seconds)

```bash
nextflow config -profile test,singularity -params-file test_data/testing.yml
```

Confirms the config tree assembles and the `ext.args` closures are syntactically reachable. Check the
resolved `process` block lists `resourceLimits` and that no `containers_*.config` is referenced
(those eight files were deleted).

### Stage 2 — stub run (~minutes)

```bash
nextflow run main.nf -profile test,singularity -params-file test_data/testing.yml -stub-run
```

Exercises every channel connection, input cardinality and output filename in the DAG without running
a single real tool. Catches wiring errors in the new POST_REFINEMENT module and the M4 `collectFile`
change for the price of a coffee, rather than after Bambu has run.

Note the BAMBU stub previously touched `bambu_output.rds`, a filename `bin/bambu.R` never writes; it
now touches `se_multiSample.rds` and `seGene_multiSample.rds` to satisfy the new named outputs.

### Stage 3 — the real run (hours)

```bash
nextflow run main.nf -profile test,singularity -params-file test_data/testing.yml
```

> [!WARNING]
> `-profile test` caps resources at `cpus: 2, memory: 12.GB` and gives `process_high_memory` 12 GB.
> That is the Bambu step, running on real ENCODE reads against human chr1. If it is OOM-killed, this
> is a test-profile sizing problem, not a defect in the branch — rerun that step with `-profile
> medium,singularity` and note it.

### Stage 4 — send for review

The full `work/` directory is large and mostly BAMs. This bundle is what is actually needed:

```bash
tar czf pulposeq-validation.tar.gz \
    .nextflow.log \
    $(find work -name '.command.sh' -o -name '.command.err' -o -name '.command.out') \
    results/pipeline_info \
    results/transcriptome_report \
    results/bambu_validated \
    results/novel_transcripts \
    results/annotated_transcripts \
    results/multiqc \
    results/bambu/*.png results/bambu/*.rds \
    results/post_refinement 2>/dev/null
```

Deliberately excluded: `work/**/*.bam`, `results/minimap2/`, `results/chopper/` — large binaries with
nothing to review in them. If a task fails, also include that task's full work directory.

---

### What the run covers

| ID | Check | Covered by |
|---|---|---|
| V1 | Strict-syntax / config resolution | Stage 0–1 |
| V3 | minimap2 preset emitted — `testing.yml` is set to `PacBio`, so expect `-ax splice:hq -uf` | Stage 3 → `grep -h "minimap2" work/*/*/.command.sh` |
| V4 | Bambu gets `--ndr` and `--stranded true`; optparse accepts the new flag | Stage 3 → `grep -h "bambu.R" work/*/*/.command.sh` |
| V5 | M4 regression — group labels correct despite changed `bamlist.txt` order | Stage 3 → `colData(readRDS(...))$group` |
| V7 | `maxlen` unset → **no** `--maxlength` in the chopper command | Stage 3 → `grep -h "chopper" work/*/*/.command.sh` |
| V9 | Publishing still works after the `outputDir` removal | Stage 3 → results/ subdirectories present |
| V10 | `BAMBU_VALIDATE` runs on `:test3` | Stage 3 → task completes |
| V11 | POST_REFINEMENT end to end | Stage 2 + 3 |

### What the run does *not* cover

One run uses one library setting and one set of valid params, so these need deliberate extra
invocations. All are cheap — none needs to run to completion.

| ID | Check | How |
|---|---|---|
| V2 | Param validation errors | Run with no `--library`, then `--library nonsense`, then `--library ONT_DRS --stranded_library false`. Each should fail or warn at initialisation, in seconds. |
| V3b | The other three presets | `-stub-run` with `--library ONT_DRS`, `--library PacBio`, `--library ONT_cDNA --stranded_library true`, then read `work/*/*/.command.sh`. Stub does not interpolate minimap2's `args`, so instead read the resolved value from `nextflow config` output, or do three short real runs killed after MINIMAP2_ALIGN starts. |
| V7b | `maxlen` set | Add `--maxlen 10000` and confirm `--maxlength 10000` appears. |
| V8 | `resourceLimits` caps retries | Force a `process_high_memory` retry under `-profile small` and confirm the second attempt requests 200 GB, not 400 GB. |
| V6 | Library assumption for the test data | `testing.yml` is set to `ONT_cDNA` / `stranded_library: false`, a conservative guess from the ENCODE files being described as CapTrap cDNA. Correct it if the real orientation is known. |

### Expected values

| `library` | `stranded_library` | minimap2 preset | Bambu |
|---|---|---|---|
| `ONT_DRS` | forced `true` | `-ax splice -uf -k14` | `--stranded true` |
| `PacBio` | forced `true` | `-ax splice:hq -uf` | `--stranded true` |
| `ONT_cDNA` | `true` | `-ax splice -uf` | `--stranded true` |
| `ONT_cDNA` | `false` | `-ax splice` | `--stranded false` |

### V11 — POST_REFINEMENT, in order of what is most likely to break

1. **RDS row-name matching.** `post_refinement.R` subsets by `rownames(object) %in% ids`, where `ids`
   is column 1 of the validated counts files. If one side carries version suffixes and the other does
   not, the intersect is empty and the script aborts with "None of the N validated IDs matched" —
   that message is deliberate, it reports disagreeing keys rather than silently plotting nothing.
   ```r
   se <- readRDS("results/bambu/se_multiSample.rds"); head(rownames(se))
   head(read.table("results/bambu_validated/BambuOutput_counts_transcript_validated.txt", header = TRUE)[[1]])
   ```
2. **`plotBambu` on a row-subset object** — confirm heatmap and PCA still render, and `colData$groupVar`
   survived the subset (it should; only rows are touched).
3. **Both plot sets reach the report** — rendered HTML should show a "Sample-level Expression
   Overview" section with four two-up comparisons. Blank panels point at `show_plots()`.
4. **No filename collisions in the RENDER_REPORT work dir** — raw plots stage as `pca.png`,
   `pca_grouped.png`, `heatmap_gene.png`, `heatmap_transcript.png`; post-refinement ones as the same
   names with a `_validated` suffix. All eight should be present, none overwritten.

### V12 — biomaRt removal and GTF enrichment

1. **Metadata parity — the decisive check.** The new `annotated_transcriptome_metadata.csv` should
   match the biomaRt-era one from the last successful run: same rows, same columns, same values.
   ```r
   old <- read.csv("<previous>/annotated_transcripts/annotated_transcriptome_metadata.csv")
   new <- read.csv("<new>/annotated_transcripts/annotated_transcriptome_metadata.csv")
   all.equal(old[order(old$ensembl_transcript_id), ], new[order(new$ensembl_transcript_id), ])
   ```
   `transcript_length` is the field most likely to diverge — biomaRt returns the **mature** length
   (sum of exon widths), not `end - start`. A systematic difference there means the exon summing in
   `read_reference_gtf()` is wrong.
2. **No network access.** `KNOWN_TRANSCRIPTS` must complete with no outbound calls; `.command.log`
   should no longer contain "Connecting to Ensembl biomaRt...".
3. **GTF attributes present:**
   ```bash
   grep -m3 "transcript_status" results/annotated_transcripts/bambu_annotated_transcriptome.gtf
   grep -m3 "class_code" results/novel_transcripts/novel_transcripts_validated.gtf
   grep -c 'gene_biotype "novel"' results/bambu_validated/BambuOutput_annotations_validated.gtf
   ```
   The last count should approximate the number of novel transcripts without a `cmp_ref_gene`
   (630 of 1498 in the chr1 test data).
4. **`ENRICH_VALIDATED_GTF` does not clobber its inputs.** The three validated GTFs are staged into
   `input/` and rewritten at the task root. Confirm the published files carry the new attributes and
   that record counts match the `SUBSET_BAMBU_GTF` outputs.
5. **Publishing moved.** `SUBSET_BAMBU_GTF` no longer publishes — its GTFs are intermediates and
   `ENRICH_VALIDATED_GTF` publishes the same three filenames to `bambu_validated/`. Confirm exactly
   one copy of each exists and it is the enriched one.
6. **Memory.** `KNOWN_TRANSCRIPTS` and `NOVEL_TRANSCRIPTS` moved from `process_single` to
   `process_medium` because they now parse the reference GTF in R. Watch actual usage; a
   whole-genome GENCODE annotation is much heavier than the chr1 test file.
7. **GENCODE run.** The point of the change is that both sources work. Run once with a GENCODE GTF
   and confirm biotypes populate via `gene_type` and `chromosome_name` reads `1`, not `chr1`. Use the
   CHR or PRI build — ALL will trip the new duplicate-identifier check by design.
8. **Params removed.** Any params file still setting `ensembl_organism_dataset` or `ensembl_version`
   now fails schema validation. `test_data/testing.yml` and `examplerun.yml` are updated in-repo.

---
## Part 2 — Applied on the `updating` branch

| Item | Change |
|---|---|
| **M4** | `transcript_reconstruction.nf` — both `collectFile` chains collapsed to single calls with `sort: true` |
| **M5** | `library` / `stranded_library` params driving the minimap2 preset and Bambu `stranded`; new `--stranded` option in `bin/bambu.R`; validation, schema, and docs |
| **P1** | `validate_counts` container `:test` → `:test3` (all 7 local modules now on one tag) |
| **P2** | Missing `\` before `$args` fixed in `subset_counts` and `validate_counts` |
| **P3** | `maxlen` wired into `CHOPPER.ext.args2` as a conditional `--maxlength`; schema description, docs row, `examplerun.yml` entry |
| **P4** | `process.resourceLimits` added to `base`, `light`, `medium`, `large`, `test` configs |
| **P5** | Inert `outputDir` / `workflow.output.mode` lines removed — committing to the `publishDir` model |
| **P6** | 8 orphaned `conf/containers_*.config` files deleted |
| **POST_REFINEMENT** | New module + `bin/post_refinement.R` regenerating the Bambu PCA/heatmaps from the validated transcriptome, between the metadata refinement steps and the report; both raw and post-refinement plots now embedded in the report |
| **Nextflow 26** | `manifest.nextflowVersion` raised to `!>=26.04.0`; README badge and strict-syntax warning updated |
| **Strict syntax** | The shared `stranded_library` / `minimap2_preset` closures in `conf/modules.config` were inlined into the `MINIMAP2_ALIGN` and `BAMBU` `ext.args` closures — the strict parser (default from 26.04) rejects top-level `def` declarations in config files, and rejects `switch`/`return` statements inside config closures |
| **biomaRt removal** | `known_transcripts.R` now reads all transcript metadata from the supplied reference annotation GTF via the new `bin/gtf_annotation_utils.R`; `ensembl_organism_dataset` and `ensembl_version` deleted; known/novel split now tested by reference membership rather than an `ENS` prefix |
| **GTF enrichment** | All custom-generated GTFs carry `transcript_status`, `gene_biotype`, `transcript_biotype`, `gene_name`, and — for novel transcripts — `class_code`, `classification` and `ref_gene_id`. New `bin/enrich_validated_gtf.R` + `ENRICH_VALIDATED_GTF` module post-processes the `subset_bambu_gtf.sh` outputs, leaving that shell script untouched |
| **test.config** | Retuned against the 2026-08-07 run's measured usage. `process_medium` 8→4 cpus and 15→12 GB, `process_high` 25→15 GB, `process_high_memory` 4→2 cpus and 40→24 GB; time limits cut except `process_high`, whose slowest task ran 84 min and now gets 3 h. `process_single` deliberately left at 5 GB — RENDER_REPORT at 68% is the tightest fit in the run. Header records the observed peaks and the reasoning |

**Deliberately left alone under P6:** the `conda.enabled = false` guards in the docker/singularity/
apptainer profiles (they actively enforce the no-conda policy), and
`modules/nf-core/multiqc/{environment.yml,.conda-lock}` — those are managed by nf-core tooling via
`modules.json`, and deleting them would break `nf-core modules update`.

**Note on P5:** `publishDir` is **not** deprecated in Nextflow 26.04 — the migration guide says
nothing about removing it, and workflow outputs merely came out of preview in 25.10. Both models are
current. The full output-DSL migration is deferred to D2 below.

---

## Part 3 — Pending, needs discussion or help

### D1. CI and nf-test — *you asked for help here*

No `.github/workflows/` directory exists. `.nf-core.yml` disables `nf_test_content` lint. No tests
for the eight local modules or the pipeline. Every local module already has a working `stub` block,
so a `-stub` CI job is cheap and would catch V1/V2 automatically. Suggested order: stub-run workflow
→ linting workflow → per-module nf-tests.

### D2. Migrate to the workflow output definition

Deferred from P5 until D1 lands, so there is something to catch regressions. Scope: every published
channel must reach the entry workflow's `publish:` section, but the subworkflows don't currently emit
most of what gets published — `QC_FILT` emits 3 channels while NanoComp alone has 11 being published;
`CLASSIFICATION_POTENTIAL_CODING` emits 5 while GffCompare and RNAmining produce 9 between them.
Roughly 50 channels to plumb through 4 subworkflows, each declared twice (`publish:` + `output {}`),
plus stripping all `publishDir` from `conf/modules.config`. Several hundred lines.

### D3. Granular step control (was H7)

You want `run only QC`, `QC+mapping`, `QC+mapping+assembly`, etc. N independent `skip_*` booleans
multiply into invalid combinations fast. Worth discussing a `--step` / `--from` / `--to` parameter
(sarek-style) instead of adding more flags.

### D4. `val` → `path` staging (was H1, downgraded)

`BAMBU` (`val bam_list`, `val sample_info`) and `RNAMINING` (`val fasta`) take files as `val`, so they
are never staged and never declared as dependencies. **Not a correctness bug** for Docker/Singularity
on a shared filesystem — ordering and `-resume` both hold. Two residual risks: `cleanup = true` /
`nextflow clean` may delete upstream work dirs that `bamlist.txt` still points at, and `-with-dag`
won't draw the BAM → BAMBU edges. Worth doing alongside other refactoring, not on its own.

### D5. Container reproducibility beyond the tag

P1 unified on `:test3`, but that is still a mutable tag in a personal Docker Hub namespace. For a
v1.0.0 carrying a DOI this wants a versioned tag or a digest pin.

---

## Part 4 — Low priority / polish

1. Strip change-log comments from source: `// <- ADDED SCRIPT PATH` (×3 in `workflows/pulposeq.nf`),
   `// FIXED: Put 104 in brackets`, `// REMOVED the unused 'def args' line`,
   `// Removed the broken CHOPPER.out.versions line entirely`, and the duplicated
   `// Running quality check in filtered reads` in `subworkflows/local/qc.nf`.
2. Dead emits: `ALIGNMENT.out.index` is always empty; `TRANSCRIPT_RECONSTRUCTION` still emits
   `bamlist`, `samp_info`, `reference`, `annotation` with no consumer. (`pca`, `pca_grouped`,
   `h_gene`, `h_transcript` are now consumed by the report via POST_REFINEMENT.) The PNGs still
   never reach MultiQC.
3. Version-collection gaps (partly intentional): `RNAMINING.out.versions` is dropped because
   `CLASSIFICATION_POTENTIAL_CODING`'s `ch_versions` is never mixed; `SUBSET_BAMBU_COUNTS.out.versions`
   and `BAMBU_VALIDATE.out.versions` are never mixed in `workflows/pulposeq.nf`.
4. `.ifEmpty(null)` on version channels in `alignment.nf` and `transcript_reconstruction.nf` — emits a
   literal `null`; deprecated pattern.
5. Reference channel in `alignment.nf` passes a String as `meta2` and a String as the path; works only
   by coercion. Replace with `channel.value([[id: file(params.reference).baseName], file(params.reference)])`.
6. Inconsistent include path — `classification_codingpotential.nf` is included with its `.nf`
   extension, the other three without.
7. Inconsistent params plumbing — `ALIGNMENT` reads `params.reference` directly; other subworkflows
   receive it via `take:`. Prefer `take:` throughout for testability.
8. Local modules use bare `path` inputs with no meta map — hard-codes one sample set per run.
9. Inconsistent script staging — R scripts are formal `path` inputs; `.sh`/`.awk` helpers rely on
   `bin/` on `PATH`, and `subset_bambu_gtf` uses `$(which subset_gtf.awk)`.
10. awk version parsing differs between sibling modules using the same container (`mawk` vs `GNU Awk`).
11. Release metadata still template state — `CHANGELOG.md` is the unfilled nf-core stub
    (`## v1.0.0 - [date]`), `manifest.doi` empty, `defaultBranch = 'main'`.
12. `nextflowVersion = '!>=25.10.4'` is a hard floor that will exclude users on institutional clusters.
13. `ext.args = {"--plot violin"}` closures in `conf/modules.config` where plain strings suffice.
14. README tells users to run `-profile test`, but `conf/test.config` carries only resource labels —
    test data must be downloaded and paths filled in first (working as designed; wording could say so).

---

## Resolved — no action

- **H2** resource profiles — `medium.config` and `large.config` set `process_high_memory { cpus = 6 }`;
  profile configs override the `base.config` baseline. Working as designed. (Minor: `light.config`
  sets only `memory`, so it inherits `cpus = 1` from base while its `process_high` gets 10.)
- **H3** test data is user-downloaded by design.
- **M1** version omissions are partly deliberate — RNAmining's version is hardcoded, and modules
  sharing a container were deduplicated. Item 3 in Part 4 lists what is dropped, if you want it back.

# pulposeq — validation status of the `updating` branch

Tracks what has been verified on this branch and what has not, ahead of merging to `dev`.

**What the branch contains:** the `library` / `stranded_library` parameters driving the minimap2
preset and Bambu strandedness; removal of the biomaRt query in favour of reading transcript metadata
from the supplied reference annotation; biotype and classification attributes written into every
generated GTF; the new `POST_REFINEMENT` and `ENRICH_VALIDATED_GTF` modules; a rewritten resource
profile scheme; the Nextflow 26.04 / strict-syntax move; and a set of report fixes. Commit history
has the detail.

---

# 1. Validations done

## Static checks

| Check | Result |
|---|---|
| **`nextflow lint .`** | **PASSED** — 42 files, 0 errors. The 3 warnings are all in nf-core template subworkflows tracked in `modules.json`; they are style-only and deliberately not fixed, since `nf-core subworkflows update` would overwrite any local edit. |
| **`nextflow config -profile test,apptainer`** | **PASSED** — config tree assembles on Apptainer/SLURM. The `MINIMAP2_ALIGN` and `BAMBU` `ext.args` ternaries survive the strict parser with the config loaded, `SUBSET_BAMBU_GTF` resolves to `publishDir { enabled = false }`, and no deleted `containers_*.config` is referenced. |
| **Label coverage** | **PASSED** — every process resolves a label in every profile, so nothing falls back silently to `base.config` defaults. Every profile's `resourceLimits` is at or above its largest request, so no allocation is clamped. |

> `nextflow config` accepts **no parameters** — its options are `-flat`, `-h`, `-profile`,
> `-properties`, `-r`, `-a`, `-sort`, `-value`. Passing `-params-file` fails with "Unknown option".
> It also cannot check the minimap2 presets: params cannot be injected, and `ext.args` is a closure
> printed unevaluated.

## Stub run

**PASSED** — 24 tasks, ~4 minutes, `-profile test,apptainer`. Confirms every channel connection,
input cardinality and output filename, including the `ENRICH_VALIDATED_GTF` and `POST_REFINEMENT`
wiring, the extra `path` inputs on `KNOWN_TRANSCRIPTS` / `NOVEL_TRANSCRIPTS`, and the corrected BAMBU
stub, which previously touched a `bambu_output.rds` that `bin/bambu.R` never writes.

> `-stub-run` is a boolean switch and takes no value. `-stub-run true` leaves `true` as a positional
> argument and warns. The same warning appears for `-profile test, apptainer` with a space — that one
> genuinely breaks, since only `test` reaches `-profile` and no container runtime is enabled.

## Real runs

Three completed: chr1 test data, a GENCODE run, and a full Ensembl 116 run. Between them:

| Check | Result |
|---|---|
| Pipeline completes end to end | **PASSED** — all tasks reach COMPLETED, output files generated as expected |
| Strict syntax at runtime | **PASSED** — parses and runs under Nextflow 26.x |
| `KNOWN_TRANSCRIPTS` without biomaRt | **PASSED** — completes with no network access |
| `POST_REFINEMENT` / `ENRICH_VALIDATED_GTF` | **PASSED** — run and produce their outputs |
| `:test3` container | **PASSED** — serves every local module including `BAMBU_VALIDATE` |
| Publishing after the `outputDir` removal | **PASSED** — results land in the expected subdirectories |
| minimap2 preset for `PacBio` | **PASSED** — `-ax splice:hq -uf` |
| minimap2 preset for unstranded `ONT_cDNA` | **PASSED** — `-ax splice`, confirmed on the 8-sample ONT run |
| **GENCODE annotation** | **PASSED** — biotypes populate via `gene_type`. Surfaced one real defect, since fixed: `chromosome_name` was being stripped of its `chr` prefix on the annotated side but not the novel side. Sequence names are now passed through untouched on both. |

## Resources

Sized from measured runs and verified in execution:

| Check | Result |
|---|---|
| `test` profile | **PASSED** — 24/24 tasks, every label between 40% and 96% utilisation |
| `large` profile | **PASSED** through CHOPPER, MINIMAP2 and NANOCOMP on 8 ONT samples of 25–30 GB |
| `--mem-per-cpu` | **PASSED** — `executor.perCpuMemAllocation` accepted at up to 12 GB/CPU |
| Label coverage | **PASSED** — every process resolves a label in every profile; no silent fallback to base defaults |

The measured figures, the scaling rules derived from them and the reporting caveats are documented
in `docs/usage.md` under *Resource requests*.

## Report

Verified against the rendered Ensembl 116 HTML:

| Check | Result |
|---|---|
| Duplicated captions on the before/after figures | **FIXED** — each caption appears once, with `(a)` / `(b)` subcaptions |
| Colliding figure numbers | **FIXED** — 25 captions numbered 1–25, no duplicates, no hand-written `**Figure N:**` left |
| Stale cross-reference | **FIXED** — `@fig-density` resolves to Figure 14; the old hardcoded "Figure 13" pointed at the wrong figure |
| Chromosome naming | **FIXED** — annotated and novel figures now agree, both using the annotation's own naming |
| Chromosome figure scaling | **FIXED** — layout derived from the number of chromosomes rather than a fixed 35×30 in |
| Tabsets, palette, legends, lollipops, log ticks | **APPLIED** — visible in the render |
| File size | 17 MB for a whole genome at 150 dpi, against 20 MB for chr1 alone at 300 dpi |

---

# 2. Validations still pending

## Priorities

Everything here should be settled before merging to `dev`.

### P1. Metadata parity — **in progress**

The decisive check for the biomaRt removal: the GTF-derived
`annotated_transcriptome_metadata.csv` must match the biomaRt-era output.

```r
old <- read.csv("<previous>/annotated_transcripts/annotated_transcriptome_metadata.csv")
new <- read.csv("<new>/annotated_transcripts/annotated_transcriptome_metadata.csv")
all.equal(old[order(old$ensembl_transcript_id), ], new[order(new$ensembl_transcript_id), ])
```

`transcript_length` is the field most likely to diverge — biomaRt reports the **mature** length, the
sum of exon widths, not `end - start`. A systematic difference there means the exon summing in
`read_reference_gtf()` is wrong.

> **The annotation release must match the baseline.** Ensembl changes content between releases: new
> transcripts, bumped versions, reassigned biotypes. A 116 GTF compared against a baseline biomaRt
> queried at 114 will differ for reasons that have nothing to do with the parsing. Either regenerate
> the baseline at the same release, or restrict the comparison to transcripts present in both.
>
> `chromosome_name` will also differ against a **GENCODE** baseline by design, since sequence names
> are no longer rewritten. Compare against Ensembl, or exclude that column.

### P2. Parameter validation — not run, needs no compute

Four deliberately bad invocations, each failing or warning at initialisation in seconds:

| Command | Expected |
|---|---|
| no `--library` | error: `--library must be provided when alignment is not skipped` |
| `--library nonsense` | rejected by the schema enum **and** `validateInputParameters()` |
| `--library ONT_DRS --stranded_library false` | warning that ONT_DRS is always stranded; run proceeds |
| `--skip_alignment`, no `--library` | no error — library is only required when alignment runs |

### P3. Remaining minimap2 presets

`PacBio` and unstranded `ONT_cDNA` are both confirmed on real runs. Still outstanding, each differing
from a confirmed case by a single flag:

| `library` | `stranded_library` | Expected preset | Bambu |
|---|---|---|---|
| `ONT_DRS` | forced `true` | `-ax splice -uf -k14` | `--stranded true` |
| `ONT_cDNA` | `true` | `-ax splice -uf` | `--stranded true` |

Read the resolved value with `grep -h "minimap2" work/*/*/.command.sh` after a real run — a stub run
does not interpolate minimap2's `args`.

### P4. Sample group labels after the `collectFile` change

The one outstanding check with real consequences if it fails. Row order in `bamlist.txt` changed
from group-ordered to path-ordered; if the group mapping broke, every group-wise figure in the report
is silently wrong and nothing else flags it.

```r
se <- readRDS("<outdir>/bambu/se_multiSample.rds")
colData(se)$group   # must match each sample's group in the samplesheet
```

`bin/bambu.R` matches sample rows to BAMs by basename rather than by position, which is why the
change should be safe — but it has not been confirmed on a multi-group dataset.

### P5. GTF attribute contents

```bash
grep -m3 "transcript_status" results/annotated_transcripts/bambu_annotated_transcriptome.gtf
grep -m3 "class_code" results/novel_transcripts/novel_transcripts_validated.gtf
grep -c 'gene_biotype "novel"' results/bambu_validated/BambuOutput_annotations_validated.gtf
```

The last count should approximate the number of novel transcripts without a `cmp_ref_gene` — 630 of
1498 in the chr1 test data. Also confirm `ENRICH_VALIDATED_GTF` did not clobber its inputs: the three
validated GTFs are staged into `input/` and rewritten at the task root, and exactly one copy of each
should exist in `bambu_validated/`, carrying the new attributes.

### P6. POST_REFINEMENT at large scale — unmeasured

The only process with no measurement at full scale. It aborted at 60 GB on the 8-sample ONT run, so
that value is known insufficient; it loads the Bambu `SummarizedExperiment`, which for that run was
built from a 640 GB assembly. Size it generously and capture the peak on the next complete large run.

### P7. Bambu sample-count ceiling — document, do not fix

Bambu grows by roughly 54 GB per additional sample (532 GB at 6, 640 GB at 8) because it holds every
sample in a single object. On a 768 GB node that caps a run at **about 10 samples**. This is a
property of the assembly step rather than a defect, and it is recorded in `docs/usage.md`. Worth
revisiting only if splitting a run and merging results becomes a requirement.

### P8. Release decisions

Not validations, but they gate a v1.0.0 with a DOI.

- **Container reproducibility.** `:test3` is a mutable tag in a personal Docker Hub namespace. If it
  is rebuilt, this branch stops being reproducible retroactively and nothing recovers it. A digest
  pin is a small change now that the image is known to work.
- **`publish_dir_mode`.** Now `symlink`, which is good for iterating on large data but produces an
  output directory that cannot be moved, archived or shared — deleting `work/` destroys it. `copy`
  is the safer default for a release, with `symlink` as a documented opt-in. Documented either way.
- **Release metadata.** `CHANGELOG.md` is still the unfilled nf-core stub (`## v1.0.0 - [date]`),
  `manifest.doi` is empty, and `defaultBranch` is `main` while work happens on `updating`.

---

## Non-urgent polishing

Nothing here affects correctness of a run.

### CI and testing

**No `.github/workflows/` exists**, `.nf-core.yml` disables `nf_test_content` lint, and there are no
tests for the local modules or the pipeline. Every local module already has a working `stub`, so a
`-stub` CI job is cheap and would catch the config and wiring classes of failure automatically.
Suggested order: stub-run workflow → linting workflow → per-module nf-tests.

This is the highest-value item in this section: it is what stops the next round of changes from
being unverified in the same way.

### Larger refactors

- **Workflow output definition.** `publishDir` is not deprecated in 26.04, so this is direction
  rather than necessity. Scope is large: every published channel must reach the entry workflow's
  `publish:` section, but the subworkflows do not currently emit most of what gets published —
  roughly 50 channels through 4 subworkflows, each declared twice, plus stripping all `publishDir`
  from `conf/modules.config`. Worth doing only after CI exists.
- **Granular step control.** "Run only QC", "QC+mapping", and so on. Independent `skip_*` booleans
  multiply into invalid combinations quickly; a `--step` / `--from` / `--to` parameter (sarek-style)
  scales better than more flags.
- **`val` → `path` staging.** `BAMBU` (`val bam_list`, `val sample_info`) and `RNAMINING`
  (`val fasta`) take files as `val`, so they are never staged and never declared as dependencies.
  **Not a correctness bug** for Docker/Singularity on a shared filesystem — ordering and `-resume`
  both hold. Two residual risks: `cleanup = true` or `nextflow clean` may remove work directories
  that `bamlist.txt` still points at, and `-with-dag` will not draw the BAM → BAMBU edges.

### Code hygiene

1. Change-log comments left in source: `// <- ADDED SCRIPT PATH`, `// FIXED: Put 104 in brackets`,
   `// REMOVED the unused 'def args' line`, `// Removed the broken CHOPPER.out.versions line
   entirely`, and a duplicated `// Running quality check in filtered reads` in `qc.nf`.
2. Dead emits: `ALIGNMENT.out.index` is always empty; `TRANSCRIPT_RECONSTRUCTION` emits `bamlist`,
   `samp_info`, `reference` and `annotation` with no consumer. Bambu's PNGs never reach MultiQC.
3. Version-collection gaps, partly deliberate: `RNAMINING.out.versions` is dropped because
   `CLASSIFICATION`'s `ch_versions` is never mixed; `SUBSET_BAMBU_COUNTS.out.versions` and
   `BAMBU_VALIDATE.out.versions` are never mixed in `workflows/pulposeq.nf`.
4. `.ifEmpty(null)` on version channels in `alignment.nf` and `transcript_reconstruction.nf` emits a
   literal `null` — a deprecated pattern.
5. The reference channel in `alignment.nf` passes a String as `meta2` and a String as the path, which
   works only by coercion. Prefer
   `channel.value([[id: file(params.reference).baseName], file(params.reference)])`.
6. Inconsistent params plumbing — `ALIGNMENT` reads `params.reference` directly while other
   subworkflows receive it via `take:`. Prefer `take:` throughout, for testability.
7. Local modules use bare `path` inputs with no meta map, which hard-codes one sample set per run.
8. Inconsistent script staging — R scripts are formal `path` inputs, while `.sh` / `.awk` helpers
   rely on `bin/` being on `PATH`, and `subset_bambu_gtf` uses `$(which subset_gtf.awk)`.
9. awk version parsing differs between sibling modules using the same container (`mawk` vs `GNU Awk`).
10. `ext.args = {"--plot violin"}` closures in `conf/modules.config` where plain strings suffice.
11. The README points users at `-profile test`, but the test data must be downloaded and paths filled
    in first. Working as designed; the wording could say so.

### Report — reviewed and deliberately not changed

- **Count bars on a log axis.** Now drawn as lollipops so position rather than bar length carries the
  value, which a log axis represents correctly. The remaining question is whether a linear axis with
  free-scale facets would communicate better — a presentation decision, not a defect.
- **`fig-format: svg`.** Would be sharper and scalable, but `showtext` draws glyphs as polygons, so
  every label and legend entry becomes vector paths. On figures with this much text that may grow the
  file rather than shrink it. `dpi: 150` was taken instead. Worth measuring if size matters more.
- **Interactivity.** Tabsets are in. Anything further — sortable tables, tooltips, linked filtering —
  needs R packages that are almost certainly not in `longnoncoder:test3`, so it is gated behind a
  container rebuild and the reproducibility decision above.

# pulposeq — pending fixes & pending validations

Carried over from the `dev`-branch review. Tracks what has been applied on the `updating` branch and
what still needs verifying before it merges to `dev`.

**Applied on `updating`:** M4, M5, P1–P6, and the new POST_REFINEMENT module.

> [!WARNING]
> None of the applied changes have been executed. They were written in an environment with no
> Nextflow, no Java and no R, so every item in Part 1 below is unverified.

---

## Part 1 — Pending validations

None of these ran on the Windows machine where the changes were written (no Java, no Nextflow, no
usable WSL distro). **Everything below is unverified.** Run from the repo root on a machine with
Nextflow and Docker/Singularity.

### V1. Config parses and closures resolve — **highest priority**

`conf/modules.config` defines two top-level `def` closures (`stranded_library`, `minimap2_preset`)
called from inside `ext.args` closures. This relies on Groovy lexical capture surviving Nextflow's
config parser. It should hold, but it is the one mechanism in the whole change never exercised.

```bash
nextflow config -profile docker
nextflow config -profile docker --library ONT_DRS
nextflow config -profile docker --library PacBio
nextflow config -profile docker --library ONT_cDNA --stranded_library true
```

Expect no `MissingMethodException` / `No such property: minimap2_preset`.

**If this fails:** inline the `switch` into `MINIMAP2_ALIGN.ext.args` and the strandedness ternary
into `BAMBU.ext.args`, accepting that the rule is then stated twice.

### V2. Parameter validation fires correctly

| Command | Expected |
|---|---|
| no `--library` | error: `--library must be provided when alignment is not skipped` |
| `--library nonsense` | rejected by schema enum **and** `validateInputParameters()` |
| `--library ONT_DRS --stranded_library false` | warning `...is always stranded...`, run proceeds |
| `--skip_alignment`, no `--library` | no error (library only required when alignment runs) |

### V3. minimap2 presets actually emitted

Needs a **real, non-stub** run — minimap2's `stub` block does not interpolate `args`.

```bash
nextflow run . -profile docker -params-file test_data/testing.yml --library ONT_DRS
grep -h "minimap2" work/*/*/.command.sh
```

| Setting | Expected preset |
|---|---|
| `ONT_DRS` | `-ax splice -uf -k14` |
| `PacBio` | `-ax splice:hq -uf` |
| `ONT_cDNA` + `stranded_library: true` | `-ax splice -uf` |
| `ONT_cDNA` + `stranded_library: false` | `-ax splice` |

### V4. Bambu receives strandedness and the new flag parses

```bash
grep -h "bambu.R" work/*/*/.command.sh
```

Expect `--ndr <value> --stranded true|false`. Confirm the run reaches `BambuOutput_*` without an
optparse error on `--stranded`, and that `--ndr null` still behaves as before — the module no longer
supplies its own default, `conf/modules.config` does.

### V5. M4 regression — group labels still correct

The specific thing the `collectFile` simplification could have broken. Row order in `bamlist.txt`
changed from group-ordered to path-ordered.

```bash
wc -l work/*/*/bamlist.txt work/*/*/sampinfo_samplesheet.tsv   # one line per sample in each
```

```r
se <- readRDS("<outdir>/bambu/se_multiSample.rds")
colData(se)$group   # must match each sample's group in the samplesheet
```

### V6. Confirm the test-data library assumption

`test_data/testing.yml` was set to `library: "ONT_cDNA"` / `stranded_library: false` based on
`download-ref-fastq.sh` describing the ENCODE files as "CapTrap cDNA". Unstranded is the
conservative default — if those CapTrap libraries come out oriented, flip to `true`.

### V7. `maxlen` reaches chopper (P3)

```bash
grep -h "chopper" work/*/*/.command.sh
```

- Default run (`maxlen` unset) → **no** `--maxlength` flag present.
- `--maxlen 10000` → `--maxlength 10000` present.
- Confirm chopper's flag really is `--maxlength` on the pinned container version (`chopper --help`).

### V8. `resourceLimits` caps retries (P4)

Force a retry on a `process_high_memory` task under `-profile small` and confirm the second attempt
requests 200 GB rather than 400 GB. Check the ceilings match your real hardware — the values were
derived from each profile's own maximum request, not from any machine:

| Profile | cpus | memory | time |
|---|---|---|---|
| base / small | 12 | 200 GB | 20 h |
| light | 10 | 30 GB | 4 h |
| medium | 20 | 100 GB | 20 h |
| large | 20 | 500 GB | 48 h |
| test | 2 | 12 GB | 12 h |

### V9. Publishing still works after the `outputDir` removal (P5)

`outputDir` / `workflow.output.mode` were deleted from `nextflow.config`. Confirm results still land
under `params.outdir` in the expected 16 subdirectories, and that `publish_dir_mode: copy` is honoured.

### V10. `BAMBU_VALIDATE` container pull (P1)

`validate_counts` moved from `:test` to `:test3`. Confirm the image pulls and the module still runs —
`:test3` was never exercised for this process before.

### V11. POST_REFINEMENT module (new feature)

The whole feature is untested. In order of what is most likely to break:

1. **RDS row-name matching.** `post_refinement.R` subsets by `rownames(object) %in% ids`, where
   `ids` is column 1 of the validated counts files. If Bambu's SE rownames carry version suffixes
   (or the counts files do) and the other doesn't, the intersect is empty and the script errors out
   with the "None of the N validated IDs matched" message. That message is deliberate — it tells you
   the keys disagree rather than silently plotting nothing. Check with:
   ```r
   se <- readRDS("<outdir>/bambu/se_multiSample.rds"); head(rownames(se))
   head(read.table("<outdir>/bambu_validated/BambuOutput_counts_transcript_validated.txt", header = TRUE)[[1]])
   ```
2. **`plotBambu` on a row-subset object** — confirm heatmap and PCA still render once rows are
   removed, and that `colData$groupVar` survived the subset (it should; only rows are touched).
3. **Both plot sets reach the report** — the rendered HTML should show a "Sample-level Expression
   Overview" section with four two-up comparisons. If the panels are blank, the `show_plots()` vector
   handling is the thing to look at.
4. **Filename collisions in the RENDER_REPORT work dir** — raw plots stage as `pca.png`,
   `pca_grouped.png`, `heatmap_gene.png`, `heatmap_transcript.png`; post-refinement ones as the same
   names with a `_validated` suffix. Confirm all eight are present in the task's work directory and
   none overwrote another.
5. **`-stub` still works** — the BAMBU stub previously touched `bambu_output.rds`, a filename
   `bin/bambu.R` never writes. It now touches `se_multiSample.rds` and `seGene_multiSample.rds` to
   satisfy the new named outputs. Verify a full `-stub` run completes.

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

**Deliberately left alone under P6:** the `conda.enabled = false` guards in the docker/singularity/
apptainer profiles (they actively enforce the no-conda policy), and
`modules/nf-core/multiqc/{environment.yml,.conda-lock}` — those are managed by nf-core tooling via
`modules.json`, and deleting them would break `nf-core modules update`.

**Note on P5:** `publishDir` is **not** deprecated in Nextflow 25.10 — the migration guide says
nothing about removing it, and workflow outputs merely came out of preview. Both models are current.
The full output-DSL migration is deferred to D1 below.

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

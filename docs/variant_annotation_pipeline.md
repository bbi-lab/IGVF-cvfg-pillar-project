# Variant-Annotation Pipeline

`scripts/variant_annotation_pipeline.sh` replaces the archived
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb` notebook.
It runs the Dockerized [variant-annotation](https://github.com/bbi-lab/variant-annotation)
pipeline (variant mapping, reverse translation, ClinVar/gnomAD/SpliceAI/ClinGen
evidence repository/VEP/MaveDB/predictor annotation) against this project's MAVE
dataset, plus this project's own `recalculate_clingen_classification` and
`flag_variants` steps and final column filtering/renaming, to produce the
condensed and expanded integrated variant effect datasets.

## Prerequisites

- Docker and Docker Compose
- A `variant-annotation` checkout -- either the vendored git submodule at
  `vendor/variant-annotation`, or your own existing checkout (see below)

## Locating variant-annotation: submodule vs. override

`vendor/variant-annotation` is a git submodule pinned to a specific commit.
Fetch it with:

```bash
git submodule update --init vendor/variant-annotation
```

For a real run, though, set `VARIANT_ANNOTATION_DIR` to point at a checkout
that already has the large reference files and caches downloaded
(AlphaMissense, REVEL, SpliceAI VCFs, MANE, dbNSFP, `clinvar_cache/`, etc.) --
these are tens of GB and shouldn't be re-fetched per clone:

```bash
export VARIANT_ANNOTATION_DIR=/path/to/your/existing/variant-annotation/checkout
```

## Data flow

```
data/raw_mave_data/                        (you populate this -- see its README)
        |  scripts/run_variant_annotation_pipeline.sh stages it in
        v
data/intermediate/variant_annotation/data/ (gitignored; mounted as the
                                             variant-annotation pipeline's
                                             VARIANT_DATA_DIR for the run)
        |  scripts/variant_annotation_pipeline.sh runs Steps 1-18 +
        |  condensed/expanded frame assembly
        v
data/mave_data/                            (final gzipped outputs, tracked)
```

`data/raw_mave_data/README.md` documents exactly which files are required
there before running the pipeline.

## Running it

```bash
scripts/run_variant_annotation_pipeline.sh
```

This:
1. Rsyncs `data/raw_mave_data/` into `data/intermediate/variant_annotation/data/`
   (non-destructive -- only overwrites the hand-provided seed files; never
   touches the `variant-annotation` checkout's own working tree).
2. Exports `VARIANT_DATA_DIR` to that staged directory and `CVFG_PROJECT_DIR`
   to this repo's root, `cd`s into the `variant-annotation` checkout, and runs
   `scripts/variant_annotation_pipeline.sh` from there.
3. Gzips the two final `integrated_variant_effect_dataset*.tsv` files into
   `data/mave_data/`.

`scripts/variant_annotation_pipeline.sh` itself is not meant to be run
directly -- it assumes the working directory and environment variables the
orchestrator sets up.

## Running a single step

Each numbered step in `scripts/variant_annotation_pipeline.sh` is defined as
a `step_N` shell function, where `step_N` reads `data/cvfg_variants.<N-1>.tsv`
and writes `data/cvfg_variants.<N>.tsv` (`N` is 1-18; see that script's header
comment for the full list, and its `case` in
`recalculate_clingen_classification`/`flag_variants` above for how 16-18 map
to this project's own steps). To run just one step -- for debugging a
specific annotation step without re-running everything before it -- pass
`--step N`:

```bash
scripts/run_variant_annotation_pipeline.sh --step 9
```

This still stages `data/raw_mave_data/` and exports `VARIANT_DATA_DIR`/
`CVFG_PROJECT_DIR` exactly as a full run does (so step 9's inputs resolve the
same way they would mid-pipeline), but runs only `step_9` and skips the final
gzip step, since a single step doesn't produce the final integrated dataset
files. Its output lands at
`data/intermediate/variant_annotation/data/cvfg_variants.9.tsv`. The step
must already have its input file present from a previous run (full or
single-step) of the step before it.

The flatten step and the condensed/expanded frame assembly at the end of the
pipeline aren't part of this numbered sequence (they don't fit the
`<N-1>.tsv -> <N>.tsv` pattern) and always run together, only as the tail of
a full (no `--step`) run.

## The `VARIANT_DATA_DIR` path-mapping subtlety

Worth understanding before editing `scripts/variant_annotation_pipeline.sh` or
its wrapper scripts. `variant-annotation`'s `src/scripts/run_*.sh` wrappers
decide, per path argument, whether to route it through the whole-repo bind
mount (`.:/usr/src/app`) or the data-dir mount
(`${VARIANT_DATA_DIR:-./data}:/work`) by checking **whether the path already
exists on the host relative to the current working directory** at the moment
the wrapper runs. If it exists, it's treated as part of the repo and mapped to
`/usr/src/app/<path>`; if not, it falls back to `/work/<path>`.

Historically, `scripts/variant_annotation_pipeline.sh`'s `cvfg_variants.*`/
`integrated_variant_effect_dataset*` paths were all written with a bare
`data/` prefix (e.g. `data/cvfg_variants.0.tsv`) and relied entirely on this
fallback: these files existed directly inside `variant-annotation`'s own
`./data/` working copy before this project existed, so every step resolved
via the bind-mount branch, and `VARIANT_DATA_DIR` was never actually
exercised. Once these files instead lived only in **our own** staged
directory, the same bare `data/...` string kept working only because it no
longer existed in whichever `variant-annotation` checkout was in use, so the
`/work` fallback happened to take over -- an implicit, easy-to-break
dependency (a `variant-annotation` checkout shared across projects could
easily have its own stray `data/cvfg_variants.0.tsv` left over from a
different run, silently routing a step at our staged data). Every such path
in `scripts/variant_annotation_pipeline.sh` is therefore now written
explicitly: as `/work/data/...` for Dockerized steps (all step_N's honor an
already-`/work`- or `/usr/src/app`-prefixed path verbatim, skipping the
host-existence check entirely), or as `"$VARIANT_DATA_DIR/data/..."` for the
handful of plain, non-Dockerized commands (step_14's and step_16's
in-script `awk`/`python3` calls, and `build_condensed_and_expanded_frames`'s
`awk`/`postprocess_integrated_variant_effect_dataset.sh` calls, all of which
run directly on the host with the `variant-annotation` checkout as `cwd`, so
a bare `data/...` there is neither Docker-remapped nor meaningful as a
container path). This is also why our staged directory needs a `data/`
subfolder of its own -- `data/intermediate/variant_annotation/data/...` --
matching the `data/...` prefix every reference still uses after `/work`:
the orchestrator script rsyncs into `data/intermediate/variant_annotation/data/`
rather than `data/intermediate/variant_annotation/` directly.

Most large external reference files the pipeline references (SpliceAI VCFs,
`data/MANE.GRCh38.v1.5.summary.txt`, `clinvar_cache/`, dbNSFP, etc.) are
untouched by any of this: they still exist on the host inside the
`variant-annotation` checkout, so they keep resolving through the bind-mount
branch regardless of `VARIANT_DATA_DIR`. That's why they're left in place
rather than duplicated into this project.

**Exception: Step 12's AlphaMissense/REVEL files.** Unlike the rest, `step_12`
passes `--alphamissense-file`/`--revel-file` as `/work/data/...` paths, so
`AlphaMissense_hg38.tsv.gz` and `revel_hg38.tsv.gz` (plus their `.tbi`
indexes) must be copied into
`data/intermediate/variant_annotation/data/` -- they don't resolve from a
`variant-annotation` checkout for this step. See [`build_training_variant_files`:
a preparatory step before Step 12](#build_training_variant_files-a-preparatory-step-before-step-12)
below for why all of Step 12's file inputs are staged this way.

**Exception: `score_sets.tsv` and `Supplementary_Data_3.xlsx` (Steps 11, 14,
and 15).** These two CVFG-specific inputs live in this project's own
`data/input/maves/` (committed to git, alongside `data/input/predictors/` --
see [`build_training_variant_files`](#build_training_variant_files-a-preparatory-step-before-step-12)
below) rather than in a `variant-annotation` checkout or
`data/raw_mave_data/`. Unlike the AlphaMissense/REVEL files above, staging
them is automatic: `scripts/run_variant_annotation_pipeline.sh` copies both
into `data/intermediate/variant_annotation/data/` on every run, right after
rsyncing `data/raw_mave_data/` in. `step_11`'s `--requested-calibrations-file`
flag, `step_14`'s `merge-columns` extra-file argument, and `step_15`'s
`merge-columns` extra-file argument all reference them as
`/work/data/score_sets.tsv` and `/work/data/Supplementary_Data_3.xlsx` rather
than bare `data/...` paths, so all three steps resolve to our staged copies
regardless of whether the `variant-annotation` checkout in use happens to
have its own files at those paths -- see `data/raw_mave_data/README.md`.

## `add_mavedb_active_calibration_columns`, `recalculate_clingen_classification`, and `flag_variants`: containerized the same way

Step 13 ("Choose the active functional classification") and the
"Recalculate ClinGen classification" and "Flag variants" steps in
`scripts/variant_annotation_pipeline.sh` call
`src/scripts/run_add_mavedb_active_calibration_columns.sh`,
`src/scripts/run_recalculate_clingen_classification.sh`, and
`src/scripts/run_flag_variants.sh` (in *this* repo, invoked by absolute path
via `$CVFG_PROJECT_DIR` since the script runs with the `variant-annotation`
checkout as cwd) instead of a bare `python3 src/<module>.py`. Those wrappers
follow the same conventions as `variant-annotation`'s own wrappers
(`docker compose --profile tools run`, host-existence-based path mapping) and
share the same `VARIANT_DATA_DIR` value, so all containers read/write the
same staged data directory. `Dockerfile` and `compose.yaml` at the repo root
build a lean image from this project's own Poetry dependencies for this
purpose.

`add_mavedb_active_calibration_columns` runs on `cvfg_variants.12.tsv`
(right after `annotate_predictors`, step 12) and writes `cvfg_variants.13.tsv`
-- see `docs/add_mavedb_active_calibration_columns.md`. It's a Dockerized
Python/click port of `add_mavedb_active_calibration_columns.sh` from the
sibling `variant-annotation` project's own `src/scripts/`; unlike that
script, whose input/output paths default to a hard-coded
`data/cvfg/v13/cvfg_variants.12.tsv` -> `cvfg_variants.13.tsv`, this version
requires them as explicit arguments, matching `flag_variants`/
`recalculate_clingen_classification`'s convention.

`recalculate_clingen_classification` runs on `cvfg_variants.16.tsv` (right
after `annotate_simplified_consequence`) and writes `cvfg_variants.17.tsv`,
which `flag_variants` then reads to produce `.18.tsv` -- see
`docs/recalculate_clingen_classification.md`.

## `build_training_variant_files`: a preparatory step before Step 12

`step_12` (REVEL and AlphaMissense annotation) first calls
`src/scripts/run_build_training_variant_files.sh` (also from the CVFG pillar
project) to regenerate `revel_training_variants.tsv` and
`mutpred2_training_variants.tsv` from this project's own upstream
training-variant sources, then passes them to `run_annotate_predictors.sh`
via `--revel-training-file`/`--mutpred2-training-file` -- see
`docs/build_training_variant_files.md`.

All five of `step_12`'s predictor-file flags (`--alphamissense-file`,
`--mutpred2-properties-file`, `--revel-file`, `--revel-training-file`,
`--mutpred2-training-file`) are written as `/work/data/...` rather than a
bare `data/...` path, so every one of them resolves against our own
`data/intermediate/variant_annotation/data/` regardless of which
`variant-annotation` checkout is in use. This matters because
`run_annotate_predictors.sh` only remaps its two *positional* input/output
arguments through its `/work`-vs-repo-bind-mount host-existence check (see
[the `VARIANT_DATA_DIR` path-mapping
subtlety](#the-variant_data_dir-path-mapping-subtlety) above); extra flags
are passed through verbatim, and the container's cwd is the
variant-annotation checkout's own bind mount (`/usr/src/app`), not `/work`.
A bare `data/...` path on one of these flags would resolve against whichever
checkout is in use (vendored submodule or a `VARIANT_ANNOTATION_DIR`
override) instead of our staged data.

This means `AlphaMissense_hg38.tsv.gz` and `revel_hg38.tsv.gz` (normally
left in the `variant-annotation` checkout for every other step -- see
[above](#the-variant_data_dir-path-mapping-subtlety)) must also be copied,
together with their `.tbi` indexes, into
`data/intermediate/variant_annotation/data/` specifically for Step 12.

Unlike `recalculate_clingen_classification`/`flag_variants`, this step reads
and writes paths that are fixed under this project's own tree (`data/input/predictors/`
and `data/intermediate/variant_annotation/data/`) rather than the staged
`cvfg_variants.*.tsv` sequence, so its wrapper script doesn't need the
`VARIANT_DATA_DIR` `/work` mount -- see the wrapper's own notes.

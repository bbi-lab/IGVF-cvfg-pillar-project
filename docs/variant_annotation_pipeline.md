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

`scripts/variant_annotation_pipeline.sh`'s paths are all written with a
`data/` prefix (e.g. `data/cvfg_variants.0.tsv`). Historically these files
existed directly inside `variant-annotation`'s own `./data/` working copy, so
every step resolved via the bind-mount branch, and `VARIANT_DATA_DIR` was
never actually exercised. To keep that same `data/...` string working when we
instead want it to resolve against **our own** staged directory, our staged
directory needs a `data/` subfolder of its own --
`data/intermediate/variant_annotation/data/...` -- so that once the path no
longer exists in `variant-annotation`'s working tree, the `/work` fallback
(`/work/data/...`) correctly lands inside it. This is why the orchestrator
script rsyncs into `data/intermediate/variant_annotation/data/` rather than
`data/intermediate/variant_annotation/` directly.

The large external reference files the pipeline references (AlphaMissense,
REVEL, SpliceAI VCFs, `data/MANE.GRCh38.v1.5.summary.txt`, `clinvar_cache/`,
etc.) are untouched by any of this: they still exist on the host inside the
`variant-annotation` checkout, so they keep resolving through the bind-mount
branch regardless of `VARIANT_DATA_DIR`. That's why they're left in place
rather than duplicated into this project.

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

Unlike `recalculate_clingen_classification`/`flag_variants`, this step reads
and writes paths that are fixed under this project's own tree (`data/input/predictors/`
and `data/intermediate/variant_annotation/data/`) rather than the staged
`cvfg_variants.*.tsv` sequence, so its wrapper script doesn't need the
`VARIANT_DATA_DIR` `/work` mount -- see the wrapper's own notes.

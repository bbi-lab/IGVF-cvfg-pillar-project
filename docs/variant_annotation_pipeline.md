# Variant-Annotation Pipeline

`scripts/variant_annotation_pipeline.sh` replaces the archived
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb` notebook.
It runs the Dockerized [variant-annotation](https://github.com/bbi-lab/variant-annotation)
pipeline (variant mapping, reverse translation, ClinVar/gnomAD/SpliceAI/ClinGen
evidence repository/VEP/MaveDB/predictor annotation) against this project's MAVE
dataset, plus this project's own `postprocess_mavedb_functional_classifications`,
`recalculate_clingen_classification`, and `flag_variants` steps and final
column filtering/renaming, to produce the condensed and expanded integrated
variant effect datasets.

## Prerequisites

- Docker and Docker Compose
- A `variant-annotation` checkout -- either the vendored git submodule at
  `vendor/variant-annotation`, or your own existing checkout (see below)
- Step 7's gnomAD Hail table cache, built or refreshed once per checkout
  (see [gnomAD Hail table cache](#gnomad-hail-table-cache-prerequisite-for-step-7)
  below) -- **not needed** if you're pointing `VARIANT_ANNOTATION_DIR` at an
  existing checkout where this has already been done

## Locating variant-annotation: submodule vs. override

`vendor/variant-annotation` is a git submodule pinned to a specific commit.
Fetch it with:

```bash
git submodule update --init vendor/variant-annotation
```

For a real run, though, set `VARIANT_ANNOTATION_DIR` to point at a checkout
that already has the large reference files and caches downloaded
(AlphaMissense, REVEL, SpliceAI VCFs, dbNSFP, `clinvar_cache/`, etc.) --
these are tens of GB and shouldn't be re-fetched per clone:

```bash
export VARIANT_ANNOTATION_DIR=/path/to/your/existing/variant-annotation/checkout
```

## gnomAD Hail table cache (prerequisite for Step 7)

Step 7 (gnomAD annotation) reads from a local Hail table cache rather than
querying gnomAD directly; that cache has to be built before Step 7 can run.
Building/refreshing it is deliberately **not** part of the normal pipeline
run -- it's a one-time (per checkout) preparation step, since a full build
can take ~6-7 hours:

```bash
scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache
```

This runs `prepare_gnomad_cache` from `scripts/variant_annotation_pipeline.sh`
-- the `--download-only --refresh-cache` invocation that used to live inline
in Step 7 -- against a local Hail table copy expected at
`$VARIANT_DATA_DIR/gnomAD/gnomad.joint.v4.1.sites.ht` (i.e.
`data/intermediate/variant_annotation/gnomAD/...` when using the default
staging directory, or wherever `VARIANT_DATA_DIR` points if you've overridden
it). It skips all of the CVFG-specific data staging (`data/raw_mave_data/`,
`score_sets.tsv`, etc.) since none of it is relevant to a cache build.

The built cache is written to the `variant-annotation-gnomad-cache` Docker
named volume, which is scoped to the `variant-annotation` checkout's Compose
project rather than to any particular `VARIANT_DATA_DIR` staging directory.
**This means you only need to run `--prepare-gnomad-cache` once per
`variant-annotation` checkout** -- if `VARIANT_ANNOTATION_DIR` already points
at a checkout where Step 7 (or this prep step) has been run before, the
cache already exists in that checkout's volume and Step 7 will use it
directly; you don't need to rebuild or refresh it again. Only re-run this
prep step if you're using a fresh checkout, switching to a different gnomAD
release/version, or need to pick up new histogram fields (`--refresh-cache`
forces a rebuild even if a cache already exists).

## Data flow

```
data/raw_mave_data/                        (you populate this -- see its README)
data/input/maves/cvfg_variants.0.tsv       (Step 1 input; see exception below)
        |  scripts/run_variant_annotation_pipeline.sh stages both in
        v
data/intermediate/variant_annotation/data/ (gitignored; mounted as the
                                             variant-annotation pipeline's
                                             VARIANT_DATA_DIR for the run)
        |  scripts/variant_annotation_pipeline.sh runs Steps 1-21
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
and writes `data/cvfg_variants.<N>.tsv` (`N` is 1-21; see that script's header
comment for the full list, and its `case` in
`recalculate_clingen_classification`/`flag_variants` above for how 17-19 map
to this project's own steps) -- except `step_21`, which flattens
`cvfg_variants.20.tsv` and writes the condensed and expanded final integrated
MAVE dataset files directly rather than a numbered `cvfg_variants.21.tsv`. To
run just one step -- for debugging a specific annotation step without
re-running everything before it -- pass `--step N`:

```bash
scripts/run_variant_annotation_pipeline.sh --step 9
```

This still stages `data/raw_mave_data/` and exports `VARIANT_DATA_DIR`/
`CVFG_PROJECT_DIR` exactly as a full run does (so step 9's inputs resolve the
same way they would mid-pipeline), but runs only `step_9` and skips the final
gzip into `data/mave_data/` -- even `--step 21` only writes its outputs to
the staged intermediate directory. Its output lands at
`data/intermediate/variant_annotation/data/cvfg_variants.9.tsv`. The step
must already have its input file present from a previous run (full or
single-step) of the step before it.

`prepare_gnomad_cache` is not part of the numbered sequence -- see [gnomAD
Hail table cache](#gnomad-hail-table-cache-prerequisite-for-step-7) above;
it's invoked with its own flag (`--prepare-gnomad-cache`) rather than running
automatically or via `--step`.

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
handful of plain, non-Dockerized commands (step_15's and step_20's in-script
`awk` calls, both of which run directly on the host with the
`variant-annotation` checkout as `cwd`, so a bare `data/...` there is neither
Docker-remapped nor meaningful as a container path). This is also why our staged directory needs a `data/`
subfolder of its own -- `data/intermediate/variant_annotation/data/...` --
matching the `data/...` prefix every reference still uses after `/work`:
the orchestrator script rsyncs into `data/intermediate/variant_annotation/data/`
rather than `data/intermediate/variant_annotation/` directly.

Most large external reference files the pipeline references (SpliceAI VCFs,
`clinvar_cache/`, dbNSFP, etc.) are untouched by any of this: they still
exist on the host inside the `variant-annotation` checkout, so they keep
resolving through the bind-mount branch regardless of `VARIANT_DATA_DIR`.
That's why they're left in place rather than duplicated into this project.

**Exception: the MANE summary file (Step 2).** Unlike the large reference
files above, the MANE summary is small enough (a few MB gzipped) to commit
to git, so it lives in this project's own `data/input/reference/` as
`MANE.GRCh38.v1.5.summary.txt.gz` rather than relying on whichever
`variant-annotation` checkout is in use. `scripts/run_variant_annotation_pipeline.sh`
copies it into `data/intermediate/variant_annotation/data/` on every run,
alongside `score_sets.tsv`/`Supplementary_Data_3.xlsx` below. `step_2`'s
`--mane-file` flag references it as
`/work/data/MANE.GRCh38.v1.5.summary.txt.gz` rather than a bare `data/...`
path, so it resolves to our staged copy regardless of whether the
`variant-annotation` checkout in use happens to have its own MANE file at
that relative path -- `run_remap_transcript_ids.sh` remaps `--mane-file`
through the same `/work`-vs-repo-bind-mount host-existence check as its
positional input/output args, so a bare path would otherwise resolve
against that checkout's copy instead (`_open_maybe_gzipped` in
`remap_transcript_ids.py` reads the `.gz` transparently) -- see
`data/raw_mave_data/README.md`.

**Exception: Step 13's AlphaMissense/REVEL/MutPred2-properties files.** Unlike
the rest, `step_13` passes `--alphamissense-file`/`--revel-file`/
`--mutpred2-properties-file` as `/work/data/...` paths, so
`AlphaMissense_hg38.tsv.gz`, `revel_hg38.tsv.gz` (plus their `.tbi` indexes),
and `data_frame_missense_variants_MP2_properties.csv.gz` must all be copied
into `data/intermediate/variant_annotation/data/` -- they don't resolve from
a `variant-annotation` checkout for this step. Unlike the first two (large,
generic, non-CVFG reference data), `data_frame_missense_variants_MP2_properties.csv.gz`
*is* CVFG-specific -- a MutPred2 scores export for this project's assayed
variants -- but at ~367MB it's still too large to commit to git the way
`data/input/predictors/`'s other files are, so it gets the same manual-copy
treatment rather than living in `data/raw_mave_data/` (which is git-tracked)
-- see `data/raw_mave_data/README.md`. See [`build_training_variant_files`:
a preparatory step before Step 13](#build_training_variant_files-a-preparatory-step-before-step-13)
below for why all of Step 13's file inputs are staged this way.

**Exception: `score_sets.tsv` and `Supplementary_Data_3.xlsx` (Steps 11, 15,
and 16).** These two CVFG-specific inputs live in this project's own
`data/input/maves/` (committed to git, alongside `data/input/predictors/` --
see [`build_training_variant_files`](#build_training_variant_files-a-preparatory-step-before-step-13)
below) rather than in a `variant-annotation` checkout or
`data/raw_mave_data/`. Unlike the AlphaMissense/REVEL files above, staging
them is automatic: `scripts/run_variant_annotation_pipeline.sh` copies both
into `data/intermediate/variant_annotation/data/` on every run, right after
rsyncing `data/raw_mave_data/` in. `step_11`'s `--requested-calibrations-file`
flag, `step_15`'s `merge-columns` extra-file argument, and `step_16`'s
`merge-columns` extra-file argument all reference them as
`/work/data/score_sets.tsv` and `/work/data/Supplementary_Data_3.xlsx` rather
than bare `data/...` paths, so all three steps resolve to our staged copies
regardless of whether the `variant-annotation` checkout in use happens to
have its own files at those paths -- see `data/raw_mave_data/README.md`.

**Exception: `cvfg_variants.0.tsv` (Step 1) lives at `data/input/maves/cvfg_variants.0.tsv`.**
Unlike `score_sets.tsv`/`Supplementary_Data_3.xlsx` above, `scripts/run_variant_annotation_pipeline.sh`
only copies it (as `cvfg_variants.0.tsv`) into the staged directory when
Step 1 is actually about to run -- a full run, or `--step 1` -- since a
later `--step N` reads an already-staged `cvfg_variants.<N-1>.tsv` and
doesn't need it. `data/raw_mave_data/` no longer carries this file.

## `postprocess_mavedb_functional_classifications`, `add_mavedb_active_calibration_columns`, `annotate_simplified_consequence`, `recalculate_clingen_classification`, and `flag_variants`: containerized the same way

Step 12 ("Fix known MaveDB functional-classification overrides"), step 14
("Choose the active functional classification"), step 17 ("Simplified
consequence"), and the "Recalculate ClinGen classification" and "Flag
variants" steps in `scripts/variant_annotation_pipeline.sh` call
`src/scripts/run_postprocess_mavedb_functional_classifications.sh`,
`src/scripts/run_add_mavedb_active_calibration_columns.sh`,
`src/scripts/run_annotate_simplified_consequence.sh`,
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

`postprocess_mavedb_functional_classifications` runs on `cvfg_variants.11.tsv`
(right after `annotate_mavedb`, step 11) and writes `cvfg_variants.12.tsv` --
see `docs/postprocess_mavedb_functional_classifications.md`. It's a
Dockerized Python/click port of `postprocess_integrated_variant_effect_dataset.sh`
from the sibling `variant-annotation` project's own `src/scripts/` (an ad
hoc, uncommitted script there); this version drops the no-longer-needed
`CHEK2_McCarthy_Leo_2024` dataset rename and re-keys its functional-class
override on MaveDB variant URN prefix rather than dataset name, since dataset
names aren't merged in until step 15.

`add_mavedb_active_calibration_columns` runs on `cvfg_variants.13.tsv`
(right after `annotate_predictors`, step 13) and writes `cvfg_variants.14.tsv`
-- see `docs/add_mavedb_active_calibration_columns.md`. It's a Dockerized
Python/click port of `add_mavedb_active_calibration_columns.sh` from the
sibling `variant-annotation` project's own `src/scripts/`; unlike that
script, whose input/output paths default to a hard-coded
`data/cvfg/v13/cvfg_variants.12.tsv` -> `cvfg_variants.13.tsv`, this version
requires them as explicit arguments, matching `flag_variants`/
`recalculate_clingen_classification`'s convention.

`annotate_simplified_consequence` runs on `cvfg_variants.16.tsv` and writes
`cvfg_variants.17.tsv` -- see `docs/annotate_simplified_consequence.md`.
Unlike the other four, it isn't a port of a `variant-annotation` script:
it's a Dockerized port of the `get_simplified_consequence` cell from the
archived `notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`
notebook, reading a small CVFG-specific VEP-term-to-SO-summary-term mapping
committed at `data/input/consequence/extended_ensembl_consequence.csv.gz`
(hence needing both the repo bind mount and the `VARIANT_DATA_DIR` `/work`
mount, like `flag_variants`'s `--filtering-dir`).

`recalculate_clingen_classification` runs on `cvfg_variants.17.tsv` (right
after `annotate_simplified_consequence`) and writes `cvfg_variants.18.tsv`,
which `flag_variants` then reads to produce `.19.tsv` -- see
`docs/recalculate_clingen_classification.md`.

## `build_training_variant_files`: a preparatory step before Step 13

`step_13` (REVEL and AlphaMissense annotation) first calls
`src/scripts/run_build_training_variant_files.sh` (also from the CVFG pillar
project) to regenerate `revel_training_variants.tsv` and
`mutpred2_training_variants.tsv` from this project's own upstream
training-variant sources, then passes them to `run_annotate_predictors.sh`
via `--revel-training-file`/`--mutpred2-training-file` -- see
`docs/build_training_variant_files.md`.

All five of `step_13`'s predictor-file flags (`--alphamissense-file`,
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
`data/intermediate/variant_annotation/data/` specifically for Step 13.

Unlike `recalculate_clingen_classification`/`flag_variants`, this step reads
and writes paths that are fixed under this project's own tree (`data/input/predictors/`
and `data/intermediate/variant_annotation/data/`) rather than the staged
`cvfg_variants.*.tsv` sequence, so its wrapper script doesn't need the
`VARIANT_DATA_DIR` `/work` mount -- see the wrapper's own notes.

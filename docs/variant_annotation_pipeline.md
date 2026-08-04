# Variant-Annotation Pipeline

`scripts/variant_annotation_pipeline.sh` replaces the archived
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb` notebook.
It runs the Dockerized [variant-annotation](https://github.com/bbi-lab/variant-annotation)
pipeline (variant mapping, reverse translation, ClinVar/gnomAD/SpliceAI/ClinGen
evidence repository/VEP/MaveDB/predictor annotation) against this project's MAVE
dataset, plus this project's own `flag_variants` step and final column
filtering/renaming, to produce the condensed and expanded integrated variant
effect datasets.

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
        |  scripts/variant_annotation_pipeline.sh runs Steps 1-17 +
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

## `flag_variants`: containerized the same way

The "Flag variants" step in `scripts/variant_annotation_pipeline.sh` calls
`src/scripts/run_flag_variants.sh` (in *this* repo, invoked by absolute path
via `$CVFG_PROJECT_DIR` since the script runs with the `variant-annotation`
checkout as cwd) instead of a bare `python3 src/flag_variants.py`. That
wrapper follows the same conventions as `variant-annotation`'s own wrappers
(`docker compose --profile tools run`, host-existence-based path mapping) and
shares the same `VARIANT_DATA_DIR` value, so both containers read/write the
same staged data directory. `Dockerfile` and `compose.yaml` at the repo root
build a lean image from this project's own Poetry dependencies for this
purpose.

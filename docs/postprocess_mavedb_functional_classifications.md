# Postprocess MaveDB Functional Classifications

`src/postprocess_mavedb_functional_classifications.py` applies known one-off
overrides to MaveDB functional-classification categories, keyed by MaveDB
variant URN prefix rather than dataset name (dataset names haven't been
merged in yet at this point in the pipeline -- that happens in a later step,
from `score_sets.tsv`).

## Current overrides

| Variant URN prefix | Dataset | Matches label | Overrides classification to |
|---|---|---|---|
| `urn:mavedb:00000674-b-1#` | KCNE1_Muhammad_2024_potassium_flux | `gain-of-function` | `not_specified` |

Each override is applied independently to all three MaveDB calibration
column groups present in the input --
`mavedb.primary_calibration.*`, `mavedb.investigator_provided_calibration.*`,
and the optional `mavedb.requested_calibration.*` -- so whichever one is
later chosen as the "active" calibration
(`add_mavedb_active_calibration_columns.py`, step 14) already carries the
corrected category.

## Source logic

Ported from `postprocess_integrated_variant_effect_dataset.sh` in the sibling
`variant-annotation` project. That script did two things, keyed by the
`Dataset` column on the fully-assembled integrated dataset:

1. Renamed `CHEK2_McCarthy_Leo_2024` to `CHEK2_McCarthy-Leo_2024`. **Not
   ported** -- the dataset is written with the hyphenated name from the
   start, so this is no longer needed.
2. Overrode `auth_reported_func_class_category` to `not_specified` for
   `KCNE1_Muhammad_2024_potassium_flux` rows with
   `auth_reported_func_class == "gain-of-function"`. Ported here, but moved
   earlier in the pipeline (step 12, right after `annotate_mavedb`, step 11)
   and re-keyed on the MaveDB variant URN prefix instead of the dataset name,
   since `Dataset`/`auth_reported_func_class*` don't exist until several
   steps later.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.postprocess_mavedb_functional_classifications \
  data/cvfg_variants.11.tsv data/cvfg_variants.12.tsv
```

Via Docker (same image as `flag_variants`/`recalculate_clingen_classification`/
`add_mavedb_active_calibration_columns`, see `compose.yaml`):

```bash
src/scripts/run_postprocess_mavedb_functional_classifications.sh \
  data/cvfg_variants.11.tsv data/cvfg_variants.12.tsv
```

Like `run_add_mavedb_active_calibration_columns.sh`, this wrapper maps its
input/output paths against the `/work` staging mount
(`${VARIANT_DATA_DIR:-./data}`), since it reads from the same staged pipeline
data those files live in rather than a committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--variant-urn-col` | `variant_urn` | Input column containing MaveDB variant URNs |

## Stdout report

After applying overrides, the script prints how many row/calibration-group
pairs were changed, broken down by URN prefix, matched label, and calibration
group -- e.g. to confirm the KCNE1 override actually found matching rows.

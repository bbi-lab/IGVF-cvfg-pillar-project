# Add MaveDB Active-Calibration Columns

`src/add_mavedb_active_calibration_columns.py` picks one MaveDB functional
classification per row and copies it into a single
`mavedb.active_calibration.*` column group:

| Column | Description |
|---|---|
| `mavedb.active_calibration.urn` | URN of the calibration that was selected |
| `mavedb.active_calibration.name` | Title of the selected calibration |
| `mavedb.active_calibration.url` | URL to the score set page on MaveDB |
| `mavedb.active_calibration.functional_class_label` | Functional classification label for this variant under the selected calibration |
| `mavedb.active_calibration.functional_classification` | The `functionalClassification` value (`normal`, `abnormal`, or `not_specified`) under the selected calibration |

## Source logic

`annotate_mavedb.py` (step 11 of the variant-annotation pipeline; see
`docs/variant_annotation_pipeline.md`) writes up to three sets of calibration
columns per row -- `mavedb.primary_calibration.*`,
`mavedb.investigator_provided_calibration.*`, and the optional
`mavedb.requested_calibration.*` (only present when that step was run with
`--requested-calibrations-file`). This script selects one per row, in
priority order:

1. `mavedb.requested_calibration.urn` is non-empty -> use the requested
   calibration.
2. otherwise, `mavedb.primary_calibration.name` is one of
   `--preferred-primary-calibration-name` (default: `Fayer calibration`,
   `Scott calibration`) -> use the primary calibration.
3. anything else (including a blank primary-calibration name) -> use the
   investigator-provided calibration.

The `mavedb.requested_calibration.*` columns are optional: if entirely absent
from the input, every row falls through to steps 2/3 as if no calibration had
been requested for it. Raises `click.ClickException` if the required
`mavedb.primary_calibration.*`/`mavedb.investigator_provided_calibration.*`
columns are missing, or if only some of the `mavedb.requested_calibration.*`
columns are present.

Ported from `add_mavedb_active_calibration_columns.sh` in the sibling
`variant-annotation` project (an awk script with the column names and
priority logic above hard-coded); this version is a Dockerized Python/click
script with those column prefixes and the preferred-primary-calibration
names all exposed as CLI options, and required (rather than defaulted)
input/output file arguments.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.add_mavedb_active_calibration_columns \
  data/cvfg_variants.13.tsv data/cvfg_variants.14.tsv
```

Via Docker (same image as `flag_variants`/`recalculate_clingen_classification`,
see `compose.yaml`):

```bash
src/scripts/run_add_mavedb_active_calibration_columns.sh \
  data/cvfg_variants.13.tsv data/cvfg_variants.14.tsv
```

Like `run_recalculate_clingen_classification.sh`, this wrapper maps its
input/output paths against the `/work` staging mount
(`${VARIANT_DATA_DIR:-./data}`), since it reads from the same staged pipeline
data those files live in rather than a committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--primary-prefix` | `mavedb.primary_calibration` | Column prefix for the primary-calibration columns |
| `--investigator-prefix` | `mavedb.investigator_provided_calibration` | Column prefix for the investigator-provided-calibration columns |
| `--requested-prefix` | `mavedb.requested_calibration` | Column prefix for the optional requested-calibration columns |
| `--output-prefix` | `mavedb.active_calibration` | Column prefix for the five columns this script writes |
| `--preferred-primary-calibration-name` | `Fayer calibration`, `Scott calibration` | Primary-calibration name preferred over the investigator-provided calibration when no calibration was requested for a row; repeatable for more than one name |

## Stdout report

After selecting, the script prints how many rows used each calibration
source (`requested`, `primary`, `investigator_provided`), letting you sanity-check,
e.g., that the requested-calibrations file is actually being picked up for
the score sets it covers.

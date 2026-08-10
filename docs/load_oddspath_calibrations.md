# Load OddsPath Calibrations

`src/load_oddspath_calibrations.py` loads `OddsPath_calculations.ipynb`'s
per-dataset likelihood-ratio output into Supplementary Data 4, so the
workbook stays in sync with the latest notebook run without hand-editing
spreadsheet rows -- the same role `load_excalibr_calibrations` plays for the
`ExCALIBR_calibrations` sheet (see
[`docs/load_excalibr_calibrations.md`](load_excalibr_calibrations.md)).

It reads the notebook's output CSV (default:
`data/output/mave_calibration/OddsPath_calibrations.csv.gz`) and writes one
row per dataset, sorted by `Dataset`, into the `OddsPath_calibrations` sheet
of a workbook (default:
`data/output/supplementary_data/Supplementary_Data_4.xlsx`),
**overwriting that sheet's existing contents in place**. All other sheets,
the sheet's tab position, and its column widths are left untouched.

## Column mapping

Unlike `load_excalibr_calibrations`, there's no reformatting step: the CSV
already has the same fourteen columns, in the same order, as the sheet
(`Dataset`, `Total Controls`, `OddsNormal`, `OddsAbnormal`,
`Pathogenic Controls`, `Benign Controls`, `Prior Probability Pathogenic`,
`Total Assay Abnormal`, `True Path in Abnormal`, `Total Assay Normal`,
`True Path in Normal`, `Pseudocount Details`, `Evidence Code Normal`,
`Evidence Code Abnormal`), so each row is copied through as-is. The script
raises `ValueError` (a `click.ClickException` from the CLI) if the CSV's
columns don't match exactly, so a schema change in
`OddsPath_calculations.ipynb` fails loudly instead of silently writing a
mismatched sheet.

Two things worth knowing about the values themselves:

- A missing value (`NaN` in the CSV, e.g. `Pseudocount Details` when no
  pseudocount adjustment was needed) is written as a blank cell, not the
  literal string `"nan"`.
- `OddsNormal`/`OddsAbnormal` can legitimately hold a string instead of a
  number -- e.g. `"No benign controls"`, `"No functionally abnormal
  controls"` -- when a dataset lacks the controls needed to compute that
  likelihood ratio (see `Odds_LR` in the notebook). These pass through
  unchanged. Because pandas reads a CSV column as all-strings the moment any
  row in it is non-numeric, the script re-parses each cell in these two
  columns individually (`_coerce_numeric_or_string`) so datasets that *do*
  have a numeric likelihood ratio are still written as numbers, not as the
  string `"192.5806"`.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.load_oddspath_calibrations \
  data/output/mave_calibration/OddsPath_calibrations.csv.gz \
  data/output/supplementary_data/Supplementary_Data_4.xlsx
```

Both positional arguments default to the paths above, so a bare invocation
works from the repo root.

Via Docker (same image as `flag_variants`/`mave_dataset_stats`/
`load_excalibr_calibrations`, see `compose.yaml`):

```bash
src/scripts/run_load_oddspath_calibrations.sh [csv] [workbook]
```

Like `run_load_excalibr_calibrations.sh`, this wrapper doesn't map paths
against a `/work` staging mount -- the script only reads/writes
`data/output/` files, which are already available at `/usr/src/app` via the
whole-repo bind mount.

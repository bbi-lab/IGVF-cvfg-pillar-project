# Derive Score-Set URN

`src/derive_score_set_urn.py` derives a `score_set_urn` column from
`variant_urn` by stripping everything from the first `#` onward (MaveDB
variant URNs are score-set URNs with a `#<index>` suffix, e.g.
`urn:mavedb:00000097-0-2#1` -> `urn:mavedb:00000097-0-2`). Step 14 of
`scripts/variant_annotation_pipeline.sh` runs this immediately before a
`merge-columns` call that joins in `dataset_name` from `score_sets.tsv` by
`score_set_urn`.

## Why this replaced an `awk` pass

Step 14 originally derived `score_set_urn` with a line-oriented `awk` script.
MaveDB occasionally persists a multi-line `HTTPStatusError` message in the
`mavedb_mapping_error` column (quoted per RFC 4180, since it's written by a
proper CSV/TSV writer upstream in `fetch_mavedb_scores.py`). Every
CSV/TSV-aware reader in the pipeline (pandas, `csv.DictReader` in
`merge_columns.py`) correctly treats such a field as part of a single
logical row -- but `awk`'s `NR`-based line splitting does not, so it treated
each such row as two records. The derived `score_set_urn` for the first
physical "half" got appended inside the still-open quoted
`mavedb_mapping_error` field instead of as a real column, so when
`merge-columns` (which is quote-aware) re-read the file and reconstituted
the true row, there was no valid `score_set_urn` value left for it to key on
-- the join to `dataset_name` silently failed for every affected row
(observed for the five `urn:mavedb:00000097-0-2` BRCA1 rows in a real run,
the only rows in that dataset with an embedded newline in any field).

This script reads and writes the TSV with pandas instead, which parses
quoted fields (including multi-line ones) the same way `merge-columns`
does, so `score_set_urn` is derived from the correct, complete `variant_urn`
value for every row regardless of what any other column contains.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.derive_score_set_urn \
  data/cvfg_variants.13.tsv data/cvfg_variants.14.temp.tsv
```

Via Docker (same image as `add_mavedb_active_calibration_columns`/
`flag_variants`, see `compose.yaml`):

```bash
src/scripts/run_derive_score_set_urn.sh \
  data/cvfg_variants.13.tsv data/cvfg_variants.14.temp.tsv
```

Like `run_add_mavedb_active_calibration_columns.sh`, this wrapper maps its
input/output paths against the `/work` staging mount
(`${VARIANT_DATA_DIR:-./data}`), since it reads from the same staged pipeline
data those files live in rather than a committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--variant-urn-column` | `variant_urn` | Column to derive `score_set_urn` from |
| `--output-column` | `score_set_urn` | Column name to write the derived URN to |

Raises `click.ClickException` if `--variant-urn-column` is missing from the
input.

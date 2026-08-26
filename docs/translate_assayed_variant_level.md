# Translate Assayed Variant Level

`src/translate_assayed_variant_level.py` translates the codes used in the
`assayed_variant_level` column: `protein` -> `aa`, `dna` -> `nt`. Values not
in that map are left unchanged. This is the code-translation half of step 20
of `scripts/variant_annotation_pipeline.sh`; see
`docs/merge_rna_score_columns.md` for the other half.

## Why this replaced an `awk` pass

Step 20 originally translated these codes with a line-oriented `awk` script.
Like step 14's old `awk` pass (see `docs/derive_score_set_urn.md`), it read
the file line-by-line, so a value spanning multiple physical lines -- e.g. a
multi-line `mavedb_mapping_error` value, quoted per RFC 4180 -- was treated
as two records instead of one, corrupting `assayed_variant_level` (and every
other column) for the affected rows.

This script reads and writes the TSV with pandas instead, which parses
quoted fields (including multi-line ones) correctly, so the translation is
applied to the right column for every row regardless of what any other
column contains.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.translate_assayed_variant_level \
  data/cvfg_variants.19.tsv data/cvfg_variants.20.codes.tsv
```

Via Docker (same image as `derive_score_set_urn`/`flag_variants`, see
`compose.yaml`):

```bash
src/scripts/run_translate_assayed_variant_level.sh \
  data/cvfg_variants.19.tsv data/cvfg_variants.20.codes.tsv
```

Like `run_derive_score_set_urn.sh`, this wrapper maps its input/output paths
against the `/work` staging mount (`${VARIANT_DATA_DIR:-./data}`), since it
reads from the same staged pipeline data those files live in rather than a
committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--column` | `assayed_variant_level` | Column to translate codes in |

Raises `click.ClickException` if `--column` is missing from the input.

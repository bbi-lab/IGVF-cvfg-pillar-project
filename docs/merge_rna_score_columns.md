# Merge RNA Score Columns

`src/merge_rna_score_columns.py` fills blank `rna_score` values from
`rna_score_d6` (the day-6 RNA score). `rna_score`, `rna_score_d6`, and
`rna_score_d20` come from `fetch_mavedb_scores.py`/`src/mavedb_scores.sql`
(see `docs/fetch_mavedb_scores.md`); some MAVE datasets only ever populate
the day-6 score, so without this merge those measurements are dropped from
the `rna_score` column that step 21 of
`scripts/variant_annotation_pipeline.sh` carries into the final integrated
dataset (`integrated_variant_effect_dataset.tsv` /
`integrated_variant_effect_dataset.condensed.tsv`).

This is the rna_score-merge half of step 20 of
`scripts/variant_annotation_pipeline.sh`, run immediately after
`translate_assayed_variant_level.py`; see
`docs/translate_assayed_variant_level.md` for the other half.

A row whose `rna_score` is already populated is left unchanged, even if
`rna_score_d6` also has a value -- `rna_score` wins whenever both are
present.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.merge_rna_score_columns \
  data/cvfg_variants.20.codes.tsv data/cvfg_variants.20.tsv
```

Via Docker (same image as `derive_score_set_urn`/`flag_variants`, see
`compose.yaml`):

```bash
src/scripts/run_merge_rna_score_columns.sh \
  data/cvfg_variants.20.codes.tsv data/cvfg_variants.20.tsv
```

Like `run_derive_score_set_urn.sh`, this wrapper maps its input/output paths
against the `/work` staging mount (`${VARIANT_DATA_DIR:-./data}`), since it
reads from the same staged pipeline data those files live in rather than a
committed repo file.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--target-column` | `rna_score` | Column to fill blanks in |
| `--source-column` | `rna_score_d6` | Column to fill blanks from |

Raises `click.ClickException` if either column is missing from the input.

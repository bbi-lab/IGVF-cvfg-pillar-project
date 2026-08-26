#!/usr/bin/env python3
"""Merge `rna_score_d6` into `rna_score` wherever `rna_score` is blank.

Step 20 of `scripts/variant_annotation_pipeline.sh`, alongside
`translate_assayed_variant_level.py`. `rna_score` and `rna_score_d6` come
from `fetch_mavedb_scores.py`/`src/mavedb_scores.sql` (see
`docs/fetch_mavedb_scores.md`); some MAVE datasets only ever populate the
day-6 RNA score, so their measurements would otherwise be dropped from the
`rna_score` column that `scripts/variant_annotation_pipeline.sh`'s step 21
carries into the final integrated dataset.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_TARGET_COLUMN = "rna_score"
DEFAULT_SOURCE_COLUMN = "rna_score_d6"


def merge_rna_score(df, target_column=DEFAULT_TARGET_COLUMN, source_column=DEFAULT_SOURCE_COLUMN):
    """Return a Series of `df[target_column]`, filled from `df[source_column]` where blank.

    Raises ValueError if either column is missing from `df`.
    """
    missing = [col for col in (target_column, source_column) if col not in df.columns]
    if missing:
        raise ValueError(f"input is missing required column(s): {', '.join(missing)}")

    return df[target_column].mask(df[target_column] == "", df[source_column])


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--target-column",
    default=DEFAULT_TARGET_COLUMN,
    show_default=True,
    help="Column to fill blanks in",
)
@click.option(
    "--source-column",
    default=DEFAULT_SOURCE_COLUMN,
    show_default=True,
    help="Column to fill blanks from",
)
def main(input, output, target_column, source_column):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        merged = merge_rna_score(df, target_column, source_column)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    filled_count = int(((df[target_column] == "") & (merged != "")).sum())
    df[target_column] = merged

    click.echo(f"Filled {target_column} from {source_column} for {filled_count} of {len(df)} row(s).")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

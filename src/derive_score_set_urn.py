#!/usr/bin/env python3
"""Derive a `score_set_urn` column from `variant_urn`.

MaveDB variant URNs are score-set URNs with a `#<index>` suffix identifying
the variant within that score set (e.g. `urn:mavedb:00000097-0-2#1`). This
script strips everything from the first `#` onward, so a later
`merge-columns` step (step 14 of `scripts/variant_annotation_pipeline.sh`)
can join in `dataset_name` from `score_sets.tsv` by `score_set_urn`.

Reads/writes the TSV with pandas rather than a line-oriented tool (the
previous `awk` implementation) so that a `variant_urn`/`mavedb_mapping_error`
value spanning multiple physical lines -- MaveDB occasionally persists a
multi-line `HTTPStatusError` message in `mavedb_mapping_error`, quoted per
RFC 4180 -- is still read as a single logical row instead of being split
across two, which silently corrupted the derived `score_set_urn` for the
affected rows (see `docs/derive_score_set_urn.md`).
"""

from pathlib import Path

import click
import pandas as pd

VARIANT_URN_COLUMN = "variant_urn"
DEFAULT_OUTPUT_COLUMN = "score_set_urn"


def derive_score_set_urn(df, variant_urn_column=VARIANT_URN_COLUMN, output_column=DEFAULT_OUTPUT_COLUMN):
    """Return a Series of `df[variant_urn_column]` with everything from the first `#` onward removed.

    Raises ValueError if `variant_urn_column` is missing from `df`.
    """
    if variant_urn_column not in df.columns:
        raise ValueError(f"input is missing required column: {variant_urn_column}")

    return df[variant_urn_column].str.split("#", n=1).str[0].rename(output_column)


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--variant-urn-column",
    default=VARIANT_URN_COLUMN,
    show_default=True,
    help="Column to derive score_set_urn from",
)
@click.option(
    "--output-column",
    default=DEFAULT_OUTPUT_COLUMN,
    show_default=True,
    help="Column name to write the derived URN to",
)
def main(input, output, variant_urn_column, output_column):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        df[output_column] = derive_score_set_urn(df, variant_urn_column, output_column)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Derived {output_column} for {len(df)} row(s).")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

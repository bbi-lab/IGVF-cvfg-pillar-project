#!/usr/bin/env python3
"""Translate `assayed_variant_level` codes ("protein" -> "aa", "dna" -> "nt").

Step 20 of `scripts/variant_annotation_pipeline.sh`.

Reads/writes the TSV with pandas rather than a line-oriented tool (the
previous `awk` implementation) so that a value spanning multiple physical
lines -- e.g. a multi-line `mavedb_mapping_error` value, quoted per RFC 4180
-- is still read as a single logical row instead of being split across two,
which silently corrupted `assayed_variant_level` (and every other column) for
the affected rows. Same underlying bug, and same fix, as
`derive_score_set_urn.py` (step 14; see `docs/derive_score_set_urn.md`).
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_COLUMN = "assayed_variant_level"
CODE_MAP = {"protein": "aa", "dna": "nt"}


def translate_assayed_variant_level(df, column=DEFAULT_COLUMN, code_map=CODE_MAP):
    """Return a Series with `code_map`'s translations applied to `df[column]`.

    Values not present in `code_map` are left unchanged. Raises ValueError if
    `column` is missing from `df`.
    """
    if column not in df.columns:
        raise ValueError(f"input is missing required column: {column}")

    return df[column].replace(code_map)


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--column",
    default=DEFAULT_COLUMN,
    show_default=True,
    help="Column to translate codes in",
)
def main(input, output, column):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        df[column] = translate_assayed_variant_level(df, column)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Translated {column} codes for {len(df)} row(s).")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

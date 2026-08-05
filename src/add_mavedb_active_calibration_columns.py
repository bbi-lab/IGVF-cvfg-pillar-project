#!/usr/bin/env python3
"""Choose the active MaveDB functional-classification calibration for each row.

`annotate_mavedb.py` (step 11 of the variant-annotation pipeline) writes up to
three sets of MaveDB calibration columns per row --
`mavedb.primary_calibration.*`, `mavedb.investigator_provided_calibration.*`,
and the optional `mavedb.requested_calibration.*` (only present when that step
was run with `--requested-calibrations-file`). This script picks one
calibration per row and copies its five fields (`urn`, `name`, `url`,
`functional_class_label`, `functional_classification`) into a single
`mavedb.active_calibration.*` column group, in this priority order:

\b
1. `mavedb.requested_calibration.urn` is non-empty -> use the requested
   calibration.
2. otherwise, `mavedb.primary_calibration.name` is one of
   `--preferred-primary-calibration-name` (default: "Fayer calibration",
   "Scott calibration") -> use the primary calibration.
3. anything else (including a blank primary-calibration name) -> use the
   investigator-provided calibration.

The `mavedb.requested_calibration.*` columns are optional: if entirely absent
from the input, every row falls through to steps 2/3 as if no calibration had
been requested for it.

Ported from `add_mavedb_active_calibration_columns.sh` in the sibling
`variant-annotation` project.
"""

from pathlib import Path

import click
import pandas as pd

FIELDS = ("urn", "name", "url", "functional_class_label", "functional_classification")

DEFAULT_PRIMARY_PREFIX = "mavedb.primary_calibration"
DEFAULT_INVESTIGATOR_PREFIX = "mavedb.investigator_provided_calibration"
DEFAULT_REQUESTED_PREFIX = "mavedb.requested_calibration"
DEFAULT_OUTPUT_PREFIX = "mavedb.active_calibration"
DEFAULT_PREFERRED_PRIMARY_CALIBRATION_NAMES = ("Fayer calibration", "Scott calibration")


def _prefixed_columns(prefix):
    return [f"{prefix}.{field}" for field in FIELDS]


def select_active_calibration(
    df,
    primary_prefix=DEFAULT_PRIMARY_PREFIX,
    investigator_prefix=DEFAULT_INVESTIGATOR_PREFIX,
    requested_prefix=DEFAULT_REQUESTED_PREFIX,
    output_prefix=DEFAULT_OUTPUT_PREFIX,
    preferred_primary_calibration_names=DEFAULT_PREFERRED_PRIMARY_CALIBRATION_NAMES,
):
    """Return `(active_calibration_df, source)` for `df`.

    `active_calibration_df` has one `{output_prefix}.{field}` column per
    `FIELDS` entry. `source` is a per-row Series of `"requested"`,
    `"primary"`, or `"investigator_provided"`, indicating which calibration
    was copied.

    Raises ValueError if any `primary_prefix.*`/`investigator_prefix.*`
    column is missing from `df`, or if only some of the `requested_prefix.*`
    columns are present.
    """
    primary_cols = _prefixed_columns(primary_prefix)
    investigator_cols = _prefixed_columns(investigator_prefix)
    requested_cols = _prefixed_columns(requested_prefix)

    missing_required = [col for col in primary_cols + investigator_cols if col not in df.columns]
    if missing_required:
        raise ValueError(f"input is missing required column(s): {', '.join(missing_required)}")

    requested_present = [col for col in requested_cols if col in df.columns]
    has_requested = len(requested_present) == len(requested_cols)
    if requested_present and not has_requested:
        missing = [col for col in requested_cols if col not in df.columns]
        raise ValueError(f"input has some but not all {requested_prefix}.* columns; missing: {', '.join(missing)}")

    if has_requested:
        use_requested = df[requested_cols[0]].str.len().gt(0)
    else:
        use_requested = pd.Series(False, index=df.index)
    use_primary = (~use_requested) & df[primary_cols[1]].isin(list(preferred_primary_calibration_names))

    source = pd.Series("investigator_provided", index=df.index)
    source.loc[use_primary] = "primary"
    source.loc[use_requested] = "requested"

    active = pd.DataFrame(index=df.index)
    for field, primary_col, investigator_col in zip(FIELDS, primary_cols, investigator_cols):
        values = df[investigator_col].where(~use_primary, df[primary_col])
        if has_requested:
            requested_col = f"{requested_prefix}.{field}"
            values = values.where(~use_requested, df[requested_col])
        active[f"{output_prefix}.{field}"] = values

    return active, source


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--primary-prefix",
    default=DEFAULT_PRIMARY_PREFIX,
    show_default=True,
    help="Column prefix for the primary-calibration columns",
)
@click.option(
    "--investigator-prefix",
    default=DEFAULT_INVESTIGATOR_PREFIX,
    show_default=True,
    help="Column prefix for the investigator-provided-calibration columns",
)
@click.option(
    "--requested-prefix",
    default=DEFAULT_REQUESTED_PREFIX,
    show_default=True,
    help="Column prefix for the optional requested-calibration columns",
)
@click.option(
    "--output-prefix",
    default=DEFAULT_OUTPUT_PREFIX,
    show_default=True,
    help="Column prefix for the five active-calibration columns this script writes",
)
@click.option(
    "--preferred-primary-calibration-name",
    "preferred_primary_calibration_names",
    multiple=True,
    default=DEFAULT_PREFERRED_PRIMARY_CALIBRATION_NAMES,
    show_default=True,
    help="Primary-calibration name preferred over the investigator-provided calibration when no "
    "calibration was requested for a row; repeatable for more than one name",
)
def main(
    input,
    output,
    primary_prefix,
    investigator_prefix,
    requested_prefix,
    output_prefix,
    preferred_primary_calibration_names,
):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        active, source = select_active_calibration(
            df,
            primary_prefix=primary_prefix,
            investigator_prefix=investigator_prefix,
            requested_prefix=requested_prefix,
            output_prefix=output_prefix,
            preferred_primary_calibration_names=preferred_primary_calibration_names,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    df[active.columns] = active

    click.echo(f"Selected an active calibration for {len(df)} row(s).")
    click.echo("By calibration source:")
    for source_name, count in source.value_counts().sort_index().items():
        click.echo(f"  {source_name}: {count}")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

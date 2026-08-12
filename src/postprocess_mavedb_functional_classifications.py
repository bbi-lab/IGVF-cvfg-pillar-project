#!/usr/bin/env python3
"""Apply known one-off overrides to MaveDB functional-classification categories.

Runs right after `annotate_mavedb.py` (step 11 of the variant-annotation
pipeline), which writes up to three sets of MaveDB calibration columns per
row -- `mavedb.primary_calibration.*`, `mavedb.investigator_provided_calibration.*`,
and the optional `mavedb.requested_calibration.*`. Dataset names haven't been
merged in yet at this point (that happens in a later step, from
`score_sets.tsv`), so overrides are keyed by MaveDB variant URN prefix
(everything up to and including the final "#") rather than dataset name.

Currently handles one override:

\b
- KCNE1_Muhammad_2024_potassium_flux (urn:mavedb:00000674-b-1#*): where a
  calibration's functional_class_label is "gain-of-function", its
  functional_classification is overridden to "not_specified".

The override is applied independently to each of the three calibration
column groups present in the input, so whichever one is later chosen as the
"active" calibration (`add_mavedb_active_calibration_columns.py`) already
carries the corrected category.

Ported from `postprocess_integrated_variant_effect_dataset.sh` in the
sibling `variant-annotation` project, which also renamed the
"CHEK2_McCarthy_Leo_2024" dataset to "CHEK2_McCarthy-Leo_2024" -- that rename
is no longer needed (the dataset is written with the hyphenated name from the
start) and isn't ported here.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_VARIANT_URN_COL = "variant_urn"

CALIBRATION_PREFIXES = (
    "mavedb.primary_calibration",
    "mavedb.investigator_provided_calibration",
    "mavedb.requested_calibration",
)

FUNCTIONAL_CLASS_OVERRIDES = [
    {
        "urn_prefix": "urn:mavedb:00000674-b-1#",
        "functional_class_label": "gain-of-function",
        "replacement_functional_classification": "not_specified",
    },
]


def apply_functional_class_overrides(
    df,
    variant_urn_col=DEFAULT_VARIANT_URN_COL,
    overrides=FUNCTIONAL_CLASS_OVERRIDES,
    calibration_prefixes=CALIBRATION_PREFIXES,
):
    """Apply `overrides` to `df` in place; return per-override row counts changed.

    Each override matches rows whose `variant_urn_col` starts with its
    `urn_prefix`, then -- for each calibration column group in
    `calibration_prefixes` that's present in `df` -- overrides
    `{prefix}.functional_classification` to `replacement_functional_classification`
    wherever `{prefix}.functional_class_label` equals `functional_class_label`.

    Raises ValueError if `variant_urn_col` isn't present in `df`.
    """
    if variant_urn_col not in df.columns:
        raise ValueError(f"Column {variant_urn_col!r} not found in input")

    counts = {}
    for override in overrides:
        urn_mask = df[variant_urn_col].str.startswith(override["urn_prefix"])
        for prefix in calibration_prefixes:
            label_col = f"{prefix}.functional_class_label"
            classification_col = f"{prefix}.functional_classification"
            if label_col not in df.columns or classification_col not in df.columns:
                continue

            change_mask = urn_mask & (df[label_col] == override["functional_class_label"])
            n_changed = int(change_mask.sum())
            if n_changed:
                df.loc[change_mask, classification_col] = override["replacement_functional_classification"]
            counts[(override["urn_prefix"], override["functional_class_label"], prefix)] = n_changed

    return counts


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--variant-urn-col",
    default=DEFAULT_VARIANT_URN_COL,
    show_default=True,
    help="Input column containing MaveDB variant URNs",
)
def main(input, output, variant_urn_col):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        counts = apply_functional_class_overrides(df, variant_urn_col=variant_urn_col)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    total = sum(counts.values())
    click.echo(f"Applied functional-classification overrides to {total} row/calibration pair(s).")
    for (urn_prefix, functional_class_label, prefix), n_changed in counts.items():
        if n_changed:
            click.echo(f"  {urn_prefix} {functional_class_label!r} ({prefix}): {n_changed}")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

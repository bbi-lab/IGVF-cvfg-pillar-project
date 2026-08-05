#!/usr/bin/env python3
"""Add a simplified_consequence column to a CVFG variants TSV.

Collapses each row's VEP consequence term(s) down to a single "SO summary
term" per DNA candidate (e.g. "missense_variant" rather than the more
granular "missense_variant^splice_region_variant"), using a curated mapping
from VEP output term to SO summary term.

Ported from the `get_simplified_consequence` function in
notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb, which
worked from a single flat `consequence` column (comma-separated VEP terms,
one row per variant) and picked the first comma-separated term with a mapping
table entry. The variant-annotation pipeline's `annotate_vep.py` instead
produces a **pipe-delimited** `vep.most_severe_mutational_consequence`
column, aligned one segment per DNA candidate (same convention as this
project's own `flag_variants.py`/`recalculate_clingen_classification.py`),
with each segment already reduced to VEP's own single most-severe term for
that candidate -- so this script maps each segment independently and
re-joins with "|", rather than picking among multiple comma-separated terms
itself.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_CONSEQUENCE_COL = "vep.most_severe_mutational_consequence"
DEFAULT_CONSEQUENCE_MAP_FILE = Path("data/input/consequence/extended_ensembl_consequence.csv.gz")
OUTPUT_COL = "simplified_consequence"

# Terms that already are their own SO summary term in some upstream sources,
# not present in the Ensembl consequence table itself. Ported verbatim from
# `extra_terms` in the notebook cited in the module docstring.
EXTRA_TERMS = ["UTR_variant", "splicing_variant", "splice_site_variant"]


def build_mapping_dict(consequence_map_file):
    """VEP output term -> SO summary term, from the consequence-map file.

    Ported from the notebook cited in the module docstring: when a VEP
    output term appears on more than one row of the mapping file, the row
    with the highest `Importance` wins.
    """
    consequence_df = pd.read_csv(consequence_map_file)
    mapping_dict = (
        consequence_df.sort_values(by="Importance", ascending=False)
        .groupby("VEP output term")[["SO summary term", "Importance"]]
        .first()
        .to_dict()["SO summary term"]
    )
    mapping_dict.update({term: term for term in EXTRA_TERMS})
    return mapping_dict


def simplify_segment(term, mapping_dict):
    """Map one VEP term, treating a missing or NaN mapping-table entry as "".

    Some VEP terms in the mapping file (e.g. `non_coding_transcript_exon_variant`)
    have no `SO summary term` at all (NaN in the source CSV) -- these are
    matched-but-unmapped, same as a term absent from the table entirely.
    """
    if not term:
        return ""
    mapped = mapping_dict.get(term)
    return mapped if isinstance(mapped, str) else ""


def simplify_consequence(value, mapping_dict):
    """Map each "|"-delimited DNA candidate's VEP term independently."""
    segments = value.split("|") if value else [""]
    return "|".join(simplify_segment(term, mapping_dict) for term in segments)


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--consequence-col",
    default=DEFAULT_CONSEQUENCE_COL,
    show_default=True,
    help="Column with pipe-delimited VEP consequence term(s) per DNA candidate",
)
@click.option(
    "--consequence-map-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CONSEQUENCE_MAP_FILE,
    show_default=True,
    help="CSV mapping 'VEP output term' to 'SO summary term' (with an 'Importance' column)",
)
def main(input, output, consequence_col, consequence_map_file):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    if consequence_col not in df.columns:
        raise click.ClickException(f"Column {consequence_col!r} not found in input")

    mapping_dict = build_mapping_dict(consequence_map_file)
    df[OUTPUT_COL] = df[consequence_col].apply(lambda v: simplify_consequence(v, mapping_dict))

    click.echo(f"Simplified consequence for {len(df)} row(s).")

    counts = df[OUTPUT_COL].apply(lambda v: v.split("|")).explode()
    counts = counts.replace("", "(no consequence)")
    for term, count in counts.value_counts().sort_index().items():
        click.echo(f"  {term}: {count}")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Add simplified_consequence, condensed_consequence, splice_variant, and
splice_var_amino columns to a CVFG variants TSV.

``simplified_consequence`` collapses each row's VEP consequence term(s) down
to a single "SO summary term" per DNA candidate (e.g. "missense_variant"
rather than the more granular "missense_variant^splice_region_variant"),
using a curated mapping from VEP output term to SO summary term. Ported from
the `get_simplified_consequence` function in
notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb, which
worked from a single flat `consequence` column (comma-separated VEP terms,
one row per variant) and picked the first comma-separated term with a
mapping table entry. The variant-annotation pipeline's `annotate_vep.py`
instead produces a **pipe-delimited** `vep.most_severe_mutational_consequence`
column, aligned one segment per DNA candidate (same convention as this
project's own `flag_variants.py`/`recalculate_clingen_classification.py`),
with each segment already reduced to VEP's own single most-severe term for
that candidate -- so this script maps each segment independently and
re-joins with "|", rather than picking among multiple comma-separated terms
itself.

``condensed_consequence`` collapses ``simplified_consequence`` to a single
value: that value if all of its pipe-delimited entries are identical, else
"conflict_consequence".

``splice_variant`` is derived, position-by-position, from the pipe-aligned
SpliceAI delta-score columns (``spliceai.ds_ag``/``ds_al``/``ds_dg``/``ds_dl``)
and ``simplified_consequence``: a position is "Yes" if any of its four delta
scores is >= 0.2 (--splice-score-threshold) or its ``simplified_consequence``
entry is "splice_site_variant", else "No".

``splice_var_amino`` collapses ``splice_variant`` to a single value: blank if
no non-blank entries, that value if all non-blank entries agree, else "Yes".

Differences from the upstream variant-annotation project's own
`src/annotate_simplified_consequence.py` (retained deliberately):

- Reads the already-VEP-reduced ``vep.most_severe_mutational_consequence``
  column (one term per DNA candidate) rather than upstream's raw
  ``vep.mutational_consequences`` column (one ``^``-delimited set of
  candidate VEP terms per DNA candidate, requiring a first-match-wins
  priority lookup within each set). Both columns come from the same
  `annotate_vep.py`; this project uses the pre-reduced one, so each segment
  maps directly through the table with no priority selection needed.
- Output column names (``simplified_consequence``, ``condensed_consequence``,
  ``splice_variant``, ``splice_var_amino``) are fixed rather than
  configurable via CLI flags, matching this project's other annotation
  scripts (e.g. `flag_variants.py`'s hard-coded ``Flag`` column) rather than
  upstream's fully-parameterized output-column-name options.
- No ``--skip``/``--limit``/``--log-level``/``--csv-field-size-limit``
  options: this script loads the whole TSV into a pandas DataFrame (like
  this project's other annotation scripts) rather than streaming rows with
  `csv.DictReader`, so those options don't apply.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_CONSEQUENCE_COL = "vep.most_severe_mutational_consequence"
DEFAULT_CONSEQUENCE_MAP_FILE = Path("data/input/reference/extended_ensembl_consequence.csv.gz")
DEFAULT_SPLICEAI_DS_AG_COL = "spliceai.ds_ag"
DEFAULT_SPLICEAI_DS_AL_COL = "spliceai.ds_al"
DEFAULT_SPLICEAI_DS_DG_COL = "spliceai.ds_dg"
DEFAULT_SPLICEAI_DS_DL_COL = "spliceai.ds_dl"
DEFAULT_SPLICE_SCORE_THRESHOLD = 0.2

SIMPLIFIED_CONSEQUENCE_COL = "simplified_consequence"
CONDENSED_CONSEQUENCE_COL = "condensed_consequence"
SPLICE_VARIANT_COL = "splice_variant"
SPLICE_VAR_AMINO_COL = "splice_var_amino"

# condensed_consequence value used when simplified_consequence's pipe-delimited
# entries disagree.
CONFLICT_CONSEQUENCE = "conflict_consequence"

# Terms that already are their own SO summary term in some upstream sources,
# not present in the Ensembl consequence table itself. Ported verbatim from
# `extra_terms` in the notebook cited in the module docstring, plus
# "no_change": annotate_vep.py's unchanged-variant fallback (used by this
# project's own step 3 reverse translation with --wt-codon-mode unambiguous,
# same as upstream) writes "no_change" directly into
# vep.most_severe_mutational_consequence.
EXTRA_TERMS = ["UTR_variant", "splicing_variant", "splice_site_variant", "no_change"]


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


def split_pipe(value):
    """Split a pipe-delimited value, treating a blank value as zero entries."""
    return value.split("|") if value else []


def compute_condensed_consequence(simplified_consequence):
    """Collapse a pipe-delimited simplified_consequence value to a single value.

    "" if blank; that value if all pipe-delimited entries are identical; else
    CONFLICT_CONSEQUENCE.
    """
    if not simplified_consequence:
        return ""
    values = simplified_consequence.split("|")
    unique = set(values)
    if len(unique) == 1:
        return next(iter(unique))
    return CONFLICT_CONSEQUENCE


def meets_splice_threshold(value, threshold):
    """Return True if *value* parses as a float >= threshold."""
    if not value:
        return False
    try:
        return float(value) >= threshold
    except ValueError:
        return False


def compute_splice_variant(ds_ag, ds_al, ds_dg, ds_dl, simplified_consequence, threshold):
    """Return a pipe-delimited Yes/No column aligned to the SpliceAI/consequence positions.

    A position is "Yes" if any of its four delta scores is >= threshold, or its
    simplified_consequence entry is "splice_site_variant"; otherwise "No".
    Returns "" if none of the inputs have any positions.
    """
    ag_list = split_pipe(ds_ag)
    al_list = split_pipe(ds_al)
    dg_list = split_pipe(ds_dg)
    dl_list = split_pipe(ds_dl)
    cons_list = split_pipe(simplified_consequence)

    n = max(len(ag_list), len(al_list), len(dg_list), len(dl_list), len(cons_list))
    if n == 0:
        return ""

    results = []
    for i in range(n):
        scores = (
            ag_list[i] if i < len(ag_list) else "",
            al_list[i] if i < len(al_list) else "",
            dg_list[i] if i < len(dg_list) else "",
            dl_list[i] if i < len(dl_list) else "",
        )
        consequence = cons_list[i] if i < len(cons_list) else ""
        is_splice = (
            any(meets_splice_threshold(score, threshold) for score in scores)
            or consequence == "splice_site_variant"
        )
        results.append("Yes" if is_splice else "No")
    return "|".join(results)


def compute_splice_var_amino(splice_variant):
    """Collapse a pipe-delimited splice_variant value to a single value.

    Blank entries are dropped; "" if none remain, that value if all remaining
    entries agree, else "Yes".
    """
    values = {v for v in splice_variant.split("|") if v} if splice_variant else set()
    if not values:
        return ""
    if len(values) == 1:
        return next(iter(values))
    return "Yes"


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
@click.option(
    "--spliceai-ds-ag-col",
    default=DEFAULT_SPLICEAI_DS_AG_COL,
    show_default=True,
    help="Column with pipe-delimited SpliceAI DS_AG scores, aligned to --consequence-col",
)
@click.option(
    "--spliceai-ds-al-col",
    default=DEFAULT_SPLICEAI_DS_AL_COL,
    show_default=True,
    help="Column with pipe-delimited SpliceAI DS_AL scores, aligned to --consequence-col",
)
@click.option(
    "--spliceai-ds-dg-col",
    default=DEFAULT_SPLICEAI_DS_DG_COL,
    show_default=True,
    help="Column with pipe-delimited SpliceAI DS_DG scores, aligned to --consequence-col",
)
@click.option(
    "--spliceai-ds-dl-col",
    default=DEFAULT_SPLICEAI_DS_DL_COL,
    show_default=True,
    help="Column with pipe-delimited SpliceAI DS_DL scores, aligned to --consequence-col",
)
@click.option(
    "--splice-score-threshold",
    type=float,
    default=DEFAULT_SPLICE_SCORE_THRESHOLD,
    show_default=True,
    help="Minimum SpliceAI delta score for a position to count as a splice variant",
)
def main(
    input,
    output,
    consequence_col,
    consequence_map_file,
    spliceai_ds_ag_col,
    spliceai_ds_al_col,
    spliceai_ds_dg_col,
    spliceai_ds_dl_col,
    splice_score_threshold,
):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    if consequence_col not in df.columns:
        raise click.ClickException(f"Column {consequence_col!r} not found in input")
    for col in (spliceai_ds_ag_col, spliceai_ds_al_col, spliceai_ds_dg_col, spliceai_ds_dl_col):
        if col not in df.columns:
            raise click.ClickException(f"Column {col!r} not found in input")

    mapping_dict = build_mapping_dict(consequence_map_file)
    df[SIMPLIFIED_CONSEQUENCE_COL] = df[consequence_col].apply(lambda v: simplify_consequence(v, mapping_dict))
    df[CONDENSED_CONSEQUENCE_COL] = df[SIMPLIFIED_CONSEQUENCE_COL].apply(compute_condensed_consequence)
    df[SPLICE_VARIANT_COL] = df.apply(
        lambda row: compute_splice_variant(
            row[spliceai_ds_ag_col],
            row[spliceai_ds_al_col],
            row[spliceai_ds_dg_col],
            row[spliceai_ds_dl_col],
            row[SIMPLIFIED_CONSEQUENCE_COL],
            splice_score_threshold,
        ),
        axis=1,
    )
    df[SPLICE_VAR_AMINO_COL] = df[SPLICE_VARIANT_COL].apply(compute_splice_var_amino)

    n_simplified = (df[SIMPLIFIED_CONSEQUENCE_COL] != "").sum()
    n_condensed = (df[CONDENSED_CONSEQUENCE_COL] != "").sum()
    n_splice_variant = (df[SPLICE_VARIANT_COL] != "").sum()
    n_splice_var_amino = (df[SPLICE_VAR_AMINO_COL] != "").sum()
    click.echo(
        f"{n_simplified}/{len(df)} rows given a non-blank {SIMPLIFIED_CONSEQUENCE_COL}, "
        f"{n_condensed} a non-blank {CONDENSED_CONSEQUENCE_COL}, "
        f"{n_splice_variant} a non-blank {SPLICE_VARIANT_COL}, "
        f"{n_splice_var_amino} a non-blank {SPLICE_VAR_AMINO_COL}."
    )

    counts = df[SIMPLIFIED_CONSEQUENCE_COL].apply(lambda v: v.split("|")).explode()
    counts = counts.replace("", "(no consequence)")
    for term, count in counts.value_counts().sort_index().items():
        click.echo(f"  {term}: {count}")

    n_conflict = (df[CONDENSED_CONSEQUENCE_COL] == CONFLICT_CONSEQUENCE).sum()
    click.echo(f"{n_conflict} row(s) with conflicting per-candidate consequences.")

    n_splice_yes = (df[SPLICE_VAR_AMINO_COL] == "Yes").sum()
    click.echo(f"{n_splice_yes} row(s) flagged as splice variants.")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

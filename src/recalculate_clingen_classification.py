#!/usr/bin/env python3
"""Recalculate ClinGen Evidence Repository classifications without functional-assay evidence.

Adds two columns to a CVFG variants TSV, recomputing each row's ClinGen
Evidence Repository (erepo) classification after discarding PS3/BS3/PP3/BP4
evidence -- the codes ClinGen assigns for functional-assay results (including
their _Supporting/_Moderate/_Strong/_Very_Strong strength variants) -- so the
result reflects clinical/population evidence only:

\b
- Updated_Classification_ClinGen_repo: ACMG/AMP classification recomputed
  from the remaining evidence codes.
- Updated_Evidence Codes_ClinGen_repo: the remaining evidence codes
  themselves (comma-separated), i.e. what the classification above was
  computed from.

Both are ported from the `filter_and_recalculate` / `classify_acmg` functions
in notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb, which
worked from a single flat erepo dump (one row per classification, evidence
codes as one comma-separated string). The variant-annotation pipeline's
`annotate_erepo.py` instead produces a pipe-delimited
`clingen_evidence_repository.Applied Evidence Codes (Met)` column, aligned
one segment per `mapped_hgvs_c` DNA candidate (same convention as this
project's own `flag_variants.py`); this script applies the ported logic
independently to each candidate's segment and re-joins the two results with
"|" in the same order. When a segment itself holds evidence from more than
one matching erepo record (`annotate_erepo.py` merges those with " | "), the
codes from all of that candidate's records are pooled before filtering and
reclassifying.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_EVIDENCE_CODES_COL = "clingen_evidence_repository.Applied Evidence Codes (Met)"
DEFAULT_ORIGINAL_CLASSIFICATION_COL = "clingen_evidence_repository.Assertion"
CLASSIFICATION_COL = "Updated_Classification_ClinGen_repo"
EVIDENCE_COL = "Updated_Evidence Codes_ClinGen_repo"

# Evidence-code -> ACMG/AMP strength. Ported verbatim from
# `acmg_evidence_codes` in the notebook cited above.
ACMG_EVIDENCE_CODES = {
    "PVS1": "Pathogenic_Very_Strong",
    "PVS1_Supporting": "Pathogenic_Supporting",
    "PVS1_Moderate": "Pathogenic_Moderate",
    "PVS1_Strong": "Pathogenic_Strong",
    "PVS1_Very": "Pathogenic_Very_Strong",
    "PS1": "Pathogenic_Strong",
    "PS1_Supporting": "Pathogenic_Supporting",
    "PS1_Moderate": "Pathogenic_Moderate",
    "PS1_Very": "Pathogenic_Very_Strong",
    "PS2": "Pathogenic_Strong",
    "PS2_Supporting": "Pathogenic_Supporting",
    "PS2_Moderate": "Pathogenic_Moderate",
    "PS2_Very": "Pathogenic_Very_Strong",
    "PS3": "Pathogenic_Strong",
    "PS3_Supporting": "Pathogenic_Supporting",
    "PS3_Moderate": "Pathogenic_Moderate",
    "PS3_Very": "Pathogenic_Very_Strong",
    "PS4": "Pathogenic_Strong",
    "PS4_Supporting": "Pathogenic_Supporting",
    "PS4_Moderate": "Pathogenic_Moderate",
    "PS4_Very": "Pathogenic_Very_Strong",
    "PM1": "Pathogenic_Moderate",
    "PM1_Supporting": "Pathogenic_Supporting",
    "PM1_Strong": "Pathogenic_Strong",
    "PM1_Very": "Pathogenic_Very_Strong",
    "PM2": "Pathogenic_Moderate",
    "PM2_Supporting": "Pathogenic_Supporting",
    "PM2_Strong": "Pathogenic_Strong",
    "PM2_Very": "Pathogenic_Very_Strong",
    "PM3": "Pathogenic_Moderate",
    "PM3_Supporting": "Pathogenic_Supporting",
    "PM3_Strong": "Pathogenic_Strong",
    "PM3_Very": "Pathogenic_Very_Strong",
    "PM4": "Pathogenic_Moderate",
    "PM4_Supporting": "Pathogenic_Supporting",
    "PM4_Strong": "Pathogenic_Strong",
    "PM4_Very": "Pathogenic_Very_Strong",
    "PM5": "Pathogenic_Moderate",
    "PM5_Supporting": "Pathogenic_Supporting",
    "PM5_Strong": "Pathogenic_Strong",
    "PM5_Very": "Pathogenic_Very_Strong",
    "PM6": "Pathogenic_Moderate",
    "PM6_Supporting": "Pathogenic_Supporting",
    "PM6_Strong": "Pathogenic_Strong",
    "PM6_Very": "Pathogenic_Very_Strong",
    "PP1": "Pathogenic_Supporting",
    "PP1_Moderate": "Pathogenic_Moderate",
    "PP1_Strong": "Pathogenic_Strong",
    "PP1_Very": "Pathogenic_Very_Strong",
    "PP2": "Pathogenic_Supporting",
    "PP2_Moderate": "Pathogenic_Moderate",
    "PP2_Strong": "Pathogenic_Strong",
    "PP2_Very": "Pathogenic_Very_Strong",
    "PP3": "Pathogenic_Supporting",
    "PP3_Moderate": "Pathogenic_Moderate",
    "PP3_Strong": "Pathogenic_Strong",
    "PP3_Very": "Pathogenic_Very_Strong",
    "PP4": "Pathogenic_Supporting",
    "PP4_Moderate": "Pathogenic_Moderate",
    "PP4_Strong": "Pathogenic_Strong",
    "PP4_Very": "Pathogenic_Very_Strong",
    "PP5": "Pathogenic_Supporting",
    "PP5_Moderate": "Pathogenic_Moderate",
    "PP5_Strong": "Pathogenic_Strong",
    "PP5_Very": "Pathogenic_Very_Strong",
    "BA1": "Benign_Standalone",
    "BS1": "Benign_Strong",
    "BS1_Supporting": "Benign_Supporting",
    "BS1_Moderate": "Benign_Moderate",
    "BS1_Very": "Benign_Very_Strong",
    "BS2": "Benign_Strong",
    "BS2_Supporting": "Benign_Supporting",
    "BS2_Moderate": "Benign_Moderate",
    "BS2_Very": "Benign_Very_Strong",
    "BS3": "Benign_Strong",
    "BS3_Supporting": "Benign_Supporting",
    "BS3_Moderate": "Benign_Moderate",
    "BS3_Very": "Benign_Very_Strong",
    "BS4": "Benign_Strong",
    "BS4_Supporting": "Benign_Supporting",
    "BS4_Moderate": "Benign_Moderate",
    "BS4_Very": "Benign_Very_Strong",
    "BP1": "Benign_Supporting",
    "BP1_Moderate": "Benign_Moderate",
    "BP1_Strong": "Benign_Strong",
    "BP1_Very": "Benign_Very_Strong",
    "BP2": "Benign_Supporting",
    "BP2_Moderate": "Benign_Moderate",
    "BP2_Strong": "Benign_Strong",
    "BP2_Very": "Benign_Very_Strong",
    "BP3": "Benign_Supporting",
    "BP3_Moderate": "Benign_Moderate",
    "BP3_Strong": "Benign_Strong",
    "BP3_Very": "Benign_Very_Strong",
    "BP4": "Benign_Supporting",
    "BP4_Moderate": "Benign_Moderate",
    "BP4_Strong": "Benign_Strong",
    "BP4_Very": "Benign_Very_Strong",
    "BP5": "Benign_Supporting",
    "BP5_Moderate": "Benign_Moderate",
    "BP5_Strong": "Benign_Strong",
    "BP5_Very": "Benign_Very_Strong",
    "BP6": "Benign_Supporting",
    "BP6_Moderate": "Benign_Moderate",
    "BP6_Strong": "Benign_Strong",
    "BP6_Very": "Benign_Very_Strong",
    "BP7": "Benign_Supporting",
    "BP7_Moderate": "Benign_Moderate",
    "BP7_Strong": "Benign_Strong",
    "BP7_Very": "Benign_Very_Strong",
}

# Functional-assay evidence codes to strip before reclassifying. Ported
# verbatim from `remove_codes` in the notebook cited above.
FUNCTIONAL_EVIDENCE_CODES = {
    "BP4",
    "BP4_Supporting",
    "BP4_Moderate",
    "BP4_Strong",
    "BP4_Very",
    "BS3",
    "BS3_Supporting",
    "BS3_Moderate",
    "BS3_Strong",
    "BS3_Very",
    "PS3",
    "PS3_Supporting",
    "PS3_Moderate",
    "PS3_Strong",
    "PS3_Very",
    "PP3",
    "PP3_Supporting",
    "PP3_Moderate",
    "PP3_Strong",
    "PP3_Very",
}


def classify_acmg(evidence_codes):
    """Combine ACMG/AMP evidence-code strengths into a classification.

    Ported verbatim from `classify_acmg` in the notebook cited in the module
    docstring, including its `"Benign_standalone"` (lowercase `s`) lookup --
    the dict maps BA1 to `"Benign_Standalone"`, so that comparison never
    matches and BA1 never actually counts toward the benign-standalone rule.
    Preserved here rather than "fixed" so this script's output matches the
    classifications already published from that notebook.
    """
    if not evidence_codes:
        return "VUS"

    strengths = [ACMG_EVIDENCE_CODES.get(e) for e in evidence_codes]

    pathogenic_very_strong = strengths.count("Pathogenic_Very_Strong")
    pathogenic_strong = strengths.count("Pathogenic_Strong")
    pathogenic_moderate = strengths.count("Pathogenic_Moderate")
    pathogenic_supporting = strengths.count("Pathogenic_Supporting")

    benign_standalone = strengths.count("Benign_standalone")
    benign_strong = strengths.count("Benign_Strong")
    benign_moderate = strengths.count("Benign_Moderate")
    benign_supporting = strengths.count("Benign_Supporting")

    if pathogenic_very_strong >= 1 and (pathogenic_strong >= 1 or pathogenic_moderate >= 2):
        return "Pathogenic"
    if pathogenic_strong >= 2:
        return "Pathogenic"
    if pathogenic_strong == 1 and pathogenic_moderate >= 2:
        return "Pathogenic"
    if pathogenic_very_strong == 1 and pathogenic_moderate >= 1:
        return "Likely Pathogenic"
    if pathogenic_strong == 1 and pathogenic_moderate == 1:
        return "Likely Pathogenic"
    if pathogenic_moderate >= 3:
        return "Likely Pathogenic"
    if pathogenic_moderate == 2 and pathogenic_supporting >= 2:
        return "Likely Pathogenic"

    if benign_standalone >= 1:
        return "Benign"
    if benign_strong >= 2:
        return "Benign"
    if benign_strong == 1 and benign_supporting >= 1:
        return "Likely Benign"
    if benign_moderate >= 2:
        return "Likely Benign"

    return "VUS"


def codes_in_segment(segment):
    """Split one DNA candidate's evidence segment into individual codes.

    A segment is normally a single comma-separated list of codes from one
    erepo record (e.g. "PM2,PP3"). When more than one expert-panel record
    matches the same candidate, `annotate_erepo.py` merges them with " | "
    (e.g. "PM2,PP3 | BS1"); split those apart first so codes from every
    matching record are pooled before filtering.
    """
    codes = []
    for record in segment.split(" | "):
        codes.extend(code.strip() for code in record.split(",") if code.strip())
    return codes


def filter_and_recalculate(evidence_segment):
    """Strip functional-assay evidence and recompute the ACMG classification.

    Ported from `filter_and_recalculate` in the notebook cited in the module
    docstring. Takes one DNA candidate's evidence-codes segment (see
    `codes_in_segment`) and returns `(updated_classification,
    updated_evidence_codes)`.
    """
    if not evidence_segment:
        return "VUS", ""

    evidence_codes = codes_in_segment(evidence_segment)
    filtered_codes = [code for code in evidence_codes if code not in FUNCTIONAL_EVIDENCE_CODES]
    updated_classification = classify_acmg(filtered_codes) if filtered_codes else "VUS"
    return updated_classification, ",".join(filtered_codes)


def recalculate_row(evidence_field):
    """Apply `filter_and_recalculate` to each "|"-delimited DNA candidate."""
    segments = evidence_field.split("|") if evidence_field else [""]
    classifications, evidence = zip(*(filter_and_recalculate(segment) for segment in segments))
    return "|".join(classifications), "|".join(evidence)


def compute_updated_classifications(df, evidence_codes_col):
    """Return (classification_series, evidence_series) for `df`.

    Raises ValueError if `evidence_codes_col` isn't present in `df`.
    """
    if evidence_codes_col not in df.columns:
        raise ValueError(f"Column {evidence_codes_col!r} not found in input")

    results = df[evidence_codes_col].apply(recalculate_row)
    classification = results.apply(lambda r: r[0])
    evidence = results.apply(lambda r: r[1])
    return classification, evidence


def classification_counts(series, empty_label=None):
    """Count per-candidate classifications in a "|"-delimited column.

    Splits each row on "|" (one value per `mapped_hgvs_c` DNA candidate,
    same convention as `recalculate_row`) and counts across all rows and
    candidates. `empty_label` replaces empty segments in the count (e.g. a
    candidate with no erepo match) with a readable placeholder.
    """
    values = series.apply(lambda v: v.split("|")).explode()
    if empty_label is not None:
        values = values.replace("", empty_label)
    return values.value_counts()


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--evidence-codes-col",
    default=DEFAULT_EVIDENCE_CODES_COL,
    show_default=True,
    help="Column with pipe-delimited erepo 'Applied Evidence Codes (Met)' values",
)
@click.option(
    "--classification-col",
    default=DEFAULT_ORIGINAL_CLASSIFICATION_COL,
    show_default=True,
    help="Column with the unmodified erepo classification, for the stdout report's before/after breakdown",
)
def main(input, output, evidence_codes_col, classification_col):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")

    try:
        df[CLASSIFICATION_COL], df[EVIDENCE_COL] = compute_updated_classifications(df, evidence_codes_col)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    if classification_col not in df.columns:
        raise click.ClickException(f"Column {classification_col!r} not found in input")

    click.echo(f"Recalculated ClinGen classification for {len(df)} row(s).")

    click.echo("By recalculated classification (functional-assay evidence excluded):")
    for classification, count in classification_counts(df[CLASSIFICATION_COL]).sort_index().items():
        click.echo(f"  {classification}: {count}")

    click.echo(f"By unmodified classification ({classification_col}):")
    for classification, count in classification_counts(
        df[classification_col], empty_label="(no erepo match)"
    ).sort_index().items():
        click.echo(f"  {classification}: {count}")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

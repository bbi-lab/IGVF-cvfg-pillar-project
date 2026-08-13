#!/usr/bin/env python3
"""Summary statistics for the integrated MAVE variant effect dataset.

Reads the condensed variant effect dataset (one row per variant effect
measurement/score; see data/output/maves/integrated_variant_effect_dataset*.tsv.gz)
together with its dataset-level metadata (Supplementary_Data_3.xlsx) and
reports, for three groupings of datasets -- IGVF-produced only, non-IGVF
("community") only, and combined -- the number of datasets, variant effect
measurements, composite scores, distinct variants assayed, and genes
represented. Each bucket's variant effect measurements are also reported as a
percentage of the combined total (100% for the combined row). Immediately
after that table, a "Genes represented" section lists which genes are
covered only by IGVF datasets, only by non-IGVF ("community") datasets, or
by both -- CALM1/CALM2/CALM3 are always merged into one "CALM1/2/3" entry in
this list regardless of `--merge-calm-genes` (see that flag below), since
listing the same calmodulin target three times isn't useful.

It additionally reports, in a set of text tables, how many variants and
variant measurements have REVEL, AlphaMissense, and MutPred2 scores, and how
many fall into each of several clinical-attribute buckets (ClinVar VUS,
ClinVar pathogenic/benign, observed in gnomAD, or none of the above). This
clinical-attribute breakdown is reported twice: once using ClinVar 2025 for
every gene, and once using ClinVar 2025 for every gene except BRCA1, PTEN,
MSH2, and TP53, which use ClinVar 2018 instead. Each is reported four ways: at
the assayed-variant (protein-resolution, from the condensed file) or
DNA-variant (from the expanded file) level, and as distinct variants or as
variant measurements (i.e. rows). Each clinical-attribute table also carries
two extra columns breaking out, for each of the four buckets, how many of the
level's SNVs are in it (or, at the assayed-variant level, "SNV-accessible"
variants: those with at least one single-nucleotide-substitution candidate
among the DNA-level changes -- `hgvs_c`/`transcript_ref`/`transcript_alt` --
that reverse-translate to that protein change), as a percentage of the
level's total SNV(-accessible) count (also given in that table's totals
line).

By default, a variant with a conflicting or ambiguous ClinVar call --
disagreement between a pathogenic-leaning and benign-leaning classification
across its measurements/DNA candidates, or ClinVar's own "Conflicting
classifications of pathogenicity" call -- is excluded from both the VUS and
pathogenic-or-benign buckets and counted in its own "ClinVar conflict" bucket
instead, mirroring the conflict handling in
Analysis/Curation_summary_V5_cleaned.ipynb. Pass `--allow-clinvar-conflicts`
to instead fold such variants into VUS/pathogenic-or-benign via any-match
(the original behavior), which also drops the conflict bucket from the
report.

A dataset counts as "IGVF-produced" if Supplementary_Data_3's Curation sheet
marks its `IGVF Produced?` column "Yes". A row counts as a direct measurement
if its dataset's `Primary Score Set or Meta-analysis?` column reads "primary
score set"; everything else (meta-analyses, trained predictors, etc.) counts
as a composite score. This is a per-dataset determination: every row of a
given dataset is classified the same way.

Note on `(hgvs_g, hgvs_p)`: the integrated dataset has no `hgvs_g` (genomic
HGVS) column -- its DNA-level identifier is `hgvs_c` (transcript-relative
HGVS, pipe-delimited when a protein-resolution measurement corresponds to more
than one underlying DNA change). This script uses `hgvs_c` as that DNA-level
key, so "distinct variants assayed" counts distinct (hgvs_c, hgvs_p) pairs in
the condensed file.

Note on pipe-delimited annotation columns (REVEL, AM_score, MutPred2,
clinvar_sig_2025, gnomad_MAF): in the condensed file these follow the same
one-part-per-underlying-DNA-change convention as `hgvs_c`, and the parts don't
always agree (e.g. one DNA-level candidate has a REVEL score and another
doesn't). An assayed variant (or measurement row) counts as having a value, or
matching a target ClinVar significance, if *any* one of its pipe-delimited
parts does. In the expanded (DNA-variant-level) file there's only ever one
part, so this reduces to a plain presence/absence check there. "Observed in
gnomAD" follows this repo's existing convention elsewhere (see
notebooks/analysis/README_OddsPath_classifications.md) of treating any non-empty
`gnomad_MAF` as observed, regardless of the allele count.

By default, CALM1, CALM2, and CALM3 are counted as three separate genes, since
that's how the underlying data labels them. Pass `--merge-calm-genes` to count
them as a single gene target instead, since they encode the same calmodulin
protein.

Two further sections cover ExCALIBR calibration coverage and reclassification
agreement, sourced from `--excalibr-calibrations-file` (default
Supplementary_Data_4.xlsx) and `--controls-file` (default
Supplementary_Data_5.xlsx) respectively:

- **ExCALIBR calibration coverage** (Extended Data Figure 2): how many genes
  have a row in Supplementary_Data_4's `ExCALIBR_calibrations` sheet, and how
  many of those genes have at least one dataset where ExCALIBR assigned at
  least one point of evidence in either direction (i.e. at least one of that
  row's `range_-8`..`range_8` columns is non-null). That sheet's `dataset`
  values are matched against this script's own dataset metadata (the same
  `Dataset Name` -> `Gene` mapping used elsewhere) to determine gene
  membership, after stripping a trailing `_clinvar_2018` suffix (this sheet's
  convention for marking a mixed-year duplicate calibration of the same
  dataset) and Unicode-normalizing (NFC) both sides, since accented
  characters can round-trip through Excel in different composed/decomposed
  forms across files. Any `dataset` value that still doesn't resolve raises
  `ValueError`, same as an unrecognized dataset in the main metadata file.

- **Reclassification agreement** (Figure 4c): for every sheet in the controls
  file whose name starts with `controls_` (one per predictor/calibration
  combination -- REVEL, AlphaMissense, MutPred2 -- which should agree with
  each other since this doesn't depend on which predictor is active), how
  often ExCALIBR's evidence assignment (`ExC_points_2025`) and the functional
  class assignment (`OP_points`) agree with the row's ClinVar
  pathogenic-or-benign control label (`clnsig_group_18_25`). A row's points
  column is treated as assigning pathogenic evidence if positive, benign
  evidence if negative, and no evidence if zero or missing -- rows with no
  evidence assigned are excluded from the agreement percentage (reported
  separately as "no point of evidence assigned") but still counted in the
  section's total.

Both file arguments are optional and default to the paths above. Output is
written as plain text (to stdout, and optionally to `--output` as well).
"""

import unicodedata
from pathlib import Path

import click
import pandas as pd

DEFAULT_CONDENSED_FILE = Path("data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz")
DEFAULT_EXPANDED_FILE = Path("data/output/maves/integrated_variant_effect_dataset.tsv.gz")
DEFAULT_METADATA_FILE = Path("data/output/supplementary_data/Supplementary_Data_3.xlsx")
DEFAULT_EXCALIBR_CALIBRATIONS_FILE = Path("data/output/supplementary_data/Supplementary_Data_4.xlsx")
DEFAULT_CONTROLS_FILE = Path("data/output/supplementary_data/Supplementary_Data_5.xlsx")

METADATA_SHEET = "Curation"
DATASET_COL = "Dataset"
GENE_COL = "Gene"
GENOMIC_VARIANT_COL = "hgvs_c"
PROTEIN_VARIANT_COL = "hgvs_p"
VARIANT_KEY_COLS = [GENOMIC_VARIANT_COL, PROTEIN_VARIANT_COL]

METADATA_DATASET_COL = "Dataset Name"
IGVF_PRODUCED_COL = "IGVF Produced?"
SCORE_SET_TYPE_COL = "Primary Score Set or Meta-analysis?"
MEASUREMENT_VALUE = "primary score set"

CALM_GENES = frozenset({"CALM1", "CALM2", "CALM3"})
CALM_MERGED_LABEL = "CALM1/2/3"

SCORE_COLUMNS = {
    "REVEL": "REVEL",
    "AlphaMissense": "AM_score",
    "MutPred2": "MutPred2",
}
CLINVAR_COL = "clinvar_sig_2025"
CLINVAR_COL_2018 = "clinvar_sig_2018"
MIXED_YEAR_GENES = frozenset({"BRCA1", "PTEN", "MSH2", "TP53"})
GNOMAD_COL = "gnomad_MAF"
TRANSCRIPT_REF_COL = "transcript_ref"
TRANSCRIPT_ALT_COL = "transcript_alt"
VUS_LABEL = "VUS (ClinVar 2025)"
PATHOGENIC_OR_BENIGN_LABEL = "Pathogenic or benign (ClinVar 2025)"
CLINVAR_CONFLICT_LABEL = "ClinVar conflict (excluded from VUS/pathogenic-or-benign)"
GNOMAD_LABEL = "Observed in gnomAD"
NO_ANNOTATION_LABEL = "No ClinVar or gnomAD annotation"
SNV_LABEL = "SNV"
SNV_ACCESSIBLE_LABEL = "SNV-accessible"
VUS_VALUES = frozenset({"Uncertain significance"})
PATHOGENIC_VALUES = frozenset({"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"})
BENIGN_VALUES = frozenset({"Benign", "Likely benign", "Benign/Likely benign"})
PATHOGENIC_OR_BENIGN_VALUES = PATHOGENIC_VALUES | BENIGN_VALUES
CONFLICTING_CLASSIFICATION_VALUE = "Conflicting classifications of pathogenicity"

EXCALIBR_CALIBRATIONS_SHEET = "ExCALIBR_calibrations"
EXCALIBR_DATASET_COL = "dataset"
EXCALIBR_RANGE_COLUMN_PREFIX = "range_"
EXCALIBR_MIXED_YEAR_DATASET_SUFFIX = "_clinvar_2018"

CONTROLS_SHEET_PREFIX = "controls_"
CONTROLS_CLINVAR_GROUP_COL = "clnsig_group_18_25"
EXCALIBR_EVIDENCE_POINTS_COL = "ExC_points_2025"
FUNCTIONAL_CLASS_POINTS_COL = "OP_points"
RECLASSIFICATION_POINTS_COLUMNS = {
    "ExCALIBR evidence": EXCALIBR_EVIDENCE_POINTS_COL,
    "Functional class": FUNCTIONAL_CLASS_POINTS_COL,
}
AGREE_LABEL = "Agrees with ClinVar"
DISAGREE_LABEL = "Disagrees with ClinVar"
NO_EVIDENCE_LABEL = "No point of evidence assigned"
RECLASSIFICATION_LABELS_ORDER = [AGREE_LABEL, DISAGREE_LABEL, NO_EVIDENCE_LABEL]


def split_genes(gene_value):
    return [g.strip() for g in gene_value.split(",") if g.strip()]


def _pipe_parts(value):
    return [p.strip() for p in value.split("|")]


def has_any_value(series):
    """Row-wise: does any '|'-delimited part of this field hold a value?

    See the module docstring's note on pipe-delimited annotation columns.
    """
    return series.apply(lambda v: any(p != "" for p in _pipe_parts(v)))


def matches_any_value(series, targets):
    """Row-wise: does any '|'-delimited part of this field match one of `targets`?"""
    return series.apply(lambda v: any(p in targets for p in _pipe_parts(v)))


def clinvar_significance_flags(clinvar_series):
    """Row-wise pathogenic/benign/VUS/literal-conflict membership (see `matches_any_value`).

    Returns (has_pathogenic, has_benign, has_vus, has_literal_conflict) boolean Series.
    Each of these is OR-able across rows sharing a distinct-variant key (via
    `.groupby(...).any()`) before being passed to `clinvar_classification_from_flags`,
    since "does the union of parts across this group contain a pathogenic value"
    is equivalent to "does any row in the group have a pathogenic value".
    """
    return (
        matches_any_value(clinvar_series, PATHOGENIC_VALUES),
        matches_any_value(clinvar_series, BENIGN_VALUES),
        matches_any_value(clinvar_series, VUS_VALUES),
        matches_any_value(clinvar_series, {CONFLICTING_CLASSIFICATION_VALUE}),
    )


def clinvar_classification_from_flags(has_pathogenic, has_benign, has_vus, has_literal_conflict):
    """Combine pathogenic/benign/VUS/literal-conflict membership into mutually
    exclusive VUS/pathogenic-or-benign/conflict boolean flags.

    A conflict is either ClinVar's own "Conflicting classifications of
    pathogenicity" call, or disagreement between a pathogenic-leaning and
    benign-leaning call -- e.g. across several reverse-translated DNA
    candidates for one protein change, or (once the inputs are grouped via
    `.any()`) across several datasets/measurements of the same variant.
    Mirrors the conflict handling in Analysis/Curation_summary_V5_cleaned.ipynb.

    Returns (vus, pathogenic_or_benign, conflict) boolean Series.
    """
    conflict = has_literal_conflict | (has_pathogenic & has_benign)
    pathogenic_or_benign = ~conflict & (has_pathogenic | has_benign)
    vus = ~conflict & ~pathogenic_or_benign & has_vus
    return vus, pathogenic_or_benign, conflict


def is_snv_accessible(ref_series, alt_series):
    """Row-wise: does any '|'-delimited (transcript_ref, transcript_alt) pair form a SNV?

    A pair is a single-nucleotide substitution if both sides are exactly one
    base. In the expanded (DNA-variant-level) file there's only ever one
    pair, so this is simply "is this variant a SNV". In the condensed
    (assayed, protein-level) file, a row's `hgvs_c` can list several
    underlying DNA-level reverse-translation candidates in lockstep with
    `transcript_ref`/`transcript_alt`; this is True if at least one candidate
    is a SNV ("SNV-accessible").
    """

    def _any_snv(ref_value, alt_value):
        return any(len(ref) == 1 and len(alt) == 1 for ref, alt in zip(_pipe_parts(ref_value), _pipe_parts(alt_value)))

    return pd.Series([_any_snv(r, a) for r, a in zip(ref_series, alt_series)], index=ref_series.index)


def load_dataset_metadata(metadata_path):
    """Return the Curation sheet indexed by dataset name.

    Raises ValueError if any dataset referenced by the condensed file is
    missing from the metadata, since every downstream stat depends on the
    IGVF/measurement classification being complete.
    """
    metadata = pd.read_excel(metadata_path, sheet_name=METADATA_SHEET)
    metadata = metadata.set_index(METADATA_DATASET_COL)
    if metadata.index.has_duplicates:
        dupes = sorted(set(metadata.index[metadata.index.duplicated()]))
        raise ValueError(f"Duplicate dataset name(s) in {METADATA_SHEET} sheet: {dupes}")
    return metadata


def merge_calm_gene_names(genes):
    """Collapse CALM1/CALM2/CALM3 -- which the underlying MAVE data treats as
    three separate gene labels but which correspond to a single gene target
    (they encode the same calmodulin protein) -- into one `CALM_MERGED_LABEL`
    entry in `genes`.
    """
    if genes & CALM_GENES:
        genes = (genes - CALM_GENES) | {CALM_MERGED_LABEL}
    return genes


def genes_in(df, merge_calm_genes=False):
    """Return the set of genes represented in `df`.

    If `merge_calm_genes` is set, CALM1/CALM2/CALM3 are collapsed into one
    `CALM_MERGED_LABEL` entry -- see `merge_calm_gene_names`.
    """
    genes = set()
    for value in df[GENE_COL].unique():
        genes.update(split_genes(value))
    if merge_calm_genes:
        genes = merge_calm_gene_names(genes)
    return genes


def compute_bucket_stats(condensed, dataset_names, measurement_datasets, merge_calm_genes=False):
    """Compute the summary stats for one bucket of datasets."""
    sub = condensed[condensed[DATASET_COL].isin(dataset_names)]
    is_measurement_row = sub[DATASET_COL].isin(measurement_datasets)
    n_measurements = int(is_measurement_row.sum())
    n_composite = int((~is_measurement_row).sum())
    n_variants = sub[VARIANT_KEY_COLS].drop_duplicates().shape[0]
    genes = genes_in(sub, merge_calm_genes=merge_calm_genes)
    return {
        "datasets": sub[DATASET_COL].nunique(),
        "variant_effect_measurements": n_measurements,
        "composite_scores": n_composite,
        "distinct_variants_assayed": n_variants,
        "genes_represented": len(genes),
    }, genes


def compute_all_stats_from_frame(condensed, metadata, merge_calm_genes=False):
    """Compute all bucket stats given an already-loaded condensed frame and metadata."""
    condensed_datasets = set(condensed[DATASET_COL].unique())
    missing = condensed_datasets - set(metadata.index)
    if missing:
        raise ValueError(f"Dataset(s) in condensed file missing from {METADATA_SHEET} metadata: {sorted(missing)}")

    is_igvf = metadata[IGVF_PRODUCED_COL].eq("Yes")
    is_measurement = metadata[SCORE_SET_TYPE_COL].eq(MEASUREMENT_VALUE)
    measurement_datasets = set(metadata.index[is_measurement])

    igvf_datasets = set(metadata.index[is_igvf])
    non_igvf_datasets = set(metadata.index[~is_igvf])
    all_datasets = igvf_datasets | non_igvf_datasets

    igvf_stats, igvf_genes = compute_bucket_stats(condensed, igvf_datasets, measurement_datasets, merge_calm_genes)
    non_igvf_stats, non_igvf_genes = compute_bucket_stats(
        condensed, non_igvf_datasets, measurement_datasets, merge_calm_genes
    )
    combined_stats, _ = compute_bucket_stats(condensed, all_datasets, measurement_datasets, merge_calm_genes)

    non_igvf_stats["genes_not_in_igvf_data"] = len(non_igvf_genes - igvf_genes)

    combined_measurements = combined_stats["variant_effect_measurements"]
    for bucket_stats in (igvf_stats, non_igvf_stats, combined_stats):
        bucket_stats["pct_variant_effect_measurements"] = (
            100 * bucket_stats["variant_effect_measurements"] / combined_measurements
            if combined_measurements
            else float("nan")
        )

    stats = {
        "IGVF": igvf_stats,
        "Community (non-IGVF)": non_igvf_stats,
        "Combined (IGVF + community)": combined_stats,
    }
    gene_breakdown = compute_gene_breakdown(igvf_genes, non_igvf_genes)
    return stats, gene_breakdown


def compute_gene_breakdown(igvf_genes, non_igvf_genes):
    """Partition genes represented in the dataset by whether they're covered
    by IGVF-produced datasets, non-IGVF ("community") datasets, or both.

    CALM1/CALM2/CALM3 are always collapsed into one `CALM_MERGED_LABEL` entry
    here (see `merge_calm_gene_names`), regardless of `--merge-calm-genes` --
    unlike the `genes_represented`/`genes_not_in_igvf_data` counts elsewhere
    in the report, this gene-name list is for readability, and listing the
    same calmodulin target three times under three different names isn't
    useful. This can make this list's counts differ by up to two from those
    counts when `--merge-calm-genes` isn't passed.

    Returns {label: sorted gene list}.
    """
    igvf_genes = merge_calm_gene_names(igvf_genes)
    non_igvf_genes = merge_calm_gene_names(non_igvf_genes)
    return {
        "IGVF only": sorted(igvf_genes - non_igvf_genes),
        "Community (non-IGVF) only": sorted(non_igvf_genes - igvf_genes),
        "Both IGVF and community (non-IGVF)": sorted(igvf_genes & non_igvf_genes),
    }


def format_gene_breakdown(gene_breakdown):
    lines = ["=== Genes represented ==="]
    for label, genes in gene_breakdown.items():
        lines.append(f"{label} ({len(genes)}): {', '.join(genes)}")
    return "\n".join(lines)


def compute_all_stats(condensed_path, metadata_path, merge_calm_genes=False):
    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)
    return compute_all_stats_from_frame(condensed, metadata, merge_calm_genes=merge_calm_genes)


def stats_to_dataframe(stats):
    table = pd.DataFrame(stats).T.reindex(
        columns=[
            "datasets",
            "variant_effect_measurements",
            "pct_variant_effect_measurements",
            "composite_scores",
            "distinct_variants_assayed",
            "genes_represented",
            "genes_not_in_igvf_data",
        ]
    )
    count_columns = [c for c in table.columns if c != "pct_variant_effect_measurements"]
    table[count_columns] = table[count_columns].astype("Int64")
    table["pct_variant_effect_measurements"] = table["pct_variant_effect_measurements"].astype(float).round(1)
    return table


def mixed_year_clinvar_series(df):
    """Per-row ClinVar significance: ClinVar 2018 for BRCA1/PTEN/MSH2/TP53, ClinVar 2025 otherwise.

    A row's `Gene` value is assumed never to mix a `MIXED_YEAR_GENES` gene
    with a non-member gene -- true today, since the only multi-gene `Gene`
    values are the CALM1/CALM2/CALM3 combination, which isn't one of the
    four.
    """
    use_2018 = df[GENE_COL].apply(lambda v: bool(MIXED_YEAR_GENES & set(split_genes(v))))
    return df[CLINVAR_COL_2018].where(use_2018, df[CLINVAR_COL])


def non_clinvar_row_flags(df, snv_label):
    """Per-row boolean flags that don't depend on ClinVar significance: score
    coverage, gnomAD observation, and the SNV(-accessible) flag.
    """
    flags = pd.DataFrame(index=df.index)
    for label, col in SCORE_COLUMNS.items():
        flags[label] = has_any_value(df[col])
    flags[GNOMAD_LABEL] = has_any_value(df[GNOMAD_COL])
    flags[snv_label] = is_snv_accessible(df[TRANSCRIPT_REF_COL], df[TRANSCRIPT_ALT_COL])
    return flags


def variant_flags(df, snv_label, clinvar_series=None, allow_clinvar_conflicts=False):
    """Per-row boolean flags for score coverage and clinical attributes.

    `snv_label` is the column name to store the SNV(-accessible) flag under
    -- `SNV_LABEL` for the DNA-level file, `SNV_ACCESSIBLE_LABEL` for the
    assayed (protein-level) file. See `is_snv_accessible`.

    `clinvar_series` overrides which ClinVar significance values are used for
    the VUS/pathogenic-or-benign/no-annotation flags -- defaults to
    `df[CLINVAR_COL]` (ClinVar 2025 for every gene). Pass
    `mixed_year_clinvar_series(df)` to use ClinVar 2018 for BRCA1/PTEN/MSH2/TP53
    instead.

    By default (`allow_clinvar_conflicts=False`), a row whose own
    pipe-delimited ClinVar values disagree (or literally report "Conflicting
    classifications of pathogenicity") is flagged `CLINVAR_CONFLICT_LABEL`
    instead of VUS/pathogenic-or-benign. Pass `allow_clinvar_conflicts=True`
    to fold such rows into VUS/pathogenic-or-benign via any-match instead
    (the original behavior), which also drops the conflict column entirely.
    """
    if clinvar_series is None:
        clinvar_series = df[CLINVAR_COL]
    flags = non_clinvar_row_flags(df, snv_label)
    if allow_clinvar_conflicts:
        flags[VUS_LABEL] = matches_any_value(clinvar_series, VUS_VALUES)
        flags[PATHOGENIC_OR_BENIGN_LABEL] = matches_any_value(clinvar_series, PATHOGENIC_OR_BENIGN_VALUES)
    else:
        vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(*clinvar_significance_flags(clinvar_series))
        flags[VUS_LABEL] = vus
        flags[PATHOGENIC_OR_BENIGN_LABEL] = pathogenic_or_benign
        flags[CLINVAR_CONFLICT_LABEL] = conflict
    flags[NO_ANNOTATION_LABEL] = ~has_any_value(clinvar_series) & ~flags[GNOMAD_LABEL]
    return flags


def distinct_variant_flags(df, snv_label, clinvar_series=None, allow_clinvar_conflicts=False):
    """Collapse per-row flags to one row per distinct (hgvs_c, hgvs_p) variant.

    Score coverage, gnomAD observation, and no-annotation flags count as
    True for a distinct variant if any of its (possibly several, across
    datasets and DNA-level candidates) measurement rows does. See
    `variant_flags` for `clinvar_series` and `allow_clinvar_conflicts`;
    when conflicts aren't allowed, the pathogenic/benign/VUS/literal-conflict
    membership used to classify VUS/pathogenic-or-benign/conflict (see
    `clinvar_classification_from_flags`) is OR'd across *all* of the
    variant's rows first, so e.g. a variant called Pathogenic by one
    measurement and Benign by another is flagged a conflict rather than
    both/either.
    """
    if clinvar_series is None:
        clinvar_series = df[CLINVAR_COL]

    non_clinvar = non_clinvar_row_flags(df, snv_label)
    keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), non_clinvar.reset_index(drop=True)], axis=1)
    distinct = keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()

    has_clinvar_value = has_any_value(clinvar_series).reset_index(drop=True)
    hcv_keyed = pd.concat(
        [df[VARIANT_KEY_COLS].reset_index(drop=True), has_clinvar_value.rename("_has_clinvar_value")], axis=1
    )
    hcv_distinct = hcv_keyed.groupby(VARIANT_KEY_COLS, as_index=False)["_has_clinvar_value"].any()
    distinct = distinct.merge(hcv_distinct, on=VARIANT_KEY_COLS)

    if allow_clinvar_conflicts:
        clinvar_flags = pd.DataFrame(
            {
                VUS_LABEL: matches_any_value(clinvar_series, VUS_VALUES),
                PATHOGENIC_OR_BENIGN_LABEL: matches_any_value(clinvar_series, PATHOGENIC_OR_BENIGN_VALUES),
            }
        ).reset_index(drop=True)
        clinvar_keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), clinvar_flags], axis=1)
        clinvar_distinct = clinvar_keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()
        distinct = distinct.merge(clinvar_distinct, on=VARIANT_KEY_COLS)
    else:
        has_pathogenic, has_benign, has_vus, has_literal_conflict = clinvar_significance_flags(clinvar_series)
        significance_flags = pd.DataFrame(
            {
                "_has_pathogenic": has_pathogenic,
                "_has_benign": has_benign,
                "_has_vus": has_vus,
                "_has_literal_conflict": has_literal_conflict,
            }
        ).reset_index(drop=True)
        significance_keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), significance_flags], axis=1)
        significance_distinct = significance_keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()
        distinct = distinct.merge(significance_distinct, on=VARIANT_KEY_COLS)

        vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(
            distinct["_has_pathogenic"], distinct["_has_benign"], distinct["_has_vus"], distinct["_has_literal_conflict"]
        )
        distinct[VUS_LABEL] = vus
        distinct[PATHOGENIC_OR_BENIGN_LABEL] = pathogenic_or_benign
        distinct[CLINVAR_CONFLICT_LABEL] = conflict
        distinct = distinct.drop(columns=["_has_pathogenic", "_has_benign", "_has_vus", "_has_literal_conflict"])

    distinct[NO_ANNOTATION_LABEL] = ~distinct["_has_clinvar_value"] & ~distinct[GNOMAD_LABEL]
    return distinct.drop(columns=VARIANT_KEY_COLS + ["_has_clinvar_value"])


def summarize_flags(flags):
    """Return (total, DataFrame[count, pct]) for a set of boolean flag columns."""
    total = len(flags)
    counts = flags.sum().astype(int)
    pct = (100 * counts / total).round(1) if total else counts.astype(float)
    return total, pd.DataFrame({"count": counts, "pct": pct})


def format_count_table(title, total, table):
    lines = [title, f"Total: {total}"]
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


def summarize_clinical_flags(flags, snv_label):
    """Like `summarize_flags(flags[clinical_labels])`, plus two columns breaking out
    how many of the level's SNV(-accessible) variants fall into each clinical-attribute
    bucket, as a percentage of that level's total SNV(-accessible) count.

    `clinical_labels` is `CLINICAL_LABELS_ORDER` filtered to columns actually present in
    `flags` -- `CLINVAR_CONFLICT_LABEL` is only there when `variant_flags`/
    `distinct_variant_flags` were called with `allow_clinvar_conflicts=False`.

    Returns (total, snv_total, table).
    """
    clinical_labels = [label for label in CLINICAL_LABELS_ORDER if label in flags.columns]
    total, table = summarize_flags(flags[clinical_labels])
    snv_total = int(flags[snv_label].sum())
    snv_counts = flags[clinical_labels].apply(lambda col: int((col & flags[snv_label]).sum()))
    pct_snv_col = f"% of {snv_label}"
    table[snv_label] = snv_counts
    table[pct_snv_col] = (100 * snv_counts / snv_total).round(1) if snv_total else snv_counts.astype(float)
    return total, snv_total, table


def format_clinical_table(title, total, snv_total, table, snv_label):
    lines = [title, f"Total: {total} ({snv_total} {snv_label})"]
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        pct_snv_col = f"% of {snv_label}"
        body[pct_snv_col] = body[pct_snv_col].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


SCORE_LABELS = list(SCORE_COLUMNS.keys())
CLINICAL_LABELS_ORDER = [
    VUS_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    CLINVAR_CONFLICT_LABEL,
    GNOMAD_LABEL,
    NO_ANNOTATION_LABEL,
]


def build_variant_level_reports(condensed, expanded, condensed_path, expanded_path, allow_clinvar_conflicts=False):
    """Build the score-coverage and clinical-attribute text sections.

    Each is reported at four levels: assayed (protein-resolution, from the
    condensed file) vs. DNA-level (from the expanded file) variants, and
    distinct variants vs. variant measurements (rows). The clinical-attribute
    breakdown is additionally reported twice per level: once using ClinVar
    2025 throughout, and once substituting ClinVar 2018 for BRCA1, PTEN,
    MSH2, and TP53 (see `mixed_year_clinvar_series`). See `variant_flags` for
    `allow_clinvar_conflicts`.
    """
    levels = [
        ("assayed variants, distinct", condensed, condensed_path, SNV_ACCESSIBLE_LABEL, distinct_variant_flags),
        ("assayed variant measurements", condensed, condensed_path, SNV_ACCESSIBLE_LABEL, variant_flags),
        ("DNA variants, distinct", expanded, expanded_path, SNV_LABEL, distinct_variant_flags),
        ("DNA variant measurements", expanded, expanded_path, SNV_LABEL, variant_flags),
    ]

    score_sections = []
    clinical_sections = []
    clinical_sections_mixed_year = []
    for label, df, source, snv_label, flags_fn in levels:
        flags = flags_fn(df, snv_label, allow_clinvar_conflicts=allow_clinvar_conflicts)
        total, table = summarize_flags(flags[SCORE_LABELS])
        score_sections.append(format_count_table(f"Score coverage -- {label} (from {source})", total, table))

        clinical_total, snv_total, clinical_table = summarize_clinical_flags(flags, snv_label)
        clinical_sections.append(
            format_clinical_table(
                f"Clinical attributes -- {label} (from {source})", clinical_total, snv_total, clinical_table, snv_label
            )
        )

        mixed_flags = flags_fn(
            df, snv_label, clinvar_series=mixed_year_clinvar_series(df), allow_clinvar_conflicts=allow_clinvar_conflicts
        )
        mixed_total, mixed_snv_total, mixed_table = summarize_clinical_flags(mixed_flags, snv_label)
        clinical_sections_mixed_year.append(
            format_clinical_table(
                f"Clinical attributes -- {label} (from {source})",
                mixed_total,
                mixed_snv_total,
                mixed_table,
                snv_label,
            )
        )

    return score_sections, clinical_sections, clinical_sections_mixed_year


def excalibr_dataset_to_gene_map(calibration_datasets, metadata):
    """Map each ExCALIBR_calibrations `dataset` value to its `Gene`, via the
    same Curation metadata (`Dataset Name` -> `Gene`) used elsewhere in this
    script.

    A `dataset` value ending in `_clinvar_2018` is a mixed-year duplicate
    calibration of the same underlying dataset -- that suffix isn't part of
    Supplementary_Data_3's `Dataset Name`, so it's stripped before lookup.
    Both sides are also Unicode-normalized (NFC) before comparison, since
    accented characters (e.g. in author names) can round-trip through Excel
    in different composed/decomposed forms across files.

    Raises ValueError listing any `dataset` value that still doesn't resolve
    to a `Dataset Name` in the metadata after these adjustments.
    """
    normalized_index = {unicodedata.normalize("NFC", name): name for name in metadata.index}

    def resolve(dataset):
        stripped = dataset.removesuffix(EXCALIBR_MIXED_YEAR_DATASET_SUFFIX)
        return normalized_index.get(unicodedata.normalize("NFC", stripped))

    mapping = {}
    unmapped = []
    for dataset in calibration_datasets:
        resolved = resolve(dataset)
        if resolved is None:
            unmapped.append(dataset)
        else:
            mapping[dataset] = metadata.loc[resolved, GENE_COL]
    if unmapped:
        raise ValueError(
            f"{EXCALIBR_CALIBRATIONS_SHEET} dataset(s) not found in {METADATA_SHEET} metadata: {sorted(unmapped)}"
        )
    return mapping


def compute_excalibr_calibration_stats(calibrations, metadata, merge_calm_genes=False):
    """How many genes have an ExCALIBR calibration, and how many of those
    genes have at least one dataset where ExCALIBR assigned at least one
    point of evidence (pathogenic or benign) -- i.e. at least one non-null
    `range_-8`..`range_8` value in that dataset's row.

    A gene can have several calibration rows (one per dataset); this counts
    distinct genes, not rows.
    """
    dataset_to_gene = excalibr_dataset_to_gene_map(calibrations[EXCALIBR_DATASET_COL], metadata)
    range_cols = [c for c in calibrations.columns if c.startswith(EXCALIBR_RANGE_COLUMN_PREFIX)]
    row_has_evidence = calibrations[range_cols].notna().any(axis=1)
    has_evidence_by_dataset = row_has_evidence.groupby(calibrations[EXCALIBR_DATASET_COL]).any()

    genes_with_calibration = set()
    genes_with_evidence = set()
    for dataset, gene_value in dataset_to_gene.items():
        genes = set(split_genes(gene_value))
        if merge_calm_genes and genes & CALM_GENES:
            genes -= CALM_GENES
            genes.add(CALM_MERGED_LABEL)
        genes_with_calibration |= genes
        if has_evidence_by_dataset.get(dataset, False):
            genes_with_evidence |= genes

    return {
        "genes_with_excalibr_calibrations": len(genes_with_calibration),
        "genes_with_evidence_assigned": len(genes_with_evidence),
    }


def format_calibration_summary(stats):
    total = stats["genes_with_excalibr_calibrations"]
    with_evidence = stats["genes_with_evidence_assigned"]
    pct = 100 * with_evidence / total if total else float("nan")
    return "\n".join(
        [
            "=== ExCALIBR calibration coverage (Extended Data Figure 2) ===",
            f"Genes with ExCALIBR calibrations: {total}",
            f"Genes with >=1 dataset assigning >=1 point of evidence (pathogenic or benign): {with_evidence} ({pct:.1f}%)",
        ]
    )


def reclassification_flags(clinvar_group, points):
    """Per-row agree/disagree/no-evidence flags for one points column against
    the row's ClinVar pathogenic-or-benign control label.

    A row is in scope only if `clinvar_group` is a pathogenic- or
    benign-leaning control label (see `PATHOGENIC_VALUES`/`BENIGN_VALUES`);
    `points` is treated as pathogenic evidence if positive, benign evidence
    if negative, and no evidence if zero or missing. Returns (flags, in_scope),
    both boolean Series/DataFrame aligned to `clinvar_group`'s index.
    """
    is_pathogenic = clinvar_group.isin(PATHOGENIC_VALUES)
    is_benign = clinvar_group.isin(BENIGN_VALUES)
    in_scope = is_pathogenic | is_benign

    assigned_pathogenic = points > 0
    assigned_benign = points < 0
    no_evidence = ~assigned_pathogenic & ~assigned_benign

    agree = (is_pathogenic & assigned_pathogenic) | (is_benign & assigned_benign)
    disagree = (is_pathogenic & assigned_benign) | (is_benign & assigned_pathogenic)

    flags = pd.DataFrame({AGREE_LABEL: agree, DISAGREE_LABEL: disagree, NO_EVIDENCE_LABEL: no_evidence})
    return flags, in_scope


def compute_reclassification_agreement(controls_df):
    """For each points column in `RECLASSIFICATION_POINTS_COLUMNS`, the
    agree/disagree/no-evidence breakdown (see `reclassification_flags`) plus
    the agreement percentage among rows with evidence assigned.

    Returns {label: (total, determinate, agreement_pct, table)}, where `table`
    is the `summarize_flags` output over the in-scope rows.
    """
    results = {}
    for label, points_col in RECLASSIFICATION_POINTS_COLUMNS.items():
        flags, in_scope = reclassification_flags(controls_df[CONTROLS_CLINVAR_GROUP_COL], controls_df[points_col])
        flags = flags[in_scope]
        total, table = summarize_flags(flags[RECLASSIFICATION_LABELS_ORDER])
        determinate = total - int(flags[NO_EVIDENCE_LABEL].sum())
        agreement_pct = 100 * int(flags[AGREE_LABEL].sum()) / determinate if determinate else float("nan")
        results[label] = (total, determinate, agreement_pct, table)
    return results


def format_reclassification_table(title, total, determinate, agreement_pct, table):
    lines = [title, f"Total control variants: {total}", f"Determinate calls (evidence assigned): {determinate}"]
    if determinate:
        lines.append(f"Agreement with ClinVar PLP/BLB (of determinate calls): {agreement_pct:.1f}%")
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


def build_reclassification_report(controls_path):
    """One section per (`controls_`-prefixed sheet, points column) pair in
    the controls file. See the module docstring's "Reclassification
    agreement" note for why every such sheet is expected to agree.
    """
    workbook = pd.ExcelFile(controls_path)
    sheets = [name for name in workbook.sheet_names if name.startswith(CONTROLS_SHEET_PREFIX)]
    sections = []
    for sheet in sheets:
        controls_df = workbook.parse(sheet)
        agreement = compute_reclassification_agreement(controls_df)
        for label, (total, determinate, agreement_pct, table) in agreement.items():
            sections.append(
                format_reclassification_table(f"{label} -- {sheet}", total, determinate, agreement_pct, table)
            )
    return sections


def build_report_text(
    table,
    gene_breakdown,
    score_sections,
    clinical_sections,
    clinical_sections_mixed_year,
    calibration_summary,
    reclassification_sections,
    allow_clinvar_conflicts=False,
):
    conflict_note = (
        "conflicting/ambiguous ClinVar calls folded in via any-match"
        if allow_clinvar_conflicts
        else "conflicting/ambiguous ClinVar calls excluded"
    )
    parts = [
        "=== Dataset summary ===",
        table.to_string(),
        format_gene_breakdown(gene_breakdown),
        "=== Score coverage (REVEL, AlphaMissense, MutPred2) ===",
        *score_sections,
        f"=== Clinical attributes (ClinVar 2025, gnomAD; {conflict_note}) ===",
        *clinical_sections,
        f"=== Clinical attributes (ClinVar 2025, except ClinVar 2018 for BRCA1/PTEN/MSH2/TP53; gnomAD; {conflict_note}) ===",
        *clinical_sections_mixed_year,
        calibration_summary,
        "=== Reclassification agreement (Figure 4c) ===",
        *reclassification_sections,
    ]
    return "\n\n".join(parts)


@click.command(help=__doc__)
@click.argument(
    "condensed_file",
    required=False,
    default=DEFAULT_CONDENSED_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "metadata_file",
    required=False,
    default=DEFAULT_METADATA_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "expanded_file",
    required=False,
    default=DEFAULT_EXPANDED_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--excalibr-calibrations-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_EXCALIBR_CALIBRATIONS_FILE,
    help=(
        f"Path to the ExCALIBR calibrations workbook (default {DEFAULT_EXCALIBR_CALIBRATIONS_FILE}), "
        f"read from its '{EXCALIBR_CALIBRATIONS_SHEET}' sheet, for the calibration-coverage section."
    ),
)
@click.option(
    "--controls-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CONTROLS_FILE,
    help=(
        f"Path to the controls workbook (default {DEFAULT_CONTROLS_FILE}), read from every sheet "
        f"whose name starts with '{CONTROLS_SHEET_PREFIX}', for the reclassification-agreement section."
    ),
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional path to also write the full report as a text file",
)
@click.option(
    "--merge-calm-genes",
    is_flag=True,
    default=False,
    help=(
        "Count CALM1, CALM2, and CALM3 as a single gene target ('CALM1/2/3') "
        "instead of three separate genes, since they encode the same calmodulin "
        "protein in the underlying data."
    ),
)
@click.option(
    "--allow-clinvar-conflicts",
    is_flag=True,
    default=False,
    help=(
        "Fold variants with a conflicting or ambiguous ClinVar call -- "
        "disagreement between a pathogenic-leaning and benign-leaning call "
        "across a variant's measurements/DNA candidates, or ClinVar's own "
        "'Conflicting classifications of pathogenicity' call -- into the "
        "VUS/pathogenic-or-benign clinical-attribute buckets via any-match. "
        "By default such variants are excluded from both buckets into their "
        "own 'ClinVar conflict' bucket instead, mirroring the conflict "
        "handling in Analysis/Curation_summary_V5_cleaned.ipynb."
    ),
)
def main(
    condensed_file,
    metadata_file,
    expanded_file,
    excalibr_calibrations_file,
    controls_file,
    output,
    merge_calm_genes,
    allow_clinvar_conflicts,
):
    condensed = pd.read_csv(condensed_file, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_file)

    try:
        stats, gene_breakdown = compute_all_stats_from_frame(condensed, metadata, merge_calm_genes=merge_calm_genes)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    table = stats_to_dataframe(stats)

    expanded = pd.read_csv(expanded_file, sep="\t", dtype=str, keep_default_na=False)
    score_sections, clinical_sections, clinical_sections_mixed_year = build_variant_level_reports(
        condensed, expanded, condensed_file, expanded_file, allow_clinvar_conflicts=allow_clinvar_conflicts
    )

    calibrations = pd.read_excel(excalibr_calibrations_file, sheet_name=EXCALIBR_CALIBRATIONS_SHEET)
    try:
        calibration_stats = compute_excalibr_calibration_stats(
            calibrations, metadata, merge_calm_genes=merge_calm_genes
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    calibration_summary = format_calibration_summary(calibration_stats)

    reclassification_sections = build_reclassification_report(controls_file)

    report = build_report_text(
        table,
        gene_breakdown,
        score_sections,
        clinical_sections,
        clinical_sections_mixed_year,
        calibration_summary,
        reclassification_sections,
        allow_clinvar_conflicts=allow_clinvar_conflicts,
    )
    click.echo(report)

    if output:
        output.write_text(report + "\n")
        click.echo(f"\nWrote report to {output}")


if __name__ == "__main__":
    main()

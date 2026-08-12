"""Shared deduplication mechanics for `notebooks/analysis/Variant_Classification_analysis.ipynb`
and `notebooks/analysis/OddsPath_classifications.ipynb`. Both pipelines have
to resolve "more than one assay/SNV scored this variant" down to a
representative record, using the same decided-approach tie-breaks -- see
`docs/variant_classification.md` for the full rationale. This module exists
because that logic used to be copy-pasted independently into each notebook,
and `OddsPath_classifications.ipynb`'s copy silently drifted out of sync
with fixes made to `Variant_Classification_analysis.ipynb` (the signed-value
`catch_mis_2` bug, the arbitrary aa-stage tie-break, the non-deterministic
`VUS`/`gnomAD`/`Unobserved` tie-break) -- sharing one implementation is meant
to prevent that happening again.

`"v1"` is preserved byte-for-byte in every function below, for historical
audit/reproducibility purposes -- it reproduces each pipeline's original,
pre-parameter behavior exactly, including known issues (the signed-value
sort in `catch_mis_2`, the arbitrary-tie-break `VariantNotes` alphabetical
sort). It is not the recommended setting for either parameter going
forward.
"""
import numpy as np
import pandas as pd

GENOMIC_KEY_COLS = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele"]

# Rank of ClinVar review-status text (e.g. `clinvar_star_18_25`). Mirrors
# ClinVar's own star-rating tiers (higher = more authoritative);
# "practice guideline" isn't observed in the data behind either pipeline but
# is included defensively. Anything unmapped (including NaN) ranks with the
# lowest tier via `.fillna(0)` in `controls_aa_sort_key`.
CLINVAR_REVIEW_STATUS_RANK = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "no classification for the single variant": 0,
    "no classification provided": 0,
    "no assertion criteria provided": 0,
}


def controls_aa_sort_key(df, strategy):
    """aa-stage tie-break sort key for `controls` (which distinct NT variant
    represents a shared amino-acid substitution). Mutates `df` in place to
    add the columns the returned sort key references, and returns
    `(sort_by, ascending)` for use with `aa_dedup_or_mark`.

    `"v1"` sorts by `assay_priority` alone (caller must have already
    computed it via `ASSAY_PRIORITY_LIST`). Otherwise, prefers the higher-
    quality ClinVar record (`clinvar_star_18_25`), then most recent review
    date (`clinvar_date_last_reviewed_18_25`, ClinVar's `"Mon DD, YYYY"`
    text -- parsed here, not ISO), then greatest `abs(Fxn_points)` and
    `Dataset` name as deterministic fallbacks (inert today since candidates
    reaching this tie-break are already magnitude-tied by construction, but
    closing off any dependency on incidental row order for a remaining
    tie). Requires `df` to already have `Fxn_points`, `Dataset`, and (for
    non-`"v1"`) `clinvar_star_18_25`/`clinvar_date_last_reviewed_18_25`.
    """
    if strategy == "v1":
        return "assay_priority", True
    df["abs_Fxn_points"] = df["Fxn_points"].abs()
    df["clinvar_review_rank"] = df["clinvar_star_18_25"].map(CLINVAR_REVIEW_STATUS_RANK).fillna(0)
    df["clinvar_review_date"] = pd.to_datetime(
        df["clinvar_date_last_reviewed_18_25"], format="%b %d, %Y", errors="coerce"
    )
    return ["clinvar_review_rank", "clinvar_review_date", "abs_Fxn_points", "Dataset"], [False, False, False, True]


def clingen_aa_sort_key(df, strategy):
    """aa-stage tie-break sort key for `ClinGen_Repo`. Mutates `df` in place,
    returns `(sort_by, ascending)` for use with `aa_dedup_or_mark`.

    Every `ClinGen_Repo` row is already a VCEP (expert-panel) assertion --
    there's no review-status tier to rank the way `controls_aa_sort_key`
    does -- so the tie-break is: prefer non-retracted, then most recent
    `Approval Date_ClinGen_repo`, then `Published Date_ClinGen_repo`, then
    `abs(Fxn_points)` and `Dataset` name as deterministic fallbacks.
    Requires `df` to already have `Fxn_points`, `Dataset`, and (for
    non-`"v1"`) `Retracted_ClinGen_repo`/`Approval Date_ClinGen_repo`/
    `Published Date_ClinGen_repo`.
    """
    if strategy == "v1":
        return "assay_priority", True
    df["abs_Fxn_points"] = df["Fxn_points"].abs()
    df["_not_retracted"] = df["Retracted_ClinGen_repo"].fillna(0) != 1
    return (
        ["_not_retracted", "Approval Date_ClinGen_repo", "Published Date_ClinGen_repo", "abs_Fxn_points", "Dataset"],
        [False, False, False, False, True],
    )


def aa_dedup_or_mark(candidates, group_cols, genomic_key_cols, sort_by, ascending, strategy):
    """
    aa-stage tie-break for one predictor's candidate rows in controls_aa/
    clingen_aa (already pre-filtered to rows tied at their group's max
    abs(Fxn_points) and max predictor points) -- resolves which SNV
    represents a shared amino-acid change across `group_cols`.

    Two stages, in order:

    1. Collapse rows that share the exact same genomic coordinates
       (`genomic_key_cols`) to a single row -- i.e. more than one
       dataset/assay scoring the *identical* physical NT variant. This
       always fully dedups, regardless of strategy: it's a different
       question (which assay's measurement of *one* variant to trust) from
       the amino-acid-level question stage 2 resolves, and was always
       meant to produce exactly one row per NT variant. Without this
       stage, two datasets scoring the same SNV would otherwise be treated
       as if they were two competing SNVs and both kept under "abs_max"/
       "nt_then_abs_max" -- see docs/variant_classification.md. Compares
       `genomic_key_cols` as strings rather than relying on
       `drop_duplicates` directly, since `Chrom` mixes int and str values
       for the same chromosome across rows in practice.
    2. Resolve which *distinct* NT variant represents the shared
       amino-acid change across `group_cols`. "v1" reproduces the original
       behavior byte-for-byte: keeps only the winner via
       drop_duplicates(keep="first") -- the losing SNV's row is not
       recoverable under "v1". "abs_max"/"nt_then_abs_max" keep every
       distinct-NT-variant candidate: the first row per group is marked
       "primary", every other one "secondary", via a new Variant_Role
       column, instead of being dropped.
    """
    sorted_candidates = candidates.sort_values(by=sort_by, ascending=ascending, na_position="last")
    genomic_key = sorted_candidates[genomic_key_cols].astype(str)
    one_row_per_nt_variant = sorted_candidates[~genomic_key.duplicated(keep="first")]
    if strategy == "v1":
        return one_row_per_nt_variant.drop_duplicates(subset=group_cols, keep="first")
    one_row_per_nt_variant = one_row_per_nt_variant.copy()
    is_primary = ~one_row_per_nt_variant.duplicated(subset=group_cols, keep="first")
    one_row_per_nt_variant["Variant_Role"] = np.where(is_primary, "primary", "secondary")
    return one_row_per_nt_variant


def catch_mis_2(df, group_cols, points_col="Fxn_points", strategy="v1"):
    """
    Handle a variant scored by more than one assay (after the nt/aa-subset
    dedup above already resolved same-type duplicates) by keeping a single
    representative row per `group_cols`. `strategy` controls which row wins
    when an nt-type survivor and an aa-type survivor collide for the same
    variant:
      - "v1": signed Fxn_points, descending (the original behavior -- a
        positive value always beats a negative one, and between two
        negatives the one closer to zero wins)
      - "abs_max": greatest absolute Fxn_points wins
      - "nt_then_abs_max": an nt-type row always wins over an aa-type row;
        ties within a type broken by greatest absolute Fxn_points
    """
    group_cols = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele"]
    df = df.copy()

    if strategy == "v1":
        sort_by, ascending = points_col, False
    elif strategy == "abs_max":
        df["_abs_points"] = df[points_col].abs()
        sort_by, ascending = "_abs_points", False
    elif strategy == "nt_then_abs_max":
        df["_is_aa"] = df["nucleotide_or_aa"] == "aa"
        df["_abs_points"] = df[points_col].abs()
        sort_by, ascending = ["_is_aa", "_abs_points"], [True, False]
    else:
        raise ValueError(f"Unknown dedup strategy: {strategy}")

    df_sorted = df.sort_values(by=sort_by, ascending=ascending, na_position="last")
    cleaned = df_sorted.drop_duplicates(subset=group_cols, keep="first")
    return cleaned.drop(columns=["_abs_points", "_is_aa"], errors="ignore")


def dedup_by_max_abs_points(df, points_col, genomic_key_cols=GENOMIC_KEY_COLS):
    """Collapse `df` to one row per DNA variant (`genomic_key_cols`), keeping
    the candidate with the greatest `abs(points_col)`; ties broken by
    `Dataset` name (ascending) for full determinism. Compares
    `genomic_key_cols` as strings rather than relying on `drop_duplicates`
    directly, since `Chrom` mixes int and str values for the same
    chromosome across rows in practice (see `aa_dedup_or_mark`).

    Unlike `catch_mis_2`/`dedup_vus_gnomad_unobserved`, this has no nt/aa
    resolution preference and no `"v1"` legacy mode -- it's used to build a
    single DNA-level dataset from scratch, not to reproduce an existing
    per-category pipeline's history.
    """
    df = df.copy()
    df["_abs_points"] = df[points_col].abs()
    sorted_df = df.sort_values(by=["_abs_points", "Dataset"], ascending=[False, True], na_position="last")
    genomic_key = sorted_df[genomic_key_cols].astype(str)
    deduped = sorted_df[~genomic_key.duplicated(keep="first")]
    return deduped.drop(columns=["_abs_points"])


def dedup_vus_gnomad_unobserved(df, group_cols, points_col="Fxn_points", strategy="v1", variant_notes_col="VariantNotes"):
    """
    Resolve a variant scored by more than one assay to a single representative
    row for VUS/gnomAD/Unobserved. Unlike `controls`/`ClinGen_Repo`, these
    categories are evaluated at DNA resolution only, so nt- and aa-type rows
    are deduped together in one pass (no nt/aa split, no double-counting
    concern). `strategy`:
      - "v1": `variant_notes_col` tag order (the original behavior -- in
        practice equivalent to greatest absolute Fxn_points within a single
        assay type, but an nt-type row always beats an aa-type row
        regardless of magnitude whenever both cover the same variant, since
        `First_max_fxn_pts` (nt) sorts ahead of `max_fxn_pts` (aa)
        alphabetically -- an accidental bias, see
        docs/assay_priority_questions.md)
      - "abs_max": greatest absolute Fxn_points wins, nt and aa candidates
        treated identically
      - "nt_then_abs_max": an nt-type row always wins over an aa-type row;
        ties within a type broken by greatest absolute Fxn_points

    "abs_max"/"nt_then_abs_max" both append `Dataset` name (ascending) as a
    final deterministic tiebreak, since a genuine tie (same abs(Fxn_points),
    same nt/aa type) would otherwise fall through to `sort_values`' stable
    sort -- i.e. whichever row happened to arrive first in `df` -- which
    isn't a documented policy and isn't robust to unrelated upstream changes
    (e.g. an upstream row-count change shifting everyone else's position).
    See docs/variant_classification.md. "v1" is untouched, preserved
    byte-for-byte for historical audit/reproducibility.

    `variant_notes_col` lets callers use their own tag column name (e.g.
    `OddsPath_classifications.ipynb`'s `VariantNotes_OP`, distinct from
    `Variant_Classification_analysis.ipynb`'s `VariantNotes`, since both
    columns can be present on the same shared checkpoint dataframe).
    """
    df = df.copy()

    if strategy == "v1":
        sort_by, ascending = variant_notes_col, True
    elif strategy == "abs_max":
        df["_abs_points"] = df[points_col].abs()
        sort_by, ascending = ["_abs_points", "Dataset"], [False, True]
    elif strategy == "nt_then_abs_max":
        df["_is_aa"] = df["nucleotide_or_aa"] == "aa"
        df["_abs_points"] = df[points_col].abs()
        sort_by, ascending = ["_is_aa", "_abs_points", "Dataset"], [True, False, True]
    else:
        raise ValueError(f"Unknown dedup strategy: {strategy}")

    df_sorted = df.sort_values(by=sort_by, ascending=ascending, na_position="last")
    cleaned = df_sorted.drop_duplicates(subset=group_cols, keep="first")
    return cleaned.drop(columns=["_abs_points", "_is_aa"], errors="ignore")

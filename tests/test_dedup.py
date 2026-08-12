import pandas as pd
import pytest

from src.lib.dedup import (
    CLINVAR_REVIEW_STATUS_RANK,
    GENOMIC_KEY_COLS,
    aa_dedup_or_mark,
    catch_mis_2,
    clingen_aa_sort_key,
    controls_aa_sort_key,
    dedup_vus_gnomad_unobserved,
)

AA_GROUP_COLS = ["Gene", "aa_pos", "aa_ref", "aa_alt"]


def _controls_aa_frame(rows):
    """`rows`: (Gene, aa_pos, aa_ref, aa_alt, Chrom, hg38_start, ref_allele,
    alt_allele, Dataset, Fxn_points, clinvar_star_18_25, clinvar_date_last_reviewed_18_25).
    """
    return pd.DataFrame(
        rows,
        columns=[
            "Gene", "aa_pos", "aa_ref", "aa_alt", "Chrom", "hg38_start", "ref_allele", "alt_allele",
            "Dataset", "Fxn_points", "clinvar_star_18_25", "clinvar_date_last_reviewed_18_25",
        ],
    )


def _clingen_aa_frame(rows):
    """`rows`: (Gene, aa_pos, aa_ref, aa_alt, Chrom, hg38_start, ref_allele,
    alt_allele, Dataset, Fxn_points, Retracted_ClinGen_repo, Approval Date_ClinGen_repo,
    Published Date_ClinGen_repo).
    """
    return pd.DataFrame(
        rows,
        columns=[
            "Gene", "aa_pos", "aa_ref", "aa_alt", "Chrom", "hg38_start", "ref_allele", "alt_allele",
            "Dataset", "Fxn_points", "Retracted_ClinGen_repo", "Approval Date_ClinGen_repo",
            "Published Date_ClinGen_repo",
        ],
    )


def _dedup_frame(rows):
    """`rows`: (Gene, hg38_start, ref_allele, alt_allele, Dataset, Fxn_points,
    nucleotide_or_aa, VariantNotes).
    """
    return pd.DataFrame(
        rows,
        columns=["Gene", "hg38_start", "ref_allele", "alt_allele", "Dataset", "Fxn_points", "nucleotide_or_aa", "VariantNotes"],
    )


# --- controls_aa_sort_key / aa_dedup_or_mark -------------------------------------------


def test_controls_aa_two_distinct_snvs_both_kept_higher_magnitude_wins():
    df = _controls_aa_frame(
        [
            ("G1", 100, "A", "V", "1", 1000, "G", "A", "Dataset_P", 5, "criteria provided, single submitter", "Jan 01, 2020"),
            ("G1", 100, "A", "V", "1", 1003, "C", "T", "Dataset_Q", 8, "criteria provided, single submitter", "Jan 01, 2020"),
        ]
    )
    sort_by, ascending = controls_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")

    assert len(out) == 2
    assert out.loc[out.Variant_Role == "primary", "Dataset"].iloc[0] == "Dataset_Q"
    assert out.loc[out.Variant_Role == "secondary", "Dataset"].iloc[0] == "Dataset_P"


def test_controls_aa_quality_breaks_tie_ahead_of_older_date():
    df = _controls_aa_frame(
        [
            ("G1", 200, "A", "V", "1", 2000, "G", "A", "Dataset_P", 5, "criteria provided, single submitter", "Jan 01, 2020"),
            ("G1", 200, "A", "V", "1", 2003, "C", "T", "Dataset_Q", 5, "reviewed by expert panel", "Jan 01, 2019"),
        ]
    )
    sort_by, ascending = controls_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")
    assert out.loc[out.Variant_Role == "primary", "Dataset"].iloc[0] == "Dataset_Q"


def test_controls_aa_same_nt_variant_different_datasets_collapses_to_one_row():
    """Regression: two datasets scoring the *identical* physical variant must
    collapse to exactly one row, not be treated as competing SNVs."""
    df = _controls_aa_frame(
        [
            ("PAX6", 10, "Q", "*", "11", 31802817, "G", "A", "PAX6_BLX_geneticin", 8, None, None),
            ("PAX6", 10, "Q", "*", "11", 31802817, "G", "A", "PAX6_LE9_geneticin", 8, None, None),
        ]
    )
    sort_by, ascending = controls_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")
    assert len(out) == 1
    assert out["Variant_Role"].iloc[0] == "primary"


def test_controls_aa_mixed_dtype_chrom_still_collapses():
    """Regression: Chrom mixing str '11' and int 11 for the same chromosome
    (an upstream merge/concat artifact, confirmed live) must not defeat the
    same-NT-variant collapse."""
    df = _controls_aa_frame(
        [
            ("PAX6", 10, "Q", "*", "11", 31802817, "G", "A", "PAX6_BLX_geneticin", 8, None, None),
            ("PAX6", 10, "Q", "*", 11, 31802817, "G", "A", "PAX6_LE9_geneticin", 8, None, None),
        ]
    )
    assert df["Chrom"].dtype == object
    sort_by, ascending = controls_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")
    assert len(out) == 1


def test_controls_aa_v1_strategy_drops_loser_byte_for_byte():
    df = _controls_aa_frame(
        [
            ("G1", 400, "A", "V", "1", 4000, "G", "A", "Dataset_Z", 5, "reviewed by expert panel", "Jan 01, 2020"),
            ("G1", 400, "A", "V", "1", 4003, "C", "T", "Dataset_A", 5, "criteria provided, single submitter", "Jan 01, 2020"),
        ]
    )
    df["assay_priority"] = [1, 0]  # Dataset_A ranks first under v1
    sort_by, ascending = controls_aa_sort_key(df, "v1")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "v1")
    assert len(out) == 1
    assert "Variant_Role" not in out.columns
    assert out["Dataset"].iloc[0] == "Dataset_A"


def test_clingen_aa_prefers_non_retracted_then_recent_approval_date():
    df = _clingen_aa_frame(
        [
            ("G1", 100, "A", "V", "1", 1000, "G", "A", "Dataset_P", 5, 1, "2020-01-01", "2020-01-01"),
            ("G1", 100, "A", "V", "1", 1003, "C", "T", "Dataset_Q", 5, 0, "2020-01-01", "2020-01-01"),
        ]
    )
    sort_by, ascending = clingen_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")
    assert out.loc[out.Variant_Role == "primary", "Dataset"].iloc[0] == "Dataset_Q"


def test_clingen_aa_recent_approval_date_breaks_tie_among_non_retracted():
    df = _clingen_aa_frame(
        [
            ("G1", 200, "A", "V", "1", 2000, "G", "A", "Dataset_P", 5, 0, "2019-01-01", "2019-01-01"),
            ("G1", 200, "A", "V", "1", 2003, "C", "T", "Dataset_Q", 5, 0, "2023-01-01", "2019-01-01"),
        ]
    )
    sort_by, ascending = clingen_aa_sort_key(df, "abs_max")
    out = aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, "abs_max")
    assert out.loc[out.Variant_Role == "primary", "Dataset"].iloc[0] == "Dataset_Q"


@pytest.mark.parametrize("strategy", ["abs_max", "nt_then_abs_max"])
def test_final_fallback_is_deterministic_dataset_name_not_row_order(strategy):
    rows_forward = [
        ("G1", 500, "A", "V", "1", 5000, "G", "A", "Dataset_Z", 5, "criteria provided, single submitter", "Jan 01, 2020", "aa"),
        ("G1", 500, "A", "V", "1", 5003, "C", "T", "Dataset_A", 5, "criteria provided, single submitter", "Jan 01, 2020", "aa"),
    ]
    cols = list(_controls_aa_frame([]).columns) + ["nucleotide_or_aa"]
    df_forward = pd.DataFrame(rows_forward, columns=cols)
    df_reversed = pd.DataFrame(list(reversed(rows_forward)), columns=cols)

    def run(df):
        sort_by, ascending = controls_aa_sort_key(df, strategy)
        if strategy == "nt_then_abs_max":
            df["_is_aa"] = df["nucleotide_or_aa"] == "aa"
            sort_by, ascending = ["_is_aa"] + sort_by, [True] + ascending
        return aa_dedup_or_mark(df, AA_GROUP_COLS, GENOMIC_KEY_COLS, sort_by, ascending, strategy)

    out_forward = run(df_forward)
    out_reversed = run(df_reversed)
    assert out_forward.loc[out_forward.Variant_Role == "primary", "Dataset"].iloc[0] == "Dataset_A"
    assert (out_forward.loc[out_forward.Variant_Role == "primary", "Dataset"].iloc[0]
            == out_reversed.loc[out_reversed.Variant_Role == "primary", "Dataset"].iloc[0])


# --- catch_mis_2 ------------------------------------------------------------------------


def test_catch_mis_2_v1_signed_value_bug_preserved_byte_for_byte():
    """v1's documented failure mode: -1 beats -5 because -1 > -5 (signed sort)."""
    df = pd.DataFrame(
        {
            "Gene": ["G1", "G1"],
            "Chrom": ["16", "16"],
            "hg38_start": [23623123, 23623123],
            "ref_allele": ["A", "A"],
            "alt_allele": ["T", "T"],
            "Fxn_points": [-1, -5],
            "nucleotide_or_aa": ["aa", "nt"],
        }
    )
    out = catch_mis_2(df, GENOMIC_KEY_COLS, strategy="v1")
    assert out["Fxn_points"].iloc[0] == -1


def test_catch_mis_2_abs_max_picks_greater_magnitude_regardless_of_sign():
    df = pd.DataFrame(
        {
            "Gene": ["G1", "G1"],
            "Chrom": ["16", "16"],
            "hg38_start": [23623123, 23623123],
            "ref_allele": ["A", "A"],
            "alt_allele": ["T", "T"],
            "Fxn_points": [-1, -5],
            "nucleotide_or_aa": ["aa", "nt"],
        }
    )
    out = catch_mis_2(df, GENOMIC_KEY_COLS, strategy="abs_max")
    assert out["Fxn_points"].iloc[0] == -5


def test_catch_mis_2_nt_then_abs_max_prefers_nt_regardless_of_magnitude():
    df = pd.DataFrame(
        {
            "Gene": ["G1", "G1"],
            "Chrom": ["17", "17"],
            "hg38_start": [43070959, 43070959],
            "ref_allele": ["A", "A"],
            "alt_allele": ["G", "G"],
            "Fxn_points": [0, -5],
            "nucleotide_or_aa": ["aa", "nt"],
        }
    )
    out = catch_mis_2(df, GENOMIC_KEY_COLS, strategy="nt_then_abs_max")
    assert out["nucleotide_or_aa"].iloc[0] == "nt"
    assert out["Fxn_points"].iloc[0] == -5


# --- dedup_vus_gnomad_unobserved --------------------------------------------------------


def test_vus_dedup_order_independent_with_custom_variant_notes_col():
    rows = [
        ("G1", 100, "A", "T", "Dataset_Z", 5, "aa", "max_fxn_pts"),
        ("G1", 100, "A", "T", "Dataset_A", 5, "aa", "max_fxn_pts"),
    ]
    df_forward = _dedup_frame(rows).rename(columns={"VariantNotes": "VariantNotes_OP"})
    df_reversed = _dedup_frame(list(reversed(rows))).rename(columns={"VariantNotes": "VariantNotes_OP"})

    out_forward = dedup_vus_gnomad_unobserved(
        df_forward, ["Gene", "hg38_start", "ref_allele", "alt_allele"], strategy="abs_max", variant_notes_col="VariantNotes_OP"
    )
    out_reversed = dedup_vus_gnomad_unobserved(
        df_reversed, ["Gene", "hg38_start", "ref_allele", "alt_allele"], strategy="abs_max", variant_notes_col="VariantNotes_OP"
    )
    assert out_forward["Dataset"].iloc[0] == "Dataset_A"
    assert out_forward["Dataset"].iloc[0] == out_reversed["Dataset"].iloc[0]


def test_vus_dedup_genuine_magnitude_still_wins_outright():
    df = _dedup_frame(
        [
            ("G1", 200, "A", "T", "Dataset_Z", 8, "aa", "max_fxn_pts"),
            ("G1", 200, "A", "T", "Dataset_A", 3, "aa", "max_fxn_pts"),
        ]
    )
    out = dedup_vus_gnomad_unobserved(df, ["Gene", "hg38_start", "ref_allele", "alt_allele"], strategy="abs_max")
    assert out["Dataset"].iloc[0] == "Dataset_Z"


def test_vus_dedup_v1_unaffected_by_dataset_name():
    df = _dedup_frame(
        [
            ("G1", 400, "A", "T", "Dataset_Z", 5, "aa", "max_fxn_pts"),
            ("G1", 400, "A", "T", "Dataset_A", 3, "nt", "First_max_fxn_pts"),
        ]
    )
    out = dedup_vus_gnomad_unobserved(df, ["Gene", "hg38_start", "ref_allele", "alt_allele"], strategy="v1")
    assert out["Dataset"].iloc[0] == "Dataset_A"  # "First_max_fxn_pts" < "max_fxn_pts" alphabetically


def test_clinvar_review_status_rank_covers_expected_tiers():
    assert CLINVAR_REVIEW_STATUS_RANK["reviewed by expert panel"] > CLINVAR_REVIEW_STATUS_RANK["criteria provided, single submitter"]
    assert CLINVAR_REVIEW_STATUS_RANK["criteria provided, single submitter"] > CLINVAR_REVIEW_STATUS_RANK["no assertion criteria provided"]

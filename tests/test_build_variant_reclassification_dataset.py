from pathlib import Path

import pandas as pd
import pytest

from src.build_variant_reclassification_dataset import (
    OUTPUT_COLUMNS,
    add_points_columns,
    apply_notebook_exclusions,
    build_reclassification_dataset,
)

CHECKPOINT_COLS = OUTPUT_COLUMNS[:-4] + [
    "VariantNotes", "ExC_points_2025", "ExC_points_2018",
    "OP_points", "Fxn_points", "Points_REVEL_GeneSpecific_GenomeWide",
]


def _checkpoint_row(**overrides):
    row = {col: None for col in CHECKPOINT_COLS}
    row.update({
        "mavedb_variant_urn": "urn:mavedb:1",
        "Dataset": "G1_Dataset",
        "Gene": "G1",
        "Chrom": "1",
        "hg38_start": 1000,
        "ref_allele": "A",
        "alt_allele": "T",
        "hgvs_p": "p.NoMatch1Xxx",
        "auth_reported_score": 99.0,
        "Flag": None,
        "VariantNotes": None,
        "splice_var_amino": "No",
        "ExC_points_2025": 5,
        "ExC_points_2018": None,
        "OP_points": None,
        "Fxn_points": 5,
        "Points_REVEL_GeneSpecific_GenomeWide": 2,
    })
    row.update(overrides)
    return row


def _checkpoint_frame(rows):
    return pd.DataFrame([_checkpoint_row(**r) for r in rows])


@pytest.fixture
def chek2_file(tmp_path):
    path = tmp_path / "CHEK2_Gebbia_2024.xlsx"
    pd.DataFrame(
        {
            "hgvs_pro": ["p.Val1Ala", "p.Val2Ala"],
            "score": [1.0, 1.0],
            "Filter_CI": [1, 0],
        }
    ).to_excel(path, index=False)
    return path


# --- apply_notebook_exclusions -----------------------------------------------------------


def test_sfpq_dropped(chek2_file):
    df = _checkpoint_frame([{"Gene": "SFPQ"}, {"Gene": "G1"}])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["Gene"]) == ["G1"]


def test_chek2_flagged_row_dropped(chek2_file):
    """hgvs_p (transcript-prefixed, as in the real checkpoint)/
    auth_reported_score matching a Filter_CI==1 CHEK2 row (unprefixed
    hgvs_pro) gets Flag='*' and is then removed by the Flag!='*' filter."""
    df = _checkpoint_frame([
        {"Gene": "CHEK2", "hgvs_p": "NP_009125.1:p.Val1Ala", "auth_reported_score": 1.0},
        {"Gene": "CHEK2", "hgvs_p": "NP_009125.1:p.Val2Ala", "auth_reported_score": 1.0},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert len(out) == 1
    assert out["hgvs_p"].iloc[0] == "NP_009125.1:p.Val2Ala"


def test_conflicting_fxn_data_dropped(chek2_file):
    df = _checkpoint_frame([
        {"VariantNotes": "conflicting_fxn_data"},
        {"VariantNotes": "max_fxn_pts"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["VariantNotes"]) == ["max_fxn_pts"]


def test_splice_var_amino_dropped(chek2_file):
    df = _checkpoint_frame([
        {"splice_var_amino": "Yes"},
        {"splice_var_amino": "No"},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert list(out["splice_var_amino"]) == ["No"]


def test_pre_existing_flag_still_removed(chek2_file):
    df = _checkpoint_frame([
        {"Flag": "*"},
        {"Flag": None},
    ])
    out = apply_notebook_exclusions(df, chek2_file)
    assert len(out) == 1
    assert pd.isna(out["Flag"].iloc[0])


# --- add_points_columns -------------------------------------------------------------------


def test_excalibr_points_is_literal_value_even_when_op_overrides_functional():
    """TP53 uses OddsPath for Functional_points, but ExCALIBR_points should
    still hold the literal ExC_points_2025 value."""
    df = _checkpoint_frame([
        {"Gene": "TP53", "ExC_points_2025": 4, "OP_points": 8, "Fxn_points": 8},
    ])
    out = add_points_columns(df)
    assert out["ExCALIBR_points"].iloc[0] == 4
    assert out["OddsPath_points"].iloc[0] == 8
    assert out["Functional_points"].iloc[0] == 8


def test_excalibr_points_vintage_override_for_brca1_pten_msh2():
    df = _checkpoint_frame([
        {"Gene": "BRCA1", "ExC_points_2025": 2, "ExC_points_2018": 6},
        {"Gene": "G1", "ExC_points_2025": 2, "ExC_points_2018": 6},
    ])
    out = add_points_columns(df)
    assert out.loc[out.Gene == "BRCA1", "ExCALIBR_points"].iloc[0] == 6
    assert out.loc[out.Gene == "G1", "ExCALIBR_points"].iloc[0] == 2


def test_excalibr_points_vintage_override_falls_back_when_2018_missing():
    df = _checkpoint_frame([{"Gene": "PTEN", "ExC_points_2025": 3, "ExC_points_2018": None}])
    out = add_points_columns(df)
    assert out["ExCALIBR_points"].iloc[0] == 3


def test_combined_points_is_functional_plus_revel_gene_specific():
    df = _checkpoint_frame([{"Fxn_points": 5, "Points_REVEL_GeneSpecific_GenomeWide": 2}])
    out = add_points_columns(df)
    assert out["Combined_points"].iloc[0] == 7


def test_combined_points_treats_missing_revel_as_zero():
    df = _checkpoint_frame([{"Fxn_points": 5, "Points_REVEL_GeneSpecific_GenomeWide": None}])
    out = add_points_columns(df)
    assert out["Combined_points"].iloc[0] == 5


# --- build_reclassification_dataset (end to end) -------------------------------------------


def test_build_reclassification_dataset_dedups_by_dna_variant(tmp_path, chek2_file):
    checkpoint = _checkpoint_frame([
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_low", "Fxn_points": 2,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
        {"Gene": "G1", "hg38_start": 5000, "Dataset": "Dataset_high", "Fxn_points": 9,
         "Points_REVEL_GeneSpecific_GenomeWide": 0},
    ])
    checkpoint_path = tmp_path / "checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")

    out = build_reclassification_dataset(checkpoint_path, chek2_file)
    assert len(out) == 1
    assert out["Dataset"].iloc[0] == "Dataset_high"
    assert list(out.columns) == OUTPUT_COLUMNS

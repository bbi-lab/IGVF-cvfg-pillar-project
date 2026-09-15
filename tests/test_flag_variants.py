"""Tests for src.flag_variants's LDLR LA-module amino-acid-range flagging.

Only the LDLR aa-range logic is covered here (`flag_variants.py`'s other
dataset kinds -- coordinate-join, string-join, SNV filter -- have no
existing test coverage to extend, and their assay-file-dependent inputs are
out of scope for this change).
"""

import pandas as pd

from src.flag_variants import compute_flags


def _ldlr_row(dataset_name, aa_pos, mapped_hgvs_g="chr19:g.1A>G"):
    return {
        "mapped_hgvs_g": mapped_hgvs_g,
        "dataset_name": dataset_name,
        "gene_symbol": "LDLR",
        "mapped_hgvs_p_start": aa_pos,
    }


def test_flags_positions_at_each_la_module_boundary():
    # Both endpoints of every LA module range should be flagged, and the
    # single-position gap between LA4 (146-186) and LA5 (195-233) should not.
    in_range_positions = [25, 65, 66, 106, 106, 145, 146, 186, 195, 233, 234, 272]
    out_of_range_positions = [1, 24, 190, 273]

    positions = in_range_positions + out_of_range_positions
    df = pd.DataFrame([_ldlr_row("LDLR_Tabet_2025_uptake", pos) for pos in positions])

    flag = compute_flags(df, "data/filtering")

    assert flag.iloc[: len(in_range_positions)].tolist() == ["*"] * len(in_range_positions)
    assert flag.iloc[len(in_range_positions) :].tolist() == [""] * len(out_of_range_positions)


def test_flags_both_remaining_ldlr_datasets():
    df = pd.DataFrame(
        [
            _ldlr_row("LDLR_Tabet_2025_uptake", 30),
            _ldlr_row("LDLR_Tabet_2025_abundance", 30),
        ]
    )
    flag = compute_flags(df, "data/filtering")
    assert flag.tolist() == ["*", "*"]


def test_does_not_flag_other_genes_or_datasets_at_the_same_position():
    df = pd.DataFrame(
        [
            {
                "mapped_hgvs_g": "chr1:g.1A>G",
                "dataset_name": "GENEA_Study_2020",
                "gene_symbol": "GENEA",
                "mapped_hgvs_p_start": "30",
            },
        ]
    )
    flag = compute_flags(df, "data/filtering")
    assert flag.tolist() == [""]


def test_flags_every_dna_candidate_in_a_multi_candidate_row():
    df = pd.DataFrame([_ldlr_row("LDLR_Tabet_2025_uptake", 30, mapped_hgvs_g="chr19:g.1A>G|chr19:g.2A>C")])
    flag = compute_flags(df, "data/filtering")
    assert flag.tolist() == ["*|*"]

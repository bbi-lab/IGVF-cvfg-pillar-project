import pandas as pd
import pytest

from src.make_extended_data_figure_7alt_op import GROUPS, PREDICTORS, build_figure, load_counts, sheet_name


def _write_workbook(path, category_counts_by_sheet, consequences_by_sheet=None):
    """category_counts_by_sheet: {sheet_name: {category: count}}. Writes one
    row per variant (Class_OP_<predictor> repeated `count` times per category),
    plus a simplified_consequence column: "missense_variant" for every row,
    unless consequences_by_sheet[sheet] gives an explicit per-row list."""
    with pd.ExcelWriter(path) as writer:
        for (group, predictor), category_counts in category_counts_by_sheet.items():
            sheet = sheet_name(group, predictor)
            values = [category for category, count in category_counts.items() for _ in range(count)]
            consequences = (consequences_by_sheet or {}).get((group, predictor)) or ["missense_variant"] * len(values)
            pd.DataFrame({f"Class_OP_{predictor}": values, "simplified_consequence": consequences}).to_excel(
                writer, sheet_name=sheet, index=False
            )


def _default_counts():
    return {
        (group, predictor): {
            "Likely Pathogenic": 2,
            "Uncertain": 3,
            "Likely Benign": 4,
            "Benign": 5,
        }
        for group in GROUPS
        for predictor in PREDICTORS
    }


def test_sheet_name_uses_mut_override_for_unobserved_mp2():
    assert sheet_name("Unobserved", "MP2") == "Unobserved_mut_OP"
    assert sheet_name("Unobserved", "REVEL") == "Unobserved_REVEL_OP"
    assert sheet_name("VUS", "MP2") == "VUS_MP2_OP"
    assert sheet_name("gnomAD", "AM") == "gnomAD_AM_OP"


def test_load_counts_computes_totals_and_percentages(tmp_path):
    workbook = tmp_path / "data.xlsx"
    _write_workbook(workbook, _default_counts())

    counts_df = load_counts(workbook)

    assert len(counts_df) == len(GROUPS) * len(PREDICTORS)
    row = counts_df[(counts_df["group"] == "VUS") & (counts_df["predictor"] == "REVEL")].iloc[0]
    assert row["total"] == 14
    assert row["Benign"] == 5
    assert row["pct_Benign"] == pytest.approx(5 / 14 * 100)


def test_load_counts_fills_missing_category_with_zero(tmp_path):
    # Every OddsPath-based VUS/gnomAD/Unobserved sheet in the real
    # Supplementary_Data_6.xlsx lacks "Pathogenic" entirely (see module
    # docstring) -- unlike make_extended_data_figure_7alt.py's stricter
    # ExCALIBR-based counterpart, that must not raise here.
    workbook = tmp_path / "data.xlsx"
    _write_workbook(workbook, _default_counts())

    counts_df = load_counts(workbook)

    row = counts_df[(counts_df["group"] == "VUS") & (counts_df["predictor"] == "REVEL")].iloc[0]
    assert row["Pathogenic"] == 0
    assert row["pct_Pathogenic"] == 0


def test_load_counts_applies_consequence_filter(tmp_path):
    workbook = tmp_path / "data.xlsx"
    sheets = _default_counts()
    # VUS/REVEL: 2 Likely Pathogenic, 3 Uncertain, 4 Likely Benign, 5 Benign
    # (14 rows total) -- one Likely Pathogenic row and one Benign row are
    # non-missense; missense-only filtering should drop just those.
    consequences_by_sheet = {
        ("VUS", "REVEL"): (
            ["missense_variant", "synonymous_variant"]
            + ["missense_variant"] * 3
            + ["missense_variant"] * 4
            + ["missense_variant"] * 4
            + ["synonymous_variant"]
        )
    }
    _write_workbook(workbook, sheets, consequences_by_sheet=consequences_by_sheet)

    all_counts = load_counts(workbook)
    missense_counts = load_counts(workbook, allowed_consequences=["missense_variant"])

    all_row = all_counts[(all_counts["group"] == "VUS") & (all_counts["predictor"] == "REVEL")].iloc[0]
    missense_row = missense_counts[
        (missense_counts["group"] == "VUS") & (missense_counts["predictor"] == "REVEL")
    ].iloc[0]

    assert all_row["total"] == 14
    assert all_row["Benign"] == 5
    assert missense_row["total"] == 12
    assert missense_row["Benign"] == 4

    # Other sheets are untouched (every row is "missense_variant" by default).
    other_row = missense_counts[(missense_counts["group"] == "VUS") & (missense_counts["predictor"] == "AM")].iloc[0]
    assert other_row["total"] == 14


def test_load_counts_rejects_unexpected_category(tmp_path):
    workbook = tmp_path / "data.xlsx"
    sheets = _default_counts()
    sheets[("VUS", "REVEL")]["Conflicting evidence"] = 1
    _write_workbook(workbook, sheets)

    with pytest.raises(ValueError, match="unexpected categories"):
        load_counts(workbook)


def test_build_figure_writes_a_pdf(tmp_path):
    workbook = tmp_path / "data.xlsx"
    _write_workbook(workbook, _default_counts())
    counts_df = load_counts(workbook)

    output_path = tmp_path / "nested" / "extended_data_figure_7alt_op.pdf"
    build_figure(counts_df, output_path)

    assert output_path.exists()
    assert output_path.read_bytes().startswith(b"%PDF")

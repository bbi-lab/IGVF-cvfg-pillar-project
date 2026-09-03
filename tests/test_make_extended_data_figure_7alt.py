import pandas as pd
import pytest

from src.make_extended_data_figure_7alt import GROUPS, PREDICTORS, build_figure, load_counts


def _write_workbook(path, category_counts_by_sheet, consequences_by_sheet=None):
    """category_counts_by_sheet: {sheet_name: {category: count}}. Writes one
    row per variant (Class_<predictor> repeated `count` times per category),
    plus a simplified_consequence column: "missense_variant" for every row,
    unless consequences_by_sheet[sheet] gives an explicit per-row list."""
    with pd.ExcelWriter(path) as writer:
        for sheet, category_counts in category_counts_by_sheet.items():
            predictor = sheet.rsplit("_", 1)[1]
            values = [category for category, count in category_counts.items() for _ in range(count)]
            consequences = (consequences_by_sheet or {}).get(sheet) or ["missense_variant"] * len(values)
            pd.DataFrame({f"Class_{predictor}": values, "simplified_consequence": consequences}).to_excel(
                writer, sheet_name=sheet, index=False
            )


def _default_counts():
    return {
        f"{group}_{predictor}": {
            "Pathogenic": 1,
            "Likely Pathogenic": 2,
            "Uncertain": 3,
            "Likely Benign": 4,
            "Benign": 5,
        }
        for group in GROUPS
        for predictor in PREDICTORS
    }


def test_load_counts_computes_totals_and_percentages(tmp_path):
    workbook = tmp_path / "data.xlsx"
    _write_workbook(workbook, _default_counts())

    counts_df = load_counts(workbook)

    assert len(counts_df) == len(GROUPS) * len(PREDICTORS)
    row = counts_df[(counts_df["group"] == "VUS") & (counts_df["predictor"] == "REVEL")].iloc[0]
    assert row["total"] == 15
    assert row["Pathogenic"] == 1
    assert row["Benign"] == 5
    assert row["pct_Pathogenic"] == pytest.approx(1 / 15 * 100)
    assert row["pct_Benign"] == pytest.approx(5 / 15 * 100)


def test_load_counts_applies_consequence_filter(tmp_path):
    workbook = tmp_path / "data.xlsx"
    sheets = _default_counts()
    # VUS_REVEL: 2 Pathogenic, 2 Likely Pathogenic, 3 Uncertain, 4 Likely
    # Benign, 5 Benign (16 rows total) -- one Pathogenic row and one Benign
    # row are non-missense; missense-only filtering should drop just those.
    sheets["VUS_REVEL"] = {
        "Pathogenic": 2,
        "Likely Pathogenic": 2,
        "Uncertain": 3,
        "Likely Benign": 4,
        "Benign": 5,
    }
    consequences_by_sheet = {
        "VUS_REVEL": (
            ["missense_variant", "synonymous_variant"]
            + ["missense_variant"] * 2
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

    assert all_row["total"] == 16
    assert all_row["Pathogenic"] == 2
    assert all_row["Benign"] == 5
    assert missense_row["total"] == 14
    assert missense_row["Pathogenic"] == 1
    assert missense_row["Benign"] == 4

    # Other sheets are untouched (every row is "missense_variant" by default).
    other_row = missense_counts[(missense_counts["group"] == "VUS") & (missense_counts["predictor"] == "AM")].iloc[0]
    assert other_row["total"] == 15


def test_load_counts_rejects_missing_category(tmp_path):
    workbook = tmp_path / "data.xlsx"
    sheets = _default_counts()
    del sheets["VUS_REVEL"]["Benign"]
    _write_workbook(workbook, sheets)

    with pytest.raises(ValueError, match="missing expected categories"):
        load_counts(workbook)


def test_load_counts_rejects_unexpected_category(tmp_path):
    workbook = tmp_path / "data.xlsx"
    sheets = _default_counts()
    sheets["VUS_REVEL"]["Conflicting evidence"] = 1
    _write_workbook(workbook, sheets)

    with pytest.raises(ValueError, match="unexpected categories"):
        load_counts(workbook)


def test_build_figure_writes_a_pdf(tmp_path):
    workbook = tmp_path / "data.xlsx"
    _write_workbook(workbook, _default_counts())
    counts_df = load_counts(workbook)

    output_path = tmp_path / "nested" / "extended_data_figure_7alt.pdf"
    build_figure(counts_df, output_path)

    assert output_path.exists()
    assert output_path.read_bytes().startswith(b"%PDF")

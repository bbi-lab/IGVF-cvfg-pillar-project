import pandas as pd
import pytest
from click.testing import CliRunner

from src.recalculate_clingen_classification import (
    CLASSIFICATION_COL,
    EVIDENCE_COL,
    classify_acmg,
    codes_in_segment,
    compute_updated_classifications,
    filter_and_recalculate,
    main,
    recalculate_row,
)


def test_classify_acmg_pathogenic_rules():
    assert classify_acmg(["PVS1", "PS1"]) == "Pathogenic"
    assert classify_acmg(["PS1", "PS2"]) == "Pathogenic"
    assert classify_acmg(["PS1", "PM1", "PM2"]) == "Pathogenic"
    assert classify_acmg(["PVS1", "PM1"]) == "Likely Pathogenic"
    assert classify_acmg(["PS1", "PM1"]) == "Likely Pathogenic"
    assert classify_acmg(["PM1", "PM2", "PM3"]) == "Likely Pathogenic"
    assert classify_acmg(["PM1", "PM2", "PP1", "PP2"]) == "Likely Pathogenic"


def test_classify_acmg_benign_rules():
    assert classify_acmg(["BS1", "BS2"]) == "Benign"
    assert classify_acmg(["BS1", "BP1"]) == "Likely Benign"
    assert classify_acmg(["BS1", "BS2", "BS3", "BS4"]) == "Benign"
    assert classify_acmg(["BM1", "BM2"]) == "VUS"  # not real codes, just unrecognized


def test_classify_acmg_ba1_never_counts_as_standalone():
    """BA1 maps to "Benign_Standalone" but the lookup checks "Benign_standalone"
    (lowercase s) -- ported bug, preserved for parity with the notebook."""
    assert classify_acmg(["BA1"]) == "VUS"


def test_classify_acmg_empty_or_unrecognized_is_vus():
    assert classify_acmg([]) == "VUS"
    assert classify_acmg(["NOT_A_CODE"]) == "VUS"


def test_codes_in_segment_splits_comma_and_merged_records():
    assert codes_in_segment("PM2,PP3") == ["PM2", "PP3"]
    assert codes_in_segment("PM2,PP3 | BS1") == ["PM2", "PP3", "BS1"]
    assert codes_in_segment("") == []


def test_filter_and_recalculate_strips_functional_evidence():
    classification, evidence = filter_and_recalculate("PM2,PS3,PP3,BP4")
    assert evidence == "PM2"
    assert classification == "VUS"


def test_filter_and_recalculate_recomputes_from_remaining_evidence():
    classification, evidence = filter_and_recalculate("PS1,PS2,PS3")
    assert evidence == "PS1,PS2"
    assert classification == "Pathogenic"


def test_filter_and_recalculate_empty_segment_is_vus():
    assert filter_and_recalculate("") == ("VUS", "")


def test_filter_and_recalculate_all_evidence_removed_is_vus():
    assert filter_and_recalculate("PS3,BP4") == ("VUS", "")


def test_recalculate_row_aligns_pipe_delimited_candidates():
    classifications, evidence = recalculate_row("PS1,PS2|PS3,BP4|")
    assert classifications == "Pathogenic|VUS|VUS"
    assert evidence == "PS1,PS2||"


def test_recalculate_row_empty_field_is_single_vus_candidate():
    assert recalculate_row("") == ("VUS", "")


def test_compute_updated_classifications_missing_column_raises_value_error():
    df = pd.DataFrame({"other_col": ["PS1"]})
    with pytest.raises(ValueError, match="not found"):
        compute_updated_classifications(df, "clingen_evidence_repository.Applied Evidence Codes (Met)")


def test_compute_updated_classifications_adds_expected_series(tmp_path):
    df = pd.DataFrame(
        {
            "clingen_evidence_repository.Applied Evidence Codes (Met)": [
                "PS1,PS2,PS3",
                "PS3,BP4",
                "",
            ]
        }
    )
    classification, evidence = compute_updated_classifications(
        df, "clingen_evidence_repository.Applied Evidence Codes (Met)"
    )
    assert list(classification) == ["Pathogenic", "VUS", "VUS"]
    assert list(evidence) == ["PS1,PS2", "", ""]


def test_main_cli_writes_expected_columns(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame(
        {
            "id": ["v1", "v2"],
            "clingen_evidence_repository.Applied Evidence Codes (Met)": [
                "PS1,PS2,PS3",
                "PS3,BP4",
            ],
            "clingen_evidence_repository.Assertion": [
                "Pathogenic",
                "Uncertain significance",
            ],
        }
    ).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out[CLASSIFICATION_COL]) == ["Pathogenic", "VUS"]
    assert list(out[EVIDENCE_COL]) == ["PS1,PS2", ""]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"id": ["v1"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "not found" in result.output

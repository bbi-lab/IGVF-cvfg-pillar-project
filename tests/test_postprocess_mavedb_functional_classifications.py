import pandas as pd
from click.testing import CliRunner

from src.postprocess_mavedb_functional_classifications import apply_functional_class_overrides, main


def _df(rows):
    return pd.DataFrame(
        {
            "variant_urn": [r.get("urn", "") for r in rows],
            "mavedb.primary_calibration.functional_class_label": [r.get("p_label", "") for r in rows],
            "mavedb.primary_calibration.functional_classification": [r.get("p_class", "") for r in rows],
            "mavedb.investigator_provided_calibration.functional_class_label": [r.get("i_label", "") for r in rows],
            "mavedb.investigator_provided_calibration.functional_classification": [r.get("i_class", "") for r in rows],
        }
    )


def test_overrides_matching_urn_prefix_and_label():
    df = _df(
        [
            {
                "urn": "urn:mavedb:00000674-b-1#1",
                "p_label": "gain-of-function",
                "p_class": "abnormal",
                "i_label": "gain-of-function",
                "i_class": "abnormal",
            }
        ]
    )
    counts = apply_functional_class_overrides(df)

    assert df["mavedb.primary_calibration.functional_classification"].iloc[0] == "not_specified"
    assert df["mavedb.investigator_provided_calibration.functional_classification"].iloc[0] == "not_specified"
    assert sum(counts.values()) == 2


def test_leaves_other_urns_untouched():
    df = _df(
        [
            {
                "urn": "urn:mavedb:00001280-a-1#1",
                "p_label": "gain-of-function",
                "p_class": "abnormal",
            }
        ]
    )
    apply_functional_class_overrides(df)

    assert df["mavedb.primary_calibration.functional_classification"].iloc[0] == "abnormal"


def test_leaves_non_matching_label_untouched():
    df = _df(
        [
            {
                "urn": "urn:mavedb:00000674-b-1#1",
                "p_label": "loss-of-function",
                "p_class": "abnormal",
            }
        ]
    )
    apply_functional_class_overrides(df)

    assert df["mavedb.primary_calibration.functional_classification"].iloc[0] == "abnormal"


def test_missing_calibration_columns_are_skipped():
    df = pd.DataFrame({"variant_urn": ["urn:mavedb:00000674-b-1#1"]})
    counts = apply_functional_class_overrides(df)

    assert sum(counts.values()) == 0


def test_missing_variant_urn_col_raises_value_error():
    df = pd.DataFrame({"other_col": ["x"]})
    try:
        apply_functional_class_overrides(df)
        raise AssertionError("expected ValueError")
    except ValueError as exc:
        assert "variant_urn" in str(exc)


def test_main_cli_writes_overridden_classification(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    _df(
        [
            {
                "urn": "urn:mavedb:00000674-b-1#1",
                "p_label": "gain-of-function",
                "p_class": "abnormal",
            },
            {
                "urn": "urn:mavedb:00001280-a-1#1",
                "p_label": "gain-of-function",
                "p_class": "abnormal",
            },
        ]
    ).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["mavedb.primary_calibration.functional_classification"]) == ["not_specified", "abnormal"]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"id": ["v1"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "not found in input" in result.output

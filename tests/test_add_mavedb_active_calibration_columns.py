import pandas as pd
import pytest
from click.testing import CliRunner

from src.add_mavedb_active_calibration_columns import main, select_active_calibration


def _df(rows, requested=True):
    base = {
        "mavedb.primary_calibration.urn": [r.get("p_urn", "") for r in rows],
        "mavedb.primary_calibration.name": [r.get("p_name", "") for r in rows],
        "mavedb.primary_calibration.url": [r.get("p_url", "") for r in rows],
        "mavedb.primary_calibration.functional_class_label": [r.get("p_label", "") for r in rows],
        "mavedb.primary_calibration.functional_classification": [r.get("p_class", "") for r in rows],
        "mavedb.investigator_provided_calibration.urn": [r.get("i_urn", "") for r in rows],
        "mavedb.investigator_provided_calibration.name": [r.get("i_name", "") for r in rows],
        "mavedb.investigator_provided_calibration.url": [r.get("i_url", "") for r in rows],
        "mavedb.investigator_provided_calibration.functional_class_label": [r.get("i_label", "") for r in rows],
        "mavedb.investigator_provided_calibration.functional_classification": [r.get("i_class", "") for r in rows],
    }
    if requested:
        base.update(
            {
                "mavedb.requested_calibration.urn": [r.get("r_urn", "") for r in rows],
                "mavedb.requested_calibration.name": [r.get("r_name", "") for r in rows],
                "mavedb.requested_calibration.url": [r.get("r_url", "") for r in rows],
                "mavedb.requested_calibration.functional_class_label": [r.get("r_label", "") for r in rows],
                "mavedb.requested_calibration.functional_classification": [r.get("r_class", "") for r in rows],
            }
        )
    return pd.DataFrame(base)


def test_uses_requested_calibration_when_present():
    df = _df([{"r_urn": "urn:req", "r_name": "Requested", "p_name": "Fayer calibration", "i_name": "Investigator"}])
    active, source = select_active_calibration(df)
    assert list(source) == ["requested"]
    assert active["mavedb.active_calibration.name"].iloc[0] == "Requested"


def test_uses_primary_when_no_requested_and_name_preferred():
    df = _df([{"p_name": "Fayer calibration", "p_urn": "urn:p", "i_name": "Investigator"}])
    active, source = select_active_calibration(df)
    assert list(source) == ["primary"]
    assert active["mavedb.active_calibration.name"].iloc[0] == "Fayer calibration"


def test_uses_investigator_when_primary_name_not_preferred():
    df = _df([{"p_name": "Some other calibration", "i_name": "Investigator", "i_urn": "urn:i"}])
    active, source = select_active_calibration(df)
    assert list(source) == ["investigator_provided"]
    assert active["mavedb.active_calibration.name"].iloc[0] == "Investigator"


def test_uses_investigator_when_primary_name_blank():
    df = _df([{"i_name": "Investigator", "i_urn": "urn:i"}])
    _active, source = select_active_calibration(df)
    assert list(source) == ["investigator_provided"]


def test_requested_columns_optional_falls_back_to_primary_investigator_logic():
    df = _df([{"p_name": "Scott calibration"}], requested=False)
    _active, source = select_active_calibration(df)
    assert list(source) == ["primary"]


def test_missing_required_column_raises_value_error():
    df = _df([{}]).drop(columns=["mavedb.primary_calibration.urn"])
    with pytest.raises(ValueError, match="missing required column"):
        select_active_calibration(df)


def test_partial_requested_columns_raises_value_error():
    df = _df([{}]).drop(columns=["mavedb.requested_calibration.name"])
    with pytest.raises(ValueError, match="some but not all"):
        select_active_calibration(df)


def test_custom_preferred_primary_calibration_names():
    df = _df([{"p_name": "Custom calibration", "i_name": "Investigator"}])
    _active, source = select_active_calibration(df, preferred_primary_calibration_names=["Custom calibration"])
    assert list(source) == ["primary"]


def test_main_cli_writes_active_calibration_columns(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    _df(
        [
            {"r_urn": "urn:req", "r_name": "Requested"},
            {"p_name": "Fayer calibration", "p_urn": "urn:p"},
            {"i_name": "Investigator", "i_urn": "urn:i"},
        ]
    ).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["mavedb.active_calibration.name"]) == ["Requested", "Fayer calibration", "Investigator"]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"id": ["v1"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "missing required column" in result.output

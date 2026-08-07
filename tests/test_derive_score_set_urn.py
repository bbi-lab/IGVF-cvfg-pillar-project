import pandas as pd
import pytest
from click.testing import CliRunner

from src.derive_score_set_urn import derive_score_set_urn, main


def test_strips_variant_index_suffix():
    df = pd.DataFrame({"variant_urn": ["urn:mavedb:00000097-0-2#1", "urn:mavedb:00000097-0-2#23"]})
    result = derive_score_set_urn(df)
    assert list(result) == ["urn:mavedb:00000097-0-2", "urn:mavedb:00000097-0-2"]
    assert result.name == "score_set_urn"


def test_leaves_urn_without_suffix_unchanged():
    df = pd.DataFrame({"variant_urn": ["urn:mavedb:00000097-0-2"]})
    assert list(derive_score_set_urn(df)) == ["urn:mavedb:00000097-0-2"]


def test_missing_variant_urn_column_raises_value_error():
    df = pd.DataFrame({"id": ["v1"]})
    with pytest.raises(ValueError, match="missing required column"):
        derive_score_set_urn(df)


def test_custom_column_names():
    df = pd.DataFrame({"vurn": ["urn:mavedb:00000097-0-2#1"]})
    result = derive_score_set_urn(df, variant_urn_column="vurn", output_column="ssurn")
    assert list(result) == ["urn:mavedb:00000097-0-2"]
    assert result.name == "ssurn"


def test_main_cli_writes_score_set_urn_column(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"variant_urn": ["urn:mavedb:00000097-0-2#1", "urn:mavedb:00000097-0-2#2"]}).to_csv(
        input_path, sep="\t", index=False
    )

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["score_set_urn"]) == ["urn:mavedb:00000097-0-2", "urn:mavedb:00000097-0-2"]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"id": ["v1"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "missing required column" in result.output


def test_main_cli_survives_multiline_quoted_field(tmp_path):
    """Regression test for the BRCA1 (urn:mavedb:00000097-0-2) rows whose
    `mavedb_mapping_error` carries an embedded newline (a multi-line
    HTTPStatusError message MaveDB persisted verbatim). The original `awk`
    implementation of this step read the file line-by-line, so it split each
    such row into two records and appended the derived score_set_urn to the
    wrong place -- inside the still-open quoted field -- instead of as a
    real column, corrupting the value a later `merge-columns` join would key
    on. Reading/writing with pandas keeps the quoted multi-line field as one
    logical row.
    """
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    multiline_error = "HTTPStatusError: Redirect response '302 Moved Temporarily'\nFor more information check: ..."
    pd.DataFrame(
        {
            "variant_urn": ["urn:mavedb:00000097-0-2#1", "urn:mavedb:00000098-1-1#1"],
            "mavedb_mapping_error": [multiline_error, ""],
        }
    ).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert len(out) == 2
    assert list(out["score_set_urn"]) == ["urn:mavedb:00000097-0-2", "urn:mavedb:00000098-1-1"]
    assert out["mavedb_mapping_error"].iloc[0] == multiline_error

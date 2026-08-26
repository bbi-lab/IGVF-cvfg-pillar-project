import pandas as pd
import pytest
from click.testing import CliRunner

from src.translate_assayed_variant_level import main, translate_assayed_variant_level


def test_translates_known_codes():
    df = pd.DataFrame({"assayed_variant_level": ["protein", "dna", "protein"]})
    result = translate_assayed_variant_level(df)
    assert list(result) == ["aa", "nt", "aa"]


def test_leaves_unknown_codes_unchanged():
    df = pd.DataFrame({"assayed_variant_level": ["protein", "", "other"]})
    assert list(translate_assayed_variant_level(df)) == ["aa", "", "other"]


def test_missing_column_raises_value_error():
    df = pd.DataFrame({"id": ["v1"]})
    with pytest.raises(ValueError, match="missing required column"):
        translate_assayed_variant_level(df)


def test_custom_column_and_code_map():
    df = pd.DataFrame({"level": ["x", "y"]})
    result = translate_assayed_variant_level(df, column="level", code_map={"x": "z"})
    assert list(result) == ["z", "y"]


def test_main_cli_translates_codes(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"assayed_variant_level": ["protein", "dna", ""]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["assayed_variant_level"]) == ["aa", "nt", ""]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"id": ["v1"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "missing required column" in result.output


def test_main_cli_survives_multiline_quoted_field(tmp_path):
    """Regression test mirroring test_derive_score_set_urn.py's: a value
    spanning multiple physical lines (RFC 4180 quoted) must not corrupt the
    row's assayed_variant_level translation, the way the original `awk`
    implementation of this step did when it read the file line-by-line.
    """
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    multiline_error = "HTTPStatusError: Redirect response '302 Moved Temporarily'\nFor more information check: ..."
    pd.DataFrame(
        {
            "assayed_variant_level": ["protein", "dna"],
            "mavedb_mapping_error": [multiline_error, ""],
        }
    ).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert len(out) == 2
    assert list(out["assayed_variant_level"]) == ["aa", "nt"]
    assert out["mavedb_mapping_error"].iloc[0] == multiline_error

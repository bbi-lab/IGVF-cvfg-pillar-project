import pandas as pd
import pytest
from click.testing import CliRunner

from src.merge_rna_score_columns import main, merge_rna_score


def test_fills_blank_target_from_source():
    df = pd.DataFrame({"rna_score": ["", "1.5", ""], "rna_score_d6": ["2.2", "9.9", ""]})
    result = merge_rna_score(df)
    assert list(result) == ["2.2", "1.5", ""]


def test_leaves_populated_target_unchanged_even_if_source_also_populated():
    df = pd.DataFrame({"rna_score": ["1.5"], "rna_score_d6": ["9.9"]})
    assert list(merge_rna_score(df)) == ["1.5"]


def test_missing_target_column_raises_value_error():
    df = pd.DataFrame({"rna_score_d6": ["1.0"]})
    with pytest.raises(ValueError, match="missing required column"):
        merge_rna_score(df)


def test_missing_source_column_raises_value_error():
    df = pd.DataFrame({"rna_score": [""]})
    with pytest.raises(ValueError, match="missing required column"):
        merge_rna_score(df)


def test_custom_column_names():
    df = pd.DataFrame({"target": [""], "source": ["3.3"]})
    result = merge_rna_score(df, target_column="target", source_column="source")
    assert list(result) == ["3.3"]


def test_main_cli_fills_blank_rna_score(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"rna_score": ["", "1.5"], "rna_score_d6": ["2.2", "9.9"]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code == 0, result.output
    assert "Filled rna_score from rna_score_d6 for 1 of 2 row(s)." in result.output
    out = pd.read_csv(output_path, sep="\t", dtype=str, keep_default_na=False)
    assert list(out["rna_score"]) == ["2.2", "1.5"]


def test_main_cli_missing_column_is_click_exception(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "output.tsv"
    pd.DataFrame({"rna_score": [""]}).to_csv(input_path, sep="\t", index=False)

    runner = CliRunner()
    result = runner.invoke(main, [str(input_path), str(output_path)])

    assert result.exit_code != 0
    assert "missing required column" in result.output

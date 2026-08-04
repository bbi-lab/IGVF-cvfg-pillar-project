import pandas as pd
import pytest
from click.testing import CliRunner

from src.mave_dataset_stats import compute_all_stats, main, split_genes


def _write_condensed(path, rows):
    columns = ["Dataset", "Gene", "hgvs_c", "hgvs_p"]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


def _write_metadata(path, rows):
    columns = ["Dataset Name", "IGVF Produced?", "Primary Score Set or Meta-analysis?"]
    df = pd.DataFrame(rows, columns=columns)
    with pd.ExcelWriter(path) as writer:
        df.to_excel(writer, sheet_name="Curation", index=False)


@pytest.fixture
def dataset_files(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_IGVF_A", "GENEA", "c1", "p1"),
            ("DS_IGVF_A", "GENEA", "c1", "p1"),  # duplicate variant, second measurement
            ("DS_IGVF_B", "GENEB, GENEC", "c2", "p2"),
            ("DS_COMM_A", "GENEC", "c3", "p3"),
            ("DS_COMM_B", "GENED", "c4", "p4"),
        ],
    )
    _write_metadata(
        metadata_path,
        [
            ("DS_IGVF_A", "Yes", "primary score set"),
            ("DS_IGVF_B", "Yes", "primary score set"),
            ("DS_COMM_A", "No", "primary score set"),
            ("DS_COMM_B", "No", "meta-analysis"),
        ],
    )
    return condensed_path, metadata_path


def test_split_genes_handles_multi_gene_datasets():
    assert split_genes("CALM1, CALM2, CALM3") == ["CALM1", "CALM2", "CALM3"]
    assert split_genes("BRCA1") == ["BRCA1"]


def test_compute_all_stats_buckets(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats = compute_all_stats(condensed_path, metadata_path)

    igvf = stats["Community (IGVF)"]
    assert igvf["datasets"] == 2
    assert igvf["variant_effect_measurements"] == 3
    assert igvf["composite_scores"] == 0
    assert igvf["distinct_variants_assayed"] == 2
    assert igvf["genes_represented"] == 3

    non_igvf = stats["Community (non-IGVF)"]
    assert non_igvf["datasets"] == 2
    assert non_igvf["variant_effect_measurements"] == 1
    assert non_igvf["composite_scores"] == 1
    assert non_igvf["distinct_variants_assayed"] == 2
    assert non_igvf["genes_represented"] == 2
    assert non_igvf["genes_not_in_igvf_data"] == 1  # GENED only; GENEC is shared with IGVF

    combined = stats["Combined (IGVF + community)"]
    assert combined["datasets"] == 4
    assert combined["variant_effect_measurements"] == 4
    assert combined["composite_scores"] == 1
    assert combined["distinct_variants_assayed"] == 4
    assert combined["genes_represented"] == 4


def test_compute_all_stats_raises_on_missing_metadata(dataset_files):
    condensed_path, metadata_path = dataset_files
    _write_condensed(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5")],
    )
    with pytest.raises(ValueError, match="DS_UNKNOWN"):
        compute_all_stats(condensed_path, metadata_path)


def test_cli_prints_table_and_writes_output(dataset_files, tmp_path):
    condensed_path, metadata_path = dataset_files
    output_path = tmp_path / "stats.csv"

    result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path), "--output", str(output_path)])

    assert result.exit_code == 0
    assert "Community (IGVF)" in result.output
    assert output_path.exists()


def test_cli_reports_missing_metadata_as_click_error(dataset_files):
    condensed_path, metadata_path = dataset_files
    _write_condensed(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5")],
    )

    result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path)])

    assert result.exit_code == 1
    assert "DS_UNKNOWN" in result.output

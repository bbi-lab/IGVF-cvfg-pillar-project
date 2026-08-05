import pandas as pd
import pytest
from click.testing import CliRunner

from src.mave_dataset_stats import (
    GNOMAD_LABEL,
    NO_ANNOTATION_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    SNV_ACCESSIBLE_LABEL,
    SNV_LABEL,
    VUS_LABEL,
    compute_all_stats,
    distinct_variant_flags,
    format_clinical_table,
    has_any_value,
    is_snv_accessible,
    main,
    matches_any_value,
    split_genes,
    stats_to_dataframe,
    summarize_clinical_flags,
    summarize_flags,
    variant_flags,
)

BASE_COLUMNS = ["Dataset", "Gene", "hgvs_c", "hgvs_p"]
ANNOTATION_COLUMNS = [
    "REVEL",
    "AM_score",
    "MutPred2",
    "clinvar_sig_2025",
    "gnomad_MAF",
    "transcript_ref",
    "transcript_alt",
]
FULL_COLUMNS = BASE_COLUMNS + ANNOTATION_COLUMNS


def _write_condensed(path, rows):
    pd.DataFrame(rows, columns=BASE_COLUMNS).to_csv(path, sep="\t", index=False)


def _write_full_variant_file(path, rows):
    pd.DataFrame(rows, columns=FULL_COLUMNS).to_csv(path, sep="\t", index=False)


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


@pytest.fixture
def full_dataset_files(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    expanded_path = tmp_path / "expanded.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    # p3 (hgvs_c "c3") has no ClinVar/gnomAD annotation and is a SNV (A>G).
    # p4's DNA-level candidates split into c4a (annotated, SNV) and c4b
    # (unannotated, a 2-base substitution -- not a SNV), so the "of which
    # SNV" sub-breakdown differs between the assayed (100%) and DNA (50%)
    # levels.
    _write_full_variant_file(
        condensed_path,
        [
            ("DS_IGVF_A", "GENEA", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "", "A", "G"),
            ("DS_IGVF_A", "GENEA", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "", "A", "G"),
            ("DS_IGVF_B", "GENEB, GENEC", "c2", "p2", "", "", "", "Uncertain significance", "0.001", "A", "T"),
            ("DS_COMM_A", "GENEC", "c3", "p3", "", "", "", "", "", "A", "G"),
            ("DS_COMM_B", "GENED", "c4a|c4b", "p4", "0.2|", "0.2|", "0.9|0.9", "Benign|", "", "AC|A", "GT|G"),
        ],
    )
    _write_full_variant_file(
        expanded_path,
        [
            ("DS_IGVF_A", "GENEA", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "", "A", "G"),
            ("DS_IGVF_B", "GENEB, GENEC", "c2", "p2", "", "", "", "Uncertain significance", "0.001", "A", "T"),
            ("DS_COMM_A", "GENEC", "c3", "p3", "", "", "", "", "", "A", "G"),
            ("DS_COMM_B", "GENED", "c4a", "p4", "0.2", "0.2", "0.9", "Benign", "", "AC", "GT"),
            ("DS_COMM_B", "GENED", "c4b", "p4", "", "", "0.9", "", "", "AC", "GT"),
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
    return condensed_path, metadata_path, expanded_path


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


def test_stats_to_dataframe_adds_pct_of_combined_measurements(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats = compute_all_stats(condensed_path, metadata_path)
    table = stats_to_dataframe(stats)

    # 3 of the combined 4 variant_effect_measurements are IGVF, 1 is non-IGVF.
    assert table.loc["Community (IGVF)", "pct_variant_effect_measurements"] == pytest.approx(75.0)
    assert table.loc["Community (non-IGVF)", "pct_variant_effect_measurements"] == pytest.approx(25.0)
    assert table.loc["Combined (IGVF + community)", "pct_variant_effect_measurements"] == pytest.approx(100.0)


def test_merge_calm_genes_collapses_calm_paralogs(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_COMM_A", "CALM1", "c1", "p1"),
            ("DS_COMM_B", "CALM2, CALM3", "c2", "p2"),
            ("DS_COMM_C", "GENED", "c3", "p3"),
        ],
    )
    _write_metadata(
        metadata_path,
        [
            ("DS_COMM_A", "No", "primary score set"),
            ("DS_COMM_B", "No", "primary score set"),
            ("DS_COMM_C", "No", "primary score set"),
        ],
    )

    default_stats = compute_all_stats(condensed_path, metadata_path)
    assert default_stats["Community (non-IGVF)"]["genes_represented"] == 4

    merged_stats = compute_all_stats(condensed_path, metadata_path, merge_calm_genes=True)
    assert merged_stats["Community (non-IGVF)"]["genes_represented"] == 2


def test_compute_all_stats_raises_on_missing_metadata(dataset_files):
    condensed_path, metadata_path = dataset_files
    _write_condensed(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5")],
    )
    with pytest.raises(ValueError, match="DS_UNKNOWN"):
        compute_all_stats(condensed_path, metadata_path)


def test_has_any_value_and_matches_any_value_handle_pipe_delimited_parts():
    values = pd.Series(["0.5", "", "0.2|", "|0.3", "|||"])
    assert list(has_any_value(values)) == [True, False, True, True, False]

    sig = pd.Series(["Pathogenic", "Benign|", "|Uncertain significance", ""])
    assert list(matches_any_value(sig, {"Pathogenic", "Benign"})) == [True, True, False, False]


def test_is_snv_accessible_checks_any_pipe_pair():
    ref = pd.Series(["A", "AC", "A|AC", ""])
    alt = pd.Series(["G", "GT", "G|GT", ""])
    # row 0: single-base pair -> SNV. row 1: two-base pair -> not a SNV.
    # row 2: mixed -- the "A"/"G" pair is a SNV, so this is SNV-accessible.
    # row 3: empty -> not a SNV.
    assert list(is_snv_accessible(ref, alt)) == [True, False, True, False]


def test_variant_flags_and_distinct_variant_flags():
    df = pd.DataFrame(
        {
            "hgvs_c": ["c1", "c1", "c2|c3"],
            "hgvs_p": ["p1", "p1", "p2"],
            "REVEL": ["0.5", "", "|0.4"],
            "AM_score": ["", "", ""],
            "MutPred2": ["", "", ""],
            "clinvar_sig_2025": ["Pathogenic", "", ""],
            "gnomad_MAF": ["", "", ""],
            "transcript_ref": ["A", "A", "AC|A"],
            "transcript_alt": ["G", "G", "GT|G"],
        }
    )

    row_flags = variant_flags(df, SNV_ACCESSIBLE_LABEL)
    assert list(row_flags["REVEL"]) == [True, False, True]
    assert list(row_flags["Pathogenic or benign (ClinVar 2025)"]) == [True, False, False]
    assert list(row_flags["No ClinVar or gnomAD annotation"]) == [False, True, True]
    # row 2's second DNA-level candidate ("A" -> "G") is a SNV, so it's SNV-accessible.
    assert list(row_flags[SNV_ACCESSIBLE_LABEL]) == [True, True, True]

    # (c1, p1) is one distinct variant measured twice; (c2|c3, p2) is a second,
    # pipe-delimited (multi-DNA-candidate) distinct variant.
    distinct = distinct_variant_flags(df, SNV_ACCESSIBLE_LABEL)
    assert len(distinct) == 2
    assert distinct["REVEL"].sum() == 2


def test_summarize_flags_reports_count_and_pct():
    flags = pd.DataFrame({"REVEL": [True, True, False, False]})
    total, table = summarize_flags(flags)
    assert total == 4
    assert table.loc["REVEL", "count"] == 2
    assert table.loc["REVEL", "pct"] == 50.0


def test_summarize_and_format_clinical_table_breaks_out_snv_per_bucket():
    flags = pd.DataFrame(
        {
            VUS_LABEL: [True, False, False, False, False],
            PATHOGENIC_OR_BENIGN_LABEL: [False, True, False, False, False],
            GNOMAD_LABEL: [False, False, True, False, False],
            NO_ANNOTATION_LABEL: [False, False, False, True, True],
            SNV_LABEL: [True, False, True, True, False],
        }
    )

    total, snv_total, table = summarize_clinical_flags(flags, SNV_LABEL)

    assert total == 5
    assert snv_total == 3
    assert table.loc[VUS_LABEL, SNV_LABEL] == 1  # the one VUS row is a SNV
    assert table.loc[PATHOGENIC_OR_BENIGN_LABEL, SNV_LABEL] == 0  # the one path/benign row isn't
    assert table.loc[NO_ANNOTATION_LABEL, SNV_LABEL] == 1  # one of the two unannotated rows is
    assert table.loc[VUS_LABEL, f"% of {SNV_LABEL}"] == pytest.approx(100 / 3, abs=0.05)

    text = format_clinical_table("Title", total, snv_total, table, SNV_LABEL)

    assert "Total: 5 (3 SNV)" in text
    assert f"% of {SNV_LABEL}" in text


def test_cli_prints_table_and_writes_output(full_dataset_files, tmp_path):
    condensed_path, metadata_path, expanded_path = full_dataset_files
    output_path = tmp_path / "report.txt"

    result = CliRunner().invoke(
        main,
        [str(condensed_path), str(metadata_path), str(expanded_path), "--output", str(output_path)],
    )

    assert result.exit_code == 0
    assert "Community (IGVF)" in result.output
    assert "Score coverage" in result.output
    assert "Clinical attributes" in result.output
    # Assayed level: all 4 distinct variants (and all 5 measurement rows) are SNV-accessible.
    assert "Total: 4 (4 SNV-accessible)" in result.output
    assert "Total: 5 (5 SNV-accessible)" in result.output
    # DNA level: 3 of the 5 distinct (and 5 measurement-row) DNA variants are SNVs
    # (c4a/c4b's 2-base substitution isn't one).
    assert "Total: 5 (3 SNV)" in result.output
    assert "% of SNV-accessible" in result.output
    assert "% of SNV" in result.output
    assert output_path.exists()
    assert "Score coverage" in output_path.read_text()


def test_cli_reports_missing_metadata_as_click_error(full_dataset_files):
    condensed_path, metadata_path, expanded_path = full_dataset_files
    _write_full_variant_file(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5", "", "", "", "", "", "A", "G")],
    )

    result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path), str(expanded_path)])

    assert result.exit_code == 1
    assert "DS_UNKNOWN" in result.output

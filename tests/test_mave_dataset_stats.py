import re

import pandas as pd
import pytest
from click.testing import CliRunner

from src.mave_dataset_stats import (
    CLINVAR_CONFLICT_LABEL,
    GNOMAD_LABEL,
    NO_ANNOTATION_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    SNV_ACCESSIBLE_LABEL,
    SNV_LABEL,
    VUS_LABEL,
    clinvar_classification_from_flags,
    clinvar_significance_flags,
    compute_all_stats,
    distinct_variant_flags,
    format_clinical_table,
    has_any_value,
    is_snv_accessible,
    main,
    matches_any_value,
    mixed_year_clinvar_series,
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
    "clinvar_sig_2018",
    "gnomad_MAF",
    "transcript_ref",
    "transcript_alt",
]
FULL_COLUMNS = BASE_COLUMNS + ANNOTATION_COLUMNS


def _write_condensed(path, rows):
    pd.DataFrame(rows, columns=BASE_COLUMNS).to_csv(path, sep="\t", index=False)


def _write_full_variant_file(path, rows):
    pd.DataFrame(rows, columns=FULL_COLUMNS).to_csv(path, sep="\t", index=False)


def _bucket_count(section_text, label):
    """Extract the `count` column for one clinical/score-coverage bucket line."""
    match = re.search(rf"^{re.escape(label)}\s+(\d+)", section_text, re.MULTILINE)
    assert match, f"no line for {label!r} in:\n{section_text}"
    return int(match.group(1))


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
    #
    # p1 is on BRCA1 (a mixed-year gene) with a 2025 call of "Pathogenic" but
    # a 2018 call of "Uncertain significance", so it should land in a
    # different clinical-attribute bucket between the two report sections.
    # p3 is on GENEC (not a mixed-year gene) and carries a 2018 call
    # ("Pathogenic") that must be ignored in favor of its empty 2025 call.
    _write_full_variant_file(
        condensed_path,
        [
            ("DS_IGVF_A", "BRCA1", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "Uncertain significance", "", "A", "G"),
            ("DS_IGVF_A", "BRCA1", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "Uncertain significance", "", "A", "G"),
            (
                "DS_IGVF_B",
                "GENEB, GENEC",
                "c2",
                "p2",
                "",
                "",
                "",
                "Uncertain significance",
                "Uncertain significance",
                "0.001",
                "A",
                "T",
            ),
            ("DS_COMM_A", "GENEC", "c3", "p3", "", "", "", "", "Pathogenic", "", "A", "G"),
            ("DS_COMM_B", "GENED", "c4a|c4b", "p4", "0.2|", "0.2|", "0.9|0.9", "Benign|", "Benign|", "", "AC|A", "GT|G"),
        ],
    )
    _write_full_variant_file(
        expanded_path,
        [
            ("DS_IGVF_A", "BRCA1", "c1", "p1", "0.5", "0.5", "0.5", "Pathogenic", "Uncertain significance", "", "A", "G"),
            (
                "DS_IGVF_B",
                "GENEB, GENEC",
                "c2",
                "p2",
                "",
                "",
                "",
                "Uncertain significance",
                "Uncertain significance",
                "0.001",
                "A",
                "T",
            ),
            ("DS_COMM_A", "GENEC", "c3", "p3", "", "", "", "", "Pathogenic", "", "A", "G"),
            ("DS_COMM_B", "GENED", "c4a", "p4", "0.2", "0.2", "0.9", "Benign", "Benign", "", "AC", "GT"),
            ("DS_COMM_B", "GENED", "c4b", "p4", "", "", "0.9", "", "", "", "AC", "GT"),
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


def test_clinvar_classification_from_flags():
    def classify(value):
        has_pathogenic, has_benign, has_vus, has_literal_conflict = clinvar_significance_flags(pd.Series([value]))
        vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(
            has_pathogenic, has_benign, has_vus, has_literal_conflict
        )
        if conflict.iloc[0]:
            return CLINVAR_CONFLICT_LABEL
        if pathogenic_or_benign.iloc[0]:
            return PATHOGENIC_OR_BENIGN_LABEL
        if vus.iloc[0]:
            return VUS_LABEL
        return None

    assert classify("") is None
    assert classify("Uncertain significance") == VUS_LABEL
    assert classify("Pathogenic") == PATHOGENIC_OR_BENIGN_LABEL
    assert classify("Benign|Likely benign") == PATHOGENIC_OR_BENIGN_LABEL
    # ClinVar's own conflict call is a conflict regardless of what else is present.
    assert classify("Conflicting classifications of pathogenicity") == CLINVAR_CONFLICT_LABEL
    assert classify("Pathogenic|Conflicting classifications of pathogenicity") == CLINVAR_CONFLICT_LABEL
    # Disagreement between a pathogenic-leaning and benign-leaning call is also a conflict.
    assert classify("Pathogenic|Benign") == CLINVAR_CONFLICT_LABEL
    # VUS alongside a definitive call isn't flagged as a conflict by this script (unlike, say, a
    # pathogenic/benign disagreement) -- the definitive call wins.
    assert classify("Uncertain significance|Pathogenic") == PATHOGENIC_OR_BENIGN_LABEL


def test_clinvar_classification_from_flags_resolves_conflicts_across_grouped_rows():
    df = pd.DataFrame(
        {
            "hgvs_c": ["c1", "c1", "c2"],
            "hgvs_p": ["p1", "p1", "p2"],
        }
    )
    # (c1, p1)'s two rows disagree (Pathogenic vs. Benign) -- a cross-row conflict that
    # neither row's own pipe-delimited parts would reveal on its own.
    clinvar_series = pd.Series(["Pathogenic", "Benign", "Uncertain significance"])

    has_pathogenic, has_benign, has_vus, has_literal_conflict = clinvar_significance_flags(clinvar_series)
    keyed = pd.concat(
        [df, has_pathogenic.rename("p"), has_benign.rename("b"), has_vus.rename("v"), has_literal_conflict.rename("c")],
        axis=1,
    )
    grouped = keyed.groupby(["hgvs_c", "hgvs_p"], as_index=False).any()
    vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(
        grouped["p"], grouped["b"], grouped["v"], grouped["c"]
    )
    grouped["vus"], grouped["pathogenic_or_benign"], grouped["conflict"] = vus, pathogenic_or_benign, conflict
    grouped = grouped.set_index(["hgvs_c", "hgvs_p"])

    assert grouped.loc[("c1", "p1"), "conflict"]
    assert not grouped.loc[("c1", "p1"), "pathogenic_or_benign"]
    assert grouped.loc[("c2", "p2"), "vus"]


def test_variant_flags_excludes_clinvar_conflicts_by_default():
    df = pd.DataFrame(
        {
            "hgvs_c": ["c1", "c2"],
            "hgvs_p": ["p1", "p2"],
            "REVEL": ["", ""],
            "AM_score": ["", ""],
            "MutPred2": ["", ""],
            # row 0's own pipe-delimited candidates disagree; row 1 has no conflict.
            "clinvar_sig_2025": ["Pathogenic|Benign", "Pathogenic"],
            "gnomad_MAF": ["", ""],
            "transcript_ref": ["A", "A"],
            "transcript_alt": ["G", "G"],
        }
    )

    default_flags = variant_flags(df, SNV_ACCESSIBLE_LABEL)
    assert list(default_flags[CLINVAR_CONFLICT_LABEL]) == [True, False]
    assert list(default_flags[VUS_LABEL]) == [False, False]
    assert list(default_flags[PATHOGENIC_OR_BENIGN_LABEL]) == [False, True]
    # A conflicting row isn't double-counted as "no annotation" either.
    assert list(default_flags[NO_ANNOTATION_LABEL]) == [False, False]

    legacy_flags = variant_flags(df, SNV_ACCESSIBLE_LABEL, allow_clinvar_conflicts=True)
    assert CLINVAR_CONFLICT_LABEL not in legacy_flags.columns
    # Any-match folds the conflicting row into pathogenic-or-benign instead.
    assert list(legacy_flags[PATHOGENIC_OR_BENIGN_LABEL]) == [True, True]


def test_distinct_variant_flags_resolves_conflicts_across_rows():
    df = pd.DataFrame(
        {
            "hgvs_c": ["c1", "c1", "c2"],
            "hgvs_p": ["p1", "p1", "p2"],
            "REVEL": ["", "", ""],
            "AM_score": ["", "", ""],
            "MutPred2": ["", "", ""],
            # (c1, p1)'s two measurement rows disagree; (c2, p2) has a single, clean call.
            "clinvar_sig_2025": ["Pathogenic", "Benign", "Benign"],
            "gnomad_MAF": ["", "", ""],
            "transcript_ref": ["A", "A", "A"],
            "transcript_alt": ["G", "G", "G"],
        }
    )

    default_distinct = distinct_variant_flags(df, SNV_ACCESSIBLE_LABEL)
    assert len(default_distinct) == 2
    conflict_row = default_distinct[default_distinct[CLINVAR_CONFLICT_LABEL]]
    assert len(conflict_row) == 1
    assert not conflict_row[PATHOGENIC_OR_BENIGN_LABEL].iloc[0]
    clean_row = default_distinct[~default_distinct[CLINVAR_CONFLICT_LABEL]]
    assert clean_row[PATHOGENIC_OR_BENIGN_LABEL].iloc[0]

    # Any-match (the legacy behavior) instead counts (c1, p1) as pathogenic-or-benign,
    # since at least one of its rows matches.
    legacy_distinct = distinct_variant_flags(df, SNV_ACCESSIBLE_LABEL, allow_clinvar_conflicts=True)
    assert CLINVAR_CONFLICT_LABEL not in legacy_distinct.columns
    assert legacy_distinct[PATHOGENIC_OR_BENIGN_LABEL].sum() == 2


def test_mixed_year_clinvar_series_swaps_only_mixed_year_genes():
    df = pd.DataFrame(
        {
            "Gene": ["BRCA1", "PTEN", "GENEC", "CALM1, CALM2, CALM3"],
            "clinvar_sig_2025": ["Pathogenic", "Benign", "", "Uncertain significance"],
            "clinvar_sig_2018": ["Uncertain significance", "Uncertain significance", "Pathogenic", "Pathogenic"],
        }
    )

    result = mixed_year_clinvar_series(df)

    # BRCA1 and PTEN are mixed-year genes -> use the 2018 call.
    assert list(result.iloc[:2]) == ["Uncertain significance", "Uncertain significance"]
    # GENEC and the CALM1/2/3 combination aren't -> use the 2025 call.
    assert list(result.iloc[2:]) == ["", "Uncertain significance"]


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

    # The mixed-year section reclassifies p1 (BRCA1) from pathogenic/benign
    # (its 2025 call) to VUS (its 2018 call), so at the assayed-variants,
    # distinct level the VUS count rises (1 -> 2) and the pathogenic/benign
    # count falls (2 -> 1) between the two clinical-attribute sections.
    clinvar_2025_section, clinvar_mixed_section = result.output.split(
        "=== Clinical attributes (ClinVar 2025, except ClinVar 2018 for BRCA1/PTEN/MSH2/TP53; "
        "gnomAD; conflicting/ambiguous ClinVar calls excluded) ==="
    )
    assayed_distinct_2025 = clinvar_2025_section.split("Clinical attributes -- assayed variants, distinct")[1]
    assayed_distinct_mixed = clinvar_mixed_section.split("Clinical attributes -- assayed variants, distinct")[1]
    assert _bucket_count(assayed_distinct_2025, VUS_LABEL) == 1
    assert _bucket_count(assayed_distinct_2025, PATHOGENIC_OR_BENIGN_LABEL) == 2
    assert _bucket_count(assayed_distinct_mixed, VUS_LABEL) == 2
    assert _bucket_count(assayed_distinct_mixed, PATHOGENIC_OR_BENIGN_LABEL) == 1

    assert output_path.exists()
    output_text = output_path.read_text()
    assert "Score coverage" in output_text
    assert "ClinVar 2025, except ClinVar 2018 for BRCA1/PTEN/MSH2/TP53" in output_text


def test_cli_allow_clinvar_conflicts_flag_toggles_conflict_handling(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    expanded_path = tmp_path / "expanded.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    # p1's two measurement rows disagree (Pathogenic vs. Benign) -- a conflict
    # that only shows up once the rows are grouped into one distinct variant.
    rows = [
        ("DS_A", "GENEA", "c1", "p1", "", "", "", "Pathogenic", "Pathogenic", "", "A", "G"),
        ("DS_A", "GENEA", "c1", "p1", "", "", "", "Benign", "Benign", "", "A", "G"),
    ]
    _write_full_variant_file(condensed_path, rows)
    _write_full_variant_file(expanded_path, rows)
    _write_metadata(metadata_path, [("DS_A", "No", "primary score set")])

    default_result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path), str(expanded_path)])
    assert default_result.exit_code == 0
    assert "conflicting/ambiguous ClinVar calls excluded" in default_result.output
    default_distinct_section = default_result.output.split("Clinical attributes -- assayed variants, distinct")[1]
    assert _bucket_count(default_distinct_section, CLINVAR_CONFLICT_LABEL) == 1
    assert _bucket_count(default_distinct_section, PATHOGENIC_OR_BENIGN_LABEL) == 0

    legacy_result = CliRunner().invoke(
        main, [str(condensed_path), str(metadata_path), str(expanded_path), "--allow-clinvar-conflicts"]
    )
    assert legacy_result.exit_code == 0
    assert "conflicting/ambiguous ClinVar calls folded in via any-match" in legacy_result.output
    assert CLINVAR_CONFLICT_LABEL not in legacy_result.output
    legacy_distinct_section = legacy_result.output.split("Clinical attributes -- assayed variants, distinct")[1]
    assert _bucket_count(legacy_distinct_section, PATHOGENIC_OR_BENIGN_LABEL) == 1


def test_cli_reports_missing_metadata_as_click_error(full_dataset_files):
    condensed_path, metadata_path, expanded_path = full_dataset_files
    _write_full_variant_file(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5", "", "", "", "", "", "", "A", "G")],
    )

    result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path), str(expanded_path)])

    assert result.exit_code == 1
    assert "DS_UNKNOWN" in result.output

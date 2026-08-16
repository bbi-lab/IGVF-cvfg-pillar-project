import re
import unicodedata

import pandas as pd
import pytest
from click.testing import CliRunner

from src.mave_dataset_stats import (
    AGREE_LABEL,
    CALM_MERGED_LABEL,
    CLINVAR_CONFLICT_LABEL,
    DISAGREE_LABEL,
    GNOMAD_LABEL,
    NO_ANNOTATION_LABEL,
    NO_EVIDENCE_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    RECLASSIFICATION_USECOLS,
    RECLASSIFICATION_VARIANT_CLASSIFICATION_TITLE,
    SNV_ACCESSIBLE_LABEL,
    SNV_LABEL,
    VARIANT_CLASSIFICATION_TITLE,
    VUS_LABEL,
    build_reclassification_report,
    clinvar_classification_from_flags,
    clinvar_significance_flags,
    compute_all_stats,
    compute_excalibr_calibration_stats,
    compute_reclassification_agreement,
    compute_variant_classification_stats,
    compute_variant_classification_stats_from_reclassification_file,
    distinct_dna_variants,
    distinct_variant_flags,
    excalibr_dataset_to_gene_map,
    format_calibration_summary,
    format_clinical_table,
    format_reclassification_table,
    format_variant_classification_summary,
    has_any_value,
    is_snv_accessible,
    load_dataset_metadata,
    main,
    matches_any_value,
    mixed_year_clinvar_series,
    points_are_pathogenic_or_benign,
    reclassification_flags,
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


def _write_metadata(path, rows, genes=None):
    """`genes`, if given, maps `Dataset Name` -> `Gene` and adds a Gene column --
    needed only by tests that exercise the ExCALIBR-calibration code path, which
    looks up genes via metadata rather than the condensed/expanded file.
    """
    columns = ["Dataset Name", "IGVF Produced?", "Primary Score Set or Meta-analysis?"]
    df = pd.DataFrame(rows, columns=columns)
    if genes is not None:
        df["Gene"] = df["Dataset Name"].map(genes)
    with pd.ExcelWriter(path) as writer:
        df.to_excel(writer, sheet_name="Curation", index=False)


def _write_excalibr_calibrations(path, rows):
    """`rows` is a list of (dataset, range_-1, range_1) tuples; `range_-1`/`range_1`
    hold a value (any non-null placeholder) when that dataset has benign/pathogenic
    evidence, or None when it doesn't.
    """
    df = pd.DataFrame(rows, columns=["dataset", "range_-1", "range_1"])
    df.to_excel(path, sheet_name="ExCALIBR_calibrations", index=False)


def _write_controls_file(path, sheets):
    """`sheets` is {sheet_name: rows}, where each row is
    (clnsig_group_18_25, ExC_points_2025, OP_points).
    """
    columns = ["clnsig_group_18_25", "ExC_points_2025", "OP_points"]
    with pd.ExcelWriter(path) as writer:
        for sheet_name, rows in sheets.items():
            pd.DataFrame(rows, columns=columns).to_excel(writer, sheet_name=sheet_name, index=False)


def _write_variant_classification_sheets(path, sheets, mode="a"):
    """`sheets` is {sheet_name: rows}, where each row is
    (Gene, Chrom, hg38_start, ref_allele, alt_allele, Class_REVEL) -- the
    columns `compute_variant_classification_stats` reads. Appends (`mode="a"`)
    to an existing workbook at `path` (e.g. one already written by
    `_write_controls_file`), since the main CLI reads both sets of sheets
    from the same `--controls-file`; pass `mode="w"` to create a fresh file.

    Also adds `clnsig_group_18_25`/`ExC_points_2025`/`OP_points` columns
    (empty/no-evidence), since the real Supplementary_Data_5 `controls_`-
    prefixed sheets carry both this function's columns and those
    `build_reclassification_report` reads -- one of `sheets`' names
    (`controls_REVEL_GeneSpecific`) matches `CONTROLS_SHEET_PREFIX`, so
    without them `build_reclassification_report` would KeyError on it.
    """
    columns = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele", "Class_REVEL"]
    with pd.ExcelWriter(path, mode=mode, engine="openpyxl") as writer:
        for sheet_name, rows in sheets.items():
            df = pd.DataFrame(rows, columns=columns)
            df["clnsig_group_18_25"] = None
            df["ExC_points_2025"] = 0
            df["OP_points"] = 0
            df.to_excel(writer, sheet_name=sheet_name, index=False)


def _write_reclassification_file(path, rows):
    """`rows` is a list of (clinvar_sig_2025, gnomad_MAF, ref_allele,
    alt_allele, Combined_points) tuples -- the columns
    `compute_variant_classification_stats_from_reclassification_file` reads.
    """
    columns = ["clinvar_sig_2025", "gnomad_MAF", "ref_allele", "alt_allele", "Combined_points"]
    pd.DataFrame(rows, columns=columns).to_csv(path, sep="\t", index=False)


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
        genes={"DS_IGVF_A": "BRCA1", "DS_IGVF_B": "GENEB, GENEC", "DS_COMM_A": "GENEC", "DS_COMM_B": "GENED"},
    )

    excalibr_path = tmp_path / "excalibr_calibrations.xlsx"
    # DS_IGVF_A -> BRCA1 has evidence; the "_clinvar_2018"-suffixed DS_IGVF_B row
    # (-> GENEB, GENEC) and DS_COMM_A (no row at all) don't; DS_COMM_B -> GENED does.
    # So 4 genes are calibrated (BRCA1, GENEB, GENEC, GENED) and 2 have evidence.
    _write_excalibr_calibrations(
        excalibr_path,
        [
            ("DS_IGVF_A", None, "0.5 1"),
            ("DS_IGVF_B_clinvar_2018", None, None),
            ("DS_COMM_B", "-1 -0.5", None),
        ],
    )

    controls_path = tmp_path / "controls.xlsx"
    _write_controls_file(
        controls_path,
        {
            "controls_TEST_GeneSpecific": [
                ("Pathogenic", 3, 2),  # agree on both
                ("Benign", -2, -1),  # agree on both
                ("Pathogenic", -1, 0),  # ExC disagrees; OP has no evidence
                ("Likely benign", 0, 1),  # ExC has no evidence; OP disagrees
                ("Benign/Likely benign", -4, -2),  # agree on both
            ]
        },
    )
    # coords (1, 100, A, G) is shared between "controls" and "gnomAD" below, to
    # exercise the combined total's cross-category dedup -- it should count
    # once, not twice, in the combined "Distinct DNA variants classified" total.
    _write_variant_classification_sheets(
        controls_path,
        {
            "controls_REVEL_GeneSpecific": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic"),
                ("GENEX", 1, 101, "A", "G", "Benign"),
            ],
            "ClinGen_Repo_REVEL_GeneSpecific": [
                ("GENEX", 1, 102, "A", "G", "Likely Pathogenic"),
            ],
            "VUS_REVEL": [
                ("GENEX", 1, 103, "A", "G", "Pathogenic"),  # resolved
                ("GENEX", 1, 104, "A", "G", "Uncertain"),  # not resolved
                ("GENEX", 1, 105, "A", "G", "Benign"),  # resolved
            ],
            "gnomAD_REVEL": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic"),  # duplicate of controls row above
            ],
            "Unobserved_REVEL": [
                ("GENEX", 1, 106, "A", "G", "Pathogenic"),  # resolved
                ("GENEX", 1, 107, "A", "G", "Uncertain"),  # not resolved
            ],
        },
    )

    # 7 rows total. Pathogenic-or-benign (points >= 6 or <= -1): rows 1, 3, 5,
    # 6, 7 (5 of 7). VUS (clinvar_sig_2025 == "Uncertain significance"): rows
    # 1-2, of which row 1 (points 8) resolves (1 of 2). Unobserved (both
    # clinvar_sig_2025 and gnomad_MAF null, SNV only): rows 3-4 (row 7 is
    # excluded -- its ref/alt aren't single-base; row 5/6 have a
    # clinvar_sig_2025/gnomad_MAF value), of which row 3 (points -8) resolves
    # (1 of 2).
    reclassification_path = tmp_path / "reclassification.tsv"
    _write_reclassification_file(
        reclassification_path,
        [
            ("Uncertain significance", None, "A", "G", 8),  # VUS, resolved
            ("Uncertain significance", None, "A", "G", 2),  # VUS, not resolved
            (None, None, "A", "G", -8),  # Unobserved SNV, resolved
            (None, None, "A", "G", 1),  # Unobserved SNV, not resolved
            ("Pathogenic", None, "A", "G", 11),  # controls-like, not VUS/Unobserved
            (None, 0.01, "A", "G", -3),  # gnomAD-like, not Unobserved (gnomad_MAF set)
            (None, None, "AC", "G", 8),  # not a SNV, excluded from Unobserved
        ],
    )

    return condensed_path, metadata_path, expanded_path, excalibr_path, controls_path, reclassification_path


def test_split_genes_handles_multi_gene_datasets():
    assert split_genes("CALM1, CALM2, CALM3") == ["CALM1", "CALM2", "CALM3"]
    assert split_genes("BRCA1") == ["BRCA1"]


def test_compute_all_stats_buckets(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats, gene_breakdown = compute_all_stats(condensed_path, metadata_path)

    igvf = stats["IGVF"]
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

    # GENEA (DS_IGVF_A) and GENEB (DS_IGVF_B) are IGVF-only; GENED (DS_COMM_B) is
    # community-only; GENEC is shared (DS_IGVF_B and DS_COMM_A).
    assert gene_breakdown["IGVF only"] == ["GENEA", "GENEB"]
    assert gene_breakdown["Community (non-IGVF) only"] == ["GENED"]
    assert gene_breakdown["Both IGVF and community (non-IGVF)"] == ["GENEC"]


def test_stats_to_dataframe_adds_pct_of_combined_measurements(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats, _ = compute_all_stats(condensed_path, metadata_path)
    table = stats_to_dataframe(stats)

    # 3 of the combined 4 variant_effect_measurements are IGVF, 1 is non-IGVF.
    assert table.loc["IGVF", "pct_variant_effect_measurements"] == pytest.approx(75.0)
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

    default_stats, default_gene_breakdown = compute_all_stats(condensed_path, metadata_path)
    assert default_stats["Community (non-IGVF)"]["genes_represented"] == 4
    # The gene breakdown list always merges CALM1/2/3, regardless of --merge-calm-genes.
    assert default_gene_breakdown["Community (non-IGVF) only"] == [CALM_MERGED_LABEL, "GENED"]

    merged_stats, merged_gene_breakdown = compute_all_stats(condensed_path, metadata_path, merge_calm_genes=True)
    assert merged_stats["Community (non-IGVF)"]["genes_represented"] == 2
    assert merged_gene_breakdown["Community (non-IGVF) only"] == [CALM_MERGED_LABEL, "GENED"]


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
    condensed_path, metadata_path, expanded_path, excalibr_path, controls_path, reclassification_path = (
        full_dataset_files
    )
    output_path = tmp_path / "report.txt"

    result = CliRunner().invoke(
        main,
        [
            str(condensed_path),
            str(metadata_path),
            str(expanded_path),
            "--excalibr-calibrations-file",
            str(excalibr_path),
            "--controls-file",
            str(controls_path),
            "--reclassification-file",
            str(reclassification_path),
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert "IGVF" in result.output
    assert "Score coverage" in result.output
    assert "Clinical attributes" in result.output

    # BRCA1 (DS_IGVF_A) and GENEB (DS_IGVF_B) are IGVF-only; GENED (DS_COMM_B) is
    # community-only; GENEC is shared (DS_IGVF_B and DS_COMM_A).
    genes_section = result.output.split("=== Genes represented ===")[1].split("=== Score coverage")[0]
    assert "IGVF only (2): BRCA1, GENEB" in genes_section
    assert "Community (non-IGVF) only (1): GENED" in genes_section
    assert "Both IGVF and community (non-IGVF) (1): GENEC" in genes_section

    # 4 genes (BRCA1, GENEB, GENEC, GENED) are calibrated; 2 (BRCA1, GENED) have evidence.
    assert "Genes with ExCALIBR calibrations: 4" in result.output
    assert "Genes with >=1 dataset assigning >=1 point of evidence (pathogenic or benign): 2 (50.0%)" in result.output

    # controls_TEST_GeneSpecific: ExC_points_2025 agrees on 3/5, disagrees on 1, no evidence on 1
    # -> 3/4 determinate calls agree; OP_points agrees on 3/5, disagrees on 1, no evidence on 1
    # -> same 3/4 determinate agreement (see full_dataset_files' controls_path rows).
    reclassification_section = result.output.split("=== Reclassification agreement (Figure 4c) ===")[1]
    excalibr_reclass, functional_reclass = reclassification_section.split(
        "ExCALIBR evidence -- controls_TEST_GeneSpecific"
    )[1].split("Functional class -- controls_TEST_GeneSpecific")
    assert "Total control variants: 5" in excalibr_reclass
    assert _bucket_count(excalibr_reclass, AGREE_LABEL) == 3
    assert _bucket_count(excalibr_reclass, DISAGREE_LABEL) == 1
    assert _bucket_count(excalibr_reclass, NO_EVIDENCE_LABEL) == 1
    assert "Agreement with ClinVar PLP/BLB (of determinate calls): 75.0%" in excalibr_reclass
    assert _bucket_count(functional_reclass, AGREE_LABEL) == 3
    assert "Agreement with ClinVar PLP/BLB (of determinate calls): 75.0%" in functional_reclass

    # 8 distinct DNA variants across the 5 REVEL sheets combined (9 rows, minus
    # the (GENEX, 1, 100, A, G) coordinate shared by controls_REVEL_GeneSpecific
    # and gnomAD_REVEL); 6 of those 8 are Pathogenic/Likely Pathogenic/Benign/
    # Likely Benign (the two "Uncertain" VUS/Unobserved rows aren't).
    variant_classification_section = result.output.split(VARIANT_CLASSIFICATION_TITLE)[1].split(
        RECLASSIFICATION_VARIANT_CLASSIFICATION_TITLE
    )[0]
    assert "Distinct DNA variants classified: 8" in variant_classification_section
    assert "Pathogenic or benign: 6 of 8 (75.0%)" in variant_classification_section
    assert "ClinVar VUS resolved (reclassified pathogenic or benign): 2 of 3 (66.7%)" in variant_classification_section
    assert (
        "Unobserved variants resolved (classified pathogenic or benign): 1 of 2 (50.0%)"
        in variant_classification_section
    )

    # See full_dataset_files' reclassification_path rows for the expected counts.
    reclassification_file_section = result.output.split(RECLASSIFICATION_VARIANT_CLASSIFICATION_TITLE)[1]
    assert "Distinct DNA variants classified: 7" in reclassification_file_section
    assert "Pathogenic or benign: 5 of 7 (71.4%)" in reclassification_file_section
    assert "ClinVar VUS resolved (reclassified pathogenic or benign): 1 of 2 (50.0%)" in reclassification_file_section
    assert (
        "Unobserved variants resolved (classified pathogenic or benign): 1 of 2 (50.0%)"
        in reclassification_file_section
    )

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
    _write_metadata(metadata_path, [("DS_A", "No", "primary score set")], genes={"DS_A": "GENEA"})

    excalibr_path = tmp_path / "excalibr_calibrations.xlsx"
    _write_excalibr_calibrations(excalibr_path, [("DS_A", None, None)])
    controls_path = tmp_path / "controls.xlsx"
    _write_controls_file(controls_path, {"controls_TEST": [("Pathogenic", 1, 1)]})
    _write_variant_classification_sheets(
        controls_path,
        {
            "controls_REVEL_GeneSpecific": [("GENEA", 1, 1, "A", "G", "Pathogenic")],
            "ClinGen_Repo_REVEL_GeneSpecific": [],
            "VUS_REVEL": [],
            "gnomAD_REVEL": [],
            "Unobserved_REVEL": [],
        },
    )
    reclassification_path = tmp_path / "reclassification.tsv"
    _write_reclassification_file(reclassification_path, [("Pathogenic", None, "A", "G", 11)])
    extra_args = [
        "--excalibr-calibrations-file",
        str(excalibr_path),
        "--controls-file",
        str(controls_path),
        "--reclassification-file",
        str(reclassification_path),
    ]

    default_result = CliRunner().invoke(
        main, [str(condensed_path), str(metadata_path), str(expanded_path), *extra_args]
    )
    assert default_result.exit_code == 0
    assert "conflicting/ambiguous ClinVar calls excluded" in default_result.output
    default_distinct_section = default_result.output.split("Clinical attributes -- assayed variants, distinct")[1]
    assert _bucket_count(default_distinct_section, CLINVAR_CONFLICT_LABEL) == 1
    assert _bucket_count(default_distinct_section, PATHOGENIC_OR_BENIGN_LABEL) == 0

    legacy_result = CliRunner().invoke(
        main, [str(condensed_path), str(metadata_path), str(expanded_path), "--allow-clinvar-conflicts", *extra_args]
    )
    assert legacy_result.exit_code == 0
    assert "conflicting/ambiguous ClinVar calls folded in via any-match" in legacy_result.output
    assert CLINVAR_CONFLICT_LABEL not in legacy_result.output
    legacy_distinct_section = legacy_result.output.split("Clinical attributes -- assayed variants, distinct")[1]
    assert _bucket_count(legacy_distinct_section, PATHOGENIC_OR_BENIGN_LABEL) == 1


def test_cli_reports_missing_metadata_as_click_error(full_dataset_files):
    condensed_path, metadata_path, expanded_path, _excalibr_path, _controls_path, _reclassification_path = (
        full_dataset_files
    )
    _write_full_variant_file(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5", "", "", "", "", "", "", "A", "G")],
    )

    result = CliRunner().invoke(main, [str(condensed_path), str(metadata_path), str(expanded_path)])

    assert result.exit_code == 1
    assert "DS_UNKNOWN" in result.output


def _write_gene_metadata(path, gene_by_dataset):
    df = pd.DataFrame({"Dataset Name": list(gene_by_dataset), "Gene": list(gene_by_dataset.values())})
    with pd.ExcelWriter(path) as writer:
        df.to_excel(writer, sheet_name="Curation", index=False)


def test_excalibr_dataset_to_gene_map_strips_clinvar_2018_suffix(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    _write_gene_metadata(metadata_path, {"BRCA1_Findlay_2018": "BRCA1"})
    metadata = load_dataset_metadata(metadata_path)

    mapping = excalibr_dataset_to_gene_map(["BRCA1_Findlay_2018_clinvar_2018"], metadata)

    assert mapping == {"BRCA1_Findlay_2018_clinvar_2018": "BRCA1"}


def test_excalibr_dataset_to_gene_map_normalizes_unicode(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    # Precomposed accented character (single codepoint), as Supplementary_Data_3 stores it.
    precomposed = unicodedata.normalize("NFC", "RAD51C_Olvera-Le\u00f3n_2024")
    _write_gene_metadata(metadata_path, {precomposed: "RAD51C"})
    metadata = load_dataset_metadata(metadata_path)

    # Decomposed base letter + combining acute accent, as the ExCALIBR_calibrations sheet stores it.
    decomposed = unicodedata.normalize("NFD", precomposed)
    assert decomposed != precomposed  # sanity check that the two forms really do differ
    mapping = excalibr_dataset_to_gene_map([decomposed], metadata)

    assert mapping == {decomposed: "RAD51C"}


def test_excalibr_dataset_to_gene_map_raises_on_unmapped_dataset(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    _write_gene_metadata(metadata_path, {"BRCA2_IGVF": "BRCA2"})
    metadata = load_dataset_metadata(metadata_path)

    with pytest.raises(ValueError, match="BRCA2_Huang_2026"):
        excalibr_dataset_to_gene_map(["BRCA2_IGVF", "BRCA2_Huang_2026"], metadata)


def test_compute_excalibr_calibration_stats_counts_genes_with_evidence(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    _write_gene_metadata(
        metadata_path,
        {"DS_A": "GENEA", "DS_B": "GENEB, GENEC", "DS_D": "GENED"},
    )
    metadata = load_dataset_metadata(metadata_path)

    calibrations = pd.DataFrame(
        {
            "dataset": ["DS_A", "DS_B_clinvar_2018", "DS_D"],
            "range_-1": [None, None, "-1 -0.5"],
            "range_1": ["0.5 1", None, None],
        }
    )

    stats = compute_excalibr_calibration_stats(calibrations, metadata)

    # 4 genes are calibrated (GENEA, GENEB, GENEC, GENED); GENEA and GENED have
    # evidence, GENEB/GENEC (from the all-null DS_B row) don't.
    assert stats["genes_with_excalibr_calibrations"] == 4
    assert stats["genes_with_evidence_assigned"] == 2


def test_compute_excalibr_calibration_stats_merges_calm_genes(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    _write_gene_metadata(metadata_path, {"DS_CALM": "CALM1, CALM2, CALM3"})
    metadata = load_dataset_metadata(metadata_path)
    calibrations = pd.DataFrame({"dataset": ["DS_CALM"], "range_-1": [None], "range_1": [None]})

    default_stats = compute_excalibr_calibration_stats(calibrations, metadata)
    assert default_stats["genes_with_excalibr_calibrations"] == 3

    merged_stats = compute_excalibr_calibration_stats(calibrations, metadata, merge_calm_genes=True)
    assert merged_stats["genes_with_excalibr_calibrations"] == 1


def test_format_calibration_summary():
    text = format_calibration_summary({"genes_with_excalibr_calibrations": 4, "genes_with_evidence_assigned": 2})
    assert "Genes with ExCALIBR calibrations: 4" in text
    assert "2 (50.0%)" in text


def test_reclassification_flags_agree_disagree_and_no_evidence():
    clinvar_group = pd.Series(["Pathogenic", "Benign", "Pathogenic", "Likely benign"])
    points = pd.Series([3, -2, -1, 0])

    flags, in_scope = reclassification_flags(clinvar_group, points)

    assert list(in_scope) == [True, True, True, True]
    assert list(flags[AGREE_LABEL]) == [True, True, False, False]
    assert list(flags[DISAGREE_LABEL]) == [False, False, True, False]
    assert list(flags[NO_EVIDENCE_LABEL]) == [False, False, False, True]


def test_reclassification_flags_excludes_out_of_scope_clinvar_calls():
    clinvar_group = pd.Series(["Uncertain significance", "Pathogenic"])
    points = pd.Series([3, 3])

    _flags, in_scope = reclassification_flags(clinvar_group, points)

    assert list(in_scope) == [False, True]


def test_compute_reclassification_agreement():
    controls_df = pd.DataFrame(
        {
            "clnsig_group_18_25": ["Pathogenic", "Benign", "Pathogenic", "Likely benign", "Benign/Likely benign"],
            "ExC_points_2025": [3, -2, -1, 0, -4],
            "OP_points": [2, -1, 0, 1, -2],
        }
    )

    results = compute_reclassification_agreement(controls_df)

    excalibr_total, excalibr_determinate, excalibr_pct, excalibr_table = results["ExCALIBR evidence"]
    assert excalibr_total == 5
    assert excalibr_determinate == 4
    assert excalibr_pct == pytest.approx(75.0)
    assert excalibr_table.loc[AGREE_LABEL, "count"] == 3

    functional_total, functional_determinate, functional_pct, functional_table = results["Functional class"]
    assert functional_total == 5
    assert functional_determinate == 4
    assert functional_pct == pytest.approx(75.0)
    assert functional_table.loc[AGREE_LABEL, "count"] == 3


def test_format_reclassification_table_reports_totals_and_agreement():
    controls_df = pd.DataFrame(
        {
            "clnsig_group_18_25": ["Pathogenic", "Benign"],
            "ExC_points_2025": [3, -2],
            "OP_points": [3, -2],
        }
    )
    total, determinate, pct, table = compute_reclassification_agreement(controls_df)["ExCALIBR evidence"]

    text = format_reclassification_table("ExCALIBR evidence -- test", total, determinate, pct, table)

    assert "Total control variants: 2" in text
    assert "Determinate calls (evidence assigned): 2" in text
    assert "Agreement with ClinVar PLP/BLB (of determinate calls): 100.0%" in text


def test_build_reclassification_report_covers_every_controls_prefixed_sheet(tmp_path):
    controls_path = tmp_path / "controls.xlsx"
    _write_controls_file(
        controls_path,
        {
            "controls_REVEL_GeneSpecific": [("Pathogenic", 1, 1), ("Benign", -1, -1)],
            "controls_AM_GeneSpecific": [("Pathogenic", 2, 2)],
            "not_a_controls_sheet": [("Pathogenic", 1, 1)],
        },
    )

    sections = build_reclassification_report(pd.ExcelFile(controls_path))

    joined = "\n\n".join(sections)
    assert "ExCALIBR evidence -- controls_REVEL_GeneSpecific" in joined
    assert "Functional class -- controls_REVEL_GeneSpecific" in joined
    assert "ExCALIBR evidence -- controls_AM_GeneSpecific" in joined
    assert "not_a_controls_sheet" not in joined
    assert len(sections) == 4  # 2 controls_ sheets x 2 points columns


def test_distinct_dna_variants_dedups_by_genomic_coordinates():
    df = pd.DataFrame(
        {
            "Gene": ["GENEA", "GENEA", "GENEB"],
            "Chrom": [1, 1, 2],
            "hg38_start": [100, 100, 200],
            "ref_allele": ["A", "A", "C"],
            "alt_allele": ["G", "G", "T"],
            "Class_REVEL": ["Pathogenic", "Benign", "Uncertain"],
        }
    )

    result = distinct_dna_variants(df)

    # The two GENEA rows share coordinates -- only the first (Pathogenic) survives.
    assert len(result) == 2
    assert result[result["Gene"] == "GENEA"]["Class_REVEL"].tolist() == ["Pathogenic"]


def test_compute_variant_classification_stats(tmp_path):
    controls_path = tmp_path / "controls.xlsx"
    # coords (1, 100, A, G) is shared between "controls" and "gnomAD", to check
    # that the combined total counts it once, not twice.
    _write_variant_classification_sheets(
        controls_path,
        {
            "controls_REVEL_GeneSpecific": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic"),
                ("GENEX", 1, 101, "A", "G", "Benign"),
            ],
            "ClinGen_Repo_REVEL_GeneSpecific": [
                ("GENEX", 1, 102, "A", "G", "Likely Pathogenic"),
            ],
            "VUS_REVEL": [
                ("GENEX", 1, 103, "A", "G", "Pathogenic"),
                ("GENEX", 1, 104, "A", "G", "Uncertain"),
                ("GENEX", 1, 105, "A", "G", "Benign"),
            ],
            "gnomAD_REVEL": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic"),
            ],
            "Unobserved_REVEL": [
                ("GENEX", 1, 106, "A", "G", "Pathogenic"),
                ("GENEX", 1, 107, "A", "G", "Uncertain"),
            ],
        },
        mode="w",
    )

    stats = compute_variant_classification_stats(pd.ExcelFile(controls_path))

    assert stats == {
        "total_classified": 8,
        "total_pathogenic_or_benign": 6,
        "vus_total": 3,
        "vus_resolved": 2,
        "unobserved_total": 2,
        "unobserved_resolved": 1,
    }


def test_format_variant_classification_summary():
    stats = {
        "total_classified": 8,
        "total_pathogenic_or_benign": 6,
        "vus_total": 3,
        "vus_resolved": 2,
        "unobserved_total": 2,
        "unobserved_resolved": 1,
    }

    text = format_variant_classification_summary(stats)

    assert text.startswith(VARIANT_CLASSIFICATION_TITLE)
    assert "Distinct DNA variants classified: 8" in text
    assert "Pathogenic or benign: 6 of 8 (75.0%)" in text
    assert "ClinVar VUS resolved (reclassified pathogenic or benign): 2 of 3 (66.7%)" in text
    assert "Unobserved variants resolved (classified pathogenic or benign): 1 of 2 (50.0%)" in text


def test_format_variant_classification_summary_uses_given_title():
    stats = {
        "total_classified": 1,
        "total_pathogenic_or_benign": 1,
        "vus_total": 0,
        "vus_resolved": 0,
        "unobserved_total": 0,
        "unobserved_resolved": 0,
    }

    text = format_variant_classification_summary(stats, title=RECLASSIFICATION_VARIANT_CLASSIFICATION_TITLE)

    assert text.startswith(RECLASSIFICATION_VARIANT_CLASSIFICATION_TITLE)


def test_format_variant_classification_summary_handles_zero_totals():
    stats = {
        "total_classified": 0,
        "total_pathogenic_or_benign": 0,
        "vus_total": 0,
        "vus_resolved": 0,
        "unobserved_total": 0,
        "unobserved_resolved": 0,
    }

    text = format_variant_classification_summary(stats)

    assert "Distinct DNA variants classified: 0" in text
    assert "(nan%)" in text


def test_points_are_pathogenic_or_benign():
    points = pd.Series([12, 6, 5, 0, -1, -6, -7, -12])

    result = points_are_pathogenic_or_benign(points)

    # >=6 (Pathogenic/Likely Pathogenic) or <=-1 (Likely Benign/Benign) is True;
    # 0-5 (Uncertain) is False.
    assert list(result) == [True, True, False, False, True, True, True, True]


def test_compute_variant_classification_stats_from_reclassification_file(tmp_path):
    reclassification_path = tmp_path / "reclassification.tsv"
    _write_reclassification_file(
        reclassification_path,
        [
            ("Uncertain significance", None, "A", "G", 8),  # VUS, resolved
            ("Uncertain significance", None, "A", "G", 2),  # VUS, not resolved
            (None, None, "A", "G", -8),  # Unobserved SNV, resolved
            (None, None, "A", "G", 1),  # Unobserved SNV, not resolved
            ("Pathogenic", None, "A", "G", 11),  # controls-like, not VUS/Unobserved
            (None, 0.01, "A", "G", -3),  # gnomAD-like, not Unobserved (gnomad_MAF set)
            (None, None, "AC", "G", 8),  # not a SNV, excluded from Unobserved
        ],
    )
    df = pd.read_csv(reclassification_path, sep="\t", usecols=RECLASSIFICATION_USECOLS)

    stats = compute_variant_classification_stats_from_reclassification_file(df)

    assert stats == {
        "total_classified": 7,
        "total_pathogenic_or_benign": 5,
        "vus_total": 2,
        "vus_resolved": 1,
        "unobserved_total": 2,
        "unobserved_resolved": 1,
    }

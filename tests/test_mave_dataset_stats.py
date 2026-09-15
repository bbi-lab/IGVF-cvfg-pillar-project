import re
import unicodedata

import pandas as pd
import pytest
from click.testing import CliRunner

from src.mave_dataset_stats import (
    AGREE_LABEL,
    BENIGN_VALUES,
    CALM_MERGED_LABEL,
    CLINGEN_EVIDENCE_REPOSITORY_TITLE,
    CLINVAR_CONFLICT_LABEL,
    COMBINED_EVIDENCE_LABEL_BY_PREDICTOR,
    COMBINED_REVEL_EVIDENCE_LABEL,
    COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR,
    COMPOSITE_SCORE_DATASETS_TITLE,
    CONCORDANT_LABEL,
    CONFLICTING_EVIDENCE_VALUE,
    CONTROL_BLB_LABEL,
    CONTROL_CONCORDANCE_EVIDENCE_LABEL,
    CONTROL_PLP_LABEL,
    CONTROL_VUS_LABEL,
    DISAGREE_LABEL,
    DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL,
    DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL,
    DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL,
    DISCORDANT_LABEL,
    DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL,
    FUNCTIONAL_POINTS_COL,
    GENE_DISCORDANCE_TITLE,
    GNOMAD_LABEL,
    IGVF_DATASET_MEASUREMENT_COUNTS_TITLE,
    MISSENSE_CONSEQUENCE_VALUE,
    MISSENSE_CONTROL_CONCORDANCE_SOURCES,
    MISSENSE_CONTROL_CONCORDANCE_TITLE,
    MUTPRED2_TRAINING_STEP_LABEL,
    NO_ANNOTATION_LABEL,
    NO_EVIDENCE_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    PATHOGENIC_VALUES,
    REVEL_TRAINING_STEP_LABEL,
    SIMPLIFIED_CONSEQUENCE_COL,
    SNV_ACCESSIBLE_LABEL,
    SNV_LABEL,
    UNIVERSAL_CALIBRATION_CATEGORY_SHEETS_BY_PREDICTOR,
    UNIVERSAL_CALIBRATION_CLASS_COL_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_CHI_SQUARED_COMPARISONS,
    VARIANT_CLASSIFICATION_CHI_SQUARED_TITLE,
    VARIANT_CLASSIFICATION_CLASS_COL_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_CONFLICTING_COL_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_POINTS_COL_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_PREDICTOR_POINTS_COL_BY_PREDICTOR,
    VARIANT_CLASSIFICATION_PREDICTORS,
    VARIANT_CLASSIFICATION_TITLE,
    VUS_LABEL,
    build_reclassification_report,
    clinvar_classification_from_flags,
    clinvar_significance_flags,
    compute_all_stats,
    compute_clingen_evidence_repository_stats,
    compute_composite_score_datasets,
    compute_control_concordance,
    compute_excalibr_calibration_stats,
    compute_gene_discordance_stats,
    compute_igvf_dataset_measurement_counts,
    compute_reclassification_agreement,
    compute_reclassification_filter_funnel,
    compute_variant_classification_chi_squared_tests,
    compute_variant_classification_stats,
    control_concordance_flags,
    distinct_dna_variants,
    distinct_variant_flags,
    excalibr_dataset_to_gene_map,
    format_calibration_summary,
    format_clinical_table,
    format_clingen_evidence_repository_summary,
    format_composite_score_datasets,
    format_control_concordance_report,
    format_gene_discordance_summary,
    format_genomic_variant_count,
    format_igvf_dataset_measurement_counts,
    format_reclassification_filter_funnel,
    format_reclassification_table,
    format_variant_classification_chi_squared_tests,
    format_variant_classification_table,
    funnel_distinct_dna_variants,
    has_any_value,
    is_snv_accessible,
    load_dataset_metadata,
    main,
    matches_any_value,
    mixed_year_clinvar_series,
    reclassification_flags,
    split_genes,
    stats_to_dataframe,
    summarize_clinical_flags,
    summarize_flags,
    variant_flags,
)

BASE_COLUMNS = ["Dataset", "Gene", "hgvs_c", "hgvs_p", "nucleotide_or_aa"]
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
# Genomic-coordinate columns compute_reclassification_filter_funnel reads from
# the expanded file (irrelevant to condensed-file-only tests, but harmless
# there -- both share this schema for simplicity, see _write_full_variant_file).
GENOMIC_COLUMNS = [
    "mavedb_variant_urn",
    "Chrom",
    "hg38_start",
    "ref_allele",
    "alt_allele",
    "aa_pos",
    "aa_ref",
    "aa_alt",
]
FULL_COLUMNS = BASE_COLUMNS + ANNOTATION_COLUMNS + GENOMIC_COLUMNS


def _write_condensed(path, rows, rna_scores=None):
    """`rna_scores`, if given, is a list of `rna_score` values aligned with
    `rows` (default: empty for every row, i.e. no RNA scores)."""
    df = pd.DataFrame(rows, columns=BASE_COLUMNS)
    df["rna_score"] = rna_scores if rna_scores is not None else ""
    df.to_csv(path, sep="\t", index=False)


def _write_full_variant_file(path, rows):
    df = pd.DataFrame(rows, columns=FULL_COLUMNS)
    df["rna_score"] = ""
    df.to_csv(path, sep="\t", index=False)


def _write_checkpoint_file(path, rows):
    """`rows` is a list of (Dataset, Gene, mavedb_variant_urn, Chrom, hg38_start,
    ref_allele, alt_allele, hgvs_p, auth_reported_score, Flag, VariantNotes,
    splice_var_amino, splice_measure, revel_train_amino, mp2_train_amino,
    Assertion_ClinGen_repo, Updated_Classification_ClinGen_repo) tuples --
    the columns compute_reclassification_filter_funnel and
    compute_clingen_evidence_repository_stats read from
    Variant_Classification_analysis.ipynb's checkpoint file. `splice_measure
    == "Yes"` (a splice-aware dataset) exempts a row from the
    `splice_var_amino` exclusion -- see SPLICE_MEASURE_COL/SPLICE_MEASURE_VALUE.
    """
    columns = [
        "Dataset",
        "Gene",
        "mavedb_variant_urn",
        "Chrom",
        "hg38_start",
        "ref_allele",
        "alt_allele",
        "hgvs_p",
        "auth_reported_score",
        "Flag",
        "VariantNotes",
        "splice_var_amino",
        "splice_measure",
        "revel_train_amino",
        "mp2_train_amino",
        "Assertion_ClinGen_repo",
        "Updated_Classification_ClinGen_repo",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)


def _write_chek2_file(path, rows=()):
    """`rows` is a list of (hgvs_pro, score, Filter_CI) tuples -- empty by
    default, since most tests don't need any CHEK2 QC matches.
    """
    columns = ["hgvs_pro", "score", "Filter_CI"]
    pd.DataFrame(rows, columns=columns).to_excel(path, index=False)


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


def _write_variant_classification_sheets_by_predictor(path, category_rows, mode="w"):
    """Writes the same `category_rows` to all three predictors' own sheets
    (REVEL/AlphaMissense/MutPred2) -- `category_rows` is {category: rows}
    where `category` is one of `VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_
    PREDICTOR`'s five keys ("controls"/"ClinGen_Repo"/"VUS"/"gnomAD"/
    "Unobserved") and each row is (Gene, Chrom, hg38_start, ref_allele,
    alt_allele, class_value, points_value, fxn_points, predictor_points) --
    written under each predictor's own sheet name and `Class_*`/
    `Total_Points_*`/predictor-points column names, plus the shared
    `Fxn_points` column, the columns `compute_variant_classification_stats`
    reads. Reusing the same rows for every predictor keeps expected numbers
    identical across REVEL/AlphaMissense/MutPred2 -- sufficient to exercise
    the per-predictor sheet-reading mechanics without three independent
    synthetic datasets. Appends (`mode="a"`) to an existing workbook (e.g.
    one already written by `_write_controls_file`); pass `mode="w"` (the
    default) to create a fresh file.

    Also adds `clnsig_group_18_25`/`ExC_points_2025`/`OP_points`/
    `Updated_Classification_ClinGen_repo` placeholder columns: the real
    Supplementary_Data_5 `controls_*` sheets carry both this function's
    columns and those `build_reclassification_report` reads, and every
    `controls_*` sheet (now three, one per predictor) gets scanned by that
    function too -- without these placeholders it would KeyError on them.

    The predictor's own `Conflicting_*` column is derived from `fxn_points`/
    `predictor_points`, mirroring `Variant_Classification_analysis.ipynb`'s
    `split_zero` (see docs/variant_classification.md#functional-vs-predictor-
    conflict-flagged-not-excluded): `"Conflicting evidence"` when the two
    sides are strictly opposite in sign, otherwise the row's own
    `points_value`.
    """
    with pd.ExcelWriter(path, mode=mode, engine="openpyxl") as writer:
        for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
            category_sheets = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]
            class_col = VARIANT_CLASSIFICATION_CLASS_COL_BY_PREDICTOR[predictor]
            points_col = VARIANT_CLASSIFICATION_POINTS_COL_BY_PREDICTOR[predictor]
            predictor_points_col = VARIANT_CLASSIFICATION_PREDICTOR_POINTS_COL_BY_PREDICTOR[predictor]
            conflicting_col = VARIANT_CLASSIFICATION_CONFLICTING_COL_BY_PREDICTOR[predictor]
            columns = [
                "Gene",
                "Chrom",
                "hg38_start",
                "ref_allele",
                "alt_allele",
                class_col,
                points_col,
                FUNCTIONAL_POINTS_COL,
                predictor_points_col,
            ]
            for category, rows in category_rows.items():
                df = pd.DataFrame(rows, columns=columns)
                df["clnsig_group_18_25"] = None
                df["ExC_points_2025"] = 0
                df["OP_points"] = 0
                df["Updated_Classification_ClinGen_repo"] = None
                opposite_signs = (df[FUNCTIONAL_POINTS_COL] * df[predictor_points_col]) < 0
                df[conflicting_col] = df[points_col].where(~opposite_signs, CONFLICTING_EVIDENCE_VALUE)
                df.to_excel(writer, sheet_name=category_sheets[category], index=False)


def _write_universal_calibration_sheets_by_predictor(path, mode="w"):
    """Minimal Supplementary Data 6-style workbook -- one ClinVar and one
    ClinGen row per predictor's own `controls_*_OP`/`ClinGen_*_OP` sheet
    (`UNIVERSAL_CALIBRATION_CATEGORY_SHEETS_BY_PREDICTOR`), concordant on
    every row (`Class_OP_<predictor>` agrees with the control label).
    Structurally valid for `compute_control_concordance`'s universal-
    calibration evidence source -- the CLI tests using this only check the
    row's presence, not its exact counts.
    """
    with pd.ExcelWriter(path, mode=mode, engine="openpyxl") as writer:
        for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
            sheets = UNIVERSAL_CALIBRATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]
            class_col = UNIVERSAL_CALIBRATION_CLASS_COL_BY_PREDICTOR[predictor]
            pd.DataFrame(
                [
                    {
                        "Gene": "GENEX",
                        "clnsig_group_18_25": "Pathogenic",
                        "Updated_Classification_ClinGen_repo": None,
                        class_col: "Pathogenic",
                        "simplified_consequence": "missense_variant",
                    }
                ]
            ).to_excel(writer, sheet_name=sheets["controls"], index=False)
            pd.DataFrame(
                [
                    {
                        "Gene": "GENEX",
                        "clnsig_group_18_25": None,
                        "Updated_Classification_ClinGen_repo": "Pathogenic",
                        class_col: "Pathogenic",
                        "simplified_consequence": "missense_variant",
                    }
                ]
            ).to_excel(writer, sheet_name=sheets["ClinGen_Repo"], index=False)


def _write_gene_discordance_sheet(path, rows, mode="w"):
    """`rows` is a list of (Gene, Chrom, hg38_start, ref_allele, alt_allele,
    clnsig_group_18_25, Class_REVEL) tuples -- the columns
    `compute_gene_discordance_stats` reads from `controls_REVEL_GeneSpecific`.
    """
    columns = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele", "clnsig_group_18_25", "Class_REVEL"]
    with pd.ExcelWriter(path, mode=mode, engine="openpyxl") as writer:
        pd.DataFrame(rows, columns=columns).to_excel(writer, sheet_name="controls_REVEL_GeneSpecific", index=False)


@pytest.fixture
def dataset_files(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_IGVF_A", "GENEA", "c1", "p1", "nt"),
            ("DS_IGVF_A", "GENEA", "c1", "p1", "nt"),  # duplicate variant, second measurement
            ("DS_IGVF_B", "GENEB, GENEC", "c2", "p2", "aa"),
            ("DS_COMM_A", "GENEC", "c3", "p3", "nt"),
            ("DS_COMM_B", "GENED", "c4", "p4", "aa"),
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
    # Trailing 8 fields on every row are GENOMIC_COLUMNS (mavedb_variant_urn,
    # Chrom, hg38_start, ref_allele, alt_allele, aa_pos, aa_ref, aa_alt) --
    # arbitrary but distinct/non-special-cased (no LDLR/F9/TP53/SFPQ genes
    # here), just enough for compute_reclassification_filter_funnel to run.
    _write_full_variant_file(
        condensed_path,
        [
            (
                "DS_IGVF_A",
                "BRCA1",
                "c1",
                "p1",
                "nt",
                "0.5",
                "0.5",
                "0.5",
                "Pathogenic",
                "Uncertain significance",
                "",
                "A",
                "G",
                "urn:mavedb:1a",
                "17",
                "100",
                "A",
                "G",
                "1",
                "I",
                "M",
            ),
            (
                # Same variant, second measurement -- distinct condensed row/urn,
                # like a real duplicate MaveDB submission (urn is 1:1 with
                # condensed rows in the real data, never repeated within it).
                "DS_IGVF_A",
                "BRCA1",
                "c1",
                "p1",
                "nt",
                "0.5",
                "0.5",
                "0.5",
                "Pathogenic",
                "Uncertain significance",
                "",
                "A",
                "G",
                "urn:mavedb:1b",
                "17",
                "100",
                "A",
                "G",
                "1",
                "I",
                "M",
            ),
            (
                "DS_IGVF_B",
                "GENEB, GENEC",
                "c2",
                "p2",
                "aa",
                "",
                "",
                "",
                "Uncertain significance",
                "Uncertain significance",
                "0.001",
                "A",
                "T",
                "urn:mavedb:2",
                "17",
                "200",
                "A",
                "T",
                "2",
                "I",
                "I",
            ),
            (
                "DS_COMM_A",
                "GENEC",
                "c3",
                "p3",
                "nt",
                "",
                "",
                "",
                "",
                "Pathogenic",
                "",
                "A",
                "G",
                "urn:mavedb:3",
                "17",
                "300",
                "A",
                "G",
                "3",
                "I",
                "I",
            ),
            (
                "DS_COMM_B",
                "GENED",
                "c4a|c4b",
                "p4",
                "aa",
                "0.2|",
                "0.2|",
                "0.9|0.9",
                "Benign|",
                "Benign|",
                "",
                "AC|A",
                "GT|G",
                "urn:mavedb:4",
                "17",
                "400",
                "AC",
                "GT",
                "4",
                "I",
                "I",
            ),
        ],
    )
    _write_full_variant_file(
        expanded_path,
        [
            (
                "DS_IGVF_A",
                "BRCA1",
                "c1",
                "p1",
                "nt",
                "0.5",
                "0.5",
                "0.5",
                "Pathogenic",
                "Uncertain significance",
                "",
                "A",
                "G",
                "urn:mavedb:1a",
                "17",
                "100",
                "A",
                "G",
                "1",
                "I",
                "M",
            ),
            (
                "DS_IGVF_B",
                "GENEB, GENEC",
                "c2",
                "p2",
                "aa",
                "",
                "",
                "",
                "Uncertain significance",
                "Uncertain significance",
                "0.001",
                "A",
                "T",
                "urn:mavedb:2",
                "17",
                "200",
                "A",
                "T",
                "2",
                "I",
                "I",
            ),
            (
                "DS_COMM_A",
                "GENEC",
                "c3",
                "p3",
                "nt",
                "",
                "",
                "",
                "",
                "Pathogenic",
                "",
                "A",
                "G",
                "urn:mavedb:3",
                "17",
                "300",
                "A",
                "G",
                "3",
                "I",
                "I",
            ),
            (
                "DS_COMM_B",
                "GENED",
                "c4a",
                "p4",
                "aa",
                "0.2",
                "0.2",
                "0.9",
                "Benign",
                "Benign",
                "",
                "AC",
                "GT",
                "urn:mavedb:4",
                "17",
                "400",
                "AC",
                "GT",
                "4",
                "I",
                "I",
            ),
            (
                "DS_COMM_B",
                "GENED",
                "c4b",
                "p4",
                "aa",
                "",
                "",
                "0.9",
                "",
                "",
                "",
                "AC",
                "GT",
                "urn:mavedb:4",
                "17",
                "401",
                "A",
                "G",
                "4",
                "I",
                "I",
            ),
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
    # Same rows written to all three predictors' sheets (REVEL/AlphaMissense/
    # MutPred2), so each predictor's row in the resulting table shows the same
    # numbers -- see _write_variant_classification_sheets_by_predictor.
    # Every row's Fxn_points/predictor-points are "don't care" (fxn = total,
    # predictor = 0) except row 105, which is single-source (fxn = 0,
    # predictor = -1) to exercise the -1-point single-source/conflicting
    # split -- see _write_variant_classification_sheets_by_predictor.
    _write_variant_classification_sheets_by_predictor(
        controls_path,
        {
            "controls": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic", 6, 6, 0),
                ("GENEX", 1, 101, "A", "G", "Benign", -2, -2, 0),
            ],
            "ClinGen_Repo": [
                ("GENEX", 1, 102, "A", "G", "Likely Pathogenic", 6, 6, 0),
            ],
            "VUS": [
                ("GENEX", 1, 103, "A", "G", "Pathogenic", 6, 6, 0),  # resolved
                ("GENEX", 1, 104, "A", "G", "Uncertain", 2, 2, 0),  # not resolved
                # resolved, at the -1-point threshold, single source (predictor only)
                ("GENEX", 1, 105, "A", "G", "Benign", -1, 0, -1),
            ],
            "gnomAD": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic", 6, 6, 0),  # duplicate of controls row above
            ],
            "Unobserved": [
                ("GENEX", 1, 106, "A", "G", "Pathogenic", 6, 6, 0),  # resolved
                ("GENEX", 1, 107, "A", "G", "Uncertain", 3, 3, 0),  # not resolved
            ],
        },
        mode="a",
    )

    # Mirrors the expanded_path rows above (same mavedb_variant_urn/genomic
    # coordinates) at the "post pre-checkpoint-filter" stage -- nothing here
    # is flagged/tagged for exclusion, so all 4 distinct DNA variants survive
    # compute_reclassification_filter_funnel's post-checkpoint steps too.
    # p1 is a ClinGen Evidence Repository control that retains a determinate
    # classification after evidence removal; p2 is one that doesn't (drops
    # to Uncertain) -- exercising compute_clingen_evidence_repository_stats'
    # pre-/post-removal split. p3/p4 carry no ClinGen data.
    checkpoint_path = tmp_path / "checkpoint.csv"
    _write_checkpoint_file(
        checkpoint_path,
        [
            (
                "DS_IGVF_A", "BRCA1", "urn:mavedb:1a", "17", 100, "A", "G", "p1", 0.5, None, None, "No", "No", "No", "No",
                "Pathogenic", "Pathogenic",
            ),
            (
                "DS_IGVF_B", "GENEB, GENEC", "urn:mavedb:2", "17", 200, "A", "T", "p2", None, None, None, "No", "No", "No", "No",
                "Benign", "Uncertain",
            ),
            (
                "DS_COMM_A", "GENEC", "urn:mavedb:3", "17", 300, "A", "G", "p3", None, None, None, "No", "No", "No", "No",
                None, None,
            ),
            (
                "DS_COMM_B", "GENED", "urn:mavedb:4", "17", 400, "AC", "GT", "p4", 0.2, None, None, "No", "No", "No", "No",
                None, None,
            ),
            (
                "DS_COMM_B", "GENED", "urn:mavedb:4", "17", 401, "A", "G", "p4", None, None, None, "No", "No", "No", "No",
                None, None,
            ),
        ],
    )
    chek2_path = tmp_path / "chek2.xlsx"
    _write_chek2_file(chek2_path)

    universal_controls_path = tmp_path / "universal_controls.xlsx"
    _write_universal_calibration_sheets_by_predictor(universal_controls_path)

    return (
        condensed_path,
        metadata_path,
        expanded_path,
        excalibr_path,
        controls_path,
        checkpoint_path,
        chek2_path,
        universal_controls_path,
    )


def test_split_genes_handles_multi_gene_datasets():
    assert split_genes("CALM1, CALM2, CALM3") == ["CALM1", "CALM2", "CALM3"]
    assert split_genes("BRCA1") == ["BRCA1"]


def test_compute_all_stats_buckets(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats, gene_breakdown = compute_all_stats(condensed_path, metadata_path)

    igvf = stats["IGVF"]
    assert igvf["datasets"] == 2
    assert igvf["variant_effect_measurements"] == 3
    assert igvf["rna_scores"] == 0
    assert igvf["composite_scores"] == 0
    # p2 (DS_IGVF_B) is the only "aa"-resolution row; c1 (DS_IGVF_A, x2) is the only "nt"-resolution variant.
    assert igvf["distinct_protein_variants_assayed"] == 1
    assert igvf["distinct_dna_variants_assayed"] == 1
    assert igvf["distinct_variants_assayed"] == 2
    assert igvf["genes_represented"] == 3

    non_igvf = stats["Community (non-IGVF)"]
    assert non_igvf["datasets"] == 2
    assert non_igvf["variant_effect_measurements"] == 1
    assert non_igvf["rna_scores"] == 0
    assert non_igvf["composite_scores"] == 1
    # p4 (DS_COMM_B) is "aa"-resolution; c3 (DS_COMM_A) is "nt"-resolution.
    assert non_igvf["distinct_protein_variants_assayed"] == 1
    assert non_igvf["distinct_dna_variants_assayed"] == 1
    assert non_igvf["distinct_variants_assayed"] == 2
    assert non_igvf["genes_represented"] == 2
    assert non_igvf["genes_not_in_igvf_data"] == 1  # GENED only; GENEC is shared with IGVF

    combined = stats["Combined (IGVF + community)"]
    assert combined["datasets"] == 4
    assert combined["variant_effect_measurements"] == 4
    assert combined["rna_scores"] == 0
    assert combined["composite_scores"] == 1
    # p2 and p4 are "aa"-resolution; c1 and c3 are "nt"-resolution.
    assert combined["distinct_protein_variants_assayed"] == 2
    assert combined["distinct_dna_variants_assayed"] == 2
    assert combined["distinct_variants_assayed"] == 4
    assert combined["genes_represented"] == 4

    # GENEA (DS_IGVF_A) and GENEB (DS_IGVF_B) are IGVF-only; GENED (DS_COMM_B) is
    # community-only; GENEC is shared (DS_IGVF_B and DS_COMM_A).
    assert gene_breakdown["IGVF only"] == ["GENEA", "GENEB"]
    assert gene_breakdown["Community (non-IGVF) only"] == ["GENED"]
    assert gene_breakdown["Both IGVF and community (non-IGVF)"] == ["GENEC"]


def test_compute_igvf_dataset_measurement_counts(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_IGVF_BIG", "GENEA", "c1", "p1", "nt"),
            ("DS_IGVF_BIG", "GENEA", "c2", "p2", "nt"),
            ("DS_IGVF_BIG", "GENEA", "c3", "p3", "nt"),
            ("DS_IGVF_SMALL", "GENEB", "c4", "p4", "nt"),
            ("DS_IGVF_META", "GENEC", "c5", "p5", "nt"),
            ("DS_COMM", "GENED", "c6", "p6", "nt"),  # not IGVF-produced, excluded
        ],
        rna_scores=["-0.1", "", "0.2", "", "", ""],
    )
    _write_metadata(
        metadata_path,
        [
            ("DS_IGVF_BIG", "Yes", "primary score set"),
            ("DS_IGVF_SMALL", "Yes", "primary score set"),
            ("DS_IGVF_META", "Yes", "meta-analysis"),
            ("DS_COMM", "No", "primary score set"),
        ],
    )

    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)

    table = compute_igvf_dataset_measurement_counts(condensed, metadata)

    # Sorted by variant_effect_measurements descending; DS_IGVF_META has 0
    # (it's a meta-analysis, so its row counts as a composite score instead).
    assert table["Dataset"].tolist() == ["DS_IGVF_BIG", "DS_IGVF_SMALL", "DS_IGVF_META"]
    assert table["variant_effect_measurements"].tolist() == [3, 1, 0]
    assert table["rna_scores"].tolist() == [2, 0, 0]
    assert table["composite_scores"].tolist() == [0, 0, 1]

    text = format_igvf_dataset_measurement_counts(table)
    assert text.startswith(IGVF_DATASET_MEASUREMENT_COUNTS_TITLE)
    assert "DS_IGVF_BIG" in text
    assert "DS_COMM" not in text


def test_compute_composite_score_datasets_uses_dataset_files_fixture(dataset_files):
    condensed_path, metadata_path = dataset_files
    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)

    table = compute_composite_score_datasets(condensed, metadata)

    # DS_COMM_B (GENED, community) is the only meta-analysis dataset in this fixture.
    assert table["Dataset"].tolist() == ["DS_COMM_B"]
    assert table["Gene"].tolist() == ["GENED"]
    assert table["IGVF / Community"].tolist() == ["Community"]
    assert table["Scores"].tolist() == [1]


def test_compute_composite_score_datasets_sorts_by_scores_descending(tmp_path):
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_META_SMALL", "GENEA", "c1", "p1", "nt"),
            ("DS_META_BIG", "GENEB, GENEC", "c2", "p2", "aa"),
            ("DS_META_BIG", "GENEB, GENEC", "c3", "p3", "aa"),
            ("DS_META_BIG", "GENEB, GENEC", "c4", "p4", "aa"),
            ("DS_PRIMARY", "GENED", "c5", "p5", "nt"),  # not a meta-analysis, excluded
        ],
    )
    _write_metadata(
        metadata_path,
        [
            ("DS_META_SMALL", "Yes", "meta-analysis"),
            ("DS_META_BIG", "No", "meta-analysis"),
            ("DS_PRIMARY", "Yes", "primary score set"),
        ],
    )

    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)

    table = compute_composite_score_datasets(condensed, metadata)

    assert table["Dataset"].tolist() == ["DS_META_BIG", "DS_META_SMALL"]
    assert table["Gene"].tolist() == ["GENEB, GENEC", "GENEA"]
    assert table["IGVF / Community"].tolist() == ["Community", "IGVF"]
    assert table["Scores"].tolist() == [3, 1]


def test_format_composite_score_datasets():
    table = pd.DataFrame(
        {
            "Gene": ["GENEB, GENEC", "GENEA"],
            "Dataset": ["DS_META_BIG", "DS_META_SMALL"],
            "IGVF / Community": ["Community", "IGVF"],
            "Scores": [3, 1],
        }
    )

    text = format_composite_score_datasets(table)

    assert text.startswith(COMPOSITE_SCORE_DATASETS_TITLE)
    assert "Total composite scores: 4" in text
    assert "DS_META_BIG" in text
    assert "DS_META_SMALL" in text


def test_format_composite_score_datasets_handles_empty_table():
    table = pd.DataFrame(columns=["Gene", "Dataset", "IGVF / Community", "Scores"])

    text = format_composite_score_datasets(table)

    assert text.startswith(COMPOSITE_SCORE_DATASETS_TITLE)
    assert "Total composite scores: 0" in text


def test_format_genomic_variant_count_counts_parsed_rows_not_lines():
    # A field containing an embedded newline (properly quoted, as pandas writes
    # it) would inflate a naive `wc -l` count relative to the true row count --
    # this must report the pandas-parsed row count (2), not the line count.
    df = pd.DataFrame({"hgvs_c": ["c.1A>G", "c.2A>G"], "note": ["single line", "line one\nline two"]})

    text = format_genomic_variant_count("some/path.tsv.gz", df)

    assert text == "Genomic variants (rows in some/path.tsv.gz): 2"


def test_stats_to_dataframe_adds_pct_of_combined_measurements(dataset_files):
    condensed_path, metadata_path = dataset_files
    stats, _ = compute_all_stats(condensed_path, metadata_path)
    table = stats_to_dataframe(stats)

    # 3 of the combined 4 variant_effect_measurements are IGVF, 1 is non-IGVF.
    assert table.loc["IGVF", "pct_variant_effect_measurements"] == pytest.approx(75.0)
    assert table.loc["Community (non-IGVF)", "pct_variant_effect_measurements"] == pytest.approx(25.0)
    assert table.loc["Combined (IGVF + community)", "pct_variant_effect_measurements"] == pytest.approx(100.0)


def test_rna_scores_is_a_breakdown_of_measurements_not_an_addition(tmp_path):
    """A row with an rna_score also has a regular auth_reported_score in the
    real data, so rna_scores is reported beside variant_effect_measurements
    as a breakdown, not counted as extra measurements on top of it."""
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_IGVF_A", "GENEA", "c1", "p1", "nt"),
            ("DS_IGVF_A", "GENEA", "c2", "p2", "nt"),
            ("DS_COMM_A", "GENEB", "c3", "p3", "nt"),
        ],
        rna_scores=["-0.5", "", ""],
    )
    _write_metadata(
        metadata_path,
        [
            ("DS_IGVF_A", "Yes", "primary score set"),
            ("DS_COMM_A", "No", "primary score set"),
        ],
    )

    stats, _ = compute_all_stats(condensed_path, metadata_path)
    table = stats_to_dataframe(stats)

    assert stats["IGVF"]["variant_effect_measurements"] == 2
    assert stats["IGVF"]["rna_scores"] == 1
    assert stats["Community (non-IGVF)"]["rna_scores"] == 0
    assert stats["Combined (IGVF + community)"]["variant_effect_measurements"] == 3
    assert stats["Combined (IGVF + community)"]["rna_scores"] == 1

    # Combined total for the RNA-inclusive percentage is 3 measurements + 1 RNA
    # score = 4; IGVF contributes 2 measurements + 1 RNA score = 3 of that.
    assert table.loc["IGVF", "pct_variant_effect_measurements_with_rna_scores"] == pytest.approx(75.0)
    assert table.loc["Community (non-IGVF)", "pct_variant_effect_measurements_with_rna_scores"] == pytest.approx(25.0)
    assert table.loc[
        "Combined (IGVF + community)", "pct_variant_effect_measurements_with_rna_scores"
    ] == pytest.approx(100.0)


def test_dataset_summary_always_merges_calm_paralogs(tmp_path):
    """CALM1/CALM2/CALM3 count as one gene in the Dataset summary's
    genes_represented/genes_not_in_igvf_data, matching the always-merged
    Genes represented list -- there's no --merge-calm-genes toggle for this
    section."""
    condensed_path = tmp_path / "condensed.tsv"
    metadata_path = tmp_path / "metadata.xlsx"

    _write_condensed(
        condensed_path,
        [
            ("DS_COMM_A", "CALM1", "c1", "p1", "nt"),
            ("DS_COMM_B", "CALM2, CALM3", "c2", "p2", "nt"),
            ("DS_COMM_C", "GENED", "c3", "p3", "nt"),
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

    stats, gene_breakdown = compute_all_stats(condensed_path, metadata_path)
    assert stats["Community (non-IGVF)"]["genes_represented"] == 2
    assert gene_breakdown["Community (non-IGVF) only"] == [CALM_MERGED_LABEL, "GENED"]


def test_compute_all_stats_raises_on_missing_metadata(dataset_files):
    condensed_path, metadata_path = dataset_files
    _write_condensed(
        condensed_path,
        [("DS_UNKNOWN", "GENEX", "c5", "p5", "nt")],
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
    (
        condensed_path,
        metadata_path,
        expanded_path,
        excalibr_path,
        controls_path,
        checkpoint_path,
        chek2_path,
        universal_controls_path,
    ) = full_dataset_files
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
            "--universal-controls-file",
            str(universal_controls_path),
            "--checkpoint-file",
            str(checkpoint_path),
            "--chek2-file",
            str(chek2_path),
            "--output",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert "IGVF" in result.output
    assert "Score coverage" in result.output
    assert "Clinical attributes" in result.output
    assert "Filtering effects on the reclassification dataset" in result.output
    assert (
        "Control concordance (ClinVar vs. ClinGen; OddsPath alone vs. combined with "
        "REVEL/AlphaMissense/MutPred2)" in result.output
    )

    # checkpoint_path's p1 (Pathogenic -> Pathogenic) and p2 (Benign -> Uncertain)
    # are the only two rows with an original determinate ClinGen classification;
    # only p1 retains one after evidence removal -- see full_dataset_files.
    clingen_section = result.output.split(CLINGEN_EVIDENCE_REPOSITORY_TITLE)[1].split(
        "=== Reclassification agreement"
    )[0]
    assert "2 distinct DNA variants across 2 genes" in clingen_section
    assert "Retained a determinate classification: 1 of 2 (50.0%) distinct DNA variants across 1 genes" in (
        clingen_section
    )
    assert "REVEL: 1 variants across 1 genes" in clingen_section
    assert "AlphaMissense: 1 variants across 1 genes" in clingen_section
    assert "MutPred2: 1 variants across 1 genes" in clingen_section

    # full_dataset_files' expanded_path has 5 rows (p1, p2, p3, and p4's c4a/c4b).
    assert f"Genomic variants (rows in {expanded_path}): 5" in result.output

    # BRCA1 (DS_IGVF_A) and GENEB (DS_IGVF_B) are IGVF-only; GENED (DS_COMM_B) is
    # community-only; GENEC is shared (DS_IGVF_B and DS_COMM_A).
    genes_section = result.output.split("=== Genes represented ===")[1].split(COMPOSITE_SCORE_DATASETS_TITLE)[0]
    assert "IGVF only (2): BRCA1, GENEB" in genes_section
    assert "Community (non-IGVF) only (1): GENED" in genes_section
    assert "Both IGVF and community (non-IGVF) (1): GENEC" in genes_section

    # DS_COMM_B is the fixture's only meta-analysis dataset (gene GENED, community,
    # 1 condensed row/composite score).
    composite_score_section = result.output.split(COMPOSITE_SCORE_DATASETS_TITLE)[1].split(
        "=== Score coverage"
    )[0]
    assert "Total composite scores: 1" in composite_score_section
    assert "GENED" in composite_score_section
    assert "DS_COMM_B" in composite_score_section
    assert "Community" in composite_score_section

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

    # 8 distinct DNA variants across the 5 category sheets combined (9 rows,
    # minus the (GENEX, 1, 100, A, G) coordinate shared by "controls" and
    # "gnomAD"); 6 of those 8 are Pathogenic/Likely Pathogenic/Benign/Likely
    # Benign (the two "Uncertain" VUS/Unobserved rows aren't). Same rows were
    # written to all three predictors' sheets, so REVEL/AlphaMissense/MutPred2
    # each show these same numbers in their own table row.
    variant_classification_section = result.output.split(VARIANT_CLASSIFICATION_TITLE)[1].split(
        VARIANT_CLASSIFICATION_CHI_SQUARED_TITLE
    )[0]
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert predictor in variant_classification_section
    assert "Distinct DNA variants classified" in variant_classification_section
    assert "ClinVar VUS unresolved:" in variant_classification_section
    assert "gnomAD variants resolved (classified pathogenic or benign):" in variant_classification_section
    assert "gnomAD variants unresolved:" in variant_classification_section
    assert "Unobserved variants unresolved:" in variant_classification_section
    assert variant_classification_section.count("6 of 8 (75.0%)") == 3
    assert variant_classification_section.count("2 of 3 (66.7%)") == 3
    # VUS: row 103 (Pathogenic, of 3 VUS) resolves pathogenic; row 105 (Benign, at the
    # -1-point threshold) resolves benign; row 104 (Uncertain, fxn 2/predictor 0) is the sole
    # unresolved VUS, experimental-only. gnomAD: its sole row (same coords as controls row 100,
    # Pathogenic, fxn 6/predictor 0) resolves pathogenic, experimental-only, leaving gnomAD with no
    # unresolved rows at all. Unobserved: row 106 (Pathogenic, of 2 Unobserved) resolves pathogenic;
    # no Unobserved row resolves benign; row 107 (Uncertain, fxn 3/predictor 0) is the sole
    # unresolved Unobserved, experimental-only. Counts below were verified against the actual
    # rendered table rather than hand-derived, given how many columns now overlap.
    assert variant_classification_section.count("1 of 2 (50.0%)") == 9
    assert variant_classification_section.count("1 of 3 (33.3%)") == 12
    assert variant_classification_section.count("0 of 2 (0.0%)") == 6
    assert variant_classification_section.count("1 of 1 (100.0%)") == 33
    assert variant_classification_section.count("0 of 1 (0.0%)") == 66
    assert variant_classification_section.count("0 of 0 (nan%)") == 51

    # Chi-squared section: gnomAD (1/1 PLP, 0/1 BLB), VUS (1/3 PLP, 1/3 BLB), and
    # Unobserved (1/2 PLP) from the same rows above -- every comparison's 2x2 table
    # is small enough that Yates' continuity correction zeroes it out (chi2=0, p=1)
    # regardless of predictor, since all three predictors' sheets carry the same rows.
    chi_squared_section = result.output.split(VARIANT_CLASSIFICATION_CHI_SQUARED_TITLE)[1].split(
        GENE_DISCORDANCE_TITLE
    )[0]
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert f"-- {predictor} --" in chi_squared_section
    assert "Pathogenic/Likely Pathogenic rate: gnomAD vs. ClinVar VUS:" in chi_squared_section
    assert "Benign/Likely Benign rate: gnomAD vs. ClinVar VUS:" in chi_squared_section
    assert "Pathogenic/Likely Pathogenic rate: Unobserved vs. ClinVar VUS:" in chi_squared_section
    assert "Pathogenic/Likely Pathogenic rate: Unobserved vs. gnomAD:" in chi_squared_section
    assert chi_squared_section.count("chi2 = 0.0000, df = 1, p = 1") == 4 * len(VARIANT_CLASSIFICATION_PREDICTORS)
    assert chi_squared_section.count("gnomAD: 1 of 1 (100.0%)") == 2 * len(VARIANT_CLASSIFICATION_PREDICTORS)
    assert chi_squared_section.count("ClinVar VUS: 1 of 3 (33.3%)") == 3 * len(VARIANT_CLASSIFICATION_PREDICTORS)
    assert chi_squared_section.count("Unobserved: 1 of 2 (50.0%)") == 2 * len(VARIANT_CLASSIFICATION_PREDICTORS)

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
    # condensed rows get distinct urns (real urns are 1:1 with condensed
    # rows); expanded's rows can share a urn, so they just reuse the first.
    rows = [
        (
            "DS_A",
            "GENEA",
            "c1",
            "p1",
            "nt",
            "",
            "",
            "",
            "Pathogenic",
            "Pathogenic",
            "",
            "A",
            "G",
            "urn:mavedb:1a",
            "17",
            "100",
            "A",
            "G",
            "1",
            "I",
            "M",
        ),
        (
            "DS_A",
            "GENEA",
            "c1",
            "p1",
            "nt",
            "",
            "",
            "",
            "Benign",
            "Benign",
            "",
            "A",
            "G",
            "urn:mavedb:1b",
            "17",
            "100",
            "A",
            "G",
            "1",
            "I",
            "M",
        ),
    ]
    _write_full_variant_file(condensed_path, rows)
    _write_full_variant_file(expanded_path, rows)
    _write_metadata(metadata_path, [("DS_A", "No", "primary score set")], genes={"DS_A": "GENEA"})

    excalibr_path = tmp_path / "excalibr_calibrations.xlsx"
    _write_excalibr_calibrations(excalibr_path, [("DS_A", None, None)])
    controls_path = tmp_path / "controls.xlsx"
    _write_controls_file(controls_path, {"controls_TEST": [("Pathogenic", 1, 1)]})
    _write_variant_classification_sheets_by_predictor(
        controls_path,
        {
            "controls": [("GENEA", 1, 1, "A", "G", "Pathogenic", 6, 6, 0)],
            "ClinGen_Repo": [],
            "VUS": [],
            "gnomAD": [],
            "Unobserved": [],
        },
        mode="a",
    )
    checkpoint_path = tmp_path / "checkpoint.csv"
    _write_checkpoint_file(
        checkpoint_path,
        [("DS_A", "GENEA", "urn:mavedb:1a", "17", 100, "A", "G", "p1", None, None, None, "No", "No", "No", "No", None, None)],
    )
    chek2_path = tmp_path / "chek2.xlsx"
    _write_chek2_file(chek2_path)
    universal_controls_path = tmp_path / "universal_controls.xlsx"
    _write_universal_calibration_sheets_by_predictor(universal_controls_path)
    extra_args = [
        "--excalibr-calibrations-file",
        str(excalibr_path),
        "--controls-file",
        str(controls_path),
        "--universal-controls-file",
        str(universal_controls_path),
        "--checkpoint-file",
        str(checkpoint_path),
        "--chek2-file",
        str(chek2_path),
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
    (
        condensed_path,
        metadata_path,
        expanded_path,
        _excalibr_path,
        _controls_path,
        _checkpoint_path,
        _chek2_path,
        _universal_controls_path,
    ) = full_dataset_files
    _write_full_variant_file(
        condensed_path,
        [
            (
                "DS_UNKNOWN",
                "GENEX",
                "c5",
                "p5",
                "nt",
                "",
                "",
                "",
                "",
                "",
                "",
                "A",
                "G",
                "urn:mavedb:5",
                "17",
                "500",
                "A",
                "G",
                "5",
                "I",
                "I",
            )
        ],
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
    # None of these genes are in EXCALIBR_EXCLUDED_GENES, so the excl. counts match.
    assert stats["genes_with_excalibr_calibrations_excl"] == 4
    assert stats["genes_with_evidence_assigned_excl"] == 2


def test_compute_excalibr_calibration_stats_excludes_f9_tp53_sfpq(tmp_path):
    metadata_path = tmp_path / "metadata.xlsx"
    _write_gene_metadata(
        metadata_path,
        {"DS_F9": "F9", "DS_TP53": "TP53", "DS_SFPQ": "SFPQ", "DS_A": "GENEA"},
    )
    metadata = load_dataset_metadata(metadata_path)

    calibrations = pd.DataFrame(
        {
            "dataset": ["DS_F9", "DS_TP53", "DS_SFPQ", "DS_A"],
            "range_-1": ["-1 -0.5", None, "-1 -0.5", "-1 -0.5"],
            "range_1": [None, "0.5 1", None, None],
        }
    )

    stats = compute_excalibr_calibration_stats(calibrations, metadata)

    # All 4 genes are calibrated and have evidence, but only GENEA survives exclusion.
    assert stats["genes_with_excalibr_calibrations"] == 4
    assert stats["genes_with_evidence_assigned"] == 4
    assert stats["genes_with_excalibr_calibrations_excl"] == 1
    assert stats["genes_with_evidence_assigned_excl"] == 1


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
    text = format_calibration_summary(
        {
            "genes_with_excalibr_calibrations": 4,
            "genes_with_evidence_assigned": 2,
            "genes_with_excalibr_calibrations_excl": 3,
            "genes_with_evidence_assigned_excl": 1,
        }
    )
    assert "Genes with ExCALIBR calibrations: 4 (3 excluding F9/TP53/SFPQ)" in text
    assert "2 (50.0%)" in text
    assert "1 (33.3%) excluding F9/TP53/SFPQ" in text


def test_funnel_distinct_dna_variants_normalizes_hg38_start_and_chrom():
    # The notebook's checkpoint file writes hg38_start/Chrom out from a
    # float-typed column (e.g. "100.0"), while the expanded file keeps them
    # as plain strings ("100") -- both must collapse to the same DNA variant.
    df = pd.DataFrame(
        {
            "Gene": ["GENEA", "GENEA"],
            "Chrom": ["17", "17.0"],
            "hg38_start": ["100", "100.0"],
            "ref_allele": ["A", "A"],
            "alt_allele": ["G", "G"],
        }
    )
    assert funnel_distinct_dna_variants(df) == 1


def _expanded_funnel_fixture():
    """One row per pre-checkpoint exclusion reason, plus survivors.

    row1/row2 share a DNA coordinate (two datasets scoring the same
    variant): the leading rows count DNA-level rows directly (5) while
    distinct_dna_variants collapses the shared coordinate (4). row3 is an
    ordinary gene, unaffected by anything at this pre-checkpoint stage.
    row4 is an F9 dataset not on the meta-analysis allowlist -- the only
    exclusion left at this stage. LDLR's LA-module exclusion no longer
    happens here: it's now a `Flag == '*'` set upstream by
    flag_variants.py, so it behaves like any other pre-existing flag,
    folded into "Other flagged variants" post-checkpoint (see
    `_checkpoint_funnel_fixture`). row5 (GENEA) survives.

    row1 also carries an `rna_score`, the fixture's sole RNA-score-bearing
    measurement -- exercises the funnel's leading "including RNA scores"/
    "- RNA scores" rows (see `compute_reclassification_filter_funnel`).
    """
    columns = [
        "Dataset",
        "Gene",
        "mavedb_variant_urn",
        "Chrom",
        "hg38_start",
        "ref_allele",
        "alt_allele",
        "rna_score",
    ]
    rows = [
        ("GENEB_Study_2020a", "GENEB", "urn:mavedb:v1", "19", 100, "A", "G", "0.5"),
        ("GENEB_Study_2020b", "GENEB", "urn:mavedb:v2", "19", 100, "A", "G", ""),
        ("GENEC_Study_2020", "GENEC", "urn:mavedb:v3", "19", 200, "C", "T", ""),
        ("F9_Popp_2025_strep_2", "F9", "urn:mavedb:v4", "X", 300, "A", "C", ""),
        ("GENEA_Study_2020", "GENEA", "urn:mavedb:v5", "1", 400, "A", "G", ""),
    ]
    return pd.DataFrame(rows, columns=columns)


def _checkpoint_funnel_fixture():
    """One row per post-checkpoint exclusion reason, plus one survivor, each
    with its own distinct DNA coordinate/urn so every step's delta is 1.

    DS_VN/DS_SPLICE_VN/DS_START_LOST each carry exactly one of the three
    `VariantNotes` tags (`conflicting_fxn_data`/`splice_variant_not_measured`/
    `start_lost_variant_not_measured`) -- each removed at its own funnel step,
    distinct from DS_SPLICE, which instead carries `splice_var_amino ==
    "Yes"` (the separate, later "at the amino-acid level" step). None of
    these datasets is splice-aware (`splice_measure == "No"` throughout), so
    DS_SPLICE is still excluded here -- see
    `test_compute_reclassification_filter_funnel_keeps_splice_aware_amino_variant`
    for the exemption itself.

    DS_REVEL (revel_train_amino="Yes", mp2_train_amino="No") and DS_SURVIVOR
    (revel_train_amino="No", mp2_train_amino="Yes") both reach the "Other
    flagged variants" branch point, then diverge: the REVEL branch excludes
    DS_REVEL and keeps DS_SURVIVOR, while the MutPred2 branch excludes
    DS_SURVIVOR and keeps DS_REVEL instead -- proving the two branches are
    independent siblings rather than one applied on top of the other.
    """
    columns = [
        "Dataset",
        "Gene",
        "mavedb_variant_urn",
        "Chrom",
        "hg38_start",
        "ref_allele",
        "alt_allele",
        "hgvs_p",
        "auth_reported_score",
        "Flag",
        "VariantNotes",
        "splice_var_amino",
        "splice_measure",
        "revel_train_amino",
        "mp2_train_amino",
    ]
    rows = [
        ("DS_SFPQ", "SFPQ", "urn:mavedb:c1", "1", 10, "A", "G", "NP_1.1:p.Ser1Ala", 0.1, None, None, "No", "No", "No", "No"),
        (
            "DS_CHEK2",
            "CHEK2",
            "urn:mavedb:c2",
            "22",
            20,
            "A",
            "G",
            "NP_2.1:p.Ser2Ala",
            0.5,
            None,
            None,
            "No",
            "No",
            "No",
            "No",
        ),
        (
            "DS_VN",
            "GENEB",
            "urn:mavedb:c3",
            "2",
            30,
            "A",
            "G",
            "NP_3.1:p.Ser3Ala",
            0.3,
            None,
            "conflicting_fxn_data",
            "No",
            "No",
            "No",
            "No",
        ),
        (
            "DS_SPLICE_VN",
            "GENEG",
            "urn:mavedb:c8",
            "7",
            80,
            "A",
            "G",
            "NP_8.1:p.Ser8Ala",
            0.8,
            None,
            "splice_variant_not_measured",
            "No",
            "No",
            "No",
            "No",
        ),
        (
            "DS_START_LOST",
            "GENEH",
            "urn:mavedb:c9",
            "8",
            90,
            "A",
            "G",
            "NP_9.1:p.Ser9Ala",
            0.9,
            None,
            "start_lost_variant_not_measured",
            "No",
            "No",
            "No",
            "No",
        ),
        (
            "DS_SPLICE",
            "GENEC",
            "urn:mavedb:c4",
            "3",
            40,
            "A",
            "G",
            "NP_4.1:p.Ser4Ala",
            0.4,
            None,
            None,
            "Yes",
            "No",
            "No",
            "No",
        ),
        (
            "DS_FLAG",
            "GENED",
            "urn:mavedb:c5",
            "4",
            50,
            "A",
            "G",
            "NP_5.1:p.Ser5Ala",
            0.5,
            "*",
            None,
            "No",
            "No",
            "No",
            "No",
        ),
        (
            "DS_REVEL",
            "GENEE",
            "urn:mavedb:c6",
            "5",
            60,
            "A",
            "G",
            "NP_6.1:p.Ser6Ala",
            0.6,
            None,
            None,
            "No",
            "No",
            "Yes",
            "No",
        ),
        (
            "DS_SURVIVOR",
            "GENEF",
            "urn:mavedb:c7",
            "6",
            70,
            "A",
            "G",
            "NP_7.1:p.Ser7Ala",
            0.7,
            None,
            None,
            "No",
            "No",
            "No",
            "Yes",
        ),
    ]
    return pd.DataFrame(rows, columns=columns)


def _chek2_funnel_fixture():
    # Matches DS_CHEK2's row above (hgvs_p prefix stripped: "p.Ser2Ala", score 0.5).
    return pd.DataFrame({"hgvs_pro": ["p.Ser2Ala"], "score": [0.5], "Filter_CI": [1]})


def _condensed_funnel_fixture():
    """One row per urn used across _expanded_funnel_fixture/_checkpoint_funnel_fixture
    (real urns are 1:1 with condensed rows), each with its own (hgvs_c, hgvs_p)
    pair -- distinct_assayed_variants maps a row's urn back to this pair, so
    with one urn per pair here the counts are the same as counting urns directly.
    """
    urns = [f"urn:mavedb:v{i}" for i in range(1, 6)] + [f"urn:mavedb:c{i}" for i in range(1, 10)]
    return pd.DataFrame(
        {
            "mavedb_variant_urn": urns,
            "hgvs_c": [f"c.{i}A>G" for i in range(len(urns))],
            "hgvs_p": [f"p.Ser{i}Ala" for i in range(len(urns))],
        }
    )


def test_compute_reclassification_filter_funnel_sequential_steps():
    expanded = _expanded_funnel_fixture()
    checkpoint = _checkpoint_funnel_fixture()
    chek2 = _chek2_funnel_fixture()
    condensed = _condensed_funnel_fixture()

    steps = compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed)
    by_label = {step["label"]: step for step in steps}

    # steps[0]: "including RNA scores" -- row1's rna_score adds 1 measurement
    # on top of the usual one-per-urn count (5), on top of the same
    # rows/dna/assayed counts as every other leading-row assertion below.
    raw = steps[0]
    assert raw["rows"] == 5
    # Every row in this fixture has its own distinct urn, so variant_effect_measurements
    # (grouped by urn) equals rows here -- see
    # test_compute_reclassification_filter_funnel_counts_measurements_by_urn_grouping
    # for a case where they differ.
    assert raw["variant_effect_measurements"] == 6
    assert raw["distinct_dna_variants"] == 4  # row1/row2 share a coordinate
    assert raw["distinct_assayed_variants"] == 5
    assert raw["rows_removed"] == 0
    assert raw["measurements_removed"] == 0
    assert raw["dna_removed"] == 0
    assert raw["assayed_removed"] == 0

    # steps[1]: "- RNA scores" -- drops back to the plain one-per-urn count (5),
    # with the same rows/dna/assayed as steps[0] (RNA scores are a column, not
    # extra rows).
    rna_scores_step = by_label["- RNA scores"]
    assert steps[1] is rna_scores_step
    assert rna_scores_step["rows"] == 5
    assert rna_scores_step["rows_removed"] == 0
    assert rna_scores_step["variant_effect_measurements"] == 5
    assert rna_scores_step["measurements_removed"] == 1
    assert rna_scores_step["distinct_dna_variants"] == 4
    assert rna_scores_step["dna_removed"] == 0
    assert rna_scores_step["distinct_assayed_variants"] == 5
    assert rna_scores_step["assayed_removed"] == 0

    f9tp53 = by_label["- F9/TP53: restricted to the meta-analysis dataset only"]
    assert f9tp53["rows"] == 4
    assert f9tp53["distinct_dna_variants"] == 3
    assert f9tp53["distinct_assayed_variants"] == 4

    # Post-checkpoint steps read `checkpoint` directly (its own 9-row fixture,
    # unrelated to the 2 pre-checkpoint survivors above).
    sfpq = by_label["- SFPQ excluded entirely (insufficient ClinVar controls)"]
    assert sfpq["rows"] == 8
    assert sfpq["distinct_dna_variants"] == 8
    assert sfpq["distinct_assayed_variants"] == 8

    # DS_VN's VariantNotes ("conflicting_fxn_data") is removed here; DS_SPLICE_VN and
    # DS_START_LOST's tags aren't this tag, so they survive this step specifically.
    conflicting = by_label["- Conflicting functional data (opposite-sign Fxn_points across assays)"]
    assert conflicting["rows"] == 7
    assert conflicting["distinct_dna_variants"] == 7
    assert conflicting["distinct_assayed_variants"] == 7

    splice_not_measured = by_label["- Splice variant not measured (no functional measurement at all)"]
    assert splice_not_measured["rows"] == 6  # DS_SPLICE_VN removed
    assert splice_not_measured["distinct_dna_variants"] == 6
    assert splice_not_measured["distinct_assayed_variants"] == 6

    start_lost = by_label["- Start-lost variant not measured (no functional measurement at all)"]
    assert start_lost["rows"] == 5  # DS_START_LOST removed
    assert start_lost["distinct_dna_variants"] == 5
    assert start_lost["distinct_assayed_variants"] == 5

    # DS_SPLICE (splice_var_amino == "Yes", a different column from VariantNotes) is
    # removed here -- distinct from the VariantNotes-tag-based steps above.
    splice = by_label[
        "- Splice variant not measured at the amino-acid level (aa-level candidates sharing a "
        "splice-affecting group, unless the assay detects splicing)"
    ]
    assert splice["rows"] == 4
    assert splice["distinct_dna_variants"] == 4
    assert splice["distinct_assayed_variants"] == 4

    chek2_step = by_label["- CHEK2 QC flag (Filter_CI == 1)"]
    assert chek2_step["rows"] == 3
    assert chek2_step["distinct_dna_variants"] == 3
    assert chek2_step["distinct_assayed_variants"] == 3

    other_flag = by_label["- Other flagged variants (pre-existing Flag == '*')"]
    assert other_flag["rows"] == 2
    assert other_flag["distinct_dna_variants"] == 2
    assert other_flag["distinct_assayed_variants"] == 2

    # REVEL and MutPred2 training exclusion are siblings baselined against
    # "Other flagged variants" (rows == 2) independently, not a chain: each
    # drops exactly one of {DS_REVEL, DS_SURVIVOR} and keeps the other, so
    # both branches end at 1 row but with rows_removed == 1 (2 - 1), not 0
    # (which is what a same-list-position default baseline would wrongly
    # compute if MutPred2 were chained after REVEL instead of after "Other
    # flagged variants").
    revel_train = by_label[REVEL_TRAINING_STEP_LABEL]
    assert revel_train["rows"] == 1
    assert revel_train["rows_removed"] == 1
    assert revel_train["variant_effect_measurements"] == 1
    assert revel_train["measurements_removed"] == 1
    assert revel_train["distinct_dna_variants"] == 1
    assert revel_train["distinct_assayed_variants"] == 1

    mp2_train = by_label[MUTPRED2_TRAINING_STEP_LABEL]
    assert mp2_train["rows"] == 1
    assert mp2_train["rows_removed"] == 1
    assert mp2_train["variant_effect_measurements"] == 1
    assert mp2_train["measurements_removed"] == 1
    assert mp2_train["distinct_dna_variants"] == 1
    assert mp2_train["distinct_assayed_variants"] == 1

    assert steps[-2] is revel_train
    assert steps[-1] is mp2_train


def test_compute_reclassification_filter_funnel_counts_measurements_by_urn_grouping():
    """`variant_effect_measurements` groups by `mavedb_variant_urn`, unlike
    `rows`, which counts DNA-level rows directly -- one measurement (urn) can
    have multiple DNA-level candidate rows (e.g. a protein-resolution
    measurement reverse-translated to several DNA candidates), so the two
    counts differ whenever that happens.
    """
    columns = [
        "Dataset",
        "Gene",
        "mavedb_variant_urn",
        "Chrom",
        "hg38_start",
        "ref_allele",
        "alt_allele",
        "aa_pos",
        "aa_ref",
        "aa_alt",
        "rna_score",
    ]
    expanded = pd.DataFrame(
        [
            # Two DNA-level candidates for the same measurement (urn:mavedb:1).
            ("DS_A", "GENEA", "urn:mavedb:1", "1", 100, "A", "G", 1, "M", "I", ""),
            ("DS_A", "GENEA", "urn:mavedb:1", "1", 100, "A", "T", 1, "M", "I", ""),
            ("DS_B", "GENEB", "urn:mavedb:2", "2", 200, "C", "G", 2, "S", "T", ""),
        ],
        columns=columns,
    )
    checkpoint = pd.DataFrame(
        columns=columns
        + [
            "hgvs_p",
            "auth_reported_score",
            "Flag",
            "VariantNotes",
            "splice_var_amino",
            "splice_measure",
            "revel_train_amino",
            "mp2_train_amino",
        ]
    )
    chek2 = pd.DataFrame(columns=["hgvs_pro", "score", "Filter_CI"])
    condensed = pd.DataFrame(
        {
            "mavedb_variant_urn": ["urn:mavedb:1", "urn:mavedb:2"],
            "hgvs_c": ["c.1A>G", "c.2C>G"],
            "hgvs_p": ["p.Met1Ile", "p.Ser2Thr"],
        }
    )

    steps = compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed)

    # No rna_score values here, so the leading "including RNA scores" row
    # (steps[0]) has no bonus over the plain one-per-urn count.
    raw = steps[0]
    assert raw["rows"] == 3
    assert raw["variant_effect_measurements"] == 2  # urn:mavedb:1's two rows collapse to one measurement


def test_compute_reclassification_filter_funnel_keeps_splice_aware_amino_variant():
    """A `splice_var_amino == 'Yes'` row from a splice-aware dataset
    (`splice_measure == 'Yes'`) is exempt from the amino-acid-level splice
    exclusion and survives, unlike an otherwise-identical row from a
    non-splice-aware dataset -- matching `src.build_variant_
    reclassification_dataset.apply_notebook_exclusions`'s own exemption.
    """
    expanded = pd.DataFrame(
        columns=["Dataset", "Gene", "mavedb_variant_urn", "Chrom", "hg38_start", "ref_allele", "alt_allele", "aa_pos", "aa_ref", "aa_alt", "rna_score"]
    )
    checkpoint = pd.DataFrame(
        [
            (
                "DS_NOT_SPLICE_AWARE", "GENEA", "urn:mavedb:1", "1", 10, "A", "G", "NP_1.1:p.Ser1Ala", 0.1,
                None, None, "Yes", "No", "No", "No",
            ),
            (
                "DS_SPLICE_AWARE", "GENEB", "urn:mavedb:2", "2", 20, "A", "G", "NP_2.1:p.Ser2Ala", 0.2,
                None, None, "Yes", "Yes", "No", "No",
            ),
        ],
        columns=[
            "Dataset", "Gene", "mavedb_variant_urn", "Chrom", "hg38_start", "ref_allele", "alt_allele", "hgvs_p",
            "auth_reported_score", "Flag", "VariantNotes", "splice_var_amino", "splice_measure",
            "revel_train_amino", "mp2_train_amino",
        ],
    )
    chek2 = pd.DataFrame(columns=["hgvs_pro", "score", "Filter_CI"])
    condensed = pd.DataFrame(
        {
            "mavedb_variant_urn": ["urn:mavedb:1", "urn:mavedb:2"],
            "hgvs_c": ["c.1A>G", "c.2A>G"],
            "hgvs_p": ["p.Ser1Ala", "p.Ser2Ala"],
        }
    )

    steps = compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed)
    by_label = {step["label"]: step for step in steps}

    splice = by_label[
        "- Splice variant not measured at the amino-acid level (aa-level candidates sharing a "
        "splice-affecting group, unless the assay detects splicing)"
    ]
    assert splice["rows"] == 1  # DS_NOT_SPLICE_AWARE removed, DS_SPLICE_AWARE kept
    assert splice["rows_removed"] == 1


def test_format_reclassification_filter_funnel_summary_lines():
    expanded = _expanded_funnel_fixture()
    checkpoint = _checkpoint_funnel_fixture()
    chek2 = _chek2_funnel_fixture()
    condensed = _condensed_funnel_fixture()
    steps = compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed)

    text = format_reclassification_filter_funnel(steps)
    assert "=== Filtering effects on the reclassification dataset ===" in text
    assert "Distinct DNA variants reclassified: 1 of 4 (25.0%)" in text
    assert "Distinct assayed variants reclassified: 1 of 5 (20.0%)" in text
    # rows/variant_effect_measurements each get their own *_removed column, same as the
    # pre-existing distinct_dna_variants/distinct_assayed_variants columns.
    assert "rows_removed" in text
    assert "variant_effect_measurements" in text
    assert "measurements_removed" in text
    # The table columns appear in this order: rows immediately followed by its own
    # removed column, then variant_effect_measurements immediately followed by its own.
    assert text.index("rows") < text.index("rows_removed") < text.index("variant_effect_measurements")
    assert text.index("variant_effect_measurements") < text.index("measurements_removed")


def test_compute_clingen_evidence_repository_stats(tmp_path):
    # s1: SFPQ, excluded entirely. s2: pre-existing Flag == '*', excluded.
    # s3: conflicting_fxn_data VariantNotes tag, excluded. s4: original
    # classification not determinate (Uncertain Significance), excluded
    # from the pre-removal population. s5/s6/s7 survive every exclusion and
    # have a determinate original classification: s5 stays Pathogenic after
    # removal (retained), s6 drops to Uncertain (not retained), s7 flips
    # direction (Likely Pathogenic -> Benign, still determinate, retained).
    # s8: splice_var_amino == "Yes" but from a splice-aware dataset
    # (splice_measure == "Yes"), so it's exempt from that exclusion and
    # survives like s5/s7 (stays Pathogenic, retained) -- proves the
    # exemption reaches this section, not just the filter-funnel one.
    # s9: original classification Uncertain Significance (excluded from
    # pre-removal, like s4) but the *recalculated* classification is
    # determinate (Likely Pathogenic) -- matches the real population gate
    # `Variant_Classification_analysis.ipynb` uses (not conditioned on the
    # original classification), so it's retained post-removal and counted
    # in post_removal_originally_vus.
    checkpoint_path = tmp_path / "checkpoint.csv"
    _write_checkpoint_file(
        checkpoint_path,
        [
            ("DS1", "SFPQ", "urn:s1", "1", 10, "A", "G", "p1", None, None, None, "No", "No", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEA", "urn:s2", "1", 20, "A", "G", "p2", None, "*", None, "No", "No", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEA", "urn:s3", "1", 30, "A", "G", "p3", None, None, "conflicting_fxn_data", "No", "No", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEA", "urn:s4", "1", 40, "A", "G", "p4", None, None, None, "No", "No", "No", "No", "Uncertain Significance", "Uncertain Significance"),
            ("DS1", "GENEA", "urn:s5", "1", 50, "A", "G", "p5", None, None, None, "No", "No", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEB", "urn:s6", "1", 60, "A", "G", "p6", None, None, None, "No", "No", "No", "No", "Likely Benign", "Uncertain"),
            ("DS1", "GENEB", "urn:s7", "1", 70, "A", "G", "p7", None, None, None, "No", "No", "No", "No", "Likely Pathogenic", "Benign"),
            ("DS1", "GENEC", "urn:s8", "1", 80, "A", "G", "p8", None, None, None, "Yes", "Yes", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEA", "urn:s9", "1", 90, "A", "G", "p9", None, None, None, "No", "No", "No", "No", "Uncertain Significance", "Likely Pathogenic"),
        ],
    )
    chek2_path = tmp_path / "chek2.xlsx"
    _write_chek2_file(chek2_path)

    controls_path = tmp_path / "controls.xlsx"
    clingen_rows = {
        "REVEL": [("GENEA", "Pathogenic"), ("GENEB", "Uncertain")],
        "AlphaMissense": [("GENEA", "Pathogenic"), ("GENEB", "Benign"), ("GENEC", "Likely Benign")],
        "MutPred2": [("GENEA", "Likely Pathogenic")],
    }
    with pd.ExcelWriter(controls_path) as writer:
        for predictor, rows in clingen_rows.items():
            sheet = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]["ClinGen_Repo"]
            pd.DataFrame(rows, columns=["Gene", "Updated_Classification_ClinGen_repo"]).to_excel(
                writer, sheet_name=sheet, index=False
            )

    checkpoint = pd.read_csv(checkpoint_path, low_memory=False)
    chek2 = pd.read_excel(chek2_path, header=0)
    controls_workbook = pd.ExcelFile(controls_path)

    stats = compute_clingen_evidence_repository_stats(checkpoint, chek2, controls_workbook)

    assert stats["pre_removal_total"] == 4
    assert stats["pre_removal_genes"] == 3
    assert stats["pre_removal_pathogenic"] == 2
    assert stats["pre_removal_likely_pathogenic"] == 1
    assert stats["pre_removal_benign"] == 0
    assert stats["pre_removal_likely_benign"] == 1

    assert stats["post_removal_total"] == 4
    assert stats["post_removal_genes"] == 3
    assert stats["post_removal_plp"] == 3
    assert stats["post_removal_blb"] == 1
    assert stats["post_removal_originally_vus"] == 1

    assert stats["REVEL_total"] == 2
    assert stats["REVEL_genes"] == 2
    assert stats["REVEL_plp"] == 1
    assert stats["REVEL_blb"] == 0

    assert stats["AlphaMissense_total"] == 3
    assert stats["AlphaMissense_genes"] == 3
    assert stats["AlphaMissense_plp"] == 1
    assert stats["AlphaMissense_blb"] == 2

    assert stats["MutPred2_total"] == 1
    assert stats["MutPred2_genes"] == 1
    assert stats["MutPred2_plp"] == 1
    assert stats["MutPred2_blb"] == 0


def test_format_clingen_evidence_repository_summary(tmp_path):
    # s10: original classification Uncertain Significance (so it's excluded
    # from the pre-removal/812-style count) but recalculates to a
    # determinate Benign -- exercises the new "Originally Uncertain
    # Significance in ClinGen" line.
    checkpoint_path = tmp_path / "checkpoint.csv"
    _write_checkpoint_file(
        checkpoint_path,
        [
            ("DS1", "GENEA", "urn:s5", "1", 50, "A", "G", "p5", None, None, None, "No", "No", "No", "No", "Pathogenic", "Pathogenic"),
            ("DS1", "GENEB", "urn:s6", "1", 60, "A", "G", "p6", None, None, None, "No", "No", "No", "No", "Likely Benign", "Uncertain"),
            ("DS1", "GENEC", "urn:s10", "1", 70, "A", "G", "p10", None, None, None, "No", "No", "No", "No", "Uncertain Significance", "Benign"),
        ],
    )
    chek2_path = tmp_path / "chek2.xlsx"
    _write_chek2_file(chek2_path)
    controls_path = tmp_path / "controls.xlsx"
    with pd.ExcelWriter(controls_path) as writer:
        for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
            sheet = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]["ClinGen_Repo"]
            pd.DataFrame([("GENEA", "Pathogenic")], columns=["Gene", "Updated_Classification_ClinGen_repo"]).to_excel(
                writer, sheet_name=sheet, index=False
            )
    checkpoint = pd.read_csv(checkpoint_path, low_memory=False)
    chek2 = pd.read_excel(chek2_path, header=0)
    controls_workbook = pd.ExcelFile(controls_path)
    stats = compute_clingen_evidence_repository_stats(checkpoint, chek2, controls_workbook)

    text = format_clingen_evidence_repository_summary(stats)

    assert "=== ClinGen Evidence Repository control set (evidence removed & reclassified) ===" in text
    assert "2 distinct DNA variants across 2 genes" in text
    assert "Pathogenic: 1" in text
    assert "Likely Benign: 1" in text
    assert "Retained a determinate classification: 2 of 2 (100.0%) distinct DNA variants across 2 genes" in text
    assert "Pathogenic or Likely Pathogenic: 1" in text
    assert "Benign or Likely Benign: 1" in text
    assert "Originally Uncertain Significance in ClinGen: 1 of 2 (50.0%)" in text
    assert (
        "Post filtering (excluding each predictor's own training variants and deduplicating "
        "amino-acid-resolution variants to a representative DNA variant), retained for analysis:"
    ) in text
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert f"{predictor}: 1 variants across 1 genes (Pathogenic or Likely Pathogenic: 1, Benign or Likely Benign: 0)" in text


def test_control_concordance_flags_concordant_discordant_vus():
    control_group = pd.Series(["Pathogenic", "Benign", "Pathogenic", "Benign", "Uncertain significance"])
    assigned_pathogenic = pd.Series([True, False, False, False, False])
    assigned_benign = pd.Series([False, True, True, False, False])

    flags, in_scope = control_concordance_flags(
        control_group, PATHOGENIC_VALUES, BENIGN_VALUES, assigned_pathogenic, assigned_benign
    )

    # row 5 (Uncertain significance) is out of scope regardless of its flags.
    assert list(in_scope) == [True, True, True, True, False]
    assert list(flags[CONCORDANT_LABEL]) == [True, True, False, False, False]
    assert list(flags[DISCORDANT_LABEL]) == [False, False, True, False, False]
    assert list(flags[CONTROL_VUS_LABEL]) == [False, False, False, True, True]
    # row 3 (control Pathogenic, assigned benign) is the sole discordant row, in the
    # PLP-to-BLB direction; no row is discordant in the opposite direction.
    assert list(flags[DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL]) == [False, False, True, False, False]
    assert list(flags[DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL]) == [False, False, False, False, False]


def _write_control_concordance_workbook(path):
    """5 ClinVar control rows + 3 ClinGen control rows per predictor, worked
    out by hand in the test docstrings below (OddsPath alone: sign of
    OP_points, read only from the REVEL sheet; combined with each predictor:
    that predictor's own Class_* category membership on its own sheet).

    Each row also carries a `Gene` -- 3 distinct genes across the 5 ClinVar
    rows (GENEA x2, GENEB x2, GENEC x1), 2 distinct across the 3 ClinGen rows
    (GENEA x2, GENEB x1) -- reused identically across all three predictors'
    sheets, since every sheet shares the same rows here and only the
    evidence columns (`OP_points`/`Class_*`) vary by predictor.
    """
    clinvar_genes = ["GENEA", "GENEA", "GENEB", "GENEB", "GENEC"]
    clingen_genes = ["GENEA", "GENEA", "GENEB"]
    clinvar_rows_revel = [
        # (clnsig_group_18_25, OP_points, Class_REVEL)
        ("Pathogenic", 2, "Pathogenic"),  # both concordant
        ("Benign", -1, "Uncertain"),  # OddsPath concordant, combined VUS
        ("Pathogenic", -1, "Likely Pathogenic"),  # OddsPath discordant, combined concordant
        ("Benign", 0, "Benign"),  # OddsPath VUS (no evidence), combined concordant
        ("Likely pathogenic", None, "Likely Benign"),  # OddsPath VUS, combined discordant
    ]
    clingen_rows_revel = [
        # (Updated_Classification_ClinGen_repo, OP_points, Class_REVEL)
        ("Pathogenic", 1, "Pathogenic"),  # both concordant
        ("Benign", 1, "Uncertain"),  # OddsPath discordant, combined VUS
        ("Likely Pathogenic", None, "Likely Benign"),  # OddsPath VUS, combined discordant
    ]
    # AM: concordant=3, discordant=1, vus=1 (ClinVar); concordant=2, discordant=0, vus=1 (ClinGen)
    clinvar_rows_am = [
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Benign"),
        ("Pathogenic", "Benign"),
        ("Benign", "Uncertain"),
        ("Likely pathogenic", "Likely Pathogenic"),
    ]
    clingen_rows_am = [
        ("Pathogenic", "Likely Pathogenic"),
        ("Benign", "Likely Benign"),
        ("Likely Pathogenic", "Uncertain"),
    ]
    # MP2: concordant=2, discordant=2, vus=1 (ClinVar); concordant=2, discordant=1, vus=0 (ClinGen)
    clinvar_rows_mp2 = [
        ("Pathogenic", "Uncertain"),
        ("Benign", "Pathogenic"),
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Benign"),
        ("Likely pathogenic", "Benign"),
    ]
    clingen_rows_mp2 = [
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Pathogenic"),
        ("Likely Pathogenic", "Likely Pathogenic"),
    ]

    def _with_genes(rows, genes):
        return [(*row, gene) for row, gene in zip(rows, genes)]

    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(
            _with_genes(clinvar_rows_revel, clinvar_genes), columns=["clnsig_group_18_25", "OP_points", "Class_REVEL", "Gene"]
        ).to_excel(writer, sheet_name="controls_REVEL_GeneSpecific", index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_revel, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "OP_points", "Class_REVEL", "Gene"],
        ).to_excel(writer, sheet_name="ClinGen_Repo_REVEL_GeneSpecific", index=False)
        pd.DataFrame(
            _with_genes(clinvar_rows_am, clinvar_genes), columns=["clnsig_group_18_25", "Class_AM", "Gene"]
        ).to_excel(writer, sheet_name="controls_AM_GeneSpecific", index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_am, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "Class_AM", "Gene"],
        ).to_excel(writer, sheet_name="ClinGen_Repo_AM_GeneSpecific", index=False)
        pd.DataFrame(
            _with_genes(clinvar_rows_mp2, clinvar_genes), columns=["clnsig_group_18_25", "Class_MP2", "Gene"]
        ).to_excel(writer, sheet_name="controls_MP2_GeneSpecific", index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_mp2, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "Class_MP2", "Gene"],
        ).to_excel(writer, sheet_name="ClinGen_Repo_MP2_GeneSpecific", index=False)


def test_compute_control_concordance_clinvar_and_clingen(tmp_path):
    path = tmp_path / "controls.xlsx"
    _write_control_concordance_workbook(path)

    concordance = compute_control_concordance(pd.ExcelFile(path))

    clinvar_oddspath_total, clinvar_oddspath_table, clinvar_oddspath_genes = concordance[
        ("ClinVar", CONTROL_CONCORDANCE_EVIDENCE_LABEL)
    ]
    assert clinvar_oddspath_total == 5
    assert clinvar_oddspath_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clinvar_oddspath_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clinvar_oddspath_table.loc[CONTROL_VUS_LABEL, "count"] == 2
    assert clinvar_oddspath_table.loc[CONCORDANT_LABEL, "pct"] == pytest.approx(40.0)
    # GENEA, GENEB, GENEC across the 5 ClinVar rows.
    assert clinvar_oddspath_genes == 3

    clinvar_combined_total, clinvar_combined_table, _clinvar_combined_genes = concordance[
        ("ClinVar", COMBINED_REVEL_EVIDENCE_LABEL)
    ]
    assert clinvar_combined_total == 5
    assert clinvar_combined_table.loc[CONCORDANT_LABEL, "count"] == 3
    assert clinvar_combined_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clinvar_combined_table.loc[CONTROL_VUS_LABEL, "count"] == 1

    clingen_oddspath_total, clingen_oddspath_table, clingen_oddspath_genes = concordance[
        ("ClinGen", CONTROL_CONCORDANCE_EVIDENCE_LABEL)
    ]
    assert clingen_oddspath_total == 3
    assert clingen_oddspath_table.loc[CONCORDANT_LABEL, "count"] == 1
    assert clingen_oddspath_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clingen_oddspath_table.loc[CONTROL_VUS_LABEL, "count"] == 1
    # GENEA (rows 1-2) and GENEB (row 3) across the 3 ClinGen rows; rows 1 and 3
    # (Pathogenic, Likely Pathogenic) are PLP, row 2 (Benign) is BLB.
    assert clingen_oddspath_genes == 2
    assert clingen_oddspath_table.loc[CONTROL_PLP_LABEL, "count"] == 2
    assert clingen_oddspath_table.loc[CONTROL_BLB_LABEL, "count"] == 1

    clingen_combined_total, clingen_combined_table, clingen_combined_genes = concordance[
        ("ClinGen", COMBINED_REVEL_EVIDENCE_LABEL)
    ]
    assert clingen_combined_total == 3
    assert clingen_combined_table.loc[CONCORDANT_LABEL, "count"] == 1
    assert clingen_combined_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clingen_combined_table.loc[CONTROL_VUS_LABEL, "count"] == 1
    assert clingen_combined_genes == 2
    assert clingen_combined_table.loc[CONTROL_PLP_LABEL, "count"] == 2
    assert clingen_combined_table.loc[CONTROL_BLB_LABEL, "count"] == 1

    am_label = COMBINED_EVIDENCE_LABEL_BY_PREDICTOR["AlphaMissense"]
    clinvar_am_total, clinvar_am_table, _clinvar_am_genes = concordance[("ClinVar", am_label)]
    assert clinvar_am_total == 5
    assert clinvar_am_table.loc[CONCORDANT_LABEL, "count"] == 3
    assert clinvar_am_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clinvar_am_table.loc[CONTROL_VUS_LABEL, "count"] == 1

    clingen_am_total, clingen_am_table, clingen_am_genes = concordance[("ClinGen", am_label)]
    assert clingen_am_total == 3
    assert clingen_am_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clingen_am_table.loc[DISCORDANT_LABEL, "count"] == 0
    assert clingen_am_table.loc[CONTROL_VUS_LABEL, "count"] == 1
    assert clingen_am_genes == 2
    assert clingen_am_table.loc[CONTROL_PLP_LABEL, "count"] == 2
    assert clingen_am_table.loc[CONTROL_BLB_LABEL, "count"] == 1

    mp2_label = COMBINED_EVIDENCE_LABEL_BY_PREDICTOR["MutPred2"]
    clinvar_mp2_total, clinvar_mp2_table, _clinvar_mp2_genes = concordance[("ClinVar", mp2_label)]
    assert clinvar_mp2_total == 5
    assert clinvar_mp2_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clinvar_mp2_table.loc[DISCORDANT_LABEL, "count"] == 2
    assert clinvar_mp2_table.loc[CONTROL_VUS_LABEL, "count"] == 1
    # Row ("Benign", "Pathogenic") is the sole BLB-to-PLP discordance; row
    # ("Likely pathogenic", "Benign") is the sole PLP-to-BLB discordance.
    assert clinvar_mp2_table.loc[DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL, "count"] == 1
    assert clinvar_mp2_table.loc[DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL, "count"] == 1

    clingen_mp2_total, clingen_mp2_table, clingen_mp2_genes = concordance[("ClinGen", mp2_label)]
    assert clingen_mp2_total == 3
    assert clingen_mp2_genes == 2
    assert clingen_mp2_table.loc[CONTROL_PLP_LABEL, "count"] == 2
    assert clingen_mp2_table.loc[CONTROL_BLB_LABEL, "count"] == 1
    assert clingen_mp2_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clingen_mp2_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clingen_mp2_table.loc[CONTROL_VUS_LABEL, "count"] == 0


def test_format_control_concordance_report(tmp_path):
    path = tmp_path / "controls.xlsx"
    _write_control_concordance_workbook(path)
    concordance = compute_control_concordance(pd.ExcelFile(path))

    text = format_control_concordance_report(concordance)

    assert (
        "=== Control concordance (ClinVar vs. ClinGen; OddsPath alone vs. combined with "
        "REVEL/AlphaMissense/MutPred2) ===" in text
    )
    assert "ClinVar controls (controls_*_GeneSpecific sheets):" in text
    assert "ClinGen controls (ClinGen_Repo_*_GeneSpecific sheets):" in text
    assert CONTROL_CONCORDANCE_EVIDENCE_LABEL in text
    assert COMBINED_REVEL_EVIDENCE_LABEL in text
    assert COMBINED_EVIDENCE_LABEL_BY_PREDICTOR["AlphaMissense"] in text
    assert COMBINED_EVIDENCE_LABEL_BY_PREDICTOR["MutPred2"] in text
    # ClinVar/OddsPath row: 2 concordant of 5 total.
    assert "2 of 5 (40.0%)" in text
    # ClinGen/AlphaMissense row: 2 concordant of 3 total.
    assert "2 of 3 (66.7%)" in text
    assert DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL.strip() in text
    assert DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL.strip() in text
    # ClinVar/MutPred2 row: 1 PLP-to-BLB and 1 BLB-to-PLP discordance, each of 5 total.
    assert text.count("1 of 5 (20.0%)") >= 2

    clinvar_section, clingen_section = text.split("ClinGen controls (ClinGen_Repo_*_GeneSpecific sheets):")
    # Genes/PLP/BLB are ClinGen-only -- ClinVar's own control set doesn't report them.
    assert "Genes" not in clinvar_section
    assert "Genes" in clingen_section
    # Every ClinGen row: 2 genes (GENEA, GENEB), PLP=2 of 3, BLB=1 of 3.
    assert clingen_section.count("2 of 3 (66.7%)") >= 4  # PLP, plus each row's own concordance count
    assert clingen_section.count("1 of 3 (33.3%)") >= 4  # BLB, plus each row's own discordance count


def _write_universal_control_concordance_workbook(path):
    """Supplementary Data 6-style counterpart to
    `_write_control_concordance_workbook`, for the universal (genome-wide)
    calibration evidence source -- 5 ClinVar + 3 ClinGen rows per predictor's
    own `UNIVERSAL_CALIBRATION_CATEGORY_SHEETS_BY_PREDICTOR` sheet/
    `UNIVERSAL_CALIBRATION_CLASS_COL_BY_PREDICTOR` column. Deliberately
    different concordance patterns from `_write_control_concordance_workbook`
    (which a test reading the wrong workbook/column would fail to reproduce).
    """
    clinvar_genes = ["GENEA", "GENEA", "GENEB", "GENEB", "GENEC"]
    clingen_genes = ["GENEA", "GENEA", "GENEB"]
    # REVEL universal: concordant=2, discordant=1, vus=2 (ClinVar); concordant=2, discordant=0, vus=1 (ClinGen)
    clinvar_rows_revel = [
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Likely Benign"),
        ("Pathogenic", "Benign"),
        ("Benign", "Uncertain"),
        ("Likely pathogenic", "Uncertain"),
    ]
    clingen_rows_revel = [
        ("Pathogenic", "Likely Pathogenic"),
        ("Benign", "Benign"),
        ("Likely Pathogenic", "Uncertain"),
    ]
    # AM universal: concordant=3, discordant=1, vus=1 (ClinVar); concordant=2, discordant=0, vus=1 (ClinGen)
    clinvar_rows_am = [
        ("Pathogenic", "Uncertain"),
        ("Benign", "Pathogenic"),
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Benign"),
        ("Likely pathogenic", "Likely Pathogenic"),
    ]
    clingen_rows_am = [
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Uncertain"),
        ("Likely Pathogenic", "Likely Pathogenic"),
    ]
    # MP2 universal: concordant=2, discordant=2, vus=1 (ClinVar); concordant=2, discordant=1, vus=0 (ClinGen)
    clinvar_rows_mp2 = [
        ("Pathogenic", "Benign"),
        ("Benign", "Pathogenic"),
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Benign"),
        ("Likely pathogenic", "Uncertain"),
    ]
    clingen_rows_mp2 = [
        ("Pathogenic", "Pathogenic"),
        ("Benign", "Pathogenic"),
        ("Likely Pathogenic", "Likely Pathogenic"),
    ]

    def _with_genes(rows, genes):
        return [(*row, gene) for row, gene in zip(rows, genes)]

    def _sheet(predictor, category):
        return UNIVERSAL_CALIBRATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor][category]

    with pd.ExcelWriter(path) as writer:
        pd.DataFrame(
            _with_genes(clinvar_rows_revel, clinvar_genes), columns=["clnsig_group_18_25", "Class_OP_REVEL", "Gene"]
        ).to_excel(writer, sheet_name=_sheet("REVEL", "controls"), index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_revel, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "Class_OP_REVEL", "Gene"],
        ).to_excel(writer, sheet_name=_sheet("REVEL", "ClinGen_Repo"), index=False)
        pd.DataFrame(
            _with_genes(clinvar_rows_am, clinvar_genes), columns=["clnsig_group_18_25", "Class_OP_AM", "Gene"]
        ).to_excel(writer, sheet_name=_sheet("AlphaMissense", "controls"), index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_am, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "Class_OP_AM", "Gene"],
        ).to_excel(writer, sheet_name=_sheet("AlphaMissense", "ClinGen_Repo"), index=False)
        pd.DataFrame(
            _with_genes(clinvar_rows_mp2, clinvar_genes), columns=["clnsig_group_18_25", "Class_OP_MP2", "Gene"]
        ).to_excel(writer, sheet_name=_sheet("MutPred2", "controls"), index=False)
        pd.DataFrame(
            _with_genes(clingen_rows_mp2, clingen_genes),
            columns=["Updated_Classification_ClinGen_repo", "Class_OP_MP2", "Gene"],
        ).to_excel(writer, sheet_name=_sheet("MutPred2", "ClinGen_Repo"), index=False)


def test_compute_control_concordance_with_universal_workbook(tmp_path):
    path = tmp_path / "controls.xlsx"
    universal_path = tmp_path / "universal_controls.xlsx"
    _write_control_concordance_workbook(path)
    _write_universal_control_concordance_workbook(universal_path)

    concordance = compute_control_concordance(pd.ExcelFile(path), universal_workbook=pd.ExcelFile(universal_path))

    revel_universal_label = COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR["REVEL"]
    clinvar_total, clinvar_table, clinvar_genes = concordance[("ClinVar", revel_universal_label)]
    assert clinvar_total == 5
    assert clinvar_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clinvar_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clinvar_table.loc[CONTROL_VUS_LABEL, "count"] == 2
    assert clinvar_genes == 3

    clingen_total, clingen_table, clingen_genes = concordance[("ClinGen", revel_universal_label)]
    assert clingen_total == 3
    assert clingen_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clingen_table.loc[DISCORDANT_LABEL, "count"] == 0
    assert clingen_table.loc[CONTROL_VUS_LABEL, "count"] == 1
    assert clingen_genes == 2

    am_label = COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR["AlphaMissense"]
    clinvar_am_total, clinvar_am_table, _clinvar_am_genes = concordance[("ClinVar", am_label)]
    assert clinvar_am_total == 5
    assert clinvar_am_table.loc[CONCORDANT_LABEL, "count"] == 3
    assert clinvar_am_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clinvar_am_table.loc[CONTROL_VUS_LABEL, "count"] == 1

    mp2_label = COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR["MutPred2"]
    clinvar_mp2_total, clinvar_mp2_table, _clinvar_mp2_genes = concordance[("ClinVar", mp2_label)]
    assert clinvar_mp2_total == 5
    assert clinvar_mp2_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clinvar_mp2_table.loc[DISCORDANT_LABEL, "count"] == 2
    assert clinvar_mp2_table.loc[DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL, "count"] == 1
    assert clinvar_mp2_table.loc[DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL, "count"] == 1

    # Without universal_workbook (the default), no universal-evidence keys are computed at all.
    gene_specific_only = compute_control_concordance(pd.ExcelFile(path))
    assert ("ClinVar", revel_universal_label) not in gene_specific_only


def test_format_control_concordance_report_includes_universal_rows(tmp_path):
    path = tmp_path / "controls.xlsx"
    universal_path = tmp_path / "universal_controls.xlsx"
    _write_control_concordance_workbook(path)
    _write_universal_control_concordance_workbook(universal_path)

    concordance = compute_control_concordance(pd.ExcelFile(path), universal_workbook=pd.ExcelFile(universal_path))
    text = format_control_concordance_report(concordance)

    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        gene_specific_label = COMBINED_EVIDENCE_LABEL_BY_PREDICTOR[predictor]
        universal_label = COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR[predictor]
        assert universal_label in text
        # Each predictor's universal row is rendered immediately after its own
        # gene-specific row.
        assert text.index(gene_specific_label) < text.index(universal_label)

    # Without a universal_workbook, compute_control_concordance never produces
    # universal-evidence keys, so format_control_concordance_report shows none
    # of these labels either.
    gene_specific_only_text = format_control_concordance_report(compute_control_concordance(pd.ExcelFile(path)))
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert COMBINED_UNIVERSAL_EVIDENCE_LABEL_BY_PREDICTOR[predictor] not in gene_specific_only_text


def _write_missense_control_concordance_workbook(path):
    """One `controls_<suffix>_GeneSpecific` (ClinVar) and one
    `ClinGen_Repo_<suffix>_GeneSpecific` (ClinGen) sheet per predictor, each
    with the same 4 rows: 3 missense_variant (2 concordant, 1 discordant --
    control Pathogenic reclassified Benign) plus 1 stop_gained row that
    should be excluded entirely once `compute_control_concordance` is called
    with `consequence_filter=MISSENSE_CONSEQUENCE_VALUE`. ClinVar's 3 missense
    rows span 2 genes (GENEA x2, GENEB x1); ClinGen's span 2 different genes
    (GENED x2, GENEE x1).
    """
    # (control classification, OP_points, Class_<predictor>, simplified_consequence, Gene)
    clinvar_rows = [
        ("Pathogenic", 5, "Pathogenic", MISSENSE_CONSEQUENCE_VALUE, "GENEA"),  # missense, concordant
        ("Benign", -1, "Likely Benign", MISSENSE_CONSEQUENCE_VALUE, "GENEA"),  # missense, concordant
        ("Pathogenic", -2, "Benign", MISSENSE_CONSEQUENCE_VALUE, "GENEB"),  # missense, discordant: PLP -> BLB
        ("Pathogenic", 5, "Pathogenic", "stop_gained", "GENEC"),  # not missense -- excluded by the filter
    ]
    clingen_rows = [
        ("Pathogenic", 5, "Pathogenic", MISSENSE_CONSEQUENCE_VALUE, "GENED"),  # missense, concordant
        ("Benign", -1, "Likely Benign", MISSENSE_CONSEQUENCE_VALUE, "GENED"),  # missense, concordant
        ("Likely Pathogenic", -2, "Benign", MISSENSE_CONSEQUENCE_VALUE, "GENEE"),  # missense, discordant: PLP -> BLB
        ("Pathogenic", 5, "Pathogenic", "stop_gained", "GENEF"),  # not missense -- excluded by the filter
    ]
    with pd.ExcelWriter(path) as writer:
        for suffix in ("REVEL", "AM", "MP2"):
            clinvar_columns = ["clnsig_group_18_25", "OP_points", f"Class_{suffix}", SIMPLIFIED_CONSEQUENCE_COL, "Gene"]
            pd.DataFrame(clinvar_rows, columns=clinvar_columns).to_excel(
                writer, sheet_name=f"controls_{suffix}_GeneSpecific", index=False
            )
            clingen_columns = [
                "Updated_Classification_ClinGen_repo",
                "OP_points",
                f"Class_{suffix}",
                SIMPLIFIED_CONSEQUENCE_COL,
                "Gene",
            ]
            pd.DataFrame(clingen_rows, columns=clingen_columns).to_excel(
                writer, sheet_name=f"ClinGen_Repo_{suffix}_GeneSpecific", index=False
            )


def test_compute_control_concordance_missense_only(tmp_path):
    path = tmp_path / "controls.xlsx"
    _write_missense_control_concordance_workbook(path)

    concordance = compute_control_concordance(
        pd.ExcelFile(path),
        control_sources=MISSENSE_CONTROL_CONCORDANCE_SOURCES,
        consequence_filter=MISSENSE_CONSEQUENCE_VALUE,
    )

    # The stop_gained row is filtered out entirely, leaving 3 in-scope rows for
    # every evidence source (OddsPath alone, reusing the REVEL sheet, and each
    # predictor's own combined evidence, since all three sheets share the same rows).
    oddspath_total, oddspath_table, _oddspath_genes = concordance[("ClinVar", CONTROL_CONCORDANCE_EVIDENCE_LABEL)]
    assert oddspath_total == 3
    assert oddspath_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert oddspath_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert oddspath_table.loc[DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL, "count"] == 1
    assert oddspath_table.loc[DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL, "count"] == 0

    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        total, table, _genes = concordance[("ClinVar", COMBINED_EVIDENCE_LABEL_BY_PREDICTOR[predictor])]
        assert total == 3
        assert table.loc[CONCORDANT_LABEL, "count"] == 2
        assert table.loc[DISCORDANT_LABEL, "count"] == 1

    # ClinGen is also part of the missense-only breakdown, with its own
    # (different) rows/genes -- the stop_gained row is excluded the same way.
    clingen_oddspath_total, clingen_oddspath_table, clingen_oddspath_genes = concordance[
        ("ClinGen", CONTROL_CONCORDANCE_EVIDENCE_LABEL)
    ]
    assert clingen_oddspath_total == 3
    assert clingen_oddspath_table.loc[CONCORDANT_LABEL, "count"] == 2
    assert clingen_oddspath_table.loc[DISCORDANT_LABEL, "count"] == 1
    assert clingen_oddspath_table.loc[DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL, "count"] == 1
    # GENED (rows 1-2) and GENEE (row 3); rows 1 and 3 (Pathogenic, Likely
    # Pathogenic) are PLP, row 2 (Benign) is BLB.
    assert clingen_oddspath_genes == 2
    assert clingen_oddspath_table.loc[CONTROL_PLP_LABEL, "count"] == 2
    assert clingen_oddspath_table.loc[CONTROL_BLB_LABEL, "count"] == 1

    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        total, table, genes = concordance[("ClinGen", COMBINED_EVIDENCE_LABEL_BY_PREDICTOR[predictor])]
        assert total == 3
        assert table.loc[CONCORDANT_LABEL, "count"] == 2
        assert table.loc[DISCORDANT_LABEL, "count"] == 1
        assert genes == 2


def test_format_control_concordance_report_missense_only(tmp_path):
    path = tmp_path / "controls.xlsx"
    _write_missense_control_concordance_workbook(path)
    concordance = compute_control_concordance(
        pd.ExcelFile(path),
        control_sources=MISSENSE_CONTROL_CONCORDANCE_SOURCES,
        consequence_filter=MISSENSE_CONSEQUENCE_VALUE,
    )

    text = format_control_concordance_report(
        concordance, control_sources=MISSENSE_CONTROL_CONCORDANCE_SOURCES, title=MISSENSE_CONTROL_CONCORDANCE_TITLE
    )

    assert text.startswith(MISSENSE_CONTROL_CONCORDANCE_TITLE)
    assert "ClinVar controls (controls_*_GeneSpecific sheets):" in text
    assert "ClinGen controls (ClinGen_Repo_*_GeneSpecific sheets):" in text

    clinvar_section, clingen_section = text.split("ClinGen controls (ClinGen_Repo_*_GeneSpecific sheets):")
    # ClinVar: 2 concordant, 1 discordant (all PLP-to-BLB), of 3 total, for every evidence source.
    assert "Genes" not in clinvar_section
    assert clinvar_section.count("2 of 3 (66.7%)") == 4
    assert clinvar_section.count("1 of 3 (33.3%)") == 4 * 2  # Discordant + its PLP-to-BLB sub-column

    # ClinGen: same shape, plus Genes=2 and a PLP/BLB population split (2 of 3 / 1 of 3).
    assert "Genes" in clingen_section
    assert clingen_section.count("2 of 3 (66.7%)") == 4 * 2  # Concordant + PLP
    assert clingen_section.count("1 of 3 (33.3%)") == 4 * 3  # Discordant + its PLP-to-BLB sub-column + BLB


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
    # that the combined total counts it once, not twice. Same rows written to
    # all three predictors' sheets, so each predictor's stats should match.
    # Row 105 is single-source (Fxn_points 0, predictor points -1); row 108 is
    # conflicting (Fxn_points +2, predictor points -3, opposite signs, still
    # summing to the -1-point threshold) -- together they cover both halves
    # of the mutually-exclusive-and-exhaustive -1-point split.
    _write_variant_classification_sheets_by_predictor(
        controls_path,
        {
            "controls": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic", 6, 6, 0),
                ("GENEX", 1, 101, "A", "G", "Benign", -2, -2, 0),
            ],
            "ClinGen_Repo": [
                ("GENEX", 1, 102, "A", "G", "Likely Pathogenic", 6, 6, 0),
            ],
            "VUS": [
                ("GENEX", 1, 103, "A", "G", "Pathogenic", 6, 6, 0),
                ("GENEX", 1, 104, "A", "G", "Uncertain", 2, 2, 0),
                ("GENEX", 1, 105, "A", "G", "Benign", -1, 0, -1),  # -1-point threshold, single source
                ("GENEX", 1, 108, "A", "G", "Benign", -1, 2, -3),  # -1-point threshold, conflicting
                ("GENEX", 1, 109, "A", "G", "Uncertain", 5, 2, 3),  # unresolved, both sources, concordant
                ("GENEX", 1, 110, "A", "G", "Uncertain", 3, 8, -5),  # unresolved, both sources, discordant
            ],
            "gnomAD": [
                ("GENEX", 1, 100, "A", "G", "Pathogenic", 6, 6, 0),
            ],
            "Unobserved": [
                ("GENEX", 1, 106, "A", "G", "Pathogenic", 6, 6, 0),
                ("GENEX", 1, 107, "A", "G", "Uncertain", 3, 3, 0),
            ],
        },
        mode="w",
    )

    stats_by_predictor = compute_variant_classification_stats(pd.ExcelFile(controls_path))

    assert set(stats_by_predictor) == set(VARIANT_CLASSIFICATION_PREDICTORS)
    expected = {
        "total_classified": 11,
        "total_pathogenic_or_benign": 7,
        "vus_total": 6,
        "vus_resolved": 3,
        # row 103 (fxn 6, predictor 0): experimental only.
        "vus_resolved_pathogenic": 1,
        "vus_resolved_pathogenic_only_experimental": 1,
        "vus_resolved_pathogenic_only_predictive": 0,
        "vus_resolved_pathogenic_both_evidence": 0,
        # row 105 (fxn 0, predictor -1): predictive only. row 108 (fxn 2, predictor -3): both.
        "vus_resolved_benign": 2,
        "vus_resolved_benign_only_experimental": 0,
        "vus_resolved_benign_only_predictive": 1,
        "vus_resolved_benign_both_evidence": 1,
        "vus_resolved_benign_at_threshold": 2,
        "vus_resolved_benign_at_threshold_single_source": 1,
        "vus_resolved_benign_at_threshold_conflicting": 1,
        # row 104 (fxn 2, predictor 0): experimental only. row 109 (fxn 2, predictor 3, both
        # same sign): concordant. row 110 (fxn 8, predictor -5, opposite signs): discordant.
        "vus_unresolved": 3,
        "vus_unresolved_concordant": 1,
        "vus_unresolved_discordant": 1,
        "vus_unresolved_only_experimental": 1,
        "vus_unresolved_only_predictive": 0,
        "vus_unresolved_neither": 0,
        "vus_unresolved_zero_or_one_source": 1,
        # row 109 (points 5) is the only VUS near-pathogenic row (row 110's points, 3, don't qualify).
        "vus_unresolved_near_pathogenic": 1,
        # the sole gnomAD row (same coords as controls row 100, fxn 6/predictor 0): resolved
        # pathogenic, experimental only.
        "gnomad_total": 1,
        "gnomad_resolved": 1,
        "gnomad_resolved_pathogenic": 1,
        "gnomad_resolved_pathogenic_only_experimental": 1,
        "gnomad_resolved_pathogenic_only_predictive": 0,
        "gnomad_resolved_pathogenic_both_evidence": 0,
        "gnomad_resolved_benign": 0,
        "gnomad_resolved_benign_only_experimental": 0,
        "gnomad_resolved_benign_only_predictive": 0,
        "gnomad_resolved_benign_both_evidence": 0,
        "gnomad_resolved_benign_at_threshold": 0,
        "gnomad_resolved_benign_at_threshold_single_source": 0,
        "gnomad_resolved_benign_at_threshold_conflicting": 0,
        "gnomad_unresolved": 0,
        "gnomad_unresolved_concordant": 0,
        "gnomad_unresolved_discordant": 0,
        "gnomad_unresolved_only_experimental": 0,
        "gnomad_unresolved_only_predictive": 0,
        "gnomad_unresolved_neither": 0,
        "gnomad_unresolved_zero_or_one_source": 0,
        "gnomad_unresolved_near_pathogenic": 0,
        "unobserved_total": 2,
        "unobserved_resolved": 1,
        # row 106 (fxn 6, predictor 0): experimental only.
        "unobserved_resolved_pathogenic": 1,
        "unobserved_resolved_pathogenic_only_experimental": 1,
        "unobserved_resolved_pathogenic_only_predictive": 0,
        "unobserved_resolved_pathogenic_both_evidence": 0,
        "unobserved_resolved_benign": 0,
        "unobserved_resolved_benign_only_experimental": 0,
        "unobserved_resolved_benign_only_predictive": 0,
        "unobserved_resolved_benign_both_evidence": 0,
        "unobserved_resolved_benign_at_threshold": 0,
        "unobserved_resolved_benign_at_threshold_single_source": 0,
        "unobserved_resolved_benign_at_threshold_conflicting": 0,
        # row 107 (fxn 3, predictor 0), the sole unresolved Unobserved: experimental only.
        "unobserved_unresolved": 1,
        "unobserved_unresolved_concordant": 0,
        "unobserved_unresolved_discordant": 0,
        "unobserved_unresolved_only_experimental": 1,
        "unobserved_unresolved_only_predictive": 0,
        "unobserved_unresolved_neither": 0,
        "unobserved_unresolved_zero_or_one_source": 1,
        "unobserved_unresolved_near_pathogenic": 0,
    }
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert stats_by_predictor[predictor] == expected


def _variant_classification_stats_by_predictor(stats):
    return {predictor: stats for predictor in VARIANT_CLASSIFICATION_PREDICTORS}


def test_format_variant_classification_table():
    stats = {
        "total_classified": 8,
        "total_pathogenic_or_benign": 6,
        "vus_total": 3,
        "vus_resolved": 2,
        # VUS P/LP (1) is experimental-only; VUS B/LB (1) is predictive-only.
        "vus_resolved_pathogenic": 1,
        "vus_resolved_pathogenic_only_experimental": 1,
        "vus_resolved_pathogenic_only_predictive": 0,
        "vus_resolved_pathogenic_both_evidence": 0,
        "vus_resolved_benign": 1,
        "vus_resolved_benign_only_experimental": 0,
        "vus_resolved_benign_only_predictive": 1,
        "vus_resolved_benign_both_evidence": 0,
        "vus_resolved_benign_at_threshold": 1,
        "vus_resolved_benign_at_threshold_single_source": 1,
        "vus_resolved_benign_at_threshold_conflicting": 0,
        # vus_total (3) - vus_resolved (2) = 1 unresolved row, experimental-only.
        "vus_unresolved": 1,
        "vus_unresolved_concordant": 0,
        "vus_unresolved_discordant": 0,
        "vus_unresolved_only_experimental": 1,
        "vus_unresolved_only_predictive": 0,
        "vus_unresolved_neither": 0,
        "vus_unresolved_zero_or_one_source": 1,
        "vus_unresolved_near_pathogenic": 0,
        # Same shape as VUS above: P/LP experimental-only, B/LB predictive-only, 1 unresolved
        # experimental-only row.
        "gnomad_total": 3,
        "gnomad_resolved": 2,
        "gnomad_resolved_pathogenic": 1,
        "gnomad_resolved_pathogenic_only_experimental": 1,
        "gnomad_resolved_pathogenic_only_predictive": 0,
        "gnomad_resolved_pathogenic_both_evidence": 0,
        "gnomad_resolved_benign": 1,
        "gnomad_resolved_benign_only_experimental": 0,
        "gnomad_resolved_benign_only_predictive": 1,
        "gnomad_resolved_benign_both_evidence": 0,
        "gnomad_resolved_benign_at_threshold": 1,
        "gnomad_resolved_benign_at_threshold_single_source": 1,
        "gnomad_resolved_benign_at_threshold_conflicting": 0,
        "gnomad_unresolved": 1,
        "gnomad_unresolved_concordant": 0,
        "gnomad_unresolved_discordant": 0,
        "gnomad_unresolved_only_experimental": 1,
        "gnomad_unresolved_only_predictive": 0,
        "gnomad_unresolved_neither": 0,
        "gnomad_unresolved_zero_or_one_source": 1,
        "gnomad_unresolved_near_pathogenic": 0,
        "unobserved_total": 2,
        "unobserved_resolved": 1,
        # Unobserved P/LP (1) is both-evidence; Unobserved has no resolved-benign row.
        "unobserved_resolved_pathogenic": 1,
        "unobserved_resolved_pathogenic_only_experimental": 0,
        "unobserved_resolved_pathogenic_only_predictive": 0,
        "unobserved_resolved_pathogenic_both_evidence": 1,
        "unobserved_resolved_benign": 0,
        "unobserved_resolved_benign_only_experimental": 0,
        "unobserved_resolved_benign_only_predictive": 0,
        "unobserved_resolved_benign_both_evidence": 0,
        "unobserved_resolved_benign_at_threshold": 0,
        "unobserved_resolved_benign_at_threshold_single_source": 0,
        "unobserved_resolved_benign_at_threshold_conflicting": 0,
        # unobserved_total (2) - unobserved_resolved (1) = 1 unresolved row, concordant.
        "unobserved_unresolved": 1,
        "unobserved_unresolved_concordant": 1,
        "unobserved_unresolved_discordant": 0,
        "unobserved_unresolved_only_experimental": 0,
        "unobserved_unresolved_only_predictive": 0,
        "unobserved_unresolved_neither": 0,
        "unobserved_unresolved_zero_or_one_source": 0,
        "unobserved_unresolved_near_pathogenic": 0,
    }

    text = format_variant_classification_table(_variant_classification_stats_by_predictor(stats))

    assert text.startswith(VARIANT_CLASSIFICATION_TITLE)
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        assert predictor in text
    assert "ClinVar VUS resolved (reclassified pathogenic or benign):" in text
    assert "ClinVar VUS unresolved:" in text
    assert "gnomAD variants resolved (classified pathogenic or benign):" in text
    assert "gnomAD variants unresolved:" in text
    assert "Unobserved variants resolved (classified pathogenic or benign):" in text
    assert "Unobserved variants unresolved:" in text
    # One row per predictor per table -- these values are identical across
    # REVEL/AlphaMissense/MutPred2 here, so each formatted count/pct string
    # appears once per predictor occurrence (x3). Counts below were verified against the
    # actual rendered table rather than hand-derived, given how many columns now overlap.
    assert text.count("6 of 8 (75.0%)") == 3
    assert text.count("2 of 3 (66.7%)") == 6  # VUS and gnomAD's "Resolved" columns
    assert text.count("1 of 3 (33.3%)") == 24
    assert text.count("1 of 2 (50.0%)") == 9
    assert text.count("0 of 2 (0.0%)") == 6
    assert text.count("1 of 1 (100.0%)") == 36
    assert text.count("0 of 1 (0.0%)") == 84
    assert text.count("0 of 0 (nan%)") == 15


def test_format_variant_classification_table_uses_given_title():
    stats = dict.fromkeys(
        [
            "total_classified",
            "total_pathogenic_or_benign",
            "vus_total",
            "vus_resolved",
            "vus_resolved_pathogenic",
            "vus_resolved_pathogenic_only_experimental",
            "vus_resolved_pathogenic_only_predictive",
            "vus_resolved_pathogenic_both_evidence",
            "vus_resolved_benign",
            "vus_resolved_benign_only_experimental",
            "vus_resolved_benign_only_predictive",
            "vus_resolved_benign_both_evidence",
            "vus_resolved_benign_at_threshold",
            "vus_resolved_benign_at_threshold_single_source",
            "vus_resolved_benign_at_threshold_conflicting",
            "vus_unresolved",
            "vus_unresolved_concordant",
            "vus_unresolved_discordant",
            "vus_unresolved_only_experimental",
            "vus_unresolved_only_predictive",
            "vus_unresolved_neither",
            "vus_unresolved_zero_or_one_source",
            "vus_unresolved_near_pathogenic",
            "gnomad_total",
            "gnomad_resolved",
            "gnomad_resolved_pathogenic",
            "gnomad_resolved_pathogenic_only_experimental",
            "gnomad_resolved_pathogenic_only_predictive",
            "gnomad_resolved_pathogenic_both_evidence",
            "gnomad_resolved_benign",
            "gnomad_resolved_benign_only_experimental",
            "gnomad_resolved_benign_only_predictive",
            "gnomad_resolved_benign_both_evidence",
            "gnomad_resolved_benign_at_threshold",
            "gnomad_resolved_benign_at_threshold_single_source",
            "gnomad_resolved_benign_at_threshold_conflicting",
            "gnomad_unresolved",
            "gnomad_unresolved_concordant",
            "gnomad_unresolved_discordant",
            "gnomad_unresolved_only_experimental",
            "gnomad_unresolved_only_predictive",
            "gnomad_unresolved_neither",
            "gnomad_unresolved_zero_or_one_source",
            "gnomad_unresolved_near_pathogenic",
            "unobserved_total",
            "unobserved_resolved",
            "unobserved_resolved_pathogenic",
            "unobserved_resolved_pathogenic_only_experimental",
            "unobserved_resolved_pathogenic_only_predictive",
            "unobserved_resolved_pathogenic_both_evidence",
            "unobserved_resolved_benign",
            "unobserved_resolved_benign_only_experimental",
            "unobserved_resolved_benign_only_predictive",
            "unobserved_resolved_benign_both_evidence",
            "unobserved_resolved_benign_at_threshold",
            "unobserved_resolved_benign_at_threshold_single_source",
            "unobserved_resolved_benign_at_threshold_conflicting",
            "unobserved_unresolved",
            "unobserved_unresolved_concordant",
            "unobserved_unresolved_discordant",
            "unobserved_unresolved_only_experimental",
            "unobserved_unresolved_only_predictive",
            "unobserved_unresolved_neither",
            "unobserved_unresolved_zero_or_one_source",
            "unobserved_unresolved_near_pathogenic",
        ],
        0,
    )

    text = format_variant_classification_table(
        _variant_classification_stats_by_predictor(stats), title="=== custom title ==="
    )

    assert text.startswith("=== custom title ===")


def test_format_variant_classification_table_handles_zero_totals():
    stats = dict.fromkeys(
        [
            "total_classified",
            "total_pathogenic_or_benign",
            "vus_total",
            "vus_resolved",
            "vus_resolved_pathogenic",
            "vus_resolved_pathogenic_only_experimental",
            "vus_resolved_pathogenic_only_predictive",
            "vus_resolved_pathogenic_both_evidence",
            "vus_resolved_benign",
            "vus_resolved_benign_only_experimental",
            "vus_resolved_benign_only_predictive",
            "vus_resolved_benign_both_evidence",
            "vus_resolved_benign_at_threshold",
            "vus_resolved_benign_at_threshold_single_source",
            "vus_resolved_benign_at_threshold_conflicting",
            "vus_unresolved",
            "vus_unresolved_concordant",
            "vus_unresolved_discordant",
            "vus_unresolved_only_experimental",
            "vus_unresolved_only_predictive",
            "vus_unresolved_neither",
            "vus_unresolved_zero_or_one_source",
            "vus_unresolved_near_pathogenic",
            "gnomad_total",
            "gnomad_resolved",
            "gnomad_resolved_pathogenic",
            "gnomad_resolved_pathogenic_only_experimental",
            "gnomad_resolved_pathogenic_only_predictive",
            "gnomad_resolved_pathogenic_both_evidence",
            "gnomad_resolved_benign",
            "gnomad_resolved_benign_only_experimental",
            "gnomad_resolved_benign_only_predictive",
            "gnomad_resolved_benign_both_evidence",
            "gnomad_resolved_benign_at_threshold",
            "gnomad_resolved_benign_at_threshold_single_source",
            "gnomad_resolved_benign_at_threshold_conflicting",
            "gnomad_unresolved",
            "gnomad_unresolved_concordant",
            "gnomad_unresolved_discordant",
            "gnomad_unresolved_only_experimental",
            "gnomad_unresolved_only_predictive",
            "gnomad_unresolved_neither",
            "gnomad_unresolved_zero_or_one_source",
            "gnomad_unresolved_near_pathogenic",
            "unobserved_total",
            "unobserved_resolved",
            "unobserved_resolved_pathogenic",
            "unobserved_resolved_pathogenic_only_experimental",
            "unobserved_resolved_pathogenic_only_predictive",
            "unobserved_resolved_pathogenic_both_evidence",
            "unobserved_resolved_benign",
            "unobserved_resolved_benign_only_experimental",
            "unobserved_resolved_benign_only_predictive",
            "unobserved_resolved_benign_both_evidence",
            "unobserved_resolved_benign_at_threshold",
            "unobserved_resolved_benign_at_threshold_single_source",
            "unobserved_resolved_benign_at_threshold_conflicting",
            "unobserved_unresolved",
            "unobserved_unresolved_concordant",
            "unobserved_unresolved_discordant",
            "unobserved_unresolved_only_experimental",
            "unobserved_unresolved_only_predictive",
            "unobserved_unresolved_neither",
            "unobserved_unresolved_zero_or_one_source",
            "unobserved_unresolved_near_pathogenic",
        ],
        0,
    )

    text = format_variant_classification_table(_variant_classification_stats_by_predictor(stats))

    assert "(nan%)" in text


def test_compute_gene_discordance_stats_counts_and_ranks_genes(tmp_path):
    controls_path = tmp_path / "controls.xlsx"
    _write_gene_discordance_sheet(
        controls_path,
        [
            ("GENEA", 1, 100, "A", "G", "Pathogenic", "Benign"),  # discordant: P/LP -> B/LB
            ("GENEA", 1, 101, "A", "G", "Likely pathogenic", "Likely Benign"),  # discordant: P/LP -> B/LB
            ("GENEA", 1, 102, "A", "G", "Pathogenic", "Pathogenic"),  # concordant
            ("GENEB", 1, 103, "A", "G", "Benign", "Pathogenic"),  # discordant: B/LB -> P/LP
            ("GENEB", 1, 104, "A", "G", "Benign", "Likely Pathogenic"),  # discordant: B/LB -> P/LP
            ("GENEC", 1, 105, "A", "G", "Pathogenic", "Benign"),  # discordant: P/LP -> B/LB
            ("GENED", 1, 106, "A", "G", "Likely benign", "Pathogenic"),  # discordant: B/LB -> P/LP
        ],
    )

    by_gene, total_discordant, total_controls = compute_gene_discordance_stats(pd.ExcelFile(controls_path))

    assert total_discordant == 6
    assert total_controls == 7
    # GENEA/GENEB tie at 2 (alphabetical), then GENEC/GENED tie at 1 (alphabetical).
    assert by_gene.index.tolist() == ["GENEA", "GENEB", "GENEC", "GENED"]
    assert by_gene["total"].tolist() == [2, 2, 1, 1]
    assert by_gene[DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL].tolist() == [2, 0, 1, 0]
    assert by_gene[DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL].tolist() == [0, 2, 0, 1]


def test_compute_gene_discordance_stats_dedups_by_dna_variant(tmp_path):
    controls_path = tmp_path / "controls.xlsx"
    _write_gene_discordance_sheet(
        controls_path,
        [
            # Same genomic coordinates -- only the first (discordant) row should survive dedup.
            ("GENEX", 1, 100, "A", "G", "Pathogenic", "Benign"),
            ("GENEX", 1, 100, "A", "G", "Pathogenic", "Pathogenic"),
        ],
    )

    by_gene, total_discordant, total_controls = compute_gene_discordance_stats(pd.ExcelFile(controls_path))

    assert total_controls == 1
    assert total_discordant == 1
    assert by_gene.loc["GENEX", "total"] == 1


def test_format_gene_discordance_summary():
    by_gene = pd.DataFrame(
        {
            "total": [3, 2, 1],
            DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL: [3, 0, 1],
            DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL: [0, 2, 0],
        },
        index=pd.Index(["GENEA", "GENEB", "GENEC"], name="Gene"),
    )

    text = format_gene_discordance_summary(by_gene, total_discordant=6, total_controls=100, top_n=2)

    assert text.startswith(GENE_DISCORDANCE_TITLE)
    assert "Total discordant control variants: 6 of 100" in text
    assert "Top 2 genes by discordant-variant count:" in text
    assert "GENEA" in text
    assert "GENEB" in text
    assert "GENEC" not in text  # truncated by top_n


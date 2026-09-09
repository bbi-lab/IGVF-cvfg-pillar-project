import matplotlib.colors as mcolors
import pandas as pd

# _gain_legend_handles/SYNERGY_CATEGORY_COLORS aren't part of the module's
# public surface, but are imported directly for a regression test (see
# test_gain_legend_handles_have_the_correct_fixed_colors) against a real bug:
# an empty gain series' auto-generated legend proxy silently fell back to
# the wrong color.
from src.ablation_variant_reclassification import (
    BOTH_ALONE_LABEL,
    CLASS_DOWNGRADED,
    CLASS_UNCHANGED,
    CLASS_UPGRADED,
    CONCORDANCE_STATUS_COLORS,
    CONCORDANT,
    CONFLICT_LABEL,
    CONTROL_SCOPES,
    DIRECTION_BENIGN,
    DIRECTION_PATHOGENIC,
    DISCORDANT,
    FUNCTIONAL_ALONE_LABEL,
    FUNCTIONAL_ARM_LABEL,
    PREDICTOR_ALONE_LABEL,
    PREDICTOR_POINTS_COLUMNS,
    SCOPE_ORDER,
    SYNERGY_CATEGORY_COLORS,
    SYNERGY_LABEL,
    UNRESOLVED,
    _gain_legend_handles,
    add_ablation_points_columns,
    build_ablation_report,
    build_arm,
    build_control_concordance_chart_data,
    build_direction_comparison_chart_data,
    build_document_chart_data,
    build_predictor_arms,
    build_redundant_gain_chart_data,
    categorize_synergy,
    class_tier,
    class_upgrade_counts,
    class_upgrade_status,
    clingen_control_missense_only_truth_direction,
    clingen_control_truth_direction,
    clinvar_control_missense_only_truth_direction,
    clinvar_control_truth_direction,
    combined_points_col,
    concordance_counts,
    concordance_status,
    control_truth_direction_series,
    describe_scope,
    is_clingen_control,
    is_clingen_control_missense_only,
    is_clinvar_control,
    is_clinvar_control_missense_only,
    is_gnomad,
    is_missense_only,
    is_missense_or_start_loss,
    is_unobserved,
    is_vus,
    mixed_year_star_series,
    predictor_control_concordance_counts,
    predictor_flags_for_direction,
    redundant_gain_counts,
    resolved_flags,
    resolved_mask,
    restrict_to_consequence,
    save_ablation_document,
    save_control_concordance_chart,
    save_document_charts_as_files,
    save_redundant_gain_chart,
    save_synergy_chart,
    save_synergy_chart_comparison,
    scope_mask,
    summarize_arm,
    synergy_counts,
)
from src.build_variant_reclassification_dataset import OUTPUT_COLUMNS

CHECKPOINT_COLS = OUTPUT_COLUMNS[:-6] + [
    "VariantNotes",
    "ExC_points_2025",
    "ExC_points_2018",
    "OP_points",
    "Fxn_points",
    "Points_REVEL_GeneSpecific_GenomeWide",
    "Points_AM_GeneSpecific_GenomeWide",
    "Points_MP2_GeneSpecific_GenomeWide",
    "Conflicting_REVEL_GeneSpecific",
    "clnsig_group_18_25",
    "clinvar_conflict_flag_18_25",
]


def _checkpoint_row(**overrides):
    row = {col: None for col in CHECKPOINT_COLS}
    row.update(
        {
            "mavedb_variant_urn": "urn:mavedb:1",
            "Dataset": "G1_Dataset",
            "Gene": "G1",
            "Chrom": "1",
            "hg38_start": 1000,
            "ref_allele": "A",
            "alt_allele": "T",
            "hgvs_p": "p.NoMatch1Xxx",
            "auth_reported_score": 99.0,
            "Flag": None,
            "VariantNotes": None,
            "splice_var_amino": "No",
            "ExC_points_2025": 5,
            "ExC_points_2018": None,
            "OP_points": None,
            "Fxn_points": 5,
            "Points_REVEL_GeneSpecific_GenomeWide": 2,
            "Points_AM_GeneSpecific_GenomeWide": None,
            "Points_MP2_GeneSpecific_GenomeWide": None,
            "clinvar_sig_2025": None,
            "gnomad_MAF": None,
        }
    )
    row.update(overrides)
    return row


def _checkpoint_frame(rows):
    return pd.DataFrame([_checkpoint_row(**r) for r in rows])


# --- add_ablation_points_columns / combined_points_col --------------------------------------


def test_revel_combined_column_is_aliased_to_existing_combined_points():
    df = _checkpoint_frame([{"Fxn_points": 5, "Points_REVEL_GeneSpecific_GenomeWide": 2}])
    out = add_ablation_points_columns(df)
    assert out[combined_points_col("REVEL")].iloc[0] == 7
    assert out["Combined_points"].iloc[0] == 7


def test_am_and_mp2_combined_columns_sum_functional_and_predictor_points():
    df = _checkpoint_frame(
        [
            {
                "Fxn_points": 5,
                "Points_AM_GeneSpecific_GenomeWide": 3,
                "Points_MP2_GeneSpecific_GenomeWide": None,
            }
        ]
    )
    out = add_ablation_points_columns(df)
    assert out[combined_points_col("AlphaMissense")].iloc[0] == 8
    assert out[combined_points_col("MutPred2")].iloc[0] == 5  # missing predictor points treated as 0


# --- build_arm: each arm dedups by its own points column --------------------------------------


def test_functional_only_arm_can_pick_a_different_winner_than_combined():
    """Row A has weaker functional but stronger REVEL evidence (wins on
    Combined_points); row B has stronger functional but weaker REVEL evidence
    (should win the functional-only arm even though it loses combined)."""
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "Fxn_points": 2, "Points_REVEL_GeneSpecific_GenomeWide": 6},
            {"Dataset": "B", "Fxn_points": 8, "Points_REVEL_GeneSpecific_GenomeWide": 0},
        ]
    )
    out = add_ablation_points_columns(df)

    functional_arm = build_arm(out, "Functional_points")
    combined_arm = build_arm(out, combined_points_col("REVEL"))

    assert functional_arm["Dataset"].iloc[0] == "B"
    assert combined_arm["Dataset"].iloc[0] == "A"


def test_predictor_only_arm_dedups_by_predictor_points_missing_treated_as_zero():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "Fxn_points": 8, "Points_REVEL_GeneSpecific_GenomeWide": None},
            {"Dataset": "B", "Fxn_points": 0, "Points_REVEL_GeneSpecific_GenomeWide": 3},
        ]
    )
    out = add_ablation_points_columns(df)
    predictor_arm = build_arm(out, PREDICTOR_POINTS_COLUMNS["REVEL"])
    assert predictor_arm["Dataset"].iloc[0] == "B"


# --- summarize_arm / is_unobserved -------------------------------------------------------------


def test_summarize_arm_counts_pathogenic_or_benign_and_vus_resolved():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "hg38_start": 1, "Fxn_points": 8, "clinvar_sig_2025": "Uncertain significance"},
            {"Dataset": "B", "hg38_start": 2, "Fxn_points": 0, "clinvar_sig_2025": "Uncertain significance"},
        ]
    )
    out = add_ablation_points_columns(df)
    arm = build_arm(out, "Functional_points")
    stats = summarize_arm(arm, "Functional_points")

    assert stats["total"] == 2
    assert stats["pathogenic_or_benign"] == 1
    assert stats["vus_total"] == 2
    assert stats["vus_resolved"] == 1


def test_is_unobserved_requires_snv_and_no_clinvar_or_gnomad():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "ref_allele": "A", "alt_allele": "T", "clinvar_sig_2025": None, "gnomad_MAF": None},
            {"Dataset": "B", "ref_allele": "AT", "alt_allele": "T", "clinvar_sig_2025": None, "gnomad_MAF": None},
            {"Dataset": "C", "ref_allele": "A", "alt_allele": "T", "clinvar_sig_2025": "Benign", "gnomad_MAF": None},
        ]
    )
    out = is_unobserved(df)
    assert list(out) == [True, False, False]


# --- cross-arm agreement -----------------------------------------------------------------------


def test_resolved_flags_align_by_genomic_key_across_arms():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "hg38_start": 1, "Fxn_points": 8, "Points_REVEL_GeneSpecific_GenomeWide": 0},
        ]
    )
    out = add_ablation_points_columns(df)
    functional_arm = build_arm(out, "Functional_points")
    combined_arm = build_arm(out, combined_points_col("REVEL"))

    functional_flags = resolved_flags(functional_arm, "Functional_points")
    combined_flags = resolved_flags(combined_arm, combined_points_col("REVEL"))

    assert functional_flags.index.equals(combined_flags.index)
    assert functional_flags.iloc[0] and combined_flags.iloc[0]


def test_is_vus_matches_clinvar_uncertain_significance_only():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "clinvar_sig_2025": "Uncertain significance"},
            {"Dataset": "B", "clinvar_sig_2025": "Pathogenic"},
            {"Dataset": "C", "clinvar_sig_2025": None},
        ]
    )
    assert list(is_vus(df)) == [True, False, False]


def test_is_gnomad_requires_gnomad_observed_regardless_of_clinvar_status():
    """`is_gnomad` reproduces the notebook's bare `gnomad_MAF.notna()` test
    (cell 104) -- a ClinVar call (row C) doesn't exclude a variant, unlike an
    earlier, incorrect version of this function that additionally required
    `clinvar_sig_2025.isna()`."""
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "clinvar_sig_2025": None, "gnomad_MAF": 0.001},
            {"Dataset": "B", "clinvar_sig_2025": None, "gnomad_MAF": None},
            {"Dataset": "C", "clinvar_sig_2025": "Benign", "gnomad_MAF": 0.001},
        ]
    )
    assert list(is_gnomad(df)) == [True, False, True]


def test_scope_mask_unions_requested_categories():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "hg38_start": 1, "clinvar_sig_2025": "Uncertain significance"},
            {"Dataset": "B", "hg38_start": 2, "clinvar_sig_2025": None, "gnomad_MAF": 0.001},
            {"Dataset": "C", "hg38_start": 3, "clinvar_sig_2025": None, "gnomad_MAF": None},
            {"Dataset": "D", "hg38_start": 4, "clinvar_sig_2025": "Pathogenic"},
        ]
    )
    out = add_ablation_points_columns(df)
    arm = build_arm(out, "Functional_points")

    assert sorted(scope_mask(arm, ("vus", "unobserved")).values) == [False, False, True, True]
    assert sorted(scope_mask(arm, ("vus", "gnomad", "unobserved")).values) == [False, True, True, True]
    assert sorted(scope_mask(arm, ("gnomad",)).values) == [False, False, False, True]


def test_describe_scope_formats_one_two_or_three_categories():
    assert describe_scope(("vus",)) == "ClinVar VUS"
    assert describe_scope(("unobserved", "vus")) == "ClinVar VUS or unobserved"
    assert describe_scope(("unobserved", "gnomad", "vus")) == "ClinVar VUS, gnomAD population, or unobserved"


# --- ClinVar/ClinGen control scopes -------------------------------------------------------------


def test_mixed_year_star_series_uses_2018_star_for_mixed_year_genes_only():
    df = _checkpoint_frame(
        [
            {"Gene": "BRCA1", "clinvar_star_2018": "reviewed by expert panel", "clinvar_star_2025": "no assertion"},
            {"Gene": "G1", "clinvar_star_2018": "reviewed by expert panel", "clinvar_star_2025": "no assertion"},
        ]
    )
    out = mixed_year_star_series(df)
    assert list(out) == ["reviewed by expert panel", "no assertion"]


def test_is_clinvar_control_requires_control_significance_no_conflict_and_one_plus_star():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
            {
                "Dataset": "B",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": "has clinvar conflict",
                "clinvar_star_2025": "reviewed by expert panel",
            },
            {
                "Dataset": "C",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "no classification for the single variant",
            },
            {
                "Dataset": "D",
                "clnsig_group_18_25": "Uncertain significance",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
        ]
    )
    assert list(is_clinvar_control(df)) == [True, False, False, False]


def test_is_clingen_control_requires_non_null_non_vus_classification():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "Updated_Classification_ClinGen_repo": "Pathogenic"},
            {"Dataset": "B", "Updated_Classification_ClinGen_repo": "VUS"},
            {"Dataset": "C", "Updated_Classification_ClinGen_repo": None},
        ]
    )
    assert list(is_clingen_control(df)) == [True, False, False]


def test_is_missense_or_start_loss_matches_only_the_two_predictor_scored_consequences():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "condensed_consequence": "missense_variant"},
            {"Dataset": "B", "condensed_consequence": "start_lost"},
            {"Dataset": "C", "condensed_consequence": "stop_gained"},
            {"Dataset": "D", "condensed_consequence": None},
        ]
    )
    assert list(is_missense_or_start_loss(df)) == [True, True, False, False]


def test_is_missense_only_excludes_start_loss_unlike_is_missense_or_start_loss():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "condensed_consequence": "missense_variant"},
            {"Dataset": "B", "condensed_consequence": "start_lost"},
            {"Dataset": "C", "condensed_consequence": "stop_gained"},
            {"Dataset": "D", "condensed_consequence": None},
        ]
    )
    assert list(is_missense_only(df)) == [True, False, False, False]


def test_restrict_to_consequence_missense_only_drops_start_loss_rows():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "condensed_consequence": "missense_variant"},
            {"Dataset": "B", "condensed_consequence": "start_lost"},
        ]
    )
    restricted = restrict_to_consequence(df, "missense_only")
    assert list(restricted["Dataset"]) == ["A"]


def test_restrict_to_consequence_all_is_a_no_op():
    df = _checkpoint_frame([{"Dataset": "A", "condensed_consequence": "missense_variant"}])
    restricted = restrict_to_consequence(df, "all")
    assert len(restricted) == len(df)


def test_is_clinvar_control_missense_only_requires_control_and_strictly_missense():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "missense_variant",
            },
            {
                # Same control status as A, but start-loss -- excluded, unlike the
                # broader is_missense_or_start_loss (--consequence missense) test.
                "Dataset": "B",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "start_lost",
            },
            {
                # Missense, but not a control variant at all.
                "Dataset": "C",
                "clnsig_group_18_25": "Uncertain significance",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "missense_variant",
            },
        ]
    )
    assert list(is_clinvar_control_missense_only(df)) == [True, False, False]


def test_is_clingen_control_missense_only_requires_control_and_strictly_missense():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "Updated_Classification_ClinGen_repo": "Pathogenic",
                "condensed_consequence": "missense_variant",
            },
            {
                # Same control status as A, but start-loss -- excluded.
                "Dataset": "B",
                "Updated_Classification_ClinGen_repo": "Pathogenic",
                "condensed_consequence": "start_lost",
            },
            {
                "Dataset": "C",
                "Updated_Classification_ClinGen_repo": "VUS",
                "condensed_consequence": "missense_variant",
            },
        ]
    )
    assert list(is_clingen_control_missense_only(df)) == [True, False, False]


def test_scope_mask_includes_clinvar_and_clingen_controls():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Updated_Classification_ClinGen_repo": "Likely Pathogenic",
            },
            {"Dataset": "C", "hg38_start": 3},
        ]
    )
    out = add_ablation_points_columns(df)
    arm = build_arm(out, "Functional_points")

    assert sorted(scope_mask(arm, ("clinvar_control",)).values) == [False, False, True]
    assert sorted(scope_mask(arm, ("clingen_control",)).values) == [False, False, True]
    assert sorted(scope_mask(arm, ("clinvar_control", "clingen_control")).values) == [False, True, True]


def test_clinvar_control_truth_direction_returns_pathogenic_benign_or_na():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "clnsig_group_18_25": "Benign",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
            {"Dataset": "C", "hg38_start": 3},  # not a control variant at all
        ]
    )
    direction = clinvar_control_truth_direction(df)
    assert direction.iloc[0] == DIRECTION_PATHOGENIC
    assert direction.iloc[1] == DIRECTION_BENIGN
    assert pd.isna(direction.iloc[2])


def test_clingen_control_truth_direction_returns_pathogenic_benign_or_na():
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "hg38_start": 1, "Updated_Classification_ClinGen_repo": "Pathogenic"},
            {"Dataset": "B", "hg38_start": 2, "Updated_Classification_ClinGen_repo": "Benign"},
            {"Dataset": "C", "hg38_start": 3, "Updated_Classification_ClinGen_repo": "VUS"},
        ]
    )
    direction = clingen_control_truth_direction(df)
    assert direction.iloc[0] == DIRECTION_PATHOGENIC
    assert direction.iloc[1] == DIRECTION_BENIGN
    assert pd.isna(direction.iloc[2])


def test_clinvar_control_missense_only_truth_direction_drops_non_missense_only():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "missense_variant",
            },
            {
                # A real ClinVar control, but a consequence no predictor scores.
                "Dataset": "B",
                "hg38_start": 2,
                "clnsig_group_18_25": "Benign",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "synonymous_variant",
            },
        ]
    )
    direction = clinvar_control_missense_only_truth_direction(df)
    assert direction.iloc[0] == DIRECTION_PATHOGENIC
    assert pd.isna(direction.iloc[1])


def test_clingen_control_missense_only_truth_direction_drops_non_missense_only():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Updated_Classification_ClinGen_repo": "Pathogenic",
                "condensed_consequence": "missense_variant",
            },
            {
                # A real ClinGen control, but start-loss -- excluded.
                "Dataset": "B",
                "hg38_start": 2,
                "Updated_Classification_ClinGen_repo": "Benign",
                "condensed_consequence": "start_lost",
            },
        ]
    )
    direction = clingen_control_missense_only_truth_direction(df)
    assert direction.iloc[0] == DIRECTION_PATHOGENIC
    assert pd.isna(direction.iloc[1])


def test_control_truth_direction_series_clinvar_wins_ties_over_clingen():
    """A variant satisfying both control scopes with conflicting truth: the
    two sources are independent in practice, but `SCOPE_ORDER` fixes a
    deterministic tiebreak (ClinVar before ClinGen)."""
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "Updated_Classification_ClinGen_repo": "Benign",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    arm = build_arm(out, "Functional_points")

    direction = control_truth_direction_series(arm, CONTROL_SCOPES)
    assert list(direction) == [DIRECTION_PATHOGENIC]


def test_control_truth_direction_series_all_na_without_a_control_scope():
    df = _checkpoint_frame([{"Dataset": "A", "hg38_start": 1, "clinvar_sig_2025": "Uncertain significance"}])
    out = add_ablation_points_columns(df)
    arm = build_arm(out, "Functional_points")

    direction = control_truth_direction_series(arm, ("vus",))
    assert direction.isna().all()


# --- direction (any/pathogenic/benign) ----------------------------------------------------------


def test_resolved_mask_pathogenic_direction_ignores_benign_points():
    points = pd.Series([8, -8, 3])
    assert list(resolved_mask(points, DIRECTION_PATHOGENIC)) == [True, False, False]


def test_resolved_mask_benign_direction_ignores_pathogenic_points():
    points = pd.Series([8, -8, 3])
    assert list(resolved_mask(points, DIRECTION_BENIGN)) == [False, True, False]


def test_direction_pathogenic_conflict_bucket_is_pathogenic_before_conflict():
    """Functional alone reaches Likely Pathogenic (8); REVEL alone doesn't
    reach either threshold (-3) but is enough, combined with functional, to
    drag the total (8 + -3 = 5) back under the Likely Pathogenic cutoff. In
    a --direction pathogenic run, this variant should land in the conflict
    bucket -- it was Likely Pathogenic by functional evidence alone, before
    combining undid that -- even though REVEL's own points never themselves
    crossed the pathogenic-direction threshold."""
    df = _checkpoint_frame(
        [{"Dataset": "A", "hg38_start": 1, "Fxn_points": 8, "Points_REVEL_GeneSpecific_GenomeWide": -3}]
    )
    out = add_ablation_points_columns(df)
    report, chart_data = build_ablation_report(out, ["REVEL"], direction=DIRECTION_PATHOGENIC)

    restricted_flags = chart_data["predictor_flags"]["REVEL"]
    category = categorize_synergy(restricted_flags, FUNCTIONAL_ARM_LABEL, "REVEL only", "Functional + REVEL")
    assert list(category) == [CONFLICT_LABEL]
    assert "Pathogenic/Likely Pathogenic agreement" in report


# --- build_ablation_report (end to end, text sanity checks) -------------------------------------


def test_build_ablation_report_includes_a_section_per_requested_predictor():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 2,
                "Points_AM_GeneSpecific_GenomeWide": 2,
            },
        ]
    )
    out = add_ablation_points_columns(df)
    report, chart_data = build_ablation_report(out, ["REVEL", "AlphaMissense"])
    assert "=== Ablation: REVEL ===" in report
    assert "=== Ablation: AlphaMissense ===" in report
    assert "=== Ablation: MutPred2 ===" not in report
    assert "Functional evidence only" in report
    assert set(chart_data["predictor_flags"]) == {"REVEL", "AlphaMissense"}


# --- categorize_synergy / synergy_counts / save_synergy_chart -----------------------------------


def _flags(functional, predictor, combined):
    return pd.DataFrame({FUNCTIONAL_ARM_LABEL: functional, "REVEL only": predictor, "Functional + REVEL": combined})


def test_categorize_synergy_labels_each_combination():
    flags = _flags(
        functional=[True, False, True, False, True, False, True],
        predictor=[False, True, True, False, False, True, True],
        combined=[True, True, True, True, False, False, False],
    )
    category = categorize_synergy(flags, FUNCTIONAL_ARM_LABEL, "REVEL only", "Functional + REVEL")
    assert list(category) == [
        FUNCTIONAL_ALONE_LABEL,
        PREDICTOR_ALONE_LABEL,
        BOTH_ALONE_LABEL,
        SYNERGY_LABEL,
        CONFLICT_LABEL,
        CONFLICT_LABEL,
        CONFLICT_LABEL,
    ]


def test_categorize_synergy_drops_variants_resolved_by_nothing():
    flags = _flags(functional=[False], predictor=[False], combined=[False])
    category = categorize_synergy(flags, FUNCTIONAL_ARM_LABEL, "REVEL only", "Functional + REVEL")
    assert category.isna().all()


def test_synergy_counts_reindexes_missing_categories_to_zero():
    flags = _flags(functional=[True], predictor=[False], combined=[True])
    counts = synergy_counts(flags, FUNCTIONAL_ARM_LABEL, "REVEL only", "Functional + REVEL")
    assert counts[FUNCTIONAL_ALONE_LABEL] == 1
    assert counts[PREDICTOR_ALONE_LABEL] == 0
    assert counts[CONFLICT_LABEL] == 0


def test_save_synergy_chart_writes_a_file(tmp_path):
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 2,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": 0,
                "Points_REVEL_GeneSpecific_GenomeWide": 8,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    _, chart_data = build_ablation_report(out, ["REVEL"])

    output_path = tmp_path / "chart.png"
    save_synergy_chart(chart_data, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# --- direction comparison (side-by-side pathogenic/benign chart) --------------------------------


def test_predictor_flags_for_direction_reuses_the_same_arms_across_directions():
    """The same pre-built arms, read under two different directions, should
    give the same answer as building fresh arms per direction -- confirming
    arm dedup (abs(points)-based) doesn't depend on direction."""
    df = _checkpoint_frame(
        [
            {"Dataset": "A", "hg38_start": 1, "Fxn_points": 8, "Points_REVEL_GeneSpecific_GenomeWide": 0},
            {"Dataset": "B", "hg38_start": 2, "Fxn_points": -8, "Points_REVEL_GeneSpecific_GenomeWide": 0},
        ]
    )
    out = add_ablation_points_columns(df)
    functional_arm = build_arm(out, "Functional_points")
    functional_flags_pathogenic = resolved_flags(functional_arm, "Functional_points", DIRECTION_PATHOGENIC)
    functional_flags_benign = resolved_flags(functional_arm, "Functional_points", DIRECTION_BENIGN)

    predictor_arm, combined_arm = build_predictor_arms("REVEL", out)

    flags_pathogenic = predictor_flags_for_direction(
        "REVEL", predictor_arm, combined_arm, functional_flags_pathogenic, DIRECTION_PATHOGENIC
    )
    flags_benign = predictor_flags_for_direction(
        "REVEL", predictor_arm, combined_arm, functional_flags_benign, DIRECTION_BENIGN
    )

    assert flags_pathogenic[FUNCTIONAL_ARM_LABEL].sum() == 1
    assert flags_benign[FUNCTIONAL_ARM_LABEL].sum() == 1
    assert not flags_pathogenic[FUNCTIONAL_ARM_LABEL].equals(flags_benign[FUNCTIONAL_ARM_LABEL])


def test_build_direction_comparison_chart_data_covers_requested_directions():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 2,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": -8,
                "Points_REVEL_GeneSpecific_GenomeWide": -2,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    chart_data_by_direction = build_direction_comparison_chart_data(out, ["REVEL"])

    assert set(chart_data_by_direction) == {DIRECTION_PATHOGENIC, DIRECTION_BENIGN}
    pathogenic = chart_data_by_direction[DIRECTION_PATHOGENIC]
    benign = chart_data_by_direction[DIRECTION_BENIGN]
    assert pathogenic["direction_label"] == "Pathogenic/Likely Pathogenic"
    assert benign["direction_label"] == "Benign/Likely Benign"
    # Same scope population either way -- direction only changes what counts as resolved.
    assert pathogenic["scope_total"] == benign["scope_total"] == 2


def test_save_synergy_chart_comparison_writes_a_file(tmp_path):
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 2,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": -8,
                "Points_REVEL_GeneSpecific_GenomeWide": -2,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    chart_data_by_direction = build_direction_comparison_chart_data(out, ["REVEL"])

    output_path = tmp_path / "comparison.png"
    save_synergy_chart_comparison(chart_data_by_direction, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# --- redundant-gain histogram (points gained within the "each alone sufficient" band) -----------


def test_redundant_gain_counts_computes_net_change_in_evidence_strength():
    functional_points = pd.Series({1: 8, 2: 5})
    predictor_points = pd.Series({1: 6, 2: 5})
    combined_points = pd.Series({1: 14, 2: 10})
    both_alone_index = pd.Index([1])

    result = redundant_gain_counts(functional_points, predictor_points, combined_points, both_alone_index)

    assert result["total"] == 1
    assert dict(result["gained_via_predictor"]) == {6: 1}  # |14| - |8| = 6, i.e. REVEL's own points
    assert dict(result["gained_via_functional"]) == {8: 1}  # |14| - |6| = 8, i.e. functional's own points


def test_redundant_gain_counts_handles_opposite_sign_erosion():
    """Both individually resolved but in opposite directions (only possible
    under --direction any); the combined score sits closer to zero than the
    stronger individual source. Gain should reflect that erosion, not just
    restate the other source's own magnitude."""
    functional_points = pd.Series({1: 8})
    predictor_points = pd.Series({1: -3})
    combined_points = pd.Series({1: 5})
    both_alone_index = pd.Index([1])

    result = redundant_gain_counts(functional_points, predictor_points, combined_points, both_alone_index)

    assert dict(result["gained_via_predictor"]) == {-3: 1}  # |5| - |8| = -3 (a loss), not |-3| = 3
    assert dict(result["gained_via_functional"]) == {2: 1}  # |5| - |-3| = 2


def test_build_redundant_gain_chart_data_only_includes_both_alone_variants():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 6,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 0,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    gain_data = build_redundant_gain_chart_data(out, ["REVEL"], direction=DIRECTION_PATHOGENIC)

    result = gain_data["REVEL"]
    assert result["total"] == 1
    assert dict(result["gained_via_predictor"]) == {6: 1}
    assert dict(result["gained_via_functional"]) == {8: 1}


def test_save_redundant_gain_chart_writes_a_file(tmp_path):
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 6,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    gain_data = build_redundant_gain_chart_data(out, ["REVEL"], direction=DIRECTION_PATHOGENIC)

    output_path = tmp_path / "gain.png"
    save_redundant_gain_chart(gain_data, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# --- --show-class-upgrade (lighter/darker split by ACMG evidence-class upgrade) -----------------


def test_class_tier_boundaries():
    points = pd.Series([12, 10, 9, 6, 5, 0, -1, -6, -7, -10])
    assert list(class_tier(points)) == [2, 2, 1, 1, 0, 0, 1, 1, 2, 2]


def test_class_upgrade_status_functional_alone_compares_against_functional_tier():
    functional_points = pd.Series({1: 8, 2: 9})  # both tier 1 (Likely Pathogenic)
    predictor_points = pd.Series({1: 2, 2: 2})  # not the baseline here -- should be ignored
    combined_points = pd.Series({1: 11, 2: 9})  # variant 1 -> tier 2 (upgrade); variant 2 stays tier 1
    category = pd.Series({1: FUNCTIONAL_ALONE_LABEL, 2: FUNCTIONAL_ALONE_LABEL})

    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    assert list(status) == [CLASS_UPGRADED, CLASS_UNCHANGED]


def test_class_upgrade_status_predictor_alone_compares_against_predictor_tier():
    functional_points = pd.Series({1: 2})  # not the baseline here -- should be ignored
    predictor_points = pd.Series({1: 8})  # tier 1
    combined_points = pd.Series({1: 11})  # tier 2 -> upgrade
    category = pd.Series({1: PREDICTOR_ALONE_LABEL})

    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    assert list(status) == [CLASS_UPGRADED]


def test_class_upgrade_status_both_alone_compares_against_the_higher_of_the_two():
    functional_points = pd.Series({1: 6})  # tier 1
    predictor_points = pd.Series({1: 9})  # tier 1 -- higher of the two tiers is still 1 (tie)
    combined_points = pd.Series({1: 15})  # tier 2 -> upgrade
    category = pd.Series({1: BOTH_ALONE_LABEL})

    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    assert list(status) == [CLASS_UPGRADED]


def test_class_upgrade_status_downgraded_when_opposing_evidence_erodes_but_stays_resolved():
    """Functional alone is tier 2 (Pathogenic); a weak opposing predictor
    score erodes the combined total to tier 1 (Likely Pathogenic) without
    dropping it out of --direction-resolved range entirely (still >= 6, so
    the category stays FUNCTIONAL_ALONE_LABEL rather than CONFLICT_LABEL --
    this is the erosion case, not a sign conflict)."""
    functional_points = pd.Series({1: 11})  # tier 2
    predictor_points = pd.Series({1: -4})  # not independently resolved as pathogenic
    combined_points = pd.Series({1: 7})  # tier 1 -- weaker than the functional-alone baseline
    category = pd.Series({1: FUNCTIONAL_ALONE_LABEL})

    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    assert list(status) == [CLASS_DOWNGRADED]


def test_class_upgrade_status_na_outside_the_three_alone_categories():
    functional_points = pd.Series({1: 8})
    predictor_points = pd.Series({1: 0})
    combined_points = pd.Series({1: 8})
    category = pd.Series({1: SYNERGY_LABEL})

    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    assert status.isna().all()


def test_class_upgrade_counts_splits_upgraded_unchanged_and_downgraded():
    functional_points = pd.Series({1: 8, 2: 9, 3: 11})
    predictor_points = pd.Series({1: 0, 2: 0, 3: -4})
    combined_points = pd.Series({1: 11, 2: 9, 3: 7})  # upgrades, unchanged, downgraded respectively
    category = pd.Series({1: FUNCTIONAL_ALONE_LABEL, 2: FUNCTIONAL_ALONE_LABEL, 3: FUNCTIONAL_ALONE_LABEL})

    counts = class_upgrade_counts(functional_points, predictor_points, combined_points, category)
    assert counts[FUNCTIONAL_ALONE_LABEL] == {"upgraded": 1, "unchanged": 1, "downgraded": 1}
    assert counts[PREDICTOR_ALONE_LABEL] == {"upgraded": 0, "unchanged": 0, "downgraded": 0}
    assert counts[BOTH_ALONE_LABEL] == {"upgraded": 0, "unchanged": 0, "downgraded": 0}


def test_build_ablation_report_show_class_upgrade_is_opt_in():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)

    _, chart_data_without = build_ablation_report(out, ["REVEL"], direction=DIRECTION_PATHOGENIC)
    assert "class_upgrade" not in chart_data_without

    _, chart_data_with = build_ablation_report(out, ["REVEL"], direction=DIRECTION_PATHOGENIC, show_class_upgrade=True)
    assert chart_data_with["class_upgrade"]["REVEL"][FUNCTIONAL_ALONE_LABEL] == {
        "upgraded": 1,
        "unchanged": 0,
        "downgraded": 0,
    }


def test_save_synergy_chart_with_class_upgrade_writes_a_file(tmp_path):
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    _, chart_data = build_ablation_report(out, ["REVEL"], direction=DIRECTION_PATHOGENIC, show_class_upgrade=True)

    output_path = tmp_path / "chart_class_upgrade.png"
    save_synergy_chart(chart_data, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_save_synergy_chart_comparison_with_class_upgrade_writes_a_file(tmp_path):
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": -8,
                "Points_REVEL_GeneSpecific_GenomeWide": -3,
                "clinvar_sig_2025": "Uncertain significance",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    chart_data_by_direction = build_direction_comparison_chart_data(out, ["REVEL"], show_class_upgrade=True)

    output_path = tmp_path / "comparison_class_upgrade.png"
    save_synergy_chart_comparison(chart_data_by_direction, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# --- control concordance (does each arm's call agree with the known ClinVar/ClinGen truth) ------


def test_concordance_status_concordant_discordant_unresolved():
    points = pd.Series({1: 8, 2: 8, 3: 2})
    truth_direction = pd.Series({1: DIRECTION_PATHOGENIC, 2: DIRECTION_BENIGN, 3: DIRECTION_PATHOGENIC})

    status = concordance_status(points, truth_direction)
    assert list(status) == [CONCORDANT, DISCORDANT, UNRESOLVED]


def test_concordance_status_na_where_truth_direction_is_na():
    points = pd.Series({1: 8})
    truth_direction = pd.Series({1: pd.NA}, dtype="object")

    status = concordance_status(points, truth_direction)
    assert status.isna().all()


def test_concordance_counts_splits_by_arm_and_status():
    functional_points = pd.Series({1: 8, 2: -8, 3: 2})  # concordant, discordant, unresolved (truth all pathogenic)
    predictor_points = pd.Series({1: 8, 2: 8, 3: 2})  # concordant, concordant, unresolved
    combined_points = pd.Series({1: 8, 2: -8, 3: 2})
    truth_direction = pd.Series({1: DIRECTION_PATHOGENIC, 2: DIRECTION_PATHOGENIC, 3: DIRECTION_PATHOGENIC})

    counts = concordance_counts(functional_points, predictor_points, combined_points, truth_direction)
    assert counts["functional"] == {"concordant": 1, "discordant": 1, "unresolved": 1}
    assert counts["predictor"] == {"concordant": 2, "discordant": 0, "unresolved": 1}
    assert counts["combined"] == {"concordant": 1, "discordant": 1, "unresolved": 1}


def test_predictor_control_concordance_counts_uses_predictor_specific_arms():
    df = _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
            },
        ]
    )
    out = add_ablation_points_columns(df)
    functional_arm = build_arm(out, "Functional_points")
    predictor_arm, combined_arm = build_predictor_arms("REVEL", out)
    truth_direction = control_truth_direction_series(functional_arm, ("clinvar_control",))

    counts = predictor_control_concordance_counts("REVEL", functional_arm, predictor_arm, combined_arm, truth_direction)
    assert counts["functional"] == {"concordant": 1, "discordant": 0, "unresolved": 0}  # Fxn_points=8, resolved
    assert counts["predictor"] == {"concordant": 0, "discordant": 0, "unresolved": 1}  # REVEL points=3, unresolved
    assert counts["combined"] == {"concordant": 1, "discordant": 0, "unresolved": 0}  # combined=11, resolved


def test_concordance_status_colors_use_the_reserved_dataviz_status_palette():
    """`--show-class-upgrade` already uses the categorical palette's
    light/dark shades on the "alone sufficient" bars, so control concordance
    -- a correctness dimension, not a category identity -- deliberately
    reuses the dataviz skill's reserved good/critical status colors instead,
    to avoid colliding with that existing meaning."""
    assert CONCORDANCE_STATUS_COLORS[CONCORDANT] == "#0ca30c"
    assert CONCORDANCE_STATUS_COLORS[DISCORDANT] == "#d03b3b"


# --- --document (one section per scope, whichever chart types requested) ------------------------


def _document_checkpoint_frame():
    """One variant in each of the seven SCOPE_ORDER categories -- functional
    alone already resolves each (Fxn_points=8, points_REVEL below threshold),
    so every section's "ablation"/"comparison" chart has a nonempty bar and
    the "gain" chart's per-predictor panel isn't entirely empty either. D and
    E are also given `condensed_consequence == "missense_variant"`, so each
    does double duty as the sole variant in its `_missense_only`-restricted
    scope too (`clinvar_control_missense_only`/`clingen_control_missense_
    only` are strict subsets of `clinvar_control`/`clingen_control`) -- no
    separate F/G rows needed. (Deliberately not `start_lost` for either: that
    consequence satisfies the broader `clinvar_control`/`clingen_control`
    scopes but not their `_missense_only` subsets, which require strictly
    `is_missense_only`.)"""
    return _checkpoint_frame(
        [
            {
                "Dataset": "A",
                "hg38_start": 1,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Uncertain significance",
            },
            {
                "Dataset": "B",
                "hg38_start": 2,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": None,
                "gnomad_MAF": 0.01,
            },
            {
                "Dataset": "C",
                "hg38_start": 3,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": None,
                "gnomad_MAF": None,
                "ref_allele": "A",
                "alt_allele": "T",
            },
            {
                "Dataset": "D",
                "hg38_start": 4,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Pathogenic",
                "clnsig_group_18_25": "Pathogenic",
                "clinvar_conflict_flag_18_25": None,
                "clinvar_star_2025": "reviewed by expert panel",
                "condensed_consequence": "missense_variant",
            },
            {
                "Dataset": "E",
                "hg38_start": 5,
                "Fxn_points": 8,
                "Points_REVEL_GeneSpecific_GenomeWide": 3,
                "clinvar_sig_2025": "Pathogenic",
                "Updated_Classification_ClinGen_repo": "Pathogenic",
                "condensed_consequence": "missense_variant",
            },
        ]
    )


def test_build_document_chart_data_has_one_section_per_scope_with_one_variant_each():
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"])

    assert set(document_data) == set(SCOPE_ORDER)
    for scope in SCOPE_ORDER:
        assert document_data[scope]["scope_total"] == 1


def test_build_document_chart_data_only_builds_requested_chart_types():
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("ablation",))

    for scope in SCOPE_ORDER:
        assert document_data[scope]["ablation"] is not None
        assert document_data[scope]["comparison"] is None
        assert document_data[scope]["gain"] is None


def test_build_document_chart_data_comparison_covers_both_directions():
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("comparison",))

    comparison = document_data["vus"]["comparison"]
    assert set(comparison) == {DIRECTION_PATHOGENIC, DIRECTION_BENIGN}
    assert document_data["vus"]["ablation"] is None
    assert document_data["vus"]["gain"] is None


def test_save_ablation_document_writes_a_file_with_all_chart_types(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"])

    output_path = tmp_path / "document.png"
    save_ablation_document(document_data, output_path, predictors=["REVEL"])
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_save_ablation_document_infers_predictors_when_omitted(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("gain",))

    output_path = tmp_path / "document_gain_only.png"
    save_ablation_document(document_data, output_path, chart_types=("gain",))
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_build_control_concordance_chart_data_only_counts_control_variants():
    """Of `_document_checkpoint_frame`'s five rows, only D (clinvar_control,
    clinvar_control_missense_only) and E (clingen_control,
    clingen_control_missense_only) carry a known truth direction --
    both Pathogenic, both with functional/combined arms resolving pathogenic
    (Fxn_points=8, combined=8+3=11) and REVEL alone unresolved (points=3).
    `build_control_concordance_chart_data`'s default `scopes` (all four
    `CONTROL_SCOPES`) doesn't double-count D/E: `control_truth_direction_
    series` fills each variant's truth direction once, from the first
    matching scope in `SCOPE_ORDER`, so a variant satisfying both a control
    scope and its `_missense_only` subset still contributes exactly one
    truth-direction row."""
    out = add_ablation_points_columns(_document_checkpoint_frame())
    concordance_data = build_control_concordance_chart_data(out, ["REVEL"])

    counts = concordance_data["REVEL"]
    assert counts["functional"] == {"concordant": 2, "discordant": 0, "unresolved": 0}
    assert counts["predictor"] == {"concordant": 0, "discordant": 0, "unresolved": 2}
    assert counts["combined"] == {"concordant": 2, "discordant": 0, "unresolved": 0}


def test_build_control_concordance_chart_data_scopes_restricts_to_one_control_source():
    out = add_ablation_points_columns(_document_checkpoint_frame())
    concordance_data = build_control_concordance_chart_data(out, ["REVEL"], scopes=("clinvar_control",))

    assert concordance_data["REVEL"]["functional"] == {"concordant": 1, "discordant": 0, "unresolved": 0}


def test_save_control_concordance_chart_writes_a_file(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    concordance_data = build_control_concordance_chart_data(out, ["REVEL"])

    output_path = tmp_path / "concordance_chart.png"
    save_control_concordance_chart(concordance_data, output_path)
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_build_document_chart_data_concordance_only_nonzero_for_control_scopes():
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("concordance",))

    for scope in SCOPE_ORDER:
        assert document_data[scope]["concordance"] is not None
        assert document_data[scope]["ablation"] is None

    assert document_data["vus"]["concordance"]["REVEL"]["functional"] == {
        "concordant": 0,
        "discordant": 0,
        "unresolved": 0,
    }
    assert document_data["clinvar_control"]["concordance"]["REVEL"]["functional"] == {
        "concordant": 1,
        "discordant": 0,
        "unresolved": 0,
    }
    assert document_data["clinvar_control_missense_only"]["concordance"]["REVEL"]["functional"] == {
        "concordant": 1,
        "discordant": 0,
        "unresolved": 0,
    }
    assert document_data["clingen_control"]["concordance"]["REVEL"]["functional"] == {
        "concordant": 1,
        "discordant": 0,
        "unresolved": 0,
    }
    assert document_data["clingen_control_missense_only"]["concordance"]["REVEL"]["functional"] == {
        "concordant": 1,
        "discordant": 0,
        "unresolved": 0,
    }


def test_save_ablation_document_with_concordance_chart_writes_a_file(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("concordance",))

    output_path = tmp_path / "document_concordance_only.png"
    save_ablation_document(document_data, output_path, chart_types=("concordance",))
    assert output_path.exists()
    assert output_path.stat().st_size > 0


# --- save_document_charts_as_files -----------------------------------------------------------


def test_save_document_charts_as_files_writes_one_file_per_block(tmp_path):
    """Of `_document_checkpoint_frame`'s seven scopes, all seven get an
    'ablation'/'comparison'/'gain' file, but only the four control scopes
    (clinvar_control/clinvar_control_missense_only/clingen_control/
    clingen_control_missense_only) carry a known truth direction and so
    get a 'concordance' file -- vus/gnomad/unobserved are skipped rather than
    writing an all-zero, uninformative chart.
    """
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"])

    output_dir = tmp_path / "split"
    written = save_document_charts_as_files(document_data, output_dir)

    assert len(written) == 7 + 7 + 4 + 7  # ablation + comparison + concordance + gain
    for path in written:
        assert path.exists()
        assert path.stat().st_size > 0

    assert (output_dir / "ablation_vus.pdf").exists()
    assert (output_dir / "concordance_clinvar_control.pdf").exists()
    assert not (output_dir / "concordance_vus.pdf").exists()
    assert not (output_dir / "concordance_gnomad.pdf").exists()
    assert not (output_dir / "concordance_unobserved.pdf").exists()


def test_save_document_charts_as_files_respects_requested_chart_types(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("gain",))

    output_dir = tmp_path / "split_gain_only"
    written = save_document_charts_as_files(document_data, output_dir, chart_types=("gain",))

    assert len(written) == len(SCOPE_ORDER)
    assert all(path.name.startswith("gain_") for path in written)


def test_save_document_charts_as_files_honors_fmt(tmp_path):
    out = add_ablation_points_columns(_document_checkpoint_frame())
    document_data = build_document_chart_data(out, ["REVEL"], chart_types=("ablation",))

    output_dir = tmp_path / "split_png"
    written = save_document_charts_as_files(document_data, output_dir, chart_types=("ablation",), fmt="png")

    assert written
    assert all(path.suffix == ".png" for path in written)


def test_gain_legend_handles_have_the_correct_fixed_colors():
    """Regression test for a real bug: `ax.bar([], [], color=...)` (an empty
    gain series, e.g. a scope/predictor with zero "both alone sufficient"
    variants) produces a `BarContainer` with no patches, and matplotlib's
    auto-generated legend proxy for it silently falls back to a default
    color instead of the one passed to `color=` -- both series would render
    with the same wrong color in the legend. `_gain_legend_handles` sidesteps
    this with explicit, fixed-color `Patch` handles instead of
    `ax.get_legend_handles_labels()`.
    """
    handles = _gain_legend_handles()
    colors_by_label = {handle.get_label(): mcolors.to_hex(handle.get_facecolor()) for handle in handles}

    assert colors_by_label["Gained via functional evidence"] == SYNERGY_CATEGORY_COLORS[FUNCTIONAL_ALONE_LABEL]
    assert colors_by_label["Gained via predictor evidence"] == SYNERGY_CATEGORY_COLORS[PREDICTOR_ALONE_LABEL]
    assert colors_by_label["Gained via functional evidence"] != colors_by_label["Gained via predictor evidence"]

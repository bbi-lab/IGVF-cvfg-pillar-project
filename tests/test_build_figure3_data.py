import gzip

import numpy as np
import pandas as pd
from click.testing import CliRunner

from src.build_figure3_data import (
    IGVF_CATEGORY_ALL,
    IGVF_CATEGORY_MIXED,
    IGVF_CATEGORY_OTHER,
    add_clinvar_snapshot_column,
    build_figure3a_gene_summary,
    build_figure3c_assay_categories,
    build_figure3c_variant_assay_combos,
    build_figure3d_counts,
    collapse_to_unique_variants,
    main,
    seen_in_gnomad,
    simplify_significance,
    summarize_clinvar_significance,
)

TESTING_REGISTRY_COLUMNS = [
    "accession_version",
    "test_type",
    "object",
    "GTR_identifier",
    "MIM_number",
    "object_name",
    "gene_or_SNOMED_CT_ID",
    "gene_symbol",
]


def test_summarize_clinvar_significance_flags_conflicts():
    assert summarize_clinvar_significance(pd.Series(["Pathogenic", "Benign"])) == "has clinvar conflict"
    assert (
        summarize_clinvar_significance(pd.Series(["Conflicting classifications of pathogenicity"]))
        == "has clinvar conflict"
    )


def test_summarize_clinvar_significance_prefers_likely_when_mixed_with_strict():
    assert summarize_clinvar_significance(pd.Series(["Pathogenic", "Likely pathogenic"])) == "Likely pathogenic"
    assert summarize_clinvar_significance(pd.Series(["Benign", "Benign/Likely benign"])) == "Likely benign"


def test_summarize_clinvar_significance_handles_vus_and_empty():
    assert summarize_clinvar_significance(pd.Series(["Uncertain significance"])) == "Uncertain significance"
    assert summarize_clinvar_significance(pd.Series([np.nan, np.nan])) == "Unseen"
    assert summarize_clinvar_significance(pd.Series(["Uncertain significance", "risk factor"])) == "VUS/conflict"


def test_seen_in_gnomad():
    assert seen_in_gnomad(pd.Series([np.nan, 0.001])) == "Seen"
    assert seen_in_gnomad(pd.Series([np.nan, np.nan])) == "Unseen"


def test_add_clinvar_snapshot_column_uses_2018_only_for_priority_genes():
    pp = pd.DataFrame(
        {
            "Gene": ["BRCA1", "PALB2"],
            "clinvar_sig_2018": ["Pathogenic", "Benign"],
            "clinvar_sig_2025": ["Likely pathogenic", "Likely benign"],
        }
    )

    result = add_clinvar_snapshot_column(pp)

    assert list(result["clinvar_18_25"]) == ["Pathogenic", "Likely benign"]


def _sample_integrated_dataset() -> pd.DataFrame:
    # None of these genes are in PRIORITY_GENES, so clinvar_18_25 always
    # follows clinvar_sig_2025 here -- the priority-gene 2018-override path is
    # covered separately by test_add_clinvar_snapshot_column_uses_2018_only_for_priority_genes.
    return pd.DataFrame(
        {
            "Gene": ["RAD51C", "RAD51C", "RAD51C", "PALB2", "PALB2"],
            "Dataset": ["TP53_A", "TP53_A", "TP53_B", "PALB2_IGVF", "PALB2_IGVF"],
            "clinvar_sig_2018": ["Pathogenic"] * 5,
            "clinvar_sig_2025": [
                "Pathogenic",
                "Pathogenic",
                "Benign",
                "Uncertain significance",
                np.nan,
            ],
            "nucleotide_or_aa": ["nt", "nt", "nt", "nt", "nt"],
            "hg38_start": [100, 100, 200, 300, 400],
            "ref_allele": ["A", "A", "C", "G", "T"],
            "alt_allele": ["T", "T", "G", "A", "C"],
            "gnomad_MAF": [np.nan, np.nan, 0.01, np.nan, np.nan],
            "aa_pos": [np.nan] * 5,
            "aa_ref": [np.nan] * 5,
            "aa_alt": [np.nan] * 5,
            "RefSeq Transcript ID": [np.nan] * 5,
        }
    )


def test_collapse_to_unique_variants_dedupes_repeated_genomic_rows():
    pp = add_clinvar_snapshot_column(_sample_integrated_dataset())

    pp_unique = collapse_to_unique_variants(pp)

    # The two identical TP53_A rows (same Gene/hg38_start/ref/alt) collapse to one.
    assert len(pp_unique) == 4
    # The nucleotide branch copies clinvar_18_25 through as-is (no
    # "Unseen"/"VUS" normalization -- that only happens for aa-level rows),
    # so a missing 2025 classification stays NaN.
    assert pp_unique["clnsig_group_18_25"].isna().sum() == 1
    assert set(pp_unique["clnsig_group_18_25"].dropna()) == {
        "Pathogenic",
        "Benign",
        "Uncertain significance",
    }


def test_collapse_to_unique_variants_dataset_scoped_counts_once_per_dataset():
    # The same RAD51C nucleotide variant (hg38_start=100, A>T) is tested by
    # both RAD51C_A and RAD51C_B; the same PALB2 protein variant
    # (p.Arg10Gly on NM_1) is tested by both PALB2_A and PALB2_B.
    pp = add_clinvar_snapshot_column(
        pd.DataFrame(
            {
                "Gene": ["RAD51C", "RAD51C", "PALB2", "PALB2"],
                "Dataset": ["RAD51C_A", "RAD51C_B", "PALB2_A", "PALB2_B"],
                "clinvar_sig_2018": ["Pathogenic"] * 4,
                "clinvar_sig_2025": ["Pathogenic", "Pathogenic", "Pathogenic", "Pathogenic"],
                "nucleotide_or_aa": ["nt", "nt", "aa", "aa"],
                "hg38_start": [100, 100, np.nan, np.nan],
                "ref_allele": ["A", "A", np.nan, np.nan],
                "alt_allele": ["T", "T", np.nan, np.nan],
                "gnomad_MAF": [np.nan, np.nan, np.nan, np.nan],
                "aa_pos": [np.nan, np.nan, 10, 10],
                "aa_ref": [np.nan, np.nan, "Arg", "Arg"],
                "aa_alt": [np.nan, np.nan, "Gly", "Gly"],
                "RefSeq Transcript ID": [np.nan, np.nan, "NM_1.2", "NM_1.2"],
            }
        )
    )

    global_unique = collapse_to_unique_variants(pp)
    dataset_scoped_unique = collapse_to_unique_variants(pp, dataset_scoped=True)

    assert len(global_unique) == 2
    assert len(dataset_scoped_unique) == 4
    assert set(dataset_scoped_unique["Dataset"]) == {"RAD51C_A", "RAD51C_B", "PALB2_A", "PALB2_B"}


def test_build_figure3d_counts_includes_controls_gnomad_and_vus():
    pp = add_clinvar_snapshot_column(_sample_integrated_dataset())
    pp_unique = collapse_to_unique_variants(pp)

    figure3d = build_figure3d_counts(pp, pp_unique)
    counts = figure3d.set_index("Group")["Count"].to_dict()

    assert counts["Pathogenic"] == 1
    assert counts["Benign"] == 1
    assert counts["VUS total"] == 1
    assert counts["gnomAD"] == 1


def test_simplify_significance_groups_and_passes_through_vus():
    assert simplify_significance("Pathogenic/Likely pathogenic") == "Likely pathogenic"
    assert simplify_significance("Benign/Likely benign") == "Likely benign"
    assert simplify_significance("Uncertain significance") == "Uncertain significance"
    assert simplify_significance("VUS/conflict") == "Uncertain significance"
    assert pd.isna(simplify_significance("has clinvar conflict"))
    assert pd.isna(simplify_significance(np.nan))


def test_build_figure3a_gene_summary_joins_gencc_uniprot_and_testing_registry():
    pp = add_clinvar_snapshot_column(_sample_integrated_dataset())
    pp_unique = collapse_to_unique_variants(pp)

    curation = pd.DataFrame({"Gene": ["RAD51C", "PALB2"], "IGVF Produced?": ["No", "Yes"]})
    gencc = pd.DataFrame(
        {
            "gene_symbol": ["RAD51C", "RAD51C", "PALB2", "OTHERGENE"],
            "classification_title": ["Definitive", "Limited", "Strong", "Definitive"],
            "disease_curie": ["MONDO:1", "MONDO:2", "MONDO:3", "MONDO:4"],
        }
    )
    uniprot = pd.DataFrame(
        {
            "Gene Names (primary)": ["RAD51C", "PALB2", "OTHERGENE"],
            "Length": [376, 1186, 100],
        }
    )
    testing_registry = pd.DataFrame(
        {
            "accession_version": ["GTR1", "GTR2", "GTR3", "GTR4"],
            "test_type": ["Clinical", "Clinical", "Research", "Clinical"],
            "object": ["gene", "condition", "gene", "gene"],
            "GTR_identifier": ["C1", "C2", "C3", "C4"],
            "MIM_number": [1, 2, 3, 4],
            "object_name": ["a", "b", "c", "d"],
            "gene_or_SNOMED_CT_ID": ["1", "2", "3", "4"],
            "gene_symbol": ["RAD51C", "RAD51C", "PALB2", "PALB2"],
        }
    )

    figure3a = build_figure3a_gene_summary(pp_unique, curation, gencc, uniprot, testing_registry)

    # RAD51C matches only its Definitive GenCC row -- the Limited-evidence
    # row is filtered out.
    rad51c_rows = figure3a.loc[figure3a["Gene"] == "RAD51C"]
    assert len(rad51c_rows) == 1
    assert rad51c_rows["possible_SNVs"].iloc[0] == 376 * 9
    # Only GTR1 is Clinical+gene for RAD51C (GTR2 is a condition row).
    assert rad51c_rows["gene_test_count"].iloc[0] == 1
    assert rad51c_rows["IGVF_produced"].iloc[0] == IGVF_CATEGORY_OTHER

    palb2_rows = figure3a.loc[figure3a["Gene"] == "PALB2"]
    assert palb2_rows["possible_SNVs"].iloc[0] == 1186 * 9
    # Only GTR4 is Clinical+gene for PALB2 (GTR3 is Research, excluded).
    assert palb2_rows["gene_test_count"].iloc[0] == 1
    assert palb2_rows["IGVF_produced"].iloc[0] == IGVF_CATEGORY_ALL

    assert "Gene Names (primary)" not in figure3a.columns
    assert "GeneSymbol" not in figure3a.columns


def test_build_figure3a_gene_summary_flags_mixed_igvf_status():
    # BRCA2 has one IGVF-produced and one community dataset -- a gene with
    # *some* but not *all* IGVF-produced datasets should be flagged as
    # mixed, not folded into either extreme.
    pp_unique = pd.DataFrame(
        {
            "Gene": ["BRCA1", "BRCA2", "BRCA2"],
            "clnsig_group_18_25": ["Pathogenic", "Pathogenic", "Benign"],
        }
    )
    curation = pd.DataFrame(
        {
            "Gene": ["BRCA1", "BRCA2", "BRCA2"],
            "IGVF Produced?": ["Yes", "Yes", "No"],
        }
    )
    empty_gencc = pd.DataFrame(columns=["gene_symbol", "classification_title", "disease_curie"])
    empty_uniprot = pd.DataFrame(columns=["Gene Names (primary)", "Length"])
    empty_registry = pd.DataFrame(columns=TESTING_REGISTRY_COLUMNS)

    figure3a = build_figure3a_gene_summary(pp_unique, curation, empty_gencc, empty_uniprot, empty_registry)
    igvf_status = figure3a.set_index("Gene")["IGVF_produced"]

    assert igvf_status["BRCA1"] == IGVF_CATEGORY_ALL
    assert igvf_status["BRCA2"] == IGVF_CATEGORY_MIXED


def test_build_figure3a_gene_summary_collapses_combined_gene_label():
    # The integrated dataset's Gene column stores CALM1/CALM2/CALM3's shared
    # dataset as one comma-joined value; it should collapse onto CALM1 alone.
    pp_unique = pd.DataFrame(
        {
            "Gene": ["CALM1, CALM2, CALM3", "CALM1, CALM2, CALM3"],
            "clnsig_group_18_25": ["Pathogenic", "Benign"],
        }
    )
    curation = pd.DataFrame({"Gene": [], "IGVF Produced?": []})
    empty_gencc = pd.DataFrame(columns=["gene_symbol", "classification_title", "disease_curie"])
    empty_uniprot = pd.DataFrame(columns=["Gene Names (primary)", "Length"])
    empty_registry = pd.DataFrame(columns=TESTING_REGISTRY_COLUMNS)

    figure3a = build_figure3a_gene_summary(pp_unique, curation, empty_gencc, empty_uniprot, empty_registry)

    assert list(figure3a["Gene"]) == ["CALM1"]
    assert figure3a.loc[0, "Pathogenic"] == 1
    assert figure3a.loc[0, "Benign"] == 1


def test_build_figure3c_assay_categories_excludes_meta_analysis_and_tags_flags():
    pp = add_clinvar_snapshot_column(_sample_integrated_dataset())
    pp_unique = collapse_to_unique_variants(pp)
    curation = pd.DataFrame(
        {
            "Dataset Name": ["TP53_A", "TP53_B", "PALB2_IGVF"],
            "Assay Name": ["SGE", "Other", "SGE"],
            "IGVF Produced?": ["No", "No", "Yes"],
            "Primary Score Set or Meta-analysis?": [
                "primary score set",
                "meta-analysis",
                "primary score set",
            ],
        }
    )

    figure3c = build_figure3c_assay_categories(pp_unique, curation).set_index("Dataset")

    assert "TP53_B" not in figure3c.index
    assert figure3c.loc["TP53_A", "SGE"] == "Yes"
    assert figure3c.loc["PALB2_IGVF", "IGVF"] == "Yes"
    assert figure3c.loc["PALB2_IGVF", "n_unique_IDs"] == 2


def test_build_figure3c_variant_assay_combos_flags_multi_assay_variants():
    # The (hg38_start=500, A>G) variant is tested by both an "Other" dataset
    # and an SGE dataset -- it should get its own "Other+SGE" combo rather
    # than being folded into just one. The (600, C>T) variant is only tested
    # by the SGE dataset. The meta-analysis dataset (excluded) also tests the
    # first variant, but shouldn't affect its combo.
    pp = pd.DataFrame(
        {
            "Gene": ["LDLRX", "LDLRX", "LDLRX", "LDLRX"],
            "Dataset": ["LDLRX_Other", "LDLRX_SGE", "LDLRX_SGE", "LDLRX_Meta"],
            "nucleotide_or_aa": ["nt", "nt", "nt", "nt"],
            "hg38_start": [500, 500, 600, 500],
            "ref_allele": ["A", "A", "C", "A"],
            "alt_allele": ["G", "G", "T", "G"],
            "aa_pos": [np.nan] * 4,
            "aa_ref": [np.nan] * 4,
            "aa_alt": [np.nan] * 4,
            "RefSeq Transcript ID": [np.nan] * 4,
        }
    )
    curation = pd.DataFrame(
        {
            "Dataset Name": ["LDLRX_Other", "LDLRX_SGE", "LDLRX_Meta"],
            "Assay Name": ["Other", "SGE", "SGE"],
            "Primary Score Set or Meta-analysis?": ["primary score set", "primary score set", "meta-analysis"],
        }
    )

    combos = build_figure3c_variant_assay_combos(pp, curation).set_index(["Gene", "assay_combo"])["n_unique_IDs"]

    assert combos[("LDLRX", "Other+SGE")] == 1
    assert combos[("LDLRX", "SGE")] == 1
    assert ("LDLRX", "Other") not in combos.index


def _brca1x_nt_aa_overlap_dataset() -> pd.DataFrame:
    # Three distinct nucleotide variants (nt1/nt2/nt3, tested by an SGE
    # dataset) all translate to the same protein change -- p.His100Phe on
    # NM_1 -- that a separate amino-acid-level variant (aa1, tested by an
    # "Other" dataset) also represents. A fourth, unrelated amino-acid
    # variant (aa2, p.Arg5Gly) has no nucleotide-level counterpart at all.
    # Matching real integrated-dataset rows, amino-acid-level rows still
    # carry their own (back-mapped) genomic coordinates -- distinct per
    # protein variant here, so the final genomic-coordinate collapse in
    # build_figure3c_variant_assay_combos doesn't fold aa1/aa2 into each
    # other.
    return pd.DataFrame(
        {
            "Gene": ["BRCA1X"] * 5,
            "Dataset": ["BRCA1X_SGE", "BRCA1X_SGE", "BRCA1X_SGE", "BRCA1X_Other", "BRCA1X_Other"],
            "nucleotide_or_aa": ["nt", "nt", "nt", "aa", "aa"],
            "hg38_start": [100, 200, 300, 1000, 2000],
            "ref_allele": ["A", "C", "G", "A", "C"],
            "alt_allele": ["T", "G", "C", "G", "A"],
            "aa_pos": [100, 100, 100, 100, 5],
            "aa_ref": ["H", "H", "H", "H", "R"],
            "aa_alt": ["F", "F", "F", "F", "G"],
            "RefSeq Transcript ID": ["NM_1.2", "NM_1.2", "NM_1.2", "NM_1.3", "NM_1.3"],
        }
    )


def _brca1x_curation() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Dataset Name": ["BRCA1X_SGE", "BRCA1X_Other"],
            "Assay Name": ["SGE", "Other"],
            "Primary Score Set or Meta-analysis?": ["primary score set", "primary score set"],
        }
    )


def test_build_figure3c_variant_assay_combos_projects_aa_onto_one_representative_nt_variant_by_default():
    combos = build_figure3c_variant_assay_combos(_brca1x_nt_aa_overlap_dataset(), _brca1x_curation()).set_index(
        ["Gene", "assay_combo"]
    )["n_unique_IDs"]

    # aa1 (p.His100Phe) is dropped entirely, but -- since it's one
    # measurement of one variant -- its "Other" category is folded onto only
    # one representative of the three nucleotide-level variants that
    # translate to the same change (nt1, the lowest hg38_start), not all
    # three: nt1 becomes "Other+SGE", nt2/nt3 stay plain "SGE". aa2
    # (p.Arg5Gly) has no nucleotide-level counterpart, so it's untouched and
    # still counted on its own under "Other". No project_aa_onto_nt argument
    # is passed here -- this is the default behavior.
    assert combos[("BRCA1X", "Other+SGE")] == 1
    assert combos[("BRCA1X", "SGE")] == 2
    assert combos[("BRCA1X", "Other")] == 1


def test_build_figure3c_variant_assay_combos_no_project_aa_onto_nt_keeps_levels_separate():
    combos = build_figure3c_variant_assay_combos(
        _brca1x_nt_aa_overlap_dataset(), _brca1x_curation(), project_aa_onto_nt=False
    ).set_index(["Gene", "assay_combo"])["n_unique_IDs"]

    # With projection explicitly disabled, nt1/nt2/nt3 keep their own
    # SGE-only combo; aa1 and aa2 are each counted separately under "Other"
    # -- the shared protein change between aa1 and nt1/nt2/nt3 has no effect.
    assert combos[("BRCA1X", "SGE")] == 3
    assert combos[("BRCA1X", "Other")] == 2


def test_main_cli_writes_all_four_outputs(tmp_path):
    integrated_dataset_path = tmp_path / "integrated.tsv.gz"
    with gzip.open(integrated_dataset_path, "wt") as f:
        _sample_integrated_dataset().to_csv(f, sep="\t", index=False)

    curation_sheet_path = tmp_path / "curation.xlsx"
    pd.DataFrame(
        {
            "Dataset Name": ["TP53_A", "TP53_B", "PALB2_IGVF"],
            "Gene": ["RAD51C", "RAD51C", "PALB2"],
            "Assay Name": ["SGE", "Other", "SGE"],
            "IGVF Produced?": ["No", "No", "Yes"],
            "Primary Score Set or Meta-analysis?": [
                "primary score set",
                "primary score set",
                "primary score set",
            ],
        }
    ).to_excel(curation_sheet_path, sheet_name="Curation", index=False)

    gencc_path = tmp_path / "gencc.csv"
    pd.DataFrame(
        {
            "gene_symbol": ["RAD51C", "PALB2"],
            "classification_title": ["Definitive", "Strong"],
            "disease_curie": ["MONDO:1", "MONDO:2"],
        }
    ).to_csv(gencc_path, index=False)

    uniprot_path = tmp_path / "uniprot.tsv.gz"
    pd.DataFrame({"Gene Names (primary)": ["RAD51C", "PALB2"], "Length": [376, 1186]}).to_csv(
        uniprot_path, sep="\t", index=False, compression="gzip"
    )

    testing_registry_path = tmp_path / "test_condition_gene.txt.gz"
    pd.DataFrame(
        {
            "accession_version": ["GTR1", "GTR2"],
            "test_type": ["Clinical", "Clinical"],
            "object": ["gene", "gene"],
            "GTR_identifier": ["C1", "C2"],
            "MIM_number": [1, 2],
            "object_name": ["a", "b"],
            "gene_or_SNOMED_CT_ID": ["1", "2"],
            "gene_symbol": ["RAD51C", "PALB2"],
        }
    ).to_csv(testing_registry_path, sep="\t", index=False, compression="gzip")

    figure3a_output_path = tmp_path / "Figure3a.csv.gz"
    figure3c_output_path = tmp_path / "Figure3c.csv.gz"
    figure3c_measurements_output_path = tmp_path / "Figure3c_measurements.csv.gz"
    figure3c_assay_combos_output_path = tmp_path / "Figure3c_assay_combos.csv.gz"
    figure3d_output_path = tmp_path / "Figure3d.csv.gz"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "--integrated-dataset",
            str(integrated_dataset_path),
            "--curation-sheet",
            str(curation_sheet_path),
            "--gencc",
            str(gencc_path),
            "--uniprot",
            str(uniprot_path),
            "--testing-registry",
            str(testing_registry_path),
            "--figure3a-output",
            str(figure3a_output_path),
            "--figure3c-output",
            str(figure3c_output_path),
            "--figure3c-measurements-output",
            str(figure3c_measurements_output_path),
            "--figure3c-assay-combos-output",
            str(figure3c_assay_combos_output_path),
            "--figure3d-output",
            str(figure3d_output_path),
        ],
    )

    assert result.exit_code == 0, result.output
    assert figure3a_output_path.exists()
    assert figure3c_output_path.exists()
    assert figure3c_measurements_output_path.exists()
    assert figure3c_assay_combos_output_path.exists()
    assert figure3d_output_path.exists()

    figure3a = pd.read_csv(figure3a_output_path)
    assert {"Gene", "possible_SNVs", "gene_test_count", "IGVF_produced"} <= set(figure3a.columns)

    expected_figure3c_columns = {
        "Gene",
        "Dataset",
        "SGE",
        "Vamp",
        "IGVF",
        "n_unique_IDs",
        "variant_dedup_scope",
    }
    figure3c = pd.read_csv(figure3c_output_path)
    assert set(figure3c.columns) == expected_figure3c_columns
    assert set(figure3c["variant_dedup_scope"]) == {"global"}

    figure3c_measurements = pd.read_csv(figure3c_measurements_output_path)
    assert set(figure3c_measurements.columns) == expected_figure3c_columns
    assert set(figure3c_measurements["variant_dedup_scope"]) == {"per-dataset"}

    figure3c_assay_combos = pd.read_csv(figure3c_assay_combos_output_path)
    assert set(figure3c_assay_combos.columns) == {"Gene", "assay_combo", "n_unique_IDs", "variant_projection"}
    # --project-aa-onto-nt defaults to on.
    assert set(figure3c_assay_combos["variant_projection"]) == {"aa_projected_onto_nt"}

    figure3d = pd.read_csv(figure3d_output_path)
    assert set(figure3d.columns) == {"Group", "Count", "Count in gnomAD"}

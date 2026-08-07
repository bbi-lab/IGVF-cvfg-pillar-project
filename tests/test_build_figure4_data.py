import pandas as pd
from click.testing import CliRunner

from src.build_figure4_data import (
    GENOME_WIDE_REVEL_THRESHOLDS,
    categorize_clinical_status,
    main,
    score_to_revel_tier,
    unique_gene_snvs_with_revel,
)


def test_score_to_revel_tier_assigns_strongest_matching_tier():
    scores = pd.Series([0.01, 0.20, 0.60, 0.90, 0.99])
    tiers = score_to_revel_tier(scores, GENOME_WIDE_REVEL_THRESHOLDS)

    assert list(tiers) == ["BP4_Strong", "BP4_Supporting", "IR", "PP3_Moderate+", "PP3_Strong"]


def test_score_to_revel_tier_handles_nan_threshold_tiers():
    # Very Strong tiers are NaN in the genome-wide table; scores must not
    # accidentally match a NaN threshold.
    scores = pd.Series([0.0, 1.0])
    tiers = score_to_revel_tier(scores, GENOME_WIDE_REVEL_THRESHOLDS)

    assert list(tiers) == ["BP4_Strong", "PP3_Strong"]


def test_unique_gene_snvs_with_revel_dedupes_by_genomic_position():
    df = pd.DataFrame(
        {
            "Gene": ["MSH2", "MSH2", "MSH2", "TP53"],
            "Chrom": ["2", "2", "2", "17"],
            "hg38_start": [100, 100, 200, 300],
            "ref_allele": ["A", "A", "C", "G"],
            "alt_allele": ["T", "T", "G", "A"],
            "REVEL": [0.5, 0.5, None, 0.8],
        }
    )

    result = unique_gene_snvs_with_revel(df, "MSH2")

    assert len(result) == 1
    assert result.iloc[0]["hg38_start"] == 100


def test_categorize_clinical_status_priority_order():
    gene_snvs = pd.DataFrame(
        {
            "clinvar_sig_2025": [
                "Pathogenic",
                "Benign",
                "Conflicting classifications of pathogenicity",
                "Uncertain significance",
                None,
                None,
            ],
            "gnomad_MAF": [None, None, None, None, 0.001, None],
        }
    )

    categories = categorize_clinical_status(gene_snvs)

    assert list(categories) == ["PLP", "BLB", "Conflicting", "VUS", "gnomAD", "allSNVs"]


def test_categorize_clinical_status_clinvar_outranks_gnomad():
    # A variant can be both ClinVar-annotated and gnomAD-observed;
    # clinical status should win over population frequency.
    gene_snvs = pd.DataFrame(
        {
            "clinvar_sig_2025": ["Pathogenic"],
            "gnomad_MAF": [0.0001],
        }
    )

    categories = categorize_clinical_status(gene_snvs)

    assert list(categories) == ["PLP"]


def test_main_cli_requires_cached_json(tmp_path):
    runner = CliRunner()
    integrated_dataset = tmp_path / "integrated.tsv.gz"
    integrated_dataset.write_text("placeholder")
    excalibr_dir = tmp_path / "excalibr"
    excalibr_dir.mkdir()
    supp_data_4 = tmp_path / "supp4.xlsx"
    supp_data_4.write_text("placeholder")

    result = runner.invoke(
        main,
        [
            "--integrated-dataset",
            str(integrated_dataset),
            "--excalibr-json-dir",
            str(excalibr_dir),
            "--supplementary-data-4",
            str(supp_data_4),
        ],
    )

    assert result.exit_code != 0
    assert "--cached-json" in result.output

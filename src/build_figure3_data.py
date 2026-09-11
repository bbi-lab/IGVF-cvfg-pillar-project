"""Rebuild data/intermediate/figures/figure_3/Figure3a.csv.gz, Figure3c.csv.gz,
and Figure3d.csv.gz from current pipeline outputs.

`curation_summary_figure3.Rmd` reads these three files as its data sources for
the bubble plot, bar/pie chart, and Euler/venn diagram panels. All three were
originally produced by an ad hoc, uncommitted notebook
(`data/output/Curation_summary_V5_cleaned.ipynb`) against a stale personal
directory layout; this script reproduces their logic against this repo's
actual current outputs, writing to the gitignored `data/intermediate/`
staging area (same convention as
`notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb`'s
`data/intermediate/figures/figure_2/` outputs) so they're regenerated on
demand rather than committed as stale snapshots.

Figure3a (bubble plot): per-gene ClinVar classification counts from the
integrated dataset, joined against three reference files under
`data/input/genes/`:
    - `gencc-submissions.csv.gz` -- GenCC gene-disease validity submissions
      (legacy UUID-based export from https://thegencc.org/download; this
      format is scheduled for removal 2026-09-30, after which re-fetching
      requires switching to GenCC's newer SGC-ID-based export and updating
      the columns this script reads), filtered to Definitive/Strong/Moderate
      classifications.
    - `uniprotkb_9606_reviewed.tsv.gz` -- reviewed (Swiss-Prot) human
      proteome from UniProt's REST API, for protein length -> possible_SNVs
      (9 possible single-nucleotide missense/nonsense changes per codon).
    - `test_condition_gene.txt.gz` -- NCBI Genetic Testing Registry's public
      bulk export (https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt),
      filtered to clinical, gene-level test records, for gene_test_count.

Figure3c (bar/pie panel): per-dataset unique-variant counts from the
integrated variant-effect dataset, tagged with SGE/Vamp-seq/IGVF flags joined
from `Supplementary_Data_3.xlsx`'s `Curation` sheet (`Assay Name` /
`IGVF Produced?` columns) rather than the original notebook's hardcoded
per-dataset sets, which had already drifted 3-4 datasets out of sync with the
curation sheet (confirmed: the curation sheet's `Assay Name`/`IGVF Produced?`/
`Primary Score Set or Meta-analysis?` columns reproduce the old hardcoded sets
exactly for every previously-known dataset). Meta-analysis datasets
(`Primary Score Set or Meta-analysis?` == "meta-analysis") are excluded, same
as the original notebook, to avoid double-counting variants already covered
by their primary score sets. Figure3a's gene-level IGVF flag is joined the
same way (curation sheet's `IGVF Produced?`, grouped up to one flag per gene)
rather than the original notebook's separate hardcoded 10-gene set --
confirmed both approaches produce the identical 10 genes.

Figure3d (Euler/venn panel): ClinVar control (Benign/Likely benign/
Pathogenic/Likely pathogenic and their combined labels), gnomAD,
unreported-SNV, and VUS counts, collapsed to one row per genomic variant (or
per protein-level variant when only amino-acid coordinates are available).
Priority genes BRCA1/PTEN/MSH2/TP53 use the dataset's 2018 ClinVar snapshot
(`clinvar_sig_2018`) rather than the current one, matching the original
notebook; every other gene uses `clinvar_sig_2025`.

Usage:
    python -m src.build_figure3_data [--integrated-dataset PATH] \\
        [--curation-sheet PATH] [--gencc PATH] [--uniprot PATH] \\
        [--testing-registry PATH] [--figure3a-output PATH] \\
        [--figure3c-output PATH] [--figure3d-output PATH]
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd

FIGURE_3_INTERMEDIATE_DIR = Path("data/intermediate/figures/figure_3")

DEFAULT_INTEGRATED_DATASET = Path("data/output/maves/integrated_variant_effect_dataset.tsv.gz")
DEFAULT_CURATION_SHEET = Path("data/input/maves/Supplementary_Data_3.xlsx")
DEFAULT_GENCC_PATH = Path("data/input/genes/gencc-submissions.csv.gz")
DEFAULT_UNIPROT_PATH = Path("data/input/genes/uniprotkb_9606_reviewed.tsv.gz")
DEFAULT_TESTING_REGISTRY_PATH = Path("data/input/genes/test_condition_gene.txt.gz")
DEFAULT_FIGURE3A_OUTPUT = FIGURE_3_INTERMEDIATE_DIR / "Figure3a.csv.gz"
DEFAULT_FIGURE3C_OUTPUT = FIGURE_3_INTERMEDIATE_DIR / "Figure3c.csv.gz"
DEFAULT_FIGURE3D_OUTPUT = FIGURE_3_INTERMEDIATE_DIR / "Figure3d.csv.gz"

# The integrated dataset's Gene column stores this dataset's three genes as
# one comma-joined value; Figure3a represents it as a single canonical gene
# (matching the previously-committed file) rather than fanning it out like
# the sunburst panel does.
GENE_NAME_OVERRIDES = {"CALM1, CALM2, CALM3": "CALM1"}

GENCC_VALIDITY_CLASSIFICATIONS = ["Definitive", "Strong", "Moderate"]

SIMPLIFIED_SIGNIFICANCE = {
    "pathogenic": "Pathogenic",
    "likely pathogenic": "Likely pathogenic",
    "pathogenic/likely pathogenic": "Likely pathogenic",
    "benign": "Benign",
    "likely benign": "Likely benign",
    "benign/likely benign": "Likely benign",
}

# NCBI's raw test_condition_gene.txt column names, renamed to match the
# columns Figure3a.csv.gz has always shipped under.
GTR_COLUMN_RENAMES = {
    "accession_version": "AccessionVersion",
    "test_type": "TestType",
    "object": "Object",
    "GTR_identifier": "GTR identifier",
    "MIM_number": "MIM number",
    "object_name": "ObjectName",
    "gene_or_SNOMED_CT_ID": "Gene/SNOMED CT ID",
    "gene_symbol": "GeneSymbol",
}

PRIORITY_GENES = {"BRCA1", "PTEN", "MSH2", "TP53"}

PATHOGENIC_SIGNIFICANCE = {
    "Pathogenic",
    "Likely pathogenic",
    "Pathogenic/Likely pathogenic",
}
BENIGN_SIGNIFICANCE = {
    "Benign",
    "Likely benign",
    "Benign/Likely benign",
}
CONFLICT_SIGNIFICANCE = {"Conflicting classifications of pathogenicity"}

CONTROL_GROUPS = [
    "Likely benign",
    "Benign",
    "Pathogenic",
    "Pathogenic/Likely pathogenic",
    "Likely pathogenic",
    "Benign/Likely benign",
]

NUCLEOTIDE_KEY = ["Gene", "hg38_start", "ref_allele", "alt_allele"]
PROTEIN_GROUP_COLUMNS = ["Gene", "aa_ref", "aa_pos", "aa_alt", "Ref_seq_transcript_ID_stripped"]


def summarize_clinvar_significance(values: pd.Series) -> str:
    significance = set(values.dropna())

    if not significance:
        return "Unseen"

    has_pathogenic = bool(significance & PATHOGENIC_SIGNIFICANCE)
    has_benign = bool(significance & BENIGN_SIGNIFICANCE)
    has_conflict = bool(significance & CONFLICT_SIGNIFICANCE)
    has_vus = "Uncertain significance" in significance

    if has_conflict or (has_pathogenic and has_benign):
        return "has clinvar conflict"

    if has_pathogenic:
        if significance & {"Likely pathogenic", "Pathogenic/Likely pathogenic"}:
            return "Likely pathogenic"
        return "Pathogenic"

    if has_benign:
        if significance & {"Likely benign", "Benign/Likely benign"}:
            return "Likely benign"
        return "Benign"

    if has_vus:
        return "Uncertain significance" if len(significance) == 1 else "VUS/conflict"

    if len(significance) == 1:
        return next(iter(significance))

    return "multiple other classifications"


def seen_in_gnomad(values: pd.Series) -> str:
    return "Seen" if values.notna().any() else "Unseen"


def add_clinvar_snapshot_column(pp: pd.DataFrame) -> pd.DataFrame:
    pp = pp.copy()
    pp["clinvar_18_25"] = np.where(
        pp["Gene"].isin(PRIORITY_GENES),
        pp["clinvar_sig_2018"],
        pp["clinvar_sig_2025"],
    )
    return pp


def collapse_to_unique_variants(pp: pd.DataFrame) -> pd.DataFrame:
    """Collapse repeated per-transcript/per-submission rows to one row per variant.

    Genomic (nucleotide-level) variants are deduped directly on their hg38
    coordinates. Protein-level (amino-acid) variants are grouped by gene/
    position/transcript first (since the same protein change can appear under
    several equivalent genomic representations) and summarized down to a
    single ClinVar/gnomAD status before being deduped.
    """
    variant_level = pp["nucleotide_or_aa"].replace({"nt": "nucleotide"})

    aa = pp.loc[variant_level.eq("aa")].copy()
    nucleotide = pp.loc[variant_level.eq("nucleotide")].copy()

    aa["Ref_seq_transcript_ID_stripped"] = (
        aa["RefSeq Transcript ID"].astype("string").str.replace(r"\.\d+$", "", regex=True)
    )
    aa["aa_pos"] = pd.to_numeric(aa["aa_pos"], errors="coerce")

    aa["clnsig_group_18_25"] = aa.groupby(PROTEIN_GROUP_COLUMNS, dropna=False)["clinvar_18_25"].transform(
        summarize_clinvar_significance
    )
    aa["gnomad_seen"] = aa.groupby(PROTEIN_GROUP_COLUMNS, dropna=False)["gnomad_MAF"].transform(seen_in_gnomad)

    nucleotide["clnsig_group_18_25"] = nucleotide["clinvar_18_25"]
    nucleotide["gnomad_seen"] = np.where(nucleotide["gnomad_MAF"].notna(), "Seen", "Unseen")

    nucleotide_unique = nucleotide.drop_duplicates(NUCLEOTIDE_KEY).copy()
    aa_unique = (
        aa.sort_values("gnomad_MAF", na_position="last").drop_duplicates(PROTEIN_GROUP_COLUMNS, keep="first").copy()
    )

    pp_unique = pd.concat([nucleotide_unique, aa_unique], ignore_index=True)
    # Amino-acid variants can still collide on genomic coordinates with each
    # other or with a nucleotide-level row; collapse once more on the final
    # coordinate-level key.
    return pp_unique.drop_duplicates(NUCLEOTIDE_KEY).copy()


def build_figure3d_counts(pp: pd.DataFrame, pp_unique: pd.DataFrame) -> pd.DataFrame:
    controls_df = pp_unique.loc[pp_unique["clnsig_group_18_25"].isin(CONTROL_GROUPS)].copy()

    final_counts = (
        controls_df.groupby("clnsig_group_18_25")
        .agg(
            Count=("clnsig_group_18_25", "size"),
            **{"Count in gnomAD": ("gnomad_MAF", lambda values: values.notna().sum())},
        )
        .reset_index()
        .rename(columns={"clnsig_group_18_25": "Group"})
    )

    pp_genomic_unique = pp.drop_duplicates(NUCLEOTIDE_KEY).copy()
    gnomad_unique_count = pp_genomic_unique["gnomad_MAF"].notna().sum()

    is_snv = pp_genomic_unique["ref_allele"].str.len().eq(1) & pp_genomic_unique["alt_allele"].str.len().eq(1)

    unreported_snvs = pp_genomic_unique.loc[
        is_snv & pp_genomic_unique["clinvar_sig_2025"].isna() & pp_genomic_unique["gnomad_MAF"].isna()
    ]
    vus = pp_genomic_unique.loc[pp_genomic_unique["clinvar_sig_2025"].eq("Uncertain significance")]

    additional_rows = pd.DataFrame(
        [
            {
                "Group": "gnomAD",
                "Count": gnomad_unique_count,
                "Count in gnomAD": gnomad_unique_count,
            },
            {
                "Group": "Unreported_SNVs",
                "Count": len(unreported_snvs),
                "Count in gnomAD": 0,
            },
            {
                "Group": "VUS total",
                "Count": len(vus),
                "Count in gnomAD": vus["gnomad_MAF"].notna().sum(),
            },
        ]
    )

    return pd.concat([final_counts, additional_rows], ignore_index=True)


def simplify_significance(value):
    if pd.isna(value):
        return pd.NA

    normalized = str(value).strip().lower()
    if normalized in SIMPLIFIED_SIGNIFICANCE:
        return SIMPLIFIED_SIGNIFICANCE[normalized]

    if "uncertain" in normalized or "vus" in normalized:
        return "Uncertain significance"

    return pd.NA


def build_figure3a_gene_summary(
    pp_unique: pd.DataFrame,
    curation: pd.DataFrame,
    gencc: pd.DataFrame,
    uniprot: pd.DataFrame,
    testing_registry: pd.DataFrame,
) -> pd.DataFrame:
    sig_group = pp_unique["clnsig_group_18_25"].map(simplify_significance)
    gene_counts = (
        pp_unique.assign(sig_group=sig_group).groupby("Gene")["sig_group"].value_counts().unstack(fill_value=0)
    )
    gene_counts.index = gene_counts.index.to_series().replace(GENE_NAME_OVERRIDES)
    # Collapse rows that now share a renamed Gene (e.g. a combined-gene label
    # folded onto a gene that also has its own separate rows elsewhere).
    gene_counts = gene_counts.groupby(level="Gene").sum().reset_index()
    gene_list = gene_counts["Gene"].dropna().unique().tolist()

    gencc_moderate_plus = (
        gencc.loc[
            gencc["gene_symbol"].isin(gene_list) & gencc["classification_title"].isin(GENCC_VALIDITY_CLASSIFICATIONS)
        ]
        .drop_duplicates(["gene_symbol", "disease_curie"])
        .copy()
    )
    gene_summary = gene_counts.merge(gencc_moderate_plus, left_on="Gene", right_on="gene_symbol", how="left")

    uniprot = (
        uniprot.loc[uniprot["Gene Names (primary)"].isin(gene_list)].drop_duplicates("Gene Names (primary)").copy()
    )
    uniprot["possible_SNVs"] = uniprot["Length"] * 9
    gene_summary = gene_summary.merge(uniprot, left_on="Gene", right_on="Gene Names (primary)", how="left")

    testing_registry = testing_registry.rename(columns=GTR_COLUMN_RENAMES)
    testing_registry_gene = testing_registry.loc[
        testing_registry["TestType"].eq("Clinical")
        & testing_registry["Object"].eq("gene")
        & testing_registry["GeneSymbol"].isin(gene_list)
    ]
    testing_registry_counts = testing_registry_gene.groupby("GeneSymbol").size().rename("gene_test_count").reset_index()
    gene_summary = gene_summary.merge(testing_registry_counts, left_on="Gene", right_on="GeneSymbol", how="left")

    igvf_genes = set(curation.loc[curation["IGVF Produced?"].astype("string").str.strip().eq("Yes"), "Gene"])
    gene_summary["IGVF_produced"] = np.where(gene_summary["Gene"].isin(igvf_genes), "Yes", "No")

    return gene_summary.drop(columns=["Gene Names (primary)", "GeneSymbol"], errors="ignore")


def build_figure3c_assay_categories(pp_unique: pd.DataFrame, curation: pd.DataFrame) -> pd.DataFrame:
    curation = curation.dropna(subset=["Dataset Name"]).copy()

    meta_analysis_datasets = set(
        curation.loc[
            curation["Primary Score Set or Meta-analysis?"].eq("meta-analysis"),
            "Dataset Name",
        ]
    )
    sge_datasets = set(curation.loc[curation["Assay Name"].eq("SGE"), "Dataset Name"])
    vamp_datasets = set(curation.loc[curation["Assay Name"].eq("Vamp-seq"), "Dataset Name"])
    igvf_datasets = set(curation.loc[curation["IGVF Produced?"].astype("string").str.strip().eq("Yes"), "Dataset Name"])

    source = pp_unique.loc[~pp_unique["Dataset"].isin(meta_analysis_datasets)].copy()
    source["SGE"] = np.where(source["Dataset"].isin(sge_datasets), "Yes", "No")
    source["Vamp"] = np.where(source["Dataset"].isin(vamp_datasets), "Yes", "No")
    source["IGVF"] = np.where(source["Dataset"].isin(igvf_datasets), "Yes", "No")

    return (
        source.groupby(["Gene", "Dataset", "SGE", "Vamp", "IGVF"], as_index=False)
        .size()
        .rename(columns={"size": "n_unique_IDs"})
    )


@click.command(help=__doc__)
@click.option(
    "--integrated-dataset",
    "integrated_dataset_path",
    default=DEFAULT_INTEGRATED_DATASET,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--curation-sheet",
    "curation_sheet_path",
    default=DEFAULT_CURATION_SHEET,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--gencc",
    "gencc_path",
    default=DEFAULT_GENCC_PATH,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--uniprot",
    "uniprot_path",
    default=DEFAULT_UNIPROT_PATH,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--testing-registry",
    "testing_registry_path",
    default=DEFAULT_TESTING_REGISTRY_PATH,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--figure3a-output",
    "figure3a_output_path",
    default=DEFAULT_FIGURE3A_OUTPUT,
    type=click.Path(path_type=Path),
)
@click.option(
    "--figure3c-output",
    "figure3c_output_path",
    default=DEFAULT_FIGURE3C_OUTPUT,
    type=click.Path(path_type=Path),
)
@click.option(
    "--figure3d-output",
    "figure3d_output_path",
    default=DEFAULT_FIGURE3D_OUTPUT,
    type=click.Path(path_type=Path),
)
def main(
    integrated_dataset_path,
    curation_sheet_path,
    gencc_path,
    uniprot_path,
    testing_registry_path,
    figure3a_output_path,
    figure3c_output_path,
    figure3d_output_path,
):
    try:
        pp = pd.read_csv(integrated_dataset_path, sep="\t", low_memory=False)
        pp = add_clinvar_snapshot_column(pp)
        pp_unique = collapse_to_unique_variants(pp)

        curation = pd.read_excel(curation_sheet_path, sheet_name="Curation")
        gencc = pd.read_csv(gencc_path)
        uniprot = pd.read_csv(uniprot_path, sep="\t")
        testing_registry = pd.read_csv(testing_registry_path, sep="\t")

        figure3a = build_figure3a_gene_summary(pp_unique, curation, gencc, uniprot, testing_registry)
        figure3c = build_figure3c_assay_categories(pp_unique, curation)
        figure3d = build_figure3d_counts(pp, pp_unique)
    except (ValueError, KeyError) as exc:
        raise click.ClickException(str(exc)) from exc

    outputs = {
        Path(figure3a_output_path): figure3a,
        Path(figure3c_output_path): figure3c,
        Path(figure3d_output_path): figure3d,
    }
    for output_path, data in outputs.items():
        output_path.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(output_path, index=False, compression="gzip")
        click.echo(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Summary statistics for the integrated MAVE variant effect dataset.

Reads the condensed variant effect dataset (one row per variant effect
measurement/score; see Data/mave_data/integrated_variant_effect_dataset*.tsv.gz)
together with its dataset-level metadata (Supplementary_Data_3.xlsx) and
reports, for three groupings of datasets -- IGVF-produced only, non-IGVF
("community") only, and combined -- the number of datasets, variant effect
measurements, composite scores, distinct variants assayed, and genes
represented.

A dataset counts as "IGVF-produced" if Supplementary_Data_3's Curation sheet
marks its `IGVF Produced?` column "Yes". A row counts as a direct measurement
if its dataset's `Primary Score Set or Meta-analysis?` column reads "primary
score set"; everything else (meta-analyses, trained predictors, etc.) counts
as a composite score. This is a per-dataset determination: every row of a
given dataset is classified the same way.

Note on `(hgvs_g, hgvs_p)`: the integrated dataset has no `hgvs_g` (genomic
HGVS) column -- its DNA-level identifier is `hgvs_c` (transcript-relative
HGVS, pipe-delimited when a protein-resolution measurement corresponds to more
than one underlying DNA change). This script uses `hgvs_c` as that DNA-level
key, so "distinct variants assayed" counts distinct (hgvs_c, hgvs_p) pairs in
the condensed file.

CLI usage:

    python -m src.mave_dataset_stats \\
        data/mave_data/integrated_variant_effect_dataset.condensed.tsv.gz \\
        data/mave_data/Supplementary_Data_3.xlsx \\
        [--output stats.csv]
"""

import argparse
from pathlib import Path

import pandas as pd

DEFAULT_CONDENSED_FILE = Path("data/mave_data/integrated_variant_effect_dataset.condensed.tsv.gz")
DEFAULT_METADATA_FILE = Path("data/mave_data/Supplementary_Data_3.xlsx")

METADATA_SHEET = "Curation"
DATASET_COL = "Dataset"
GENE_COL = "Gene"
GENOMIC_VARIANT_COL = "hgvs_c"
PROTEIN_VARIANT_COL = "hgvs_p"

METADATA_DATASET_COL = "Dataset Name"
IGVF_PRODUCED_COL = "IGVF Produced?"
SCORE_SET_TYPE_COL = "Primary Score Set or Meta-analysis?"
MEASUREMENT_VALUE = "primary score set"


def split_genes(gene_value):
    return [g.strip() for g in gene_value.split(",") if g.strip()]


def load_dataset_metadata(metadata_path):
    """Return the Curation sheet indexed by dataset name.

    Raises ValueError if any dataset referenced by the condensed file is
    missing from the metadata, since every downstream stat depends on the
    IGVF/measurement classification being complete.
    """
    metadata = pd.read_excel(metadata_path, sheet_name=METADATA_SHEET)
    metadata = metadata.set_index(METADATA_DATASET_COL)
    if metadata.index.has_duplicates:
        dupes = sorted(set(metadata.index[metadata.index.duplicated()]))
        raise ValueError(f"Duplicate dataset name(s) in {METADATA_SHEET} sheet: {dupes}")
    return metadata


def genes_in(df):
    genes = set()
    for value in df[GENE_COL].unique():
        genes.update(split_genes(value))
    return genes


def compute_bucket_stats(condensed, dataset_names, measurement_datasets):
    """Compute the five (or six) summary stats for one bucket of datasets."""
    sub = condensed[condensed[DATASET_COL].isin(dataset_names)]
    is_measurement_row = sub[DATASET_COL].isin(measurement_datasets)
    n_measurements = int(is_measurement_row.sum())
    n_composite = int((~is_measurement_row).sum())
    n_variants = sub[[GENOMIC_VARIANT_COL, PROTEIN_VARIANT_COL]].drop_duplicates().shape[0]
    genes = genes_in(sub)
    return {
        "datasets": sub[DATASET_COL].nunique(),
        "variant_effect_measurements": n_measurements,
        "composite_scores": n_composite,
        "distinct_variants_assayed": n_variants,
        "genes_represented": len(genes),
    }, genes


def compute_all_stats(condensed_path, metadata_path):
    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)

    condensed_datasets = set(condensed[DATASET_COL].unique())
    missing = condensed_datasets - set(metadata.index)
    if missing:
        raise ValueError(f"Dataset(s) in condensed file missing from {METADATA_SHEET} metadata: {sorted(missing)}")

    is_igvf = metadata[IGVF_PRODUCED_COL].eq("Yes")
    is_measurement = metadata[SCORE_SET_TYPE_COL].eq(MEASUREMENT_VALUE)
    measurement_datasets = set(metadata.index[is_measurement])

    igvf_datasets = set(metadata.index[is_igvf])
    non_igvf_datasets = set(metadata.index[~is_igvf])
    all_datasets = igvf_datasets | non_igvf_datasets

    igvf_stats, igvf_genes = compute_bucket_stats(condensed, igvf_datasets, measurement_datasets)
    non_igvf_stats, non_igvf_genes = compute_bucket_stats(condensed, non_igvf_datasets, measurement_datasets)
    combined_stats, _ = compute_bucket_stats(condensed, all_datasets, measurement_datasets)

    non_igvf_stats["genes_not_in_igvf_data"] = len(non_igvf_genes - igvf_genes)

    return {
        "Community (IGVF)": igvf_stats,
        "Community (non-IGVF)": non_igvf_stats,
        "Combined (IGVF + community)": combined_stats,
    }


def stats_to_dataframe(stats):
    table = pd.DataFrame(stats).T.reindex(
        columns=[
            "datasets",
            "variant_effect_measurements",
            "composite_scores",
            "distinct_variants_assayed",
            "genes_represented",
            "genes_not_in_igvf_data",
        ]
    )
    return table.astype("Int64")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "condensed_file",
        type=Path,
        nargs="?",
        default=DEFAULT_CONDENSED_FILE,
        help=f"Condensed variant effect dataset TSV (default: {DEFAULT_CONDENSED_FILE})",
    )
    parser.add_argument(
        "metadata_file",
        type=Path,
        nargs="?",
        default=DEFAULT_METADATA_FILE,
        help=f"Supplementary Data 3 metadata workbook (default: {DEFAULT_METADATA_FILE})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional path to also write the summary table as CSV",
    )
    args = parser.parse_args()

    try:
        stats = compute_all_stats(args.condensed_file, args.metadata_file)
    except ValueError as exc:
        parser.error(str(exc))

    table = stats_to_dataframe(stats)
    print(table.to_string())

    if args.output:
        table.to_csv(args.output)
        print(f"\nWrote summary table to {args.output}")


if __name__ == "__main__":
    main()

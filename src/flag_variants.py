#!/usr/bin/env python3
"""Add a Flag column to a CVFG variants TSV.

Flag marks DNA-level variants that should be excluded from downstream
analysis. It is a pipe-delimited list of "" / "*" aligned to the other
DNA-resolution columns (e.g. mapped_hgvs_g), since a row may describe a
single DNA variant or a protein variant with several reverse-translation
candidates.

Three kinds of datasets are handled, selected by the `dataset_name` column:

\b
1. MAVE-join datasets (BAP1_Waters_2024, CHEK2_Gebbia_2024, OTC_Lo_2023,
   RAD51C_Olvera-Leon_2024): rows are joined against the investigator's
   original assay file in
   data/filtering/<dataset_name>.tsv. If any DNA candidate in a row
   matches an assay record that meets the dataset's flagging rule, every
   position in the row's Flag list is set to "*".
2. KCNQ4_Zheng_2022_current_homozygous / KCNQ4_Zheng_2022_v12_homozygous:
   flagged per DNA candidate, marking anything that is not a SNV.
3. Everything else: Flag is a list of empty strings.
"""

from pathlib import Path

import click
import pandas as pd

DEFAULT_FILTERING_DIR = Path("data/filtering")

# Datasets joined to their assay file by genomic coordinates (chromosome,
# position, ref, alt), optionally OR'd with an HGVSc string match.
COORD_JOIN_DATASETS = {
    "BAP1_Waters_2024": {
        "assay_file": "BAP1_Waters_2024.tsv",
        "chrom_col": "ref_chr",
        "pos_col": "pos",
        "ref_col": "ref",
        "alt_col": "alt",
        "hgvsc_col": "HGVSc",
        "mark": lambda assay: assay["pam_flag"] == "Y",
    },
    "RAD51C_Olvera-León_2024": {
        "assay_file": "RAD51C_Olvera-León_2024.tsv",
        "chrom_col": "ref_chr",
        "pos_col": "pos",
        "ref_col": "ref",
        "alt_col": "alt",
        "hgvsc_col": None,
        "mark": lambda assay: assay["pam_codon"] == "Y",
    },
}

# Datasets joined to their assay file by a single protein-level HGVS string
# on the main dataframe (raw_hgvs_pro).
STRING_JOIN_DATASETS = {
    "CHEK2_Gebbia_2024": {
        "assay_file": "CHEK2_Gebbia_2024.tsv",
        "main_key_col": "raw_hgvs_pro",
        "to_assay_key": lambda assay: assay["hgvs_pro"].str.replace("*", "Ter", regex=False),
        "mark": lambda assay: (assay["Filter_CI"] == "TRUE")
        | (assay["Filter_Hypercomplement"] == "TRUE"),
    },
    "OTC_Lo_2023": {
        "assay_file": "OTC_Lo_2023.tsv",
        "main_key_col": "raw_hgvs_pro",
        "to_assay_key": lambda assay: "p." + assay["Human_Variant"],
        "mark": lambda assay: assay["Functional_Class"] == "SMG Loop",
    },
}

# Datasets flagged per-DNA-candidate for being a non-SNV.
SNV_FILTER_DATASETS = {
    "KCNQ4_Zheng_2022_current_homozygous",
    "KCNQ4_Zheng_2022_v12_homozygous",
}


def split_pipe(value):
    return value.split("|") if value else []


def dna_variant_count(mapped_hgvs_g):
    return len(split_pipe(mapped_hgvs_g))


def empty_flag(count):
    return "|".join([""] * count)


def starred_flag(count):
    return "|".join(["*"] * count)


def explode_dna_lists(sub, cols):
    """Split pipe-delimited `cols` in lockstep, one row per DNA candidate.

    Keeps the original dataframe index as a `row_id` column so matches can
    be mapped back to their source row.
    """
    work = sub[cols].copy()
    for col in cols:
        work[col] = work[col].apply(split_pipe)
    work = work.explode(cols, ignore_index=False)
    work.index.name = "row_id"
    return work.reset_index()


def load_assay(filtering_dir, filename):
    path = Path(filtering_dir) / filename
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def flagged_rows_for_coord_join(df, dataset_name, cfg, filtering_dir):
    """Return the set of row indices to flag for a coordinate-join dataset."""
    sub = df[df["dataset_name"] == dataset_name]
    if sub.empty:
        return set()

    cols = [
        "mapped_hgvs_g_chromosome",
        "mapped_hgvs_g_start",
        "mapped_hgvs_g_ref",
        "mapped_hgvs_g_alt",
    ]
    if cfg["hgvsc_col"]:
        cols = cols + ["raw_hgvs_nt"]
    long_df = explode_dna_lists(sub, cols)

    assay = load_assay(filtering_dir, cfg["assay_file"])
    mark = cfg["mark"](assay)

    assay_coord = pd.DataFrame(
        {
            "a_chrom": assay[cfg["chrom_col"]],
            "a_pos": pd.to_numeric(assay[cfg["pos_col"]], errors="coerce"),
            "a_ref": assay[cfg["ref_col"]],
            "a_alt": assay[cfg["alt_col"]],
            "_mark": mark,
        }
    )

    long_df["_chrom_key"] = "chr" + long_df["mapped_hgvs_g_chromosome"]
    long_df["_pos_key"] = pd.to_numeric(long_df["mapped_hgvs_g_start"], errors="coerce")

    coord_matches = long_df.merge(
        assay_coord,
        left_on=["_chrom_key", "_pos_key", "mapped_hgvs_g_ref", "mapped_hgvs_g_alt"],
        right_on=["a_chrom", "a_pos", "a_ref", "a_alt"],
        how="inner",
    )
    flagged_row_ids = set(coord_matches.loc[coord_matches["_mark"], "row_id"])

    if cfg["hgvsc_col"]:
        assay_hgvsc = pd.DataFrame({"a_hgvsc": assay[cfg["hgvsc_col"]], "_mark": mark})
        hgvsc_matches = long_df.merge(
            assay_hgvsc, left_on="raw_hgvs_nt", right_on="a_hgvsc", how="inner"
        )
        flagged_row_ids |= set(hgvsc_matches.loc[hgvsc_matches["_mark"], "row_id"])

    return flagged_row_ids


def flagged_rows_for_string_join(df, dataset_name, cfg, filtering_dir):
    """Return the set of row indices to flag for a single-string-join dataset."""
    sub = df[df["dataset_name"] == dataset_name]
    if sub.empty:
        return set()

    assay = load_assay(filtering_dir, cfg["assay_file"])
    mark = cfg["mark"](assay)
    assay_key = cfg["to_assay_key"](assay)
    flagged_keys = set(assay_key[mark])

    matches = sub[cfg["main_key_col"]].isin(flagged_keys)
    return set(sub.index[matches])


def flag_snv_filter_row(row):
    starts = split_pipe(row["mapped_hgvs_g_start"])
    stops = split_pipe(row["mapped_hgvs_g_stop"])
    alts = split_pipe(row["mapped_hgvs_g_alt"])
    flags = ["" if start == stop and len(alt) == 1 else "*" for start, stop, alt in zip(starts, stops, alts)]
    return "|".join(flags)


def compute_flags(df, filtering_dir):
    counts = df["mapped_hgvs_g"].apply(dna_variant_count)
    flag = counts.apply(empty_flag)

    snv_mask = df["dataset_name"].isin(SNV_FILTER_DATASETS)
    if snv_mask.any():
        flag.loc[snv_mask] = df.loc[snv_mask].apply(flag_snv_filter_row, axis=1)

    for dataset_name, cfg in COORD_JOIN_DATASETS.items():
        flagged_ids = flagged_rows_for_coord_join(df, dataset_name, cfg, filtering_dir)
        if flagged_ids:
            idx = list(flagged_ids)
            flag.loc[idx] = counts.loc[idx].apply(starred_flag)

    for dataset_name, cfg in STRING_JOIN_DATASETS.items():
        flagged_ids = flagged_rows_for_string_join(df, dataset_name, cfg, filtering_dir)
        if flagged_ids:
            idx = list(flagged_ids)
            flag.loc[idx] = counts.loc[idx].apply(starred_flag)

    return flag


@click.command(help=__doc__)
@click.argument("input", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.argument("output", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--filtering-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_FILTERING_DIR,
    show_default=True,
    help="Directory containing per-dataset assay files",
)
def main(input, output, filtering_dir):
    df = pd.read_csv(input, sep="\t", dtype=str, keep_default_na=False, engine="c")
    df["Flag"] = compute_flags(df, filtering_dir)

    n_flagged_variants = df["Flag"].apply(lambda f: f.count("*")).sum()
    n_flagged_rows = (df["Flag"].apply(lambda f: "*" in f)).sum()
    click.echo(f"Flagged {n_flagged_variants} DNA variant(s) across {n_flagged_rows} row(s) of {len(df)}.")

    per_dataset = pd.DataFrame(
        {
            "flagged_rows": df["Flag"].apply(lambda f: "*" in f),
            "flagged_variants": df["Flag"].apply(lambda f: f.count("*")),
        }
    ).groupby(df["dataset_name"]).sum()
    per_dataset = per_dataset[per_dataset["flagged_variants"] > 0].sort_index()
    for dataset_name, row in per_dataset.iterrows():
        click.echo(f"  {dataset_name}: {int(row['flagged_rows'])} row(s), {int(row['flagged_variants'])} variant(s)")

    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    main()

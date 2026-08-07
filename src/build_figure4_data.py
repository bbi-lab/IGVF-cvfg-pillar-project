"""Rebuild Main_Figures/Figure_4/figure4_data.json.gz from current pipeline outputs.

`figure4.ipynb` only *loads* that cache; nothing in the repo previously built it.
This script reconstructs the parts that a current pipeline run can actually
support, and carries the rest forward unchanged from a prior cache via
``--cached-json``.

Reconstructed from scratch (panels a/b/e/f):
    - `scoreset_2018` / `scoreset`: `Scoreset` objects for the MSH2_Jia_2021
      dataset (2018 vs. 2025 ClinVar release), from
      `data/output/maves/integrated_variant_effect_dataset.tsv.gz`.
    - `n_c`, `scoreset_flipped`, and `indv_summary`'s `prior`/`point_ranges`:
      from the per-dataset exCALIBR calibration JSON at
      `data/input/mave_calibration/excalibr/json/MSH2_Jia_2021_clinvar_2018.json`.
    - `gene_4e`, `dist_4e`, `labdat_4e`, `snvdf_4e`, `sorted_thresholds_4e`,
      `oldsorted_thresholds_4e`, `dist_4f`, `finalout_4f`: MSH2 (and, for panel
      f, a fixed 6-gene comparison set) REVEL gene-specific-vs-genome-wide
      calibration, built from the integrated dataset plus the
      `REVEL_gene_specific_calibration` sheet of `Supplementary_Data_4.xlsx`.

NOT reconstructible from anything currently in this repo -- must be carried
forward from a prior `figure4_data.json.gz` via `--cached-json`:
    - `fits`, and `indv_summary`'s `priors`/`log_lr_plus`/`score_range`/`C`,
      and the top-level `score_range` (panel a): the smoothed density-band fit
      requires exCALIBR's raw per-bootstrap fit output. The only exCALIBR
      artifact checked into `data/input/mave_calibration/` is the *aggregated*
      per-dataset summary (prior/point_ranges/etc.) -- the same fields
      `src/load_excalibr_calibrations.py` reads. The 1000-entry raw bootstrap
      ensemble isn't retained anywhere on disk.
    - `danzs_oob`, `auths_oob`, `datasets` (panel c): out-of-bag confusion
      matrices per MAVE dataset. This requires exCALIBR's out-of-bag
      evidence-direction predictions, which aren't in any current pipeline
      output either.
    - `prior`, `Post_p`, `Post_b`, `p_data`, `b_data` (panel d): this panel is
      an illustrative cartoon (the posterior curves it plots,
      `1 - x**2` and `exp(3x)/exp(3)`, are hardcoded directly in
      `plot_utils.plot_panel_d`, not derived from `p_data`/`b_data`). There is
      no recorded recipe for how the original synthetic `p_data`/`b_data`/
      `Post_p`/`Post_b` values were generated, so they can only be carried
      forward, not regenerated.

Two notes on drift from the original figure:
    - The integrated dataset no longer has a plain `ID` column (variant
      grouping key `Scoreset`/`Variant` expect); it now ships
      `mavedb_variant_urn` instead. This script renames it before constructing
      `Scoreset` objects -- confirmed to reproduce the original panel b legend
      counts almost exactly (e.g. MSH2 ClinVar BLB: 230 -> 229 today).
    - `FIGURE_4F_GENES` (BRCA1/BRCA2/F9/MSH2/TP53/TSC2) is the fixed gene list
      the *original* figure compared in panel f. It is not recoverable from
      any single data-driven rule (e.g. "genes with a fully-populated
      gene-specific REVEL calibration" matches 30 genes, not these 6, and
      doesn't even contain all of these 6) -- treat it as an editorial choice
      and confirm with the paper authors before changing it.

Usage:
    python -m src.build_figure4_data --cached-json Main_Figures/Figure_4/old_figure4_data.json.gz \\
        [--integrated-dataset PATH] [--excalibr-json-dir PATH] [--supplementary-data-4 PATH] \\
        [--output PATH]

Omitting `--cached-json` raises immediately, listing the fields above that
have no other source.
"""

import gzip
import json
import sys
from pathlib import Path

import click
import numpy as np
import pandas as pd

FIGURE_4_DIR = Path(__file__).resolve().parent.parent / "Main_Figures" / "Figure_4"
sys.path.insert(0, str(FIGURE_4_DIR))
from data_utils import Scoreset
from plot_utils import fig_json_encoder, fig_json_hook

DEFAULT_INTEGRATED_DATASET = Path("data/output/maves/integrated_variant_effect_dataset.tsv.gz")
DEFAULT_EXCALIBR_JSON_DIR = Path("data/input/mave_calibration/excalibr/json")
DEFAULT_SUPPLEMENTARY_DATA_4 = Path("data/output/supplementary_data/Supplementary_Data_4.xlsx")
DEFAULT_OUTPUT = FIGURE_4_DIR / "figure4_data.json.gz"

MSH2_GENE = "MSH2"
MSH2_DATASET = "MSH2_Jia_2021"
MSH2_EXCALIBR_DATASET = "MSH2_Jia_2021_clinvar_2018"

# Fields that only exist in a prior figure4_data.json.gz -- see module docstring.
CARRYOVER_TOP_LEVEL_FIELDS = [
    "score_range",
    "fits",
    "danzs_oob",
    "auths_oob",
    "datasets",
    "prior",
    "Post_p",
    "Post_b",
    "p_data",
    "b_data",
]
CARRYOVER_INDV_SUMMARY_FIELDS = ["priors", "log_lr_plus", "C"]

# Bergquist et al. genome-wide REVEL evidence-strength thresholds, as used in
# notebooks/analysis/Variant_Classification_analysis.ipynb (the cell defining
# `REVEL_conditions`/`REVEL_values`).
GENOME_WIDE_REVEL_THRESHOLDS = pd.Series(
    {
        "BP4_Very Strong": np.nan,
        "BP4_Strong": 0.016,
        "BP4_Moderate+": 0.052,
        "BP4_Moderate": 0.183,
        "BP4_Supporting": 0.290,
        "PP3_Supporting": 0.644,
        "PP3_Moderate": 0.773,
        "PP3_Moderate+": 0.879,
        "PP3_Strong": 0.932,
        "PP3_Very Strong": np.nan,
    }
)

REVEL_TIER_ORDER = list(GENOME_WIDE_REVEL_THRESHOLDS.index)

# Mirrors the local `strenth_to_point` dict inside plot_utils.plot_panel_e.
REVEL_TIER_TO_POINT = {
    "BP4_Very Strong": -8,
    "BP4_Strong": -4,
    "BP4_Moderate+": -3,
    "BP4_Moderate": -2,
    "BP4_Supporting": -1,
    "IR": 0,
    "PP3_Supporting": 1,
    "PP3_Moderate": 2,
    "PP3_Moderate+": 3,
    "PP3_Strong": 4,
    "PP3_Very Strong": 8,
}

# Fixed panel-4f gene comparison set -- see module docstring.
FIGURE_4F_GENES = ["BRCA1", "BRCA2", "F9", "MSH2", "TP53", "TSC2"]

BLB_SIGS = {"Benign", "Likely benign", "Benign/Likely benign"}
PLP_SIGS = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
CONFLICTING_SIG = "Conflicting classifications of pathogenicity"
VUS_SIG = "Uncertain significance"


def load_integrated_dataset(path):
    return pd.read_csv(path, sep="\t", low_memory=False)


def unique_gene_snvs_with_revel(df, gene):
    """Deduplicated (by genomic position) SNVs for `gene` that have a REVEL score."""
    gene_df = df[df["Gene"] == gene]
    gene_df = gene_df.drop_duplicates(subset=["Chrom", "hg38_start", "ref_allele", "alt_allele"])
    return gene_df[gene_df["REVEL"].notna()].copy()


def build_msh2_scoresets(integrated_df, **scoreset_kwargs):
    """Build the panel a/b `Scoreset` objects for MSH2_Jia_2021.

    `Scoreset`/`Variant` (data_utils.py) group rows by a plain `ID` column that
    no longer exists in `integrated_variant_effect_dataset.tsv.gz` -- it now
    ships `mavedb_variant_urn` instead. Renaming it here reproduces the
    original panel b legend counts closely (e.g. MSH2 ClinVar BLB: 230 in the
    published figure vs. 229 today), the remaining difference being ClinVar's
    ongoing updates since the figure was made, not a reconstruction error.
    """
    msh2_df = integrated_df[integrated_df["Dataset"] == MSH2_DATASET].rename(columns={"mavedb_variant_urn": "ID"})
    if not len(msh2_df):
        raise ValueError(f"No rows found for Dataset == {MSH2_DATASET!r} in the integrated dataset")

    scoreset_2018 = Scoreset(msh2_df, clinvar_release="2018", **scoreset_kwargs)
    scoreset = Scoreset(msh2_df, clinvar_release="2025", **scoreset_kwargs)
    return scoreset_2018, scoreset


def load_excalibr_summary(json_dir, dataset=MSH2_EXCALIBR_DATASET):
    """Load the fields of a per-dataset exCALIBR calibration JSON that this
    script can actually reconstruct: `prior`, `point_ranges`, `n_c`, and
    `scoreset_flipped`. See module docstring for what's *not* in this file.
    """
    path = Path(json_dir) / f"{dataset}.json"
    if not path.exists():
        raise FileNotFoundError(f"No exCALIBR calibration JSON found at {path}")
    data = json.loads(path.read_text())
    return {
        "prior": data["prior"],
        "point_ranges": data["point_ranges"],
        "n_c": data["n_c"],
        "scoreset_flipped": bool(data["scoreset_flipped"]),
    }


def score_to_revel_tier(scores, thresholds=GENOME_WIDE_REVEL_THRESHOLDS):
    """Map REVEL scores to evidence-tier labels ('IR' for indeterminate) given
    an ascending Series of tier thresholds (NaN entries are skipped)."""
    populated = thresholds.dropna().sort_values()
    tiers = pd.Series("IR", index=scores.index, dtype=object)
    benign_tiers = [t for t in populated.index if t.startswith("BP4")]
    pathogenic_tiers = [t for t in populated.index if t.startswith("PP3")]

    # Benign tiers: score <= threshold, strongest (smallest threshold) wins.
    for tier in reversed(benign_tiers):
        tiers[scores <= populated[tier]] = tier
    # Pathogenic tiers: score >= threshold, strongest (largest threshold) wins.
    for tier in pathogenic_tiers:
        tiers[scores >= populated[tier]] = tier
    return tiers


def categorize_clinical_status(gene_snvs, clinvar_release="2025"):
    """Assign each unique SNV to exactly one of BLB/PLP/Conflicting/VUS/gnomAD/
    allSNVs, in that priority order (ClinVar status outranks population
    frequency, with everything else -- including truly unobserved variants --
    falling back to 'allSNVs'). This mirrors the exhaustive partition seen in
    the original `finalout_4f` (its category counts sum exactly to the unique
    SNV count), but the priority order itself is inferred, not verified
    against original source code.
    """
    sig_col = f"clinvar_sig_{clinvar_release}"
    sig = gene_snvs[sig_col]
    category = pd.Series("allSNVs", index=gene_snvs.index, dtype=object)
    category[gene_snvs["gnomad_MAF"].notna()] = "gnomAD"
    category[sig == VUS_SIG] = "VUS"
    category[sig == CONFLICTING_SIG] = "Conflicting"
    category[sig.isin(BLB_SIGS)] = "BLB"
    category[sig.isin(PLP_SIGS)] = "PLP"
    return category


def load_gene_specific_revel_thresholds(supplementary_data_4_path, gene):
    df = pd.read_excel(supplementary_data_4_path, sheet_name="REVEL_gene_specific_calibration")
    df = df.set_index("Gene")
    if gene not in df.index:
        raise ValueError(f"{gene!r} has no row in the REVEL_gene_specific_calibration sheet")
    return df.loc[gene, REVEL_TIER_ORDER]


def build_panel_e_data(integrated_df, supplementary_data_4_path, gene=MSH2_GENE):
    gene_snvs = unique_gene_snvs_with_revel(integrated_df, gene)
    sig = gene_snvs["clinvar_sig_2025"]

    labdat = pd.DataFrame(
        {
            0: gene_snvs.loc[sig.isin(BLB_SIGS) | sig.isin(PLP_SIGS), "REVEL"],
            1: sig[sig.isin(BLB_SIGS) | sig.isin(PLP_SIGS)].isin(PLP_SIGS).astype(int),
        }
    ).reset_index(drop=True)

    snvdf = pd.DataFrame(
        {
            "REVEL": gene_snvs["REVEL"].values,
            "merg_clinvar_sig": "allSNVs",
            "GeneSymbol": gene,
        }
    )

    sorted_thresholds = load_gene_specific_revel_thresholds(supplementary_data_4_path, gene)
    return {
        "gene_4e": gene,
        "dist_4e": "REVEL",
        "labdat_4e": labdat,
        "snvdf_4e": snvdf,
        "sorted_thresholds_4e": sorted_thresholds,
        "oldsorted_thresholds_4e": GENOME_WIDE_REVEL_THRESHOLDS.copy(),
    }


def build_panel_f_data(integrated_df, supplementary_data_4_path, genes=FIGURE_4F_GENES):
    rows = []
    for gene in genes:
        gene_snvs = unique_gene_snvs_with_revel(integrated_df, gene)
        gene_specific_thresholds = load_gene_specific_revel_thresholds(supplementary_data_4_path, gene)

        old_tier = score_to_revel_tier(gene_snvs["REVEL"], GENOME_WIDE_REVEL_THRESHOLDS)
        new_tier = score_to_revel_tier(gene_snvs["REVEL"], gene_specific_thresholds)

        rows.append(
            pd.DataFrame(
                {
                    "GeneSymbol": gene,
                    "merg_clinvar_sig": categorize_clinical_status(gene_snvs),
                    "REVEL": gene_snvs["REVEL"].values,
                    "Old Thresh": old_tier.values,
                    "New Thresh": new_tier.values,
                    "Old_scr": old_tier.map(REVEL_TIER_TO_POINT).values,
                    "New_scr": new_tier.map(REVEL_TIER_TO_POINT).values,
                }
            )
        )
    finalout = pd.concat(rows, ignore_index=True).reset_index().rename(columns={"index": "Unnamed: 0"})
    return {"dist_4f": "REVEL", "finalout_4f": finalout}


def load_carryover_fields(cached_json_path):
    """Decode a prior figure4_data.json.gz and pull out only the fields that
    cannot currently be rebuilt from pipeline outputs (see module docstring).
    """
    with gzip.open(cached_json_path, "rt", encoding="utf-8") as f:
        cached = json.load(f, object_hook=fig_json_hook)

    top_level = {key: cached[key] for key in CARRYOVER_TOP_LEVEL_FIELDS}
    indv_summary = {key: cached["indv_summary"][key] for key in CARRYOVER_INDV_SUMMARY_FIELDS}
    return top_level, indv_summary


def build_figure4_data(
    integrated_dataset_path=DEFAULT_INTEGRATED_DATASET,
    excalibr_json_dir=DEFAULT_EXCALIBR_JSON_DIR,
    supplementary_data_4_path=DEFAULT_SUPPLEMENTARY_DATA_4,
    cached_json_path=None,
):
    if cached_json_path is None:
        raise ValueError(
            "No --cached-json supplied. Panel a's exCALIBR bootstrap density fit "
            "(fits/priors/log_lr_plus/score_range), panel c's out-of-bag confusion "
            "matrices (danzs_oob/auths_oob/datasets), and panel d's cartoon "
            "(prior/Post_p/Post_b/p_data/b_data) have no source anywhere in the "
            "current pipeline outputs -- see this module's docstring. Pass "
            "--cached-json pointing at a prior figure4_data.json.gz (or "
            "old_figure4_data.json.gz) to carry those fields forward unchanged."
        )

    integrated_df = load_integrated_dataset(integrated_dataset_path)
    scoreset_2018, scoreset = build_msh2_scoresets(integrated_df)
    excalibr_summary = load_excalibr_summary(excalibr_json_dir)
    panel_e = build_panel_e_data(integrated_df, supplementary_data_4_path)
    panel_f = build_panel_f_data(integrated_df, supplementary_data_4_path)
    carryover_top_level, carryover_indv_summary = load_carryover_fields(cached_json_path)

    n_samples = int((scoreset.sample_counts > 0).sum())

    data = {
        "scoreset_2018": scoreset_2018,
        "scoreset": scoreset,
        "indv_summary": {
            "prior": excalibr_summary["prior"],
            "point_ranges": excalibr_summary["point_ranges"],
            **carryover_indv_summary,
        },
        "n_c": excalibr_summary["n_c"],
        "n_samples": n_samples,
        "scoreset_flipped": excalibr_summary["scoreset_flipped"],
        **panel_e,
        **panel_f,
        **carryover_top_level,
    }
    return data


@click.command(help=__doc__)
@click.option(
    "--integrated-dataset",
    "integrated_dataset_path",
    default=DEFAULT_INTEGRATED_DATASET,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--excalibr-json-dir",
    "excalibr_json_dir",
    default=DEFAULT_EXCALIBR_JSON_DIR,
    type=click.Path(exists=True, file_okay=False, path_type=Path),
)
@click.option(
    "--supplementary-data-4",
    "supplementary_data_4_path",
    default=DEFAULT_SUPPLEMENTARY_DATA_4,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--cached-json",
    "cached_json_path",
    default=None,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option("--output", "output_path", default=DEFAULT_OUTPUT, type=click.Path(path_type=Path))
def main(integrated_dataset_path, excalibr_json_dir, supplementary_data_4_path, cached_json_path, output_path):
    try:
        data = build_figure4_data(
            integrated_dataset_path=integrated_dataset_path,
            excalibr_json_dir=excalibr_json_dir,
            supplementary_data_4_path=supplementary_data_4_path,
            cached_json_path=cached_json_path,
        )
    except (ValueError, FileNotFoundError) as exc:
        raise click.ClickException(str(exc)) from exc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", encoding="utf-8") as f:
        json.dump(data, f, cls=fig_json_encoder)

    click.echo(f"Wrote {output_path}")


if __name__ == "__main__":
    main()

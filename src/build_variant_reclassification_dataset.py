#!/usr/bin/env python3
"""Build a deduplicated, DNA-level export of every variant that survives
`Variant_Classification_analysis.ipynb`'s exclusion rules.

Reads the notebook's own intermediate checkpoint
(`data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`,
written at its cell 69, *before* category split), which already has the
LDLR LA-module-1 exclusion and the F9/TP53 restricted-dataset filter baked
in. Three further exclusions the notebook applies later, downstream of that
checkpoint, are re-applied here:

- `SFPQ` is dropped entirely (insufficient ClinVar controls).
- The CHEK2 QC flag (`CHEK2_Gebbia_2024.xlsx`, `Filter_CI == 1`, merged on
  (`hgvs_p` with its transcript prefix stripped, `auth_reported_score`) /
  (`hgvs_pro`, `score`)) is applied and then filtered out, along with every
  other `Flag == '*'` row. The prefix-stripping is a deliberate deviation
  from the notebook's own merge, which joins unstripped `hgvs_p` (e.g.
  `"NP_009125.1:p.Ser2Ala"`) directly against `hgvs_pro` (e.g.
  `"p.Ser2Ala"`, no prefix) -- confirmed against real data to never
  actually match, an apparently unintentional no-op in production that's
  out of scope to fix there.
- Rows tagged `conflicting_fxn_data`, `splice_variant_not_measured`, or
  `start_lost_variant_not_measured` in `VariantNotes`, or `splice_var_amino
  == 'Yes'`, are dropped.

Four points columns are added:

- `ExCALIBR_points`: the literal ExCALIBR score-interval calibration value
  (`ExC_points_2018`, falling back to `ExC_points_2025`, for `BRCA1`/
  `PTEN`/`MSH2`; `ExC_points_2025` for every other gene) -- independent of
  whether ExCALIBR is what's actually used for a given gene's evidence.
- `OddsPath_points`: `OP_points` verbatim.
- `Functional_points`: `Fxn_points` verbatim -- the pipeline's own choice of
  `ExCALIBR_points` or `OddsPath_points` per gene (`F9`/`TP53` use
  `OddsPath_points`; every other gene uses `ExCALIBR_points`).
- `Combined_points`: `Functional_points + Points_REVEL_GeneSpecific_GenomeWide`
  (gene-specific REVEL, falling back to genome-wide) -- matches the
  notebook's own `Total_Points_GeneSpecific_REVEL`.

Finally, every surviving row is collapsed to one per DNA variant
(`src.lib.dedup.dedup_by_max_abs_points`, keyed on `Gene`/`Chrom`/
`hg38_start`/`ref_allele`/`alt_allele`): the candidate with the greatest
`abs(Combined_points)` wins, ties broken by `Dataset` name for
determinism. This also reconciles amino-acid-resolution assay rows onto
their DNA coordinate, since they share the same genomic key columns as any
nt-resolution row for the same physical variant.

Output columns match `integrated_variant_effect_dataset.tsv`'s schema,
except `ID` is replaced by `mavedb_variant_urn` (a different identifier
scheme -- the checkpoint doesn't carry the raw file's `ID` values) and
`hg19_pos`/`auth_reported_rep_score` are dropped (not present in the
checkpoint), plus the four new points columns appended at the end.
"""

from pathlib import Path

import click
import pandas as pd

from src.lib.dedup import GENOMIC_KEY_COLS, dedup_by_max_abs_points

DEFAULT_CHECKPOINT_FILE = Path("data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz")
DEFAULT_CHEK2_FILE = Path("data/input/maves/CHEK2_Gebbia_2024.xlsx")
DEFAULT_OUTPUT_FILE = Path("data/output/reclassification/integrated_variant_effect_reclassification.tsv.gz")

# ExC_points vintage override for ExCALIBR_points. Deliberately excludes
# TP53 -- unlike BRCA1/PTEN/MSH2, TP53's Functional_points come entirely
# from OddsPath, not a 2018-vs-2025 ExCALIBR choice (see
# Variant_Classification_analysis.ipynb cell 19).
EXCALIBR_VINTAGE_OVERRIDE_GENES = frozenset({"BRCA1", "PTEN", "MSH2"})

DISALLOWED_VARIANT_NOTES = frozenset({
    "conflicting_fxn_data",
    "splice_variant_not_measured",
    "splice_variant_not_measured;conflicting_fxn_data",
    "start_lost_variant_not_measured",
})

OUTPUT_COLUMNS = [
    "mavedb_variant_urn", "Dataset", "Gene", "HGNC_id", "Chrom", "Strand", "hg38_start", "hg38_end",
    "ref_allele", "alt_allele", "auth_transcript_id", "transcript_pos", "transcript_ref", "transcript_alt",
    "aa_pos", "aa_ref", "aa_alt", "hgvs_c", "hgvs_p", "consequence", "simplified_consequence",
    "auth_reported_score", "auth_reported_func_class", "splice_measure", "gnomad_MAF",
    "clinvar_sig_2025", "clinvar_star_2025", "clinvar_date_last_reviewed_2025",
    "clinvar_sig_2018", "clinvar_star_2018", "clinvar_date_last_reviewed_2018",
    "nucleotide_or_aa", "Ensembl Transcript ID", "RefSeq Transcript ID",
    "Interval 1 Name", "Interval 1 Range", "Interval 1 Class",
    "Interval 2 Name", "Interval 2 Range", "Interval 2 Class",
    "Interval 3 Name", "Interval 3 Range", "Interval 3 Class",
    "Interval 4 Name", "Interval 4 Range", "Interval 4 Class",
    "Interval 5 Name", "Interval 5 Range", "Interval 5 Class",
    "Interval 6 Name", "Interval 6 Range", "Interval 6 Class",
    "Flag", "REVEL", "REVEL_train", "AM_score", "AM_class", "MutPred2", "MP2_train",
    "spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL",
    "spliceAI_DP_AG", "spliceAI_DP_AL", "spliceAI_DP_DG", "spliceAI_DP_DL",
    "ClinVar Variation Id_ClinGen_repo", "Allele Registry Id_ClinGen_repo", "Disease_ClinGen_repo",
    "Mondo Id_ClinGen_repo", "Mode of Inheritance_ClinGen_repo", "Assertion_ClinGen_repo",
    "Applied Evidence Codes (Met)_ClinGen_repo", "Applied Evidence Codes (Not Met)_ClinGen_repo",
    "Summary of interpretation_ClinGen_repo", "PubMed Articles_ClinGen_repo", "Expert Panel_ClinGen_repo",
    "Guideline_ClinGen_repo", "Approval Date_ClinGen_repo", "Published Date_ClinGen_repo",
    "Retracted_ClinGen_repo", "Evidence Repo Link_ClinGen_repo", "Uuid_ClinGen_repo",
    "Updated_Classification_ClinGen_repo", "Updated_Evidence Codes_ClinGen_repo",
    "ExCALIBR_points", "OddsPath_points", "Functional_points", "Combined_points",
]


def apply_notebook_exclusions(df: pd.DataFrame, chek2_file: Path) -> pd.DataFrame:
    """Re-apply the checkpoint-downstream exclusions from
    `Variant_Classification_analysis.ipynb` cells 71-77: `SFPQ`, the CHEK2
    QC flag, conflicting/unmeasured-splice `VariantNotes` tags, and any
    other `Flag == '*'` row.

    The CHEK2 merge key is normalized (`hgvs_p`'s transcript prefix, e.g.
    `"NP_009125.1:"`, stripped before matching against `hgvs_pro`) rather
    than reproducing the notebook's own merge verbatim: the notebook joins
    unstripped `hgvs_p` directly against `hgvs_pro`, which never actually
    matches (`hgvs_pro` carries no transcript prefix) -- a pre-existing,
    apparently unintentional no-op in production, flagged separately rather
    than fixed there.
    """
    df = df[df["Gene"] != "SFPQ"].copy()

    chek2 = pd.read_excel(chek2_file, header=0)
    df["_hgvs_p_no_transcript"] = df["hgvs_p"].str.replace(r"^[^:]+:", "", regex=True)
    df = df.merge(
        chek2[["hgvs_pro", "score", "Filter_CI"]],
        left_on=["_hgvs_p_no_transcript", "auth_reported_score"],
        right_on=["hgvs_pro", "score"],
        how="left",
    )
    df["Flag"] = df["Flag"].where(df["Filter_CI"] != 1, "*")
    df = df.drop(columns=["_hgvs_p_no_transcript", "hgvs_pro", "score", "Filter_CI"])

    df = df[
        ~df["VariantNotes"].isin(DISALLOWED_VARIANT_NOTES)
        & (df["splice_var_amino"] != "Yes")
    ]
    df = df[df["Flag"] != "*"]
    return df


def add_points_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add `ExCALIBR_points`, `OddsPath_points`, `Functional_points`, and
    `Combined_points` -- see module docstring for the exact definitions.
    """
    df = df.copy()

    df["ExCALIBR_points"] = df["ExC_points_2025"]
    vintage_mask = df["Gene"].isin(EXCALIBR_VINTAGE_OVERRIDE_GENES)
    excalibr_2018 = pd.to_numeric(df.loc[vintage_mask, "ExC_points_2018"], errors="coerce")
    excalibr_2025 = pd.to_numeric(df.loc[vintage_mask, "ExC_points_2025"], errors="coerce")
    df.loc[vintage_mask, "ExCALIBR_points"] = excalibr_2018.fillna(excalibr_2025)

    df["OddsPath_points"] = df["OP_points"]
    df["Functional_points"] = df["Fxn_points"]
    df["Combined_points"] = (
        pd.to_numeric(df["Functional_points"], errors="coerce").fillna(0)
        + pd.to_numeric(df["Points_REVEL_GeneSpecific_GenomeWide"], errors="coerce").fillna(0)
    )
    return df


def build_reclassification_dataset(checkpoint_file: Path, chek2_file: Path) -> pd.DataFrame:
    df = pd.read_csv(checkpoint_file)
    df = apply_notebook_exclusions(df, chek2_file)
    df = add_points_columns(df)
    df = dedup_by_max_abs_points(df, points_col="Combined_points", genomic_key_cols=GENOMIC_KEY_COLS)
    return df[OUTPUT_COLUMNS]


@click.command(help=__doc__)
@click.argument(
    "checkpoint_file",
    required=False,
    default=DEFAULT_CHECKPOINT_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--chek2-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CHEK2_FILE,
    help=f"Path to the CHEK2 QC workbook (default {DEFAULT_CHEK2_FILE}).",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_FILE,
    help=f"Output path (default {DEFAULT_OUTPUT_FILE}), written as gzip-compressed TSV.",
)
def main(checkpoint_file: Path, chek2_file: Path, output: Path) -> None:
    result = build_reclassification_dataset(checkpoint_file, chek2_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, sep="\t", index=False, compression="gzip")
    click.echo(f"Wrote {len(result)} rows to {output}")


if __name__ == "__main__":
    main()

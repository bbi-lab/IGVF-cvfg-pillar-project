# MAVE Dataset Stats

`src/mave_dataset_stats.py` reports summary statistics for the integrated MAVE
variant effect dataset (`data/output/maves/integrated_variant_effect_dataset*.tsv.gz`),
split into three groupings of datasets:

- **IGVF** -- datasets Supplementary Data 3 marks `IGVF Produced? = Yes`
- **Community (non-IGVF)** -- everything else
- **Combined (IGVF + community)** -- all datasets together

For each grouping it reports:

- Number of datasets
- Number of variant effect measurements (condensed-file rows belonging to a
  "primary score set" dataset, per Supplementary Data 3)
- Number of composite scores (condensed-file rows belonging to a
  meta-analysis/trained-predictor dataset, per Supplementary Data 3)
- Number of distinct variants assayed (distinct `(hgvs_c, hgvs_p)` pairs in the
  condensed file -- see note below)
- Number of genes represented

The non-IGVF grouping additionally reports the number of genes not also
represented in IGVF data.

Immediately after that table, a **Genes represented** section lists which
genes are covered only by IGVF datasets, only by non-IGVF ("community")
datasets, or by both. CALM1/CALM2/CALM3 are always merged into one
"CALM1/2/3" entry in this list regardless of `--merge-calm-genes` (see
below), since listing the same calmodulin target three times isn't useful;
this can make this list's per-category counts differ by up to two from the
`genes_represented`/`genes_not_in_igvf_data` counts in the table above when
that flag isn't passed.

It also reports score coverage (REVEL, AlphaMissense, MutPred2) and clinical
attributes (VUS, pathogenic/benign, observed in gnomAD) across the assayed
and DNA-level variants. The clinical-attribute breakdown is reported twice:
once using ClinVar 2025 for every gene, and once using ClinVar 2025 for every
gene except BRCA1, PTEN, MSH2, and TP53, which use ClinVar 2018 instead (see
`MIXED_YEAR_GENES` / `mixed_year_clinvar_series` in the script).

By default, a variant with a conflicting or ambiguous ClinVar call --
disagreement between a pathogenic-leaning and benign-leaning classification
across its measurements/DNA candidates, or ClinVar's own "Conflicting
classifications of pathogenicity" call -- is excluded from both the VUS and
pathogenic-or-benign buckets and counted in its own "ClinVar conflict" bucket
instead. This mirrors the conflict handling in
`Analysis/Curation_summary_V5_cleaned.ipynb`, which is why the two scripts'
pathogenic/benign counts can otherwise disagree (that notebook resolves each
protein variant to a single label and drops ambiguous ones, while this
script's older any-match behavior folded them into whichever bucket matched).
Pass `--allow-clinvar-conflicts` to restore that any-match behavior instead,
which also drops the conflict bucket from the report.

It also reports two further sections, sourced from separate input files:

- **ExCALIBR calibration coverage** (Extended Data Figure 2): how many genes
  have a row in `--excalibr-calibrations-file`'s `ExCALIBR_calibrations` sheet
  (default `data/output/supplementary_data/Supplementary_Data_4.xlsx`), and how many
  of those genes have at least one dataset where ExCALIBR assigned at least
  one point of evidence in either direction. That sheet's `dataset` values are
  matched to a gene via this script's own dataset metadata (the same
  `Dataset Name` -> `Gene` mapping used for the summary above), after
  stripping a trailing `_clinvar_2018` suffix and Unicode-normalizing both
  sides -- see `excalibr_dataset_to_gene_map` in the script for why.
- **Reclassification agreement** (Figure 4c): for every sheet in
  `--controls-file` (default `data/output/supplementary_data/Supplementary_Data_5.xlsx`)
  whose name starts with `controls_` -- one per predictor/calibration
  combination, which should agree with each other since this doesn't depend
  on which predictor is active -- how often ExCALIBR's evidence assignment
  (`ExC_points_2025`) and the functional class assignment (`OP_points`) agree
  with the row's ClinVar pathogenic-or-benign control label
  (`clnsig_group_18_25`).

## Inputs

- **Condensed variant effect dataset**
  (default `data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz`):
  one row per variant effect measurement/score. For protein-resolution
  datasets, one row can correspond to several DNA-level changes (pipe-delimited
  within a cell) -- the *expanded* sibling file
  (`integrated_variant_effect_dataset.tsv.gz`) explodes those into one row per
  DNA variant, but this script intentionally works from the condensed file
  since the requested counts are measurement-level, not DNA-variant-level.
- **Dataset metadata** (default `data/output/supplementary_data/Supplementary_Data_3.xlsx`,
  `Curation` sheet): one row per dataset, keyed by `Dataset Name` (joined
  against the condensed file's `Dataset` column). Used for the
  `IGVF Produced?` and `Primary Score Set or Meta-analysis?` classifications,
  and (via its `Gene` column) to resolve genes for the ExCALIBR calibration
  section.
- **ExCALIBR calibrations** (`--excalibr-calibrations-file`, default
  `data/output/supplementary_data/Supplementary_Data_4.xlsx`, `ExCALIBR_calibrations`
  sheet): one row per dataset. A row counts as assigning a point of evidence
  if any of its `range_-8`..`range_8` columns is non-null.
- **Controls** (`--controls-file`, default
  `data/output/supplementary_data/Supplementary_Data_5.xlsx`): every sheet whose name
  starts with `controls_` (e.g. `controls_REVEL_GeneSpecific`), one row per
  control variant, using its `clnsig_group_18_25`, `ExC_points_2025`, and
  `OP_points` columns.

### Note on `(hgvs_g, hgvs_p)`

The integrated dataset has no `hgvs_g` (genomic HGVS) column. Its DNA-level
identifier is `hgvs_c` (transcript-relative HGVS, pipe-delimited per row when a
protein-resolution measurement corresponds to more than one underlying DNA
change). This script uses `hgvs_c` as that DNA-level key, so "distinct
variants assayed" counts distinct `(hgvs_c, hgvs_p)` pairs in the condensed
file. If a true genomic (`NC_...:g.`) identifier is later added to the
integrated dataset, swap `GENOMIC_VARIANT_COL` in the script to point at it.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.mave_dataset_stats \
  data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz \
  data/output/supplementary_data/Supplementary_Data_3.xlsx \
  data/output/maves/integrated_variant_effect_dataset.tsv.gz \
  [--excalibr-calibrations-file data/output/supplementary_data/Supplementary_Data_4.xlsx] \
  [--controls-file data/output/supplementary_data/Supplementary_Data_5.xlsx] \
  [--output stats.txt]
```

All arguments default to the paths above, so a bare invocation works from the
repo root.

Via Docker (same image as `flag_variants`, see `compose.yaml`):

```bash
src/scripts/run_mave_dataset_stats.sh [condensed-file] [metadata-file] [--output path]
```

Unlike `run_flag_variants.sh`, this wrapper doesn't map paths against a
`/work` staging mount -- the script only reads locally-generated
`data/output/maves/`/`data/output/supplementary_data/` files, which are
already available at `/usr/src/app` via the whole-repo bind mount.

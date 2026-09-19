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
- Number of distinct protein variants assayed directly (distinct `hgvs_p`
  among condensed-file rows reported at protein resolution,
  `nucleotide_or_aa == "aa"` -- see note below)
- Number of distinct DNA variants assayed directly (distinct `hgvs_c` among
  condensed-file rows reported at DNA resolution, `nucleotide_or_aa == "nt"`)
- Number of distinct variants assayed, total (distinct `(hgvs_c, hgvs_p)`
  pairs in the condensed file, across all rows regardless of resolution)
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
and DNA-level variants. The clinical-attribute breakdown uses ClinVar 2025
for every gene except BRCA1, PTEN, MSH2, and TP53, which use ClinVar 2018
instead (see `MIXED_YEAR_GENES` / `mixed_year_clinvar_series` in the
script).

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
  sides -- see `excalibr_dataset_to_gene_map` in the script for why. Each line
  reports a second number in parentheses, excluding F9, TP53, and SFPQ
  (`EXCALIBR_EXCLUDED_GENES`) -- F9 and TP53 use OddsPath rather than ExCALIBR
  calibration (see `docs/variant_classification.md`).
- **Filtering effects on the reclassification dataset**: a sequential funnel
  from every DNA-level measurement row in the expanded file down to the
  reclassification export, applying -- in the pipeline's actual order --
  every exclusion `Variant_Classification_analysis.ipynb` and
  `src/build_variant_reclassification_dataset.py` apply first. The funnel
  starts with two rows over that same unfiltered file: "including RNA
  scores" counts each `rna_score`-bearing measurement as an *additional*
  measurement on top of the usual one-per-`mavedb_variant_urn` count, and
  "- RNA scores" drops back to that plain count (the baseline every later
  step uses). Each step's
  effect is reported two ways: distinct DNA variants (the reclassification
  pipeline's own dedup key, `Gene`/`Chrom`/`hg38_start`/`ref_allele`/
  `alt_allele`) and distinct assayed variants -- *not* distinct
  `mavedb_variant_urn` (that identifies one MaveDB score-set record, one per
  dataset, so a variant assayed by two datasets would count twice); instead
  each row's urn is mapped back to its parent condensed-file row's
  `(hgvs_c, hgvs_p)` pair, the same key `distinct_variants_assayed` uses, so
  this matches that figure exactly over the unfiltered file. An assayed
  variant is only counted as filtered out once every one of its DNA-level
  candidates has been. Two summary lines at the end report the final count
  reclassified out of the starting count, at each resolution. See
  `compute_reclassification_filter_funnel`'s docstring in the script for the
  exact steps, and the two extra inputs this section needs:
  - `--checkpoint-file` (default
    `data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`):
    the notebook's own checkpoint, saved just before its category split --
    already has the F9/TP53 restricted-dataset filter baked in, which this
    section reproduces separately (against the expanded file) purely to
    report its individual effect. LDLR's LA-module exclusion is a
    `Flag == '*'` set upstream by `flag_variants.py`, so it isn't baked into
    the checkpoint as a row-drop -- it survives into the checkpoint (still
    present, just flagged) and is folded into the "Other flagged variants"
    step below like any other pre-existing flag.
  - `--chek2-file` (default `data/input/maves/CHEK2_Gebbia_2024.xlsx`): same
    CHEK2 QC workbook `build_variant_reclassification_dataset.py` uses.
- **Reclassification agreement** (Figure 4c): for every sheet in
  `--controls-file` (default `data/output/supplementary_data/Supplementary_Data_5.xlsx`)
  whose name starts with `controls_` -- one per predictor/calibration
  combination, which should agree with each other since this doesn't depend
  on which predictor is active -- how often ExCALIBR's evidence assignment
  (`ExC_points_2025`) and the functional class assignment (`OP_points`) agree
  with the row's ClinVar pathogenic-or-benign control label
  (`clnsig_group_18_25`).
- **Control concordance** (ClinVar vs. ClinGen; OddsPath alone vs. combined
  with each of REVEL/AlphaMissense/MutPred2): for the ClinVar
  (`controls_*_GeneSpecific`) and ClinGen Evidence Repository
  (`ClinGen_Repo_*_GeneSpecific`) control sets separately, how many variants
  (and what percent) are **concordant** (evidence and the control
  classification agree on pathogenic vs. benign), **discordant** (evidence
  and the control classification disagree), or classified **VUS** (no
  determinate call from that evidence source) -- reported once for OddsPath
  calibration evidence alone (`OP_points` sign, read from the REVEL sheet)
  and once each for the combined ExCALIBR/OddsPath + REVEL/AlphaMissense/
  MutPred2 gene-specific evidence (each predictor's own precomputed
  `Class_REVEL`/`Class_AM`/`Class_MP2` column), rendered as one table per
  control source with a row per evidence source so all four can be compared
  at a glance. **Discordant** is further split into its two directions --
  control Pathogenic/Likely Pathogenic reclassified Benign/Likely Benign by
  the evidence source, and the reverse -- which sum back to the Discordant
  count. ClinVar's control label is `clnsig_group_18_25`; ClinGen's is
  `Updated_Classification_ClinGen_repo` (its own P/LP/B/LB assertion --
  `clnsig_group_18_25` isn't a clean ClinVar label for these rows, since a
  ClinGen Evidence Repository control need not have an unambiguous ClinVar
  entry of its own). The ClinGen table additionally reports each row's
  `Genes` (distinct genes among that row's in-scope ClinGen control
  variants) and its own `PLP`/`BLB` population split (ClinGen's
  classification alone, independent of the evidence source, summing to
  `Total`) -- not shown for ClinVar, whose much larger, more stable control
  set doesn't need this called out per row. See
  `compute_control_concordance`'s docstring in the script.

  Immediately after that table, a second, missense-only table repeats the
  same breakdown -- ClinVar and ClinGen subsections alike, ClinGen's
  `Genes`/`PLP`/`BLB` columns included -- with every sheet first restricted
  to `simplified_consequence == "missense_variant"` rows before scoring.
  See `MISSENSE_CONTROL_CONCORDANCE_SOURCES`/
  `MISSENSE_CONSEQUENCE_VALUE` in the script.
- **Variant classification**: how many distinct DNA variants have a
  classification, how many of those are pathogenic or benign, and how many
  ClinVar VUS / unobserved variants are "resolved" -- reclassified/classified
  pathogenic or benign -- and what percent of that category that is:
  (re)classification of DNA variants using functional evidence points
  (ExCALIBR for most genes, OddsPath for F9/TP53 -- see
  `docs/variant_classification.md`) plus REVEL predictor evidence,
  gene-specific calibration falling back to genome-wide. Reads the
  `--controls-file`'s `controls_REVEL_GeneSpecific`,
  `ClinGen_Repo_REVEL_GeneSpecific`, `VUS_REVEL`, `gnomAD_REVEL`, and
  `Unobserved_REVEL` sheets' precomputed `Class_REVEL` column, deduplicated to
  one row per distinct DNA variant (by `Gene`/`Chrom`/`hg38_start`/
  `ref_allele`/`alt_allele`) -- restricted to variants falling into one of
  those five categories, and (for `controls`/`ClinGen_Repo`) deduplicated by
  that category's own DNA-resolution-preferred rule (see
  `docs/variant_classification.md#decided-approach`).
- **Chi-squared tests**: for each predictor, whether the Pathogenic/Likely
  Pathogenic and Benign/Likely Benign rates cited in the manuscript's Fig.
  6d/e paragraphs actually differ between gnomAD, ClinVar VUS, and
  Unobserved variants: gnomAD vs. ClinVar VUS (both rates), and Unobserved
  vs. ClinVar VUS and vs. gnomAD (Pathogenic/Likely Pathogenic rate only --
  the two comparisons the manuscript draws for the Unobserved category).
  Each comparison is a 2x2 Pearson's chi-squared test of independence with
  Yates' continuity correction (R's `chisq.test()` default for a 2x2 table,
  equivalent to `prop.test(..., correct = TRUE)`) run on the two groups'
  raw counts/totals from the variant classification section's own stats --
  the rate over the full category, not just the "resolved"
  (pathogenic-or-benign) subset. See
  `compute_variant_classification_chi_squared_tests` in the script.
- **Gene-level discordance**: the top 5 genes by number of control variants
  where the `controls_REVEL_GeneSpecific` sheet's `Class_REVEL`
  classification disagrees with the row's ClinVar pathogenic-or-benign
  control label (`clnsig_group_18_25`) -- i.e. a ClinVar Pathogenic/Likely
  pathogenic control classified Benign/Likely benign, or a ClinVar
  Benign/Likely benign control classified Pathogenic/Likely Pathogenic.
  Deduplicated to one row per distinct DNA variant the same way as the
  variant classification section above. Each of the top 5 genes is broken
  down by direction of discordance (ClinVar P/LP reclassified B/LB, vs.
  ClinVar B/LB reclassified P/LP). See `compute_gene_discordance_stats` in
  the script.

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
  `OP_points` columns -- plus, for the Supplementary Data 5 variant
  classification section, the `controls_REVEL_GeneSpecific`,
  `ClinGen_Repo_REVEL_GeneSpecific`, `VUS_REVEL`, `gnomAD_REVEL`, and
  `Unobserved_REVEL` sheets, using their `Gene`, `Chrom`, `hg38_start`,
  `ref_allele`, `alt_allele`, and `Class_REVEL` columns.
- **Checkpoint file** (`--checkpoint-file`, default
  `data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`):
  `Variant_Classification_analysis.ipynb`'s pre-category-split checkpoint, for
  the filter-funnel section's post-checkpoint steps -- uses its `Gene`,
  `Chrom`, `hg38_start`, `ref_allele`, `alt_allele`, `mavedb_variant_urn`,
  `hgvs_p`, `auth_reported_score`, `Flag`, `VariantNotes`, `splice_var_amino`,
  `splice_measure` (a splice-aware dataset, `splice_measure == "Yes"`, is
  exempt from the `splice_var_amino` exclusion), and `revel_train_amino`
  columns.
- **CHEK2 QC workbook** (`--chek2-file`, default
  `data/input/maves/CHEK2_Gebbia_2024.xlsx`): same file
  `build_variant_reclassification_dataset.py` uses, for the filter-funnel
  section's CHEK2 QC-flag step -- uses its `hgvs_pro`, `score`, and
  `Filter_CI` columns.

### Note on `(hgvs_g, hgvs_p)`

The integrated dataset has no `hgvs_g` (genomic HGVS) column. Its DNA-level
identifier is `hgvs_c` (transcript-relative HGVS, pipe-delimited per row when a
protein-resolution measurement corresponds to more than one underlying DNA
change). This script uses `hgvs_c` as that DNA-level key, so "distinct
variants assayed, total" counts distinct `(hgvs_c, hgvs_p)` pairs in the
condensed file. If a true genomic (`NC_...:g.`) identifier is later added to
the integrated dataset, swap `GENOMIC_VARIANT_COL` in the script to point at
it.

### Note on `nucleotide_or_aa`

`nucleotide_or_aa` records each row's original assay resolution, from the
variant-annotation pipeline's reverse-translation step (see
`docs/variant_annotation_pipeline.md`): `"nt"` when the source assay reported
a DNA-level variant directly, `"aa"` when it only reported a protein-level
variant and one or more DNA candidates were reverse-translated into `hgvs_c`
(the pipe-delimited case above). This is why the protein- and DNA-level
distinct-variant counts don't sum to the total distinct `(hgvs_c, hgvs_p)`
count above -- they're a resolution-based partition of the same rows, not an
independent count.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.mave_dataset_stats \
  data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz \
  data/output/supplementary_data/Supplementary_Data_3.xlsx \
  data/output/maves/integrated_variant_effect_dataset.tsv.gz \
  [--excalibr-calibrations-file data/output/supplementary_data/Supplementary_Data_4.xlsx] \
  [--controls-file data/output/supplementary_data/Supplementary_Data_5.xlsx] \
  [--checkpoint-file data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz] \
  [--chek2-file data/input/maves/CHEK2_Gebbia_2024.xlsx] \
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

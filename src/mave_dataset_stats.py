#!/usr/bin/env python3
"""Summary statistics for the integrated MAVE variant effect dataset.

Reads the condensed variant effect dataset (one row per variant effect
measurement/score; see data/output/maves/integrated_variant_effect_dataset*.tsv.gz)
together with its dataset-level metadata (Supplementary_Data_3.xlsx) and
reports, for three groupings of datasets -- IGVF-produced only, non-IGVF
("community") only, and combined -- the number of datasets, variant effect
measurements, RNA scores, composite scores, distinct variants assayed, and
genes represented. `rna_scores` is a breakdown of `variant_effect_measurements`
(every row with a non-empty `rna_score` also carries a regular
`auth_reported_score`), not an addition to it, reported beside it for
visibility. Distinct variants assayed is reported three ways: distinct
protein variants assayed directly (distinct `hgvs_p` among rows reported at
protein resolution, `nucleotide_or_aa == "aa"`), distinct DNA variants
assayed directly (distinct `hgvs_c` among rows reported at DNA resolution,
`nucleotide_or_aa == "nt"`), and the total distinct (`hgvs_c`, `hgvs_p`)
pairs across all rows regardless of resolution -- see `compute_bucket_stats`.
Each bucket's variant effect measurements are also reported as a percentage
of the combined total (100% for the combined row), both on their own
(`pct_variant_effect_measurements`) and combined with RNA scores
(`pct_variant_effect_measurements_with_rna_scores`, out of the combined total
of measurements plus RNA scores). Immediately after that table, a "Genes
represented" section lists which genes are covered only by IGVF datasets,
only by non-IGVF ("community") datasets, or by both. CALM1/CALM2/CALM3 are
always counted and listed as a single "CALM1/2/3" gene, both in this table's
`genes_represented`/`genes_not_in_igvf_data` counts and in the "Genes
represented" list below, since treating the same calmodulin target as three
separate genes isn't useful.

Next, an "IGVF-produced datasets" table breaks the Dataset summary's IGVF
`variant_effect_measurements`, `rna_scores`, and `composite_scores` counts out
per dataset -- one row per IGVF-produced dataset, sorted by
`variant_effect_measurements` descending -- see
`compute_igvf_dataset_measurement_counts`.

It additionally reports, in a set of text tables, how many variants and
variant measurements have REVEL, AlphaMissense, and MutPred2 scores, and how
many fall into each of several clinical-attribute buckets (ClinVar VUS,
ClinVar pathogenic/benign, observed in gnomAD, or none of the above). This
clinical-attribute breakdown is reported twice: once using ClinVar 2025 for
every gene, and once using ClinVar 2025 for every gene except BRCA1, PTEN,
MSH2, and TP53, which use ClinVar 2018 instead. Each is reported four ways: at
the assayed-variant (protein-resolution, from the condensed file) or
DNA-variant (from the expanded file) level, and as distinct variants or as
variant measurements (i.e. rows). Each clinical-attribute table also carries
two extra columns breaking out, for each of the four buckets, how many of the
level's SNVs are in it (or, at the assayed-variant level, "SNV-accessible"
variants: those with at least one single-nucleotide-substitution candidate
among the DNA-level changes -- `hgvs_c`/`transcript_ref`/`transcript_alt` --
that reverse-translate to that protein change), as a percentage of the
level's total SNV(-accessible) count (also given in that table's totals
line).

By default, a variant with a conflicting or ambiguous ClinVar call --
disagreement between a pathogenic-leaning and benign-leaning classification
across its measurements/DNA candidates, or ClinVar's own "Conflicting
classifications of pathogenicity" call -- is excluded from both the VUS and
pathogenic-or-benign buckets and counted in its own "ClinVar conflict" bucket
instead, mirroring the conflict handling in
Analysis/Curation_summary_V5_cleaned.ipynb. Pass `--allow-clinvar-conflicts`
to instead fold such variants into VUS/pathogenic-or-benign via any-match
(the original behavior), which also drops the conflict bucket from the
report.

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

Note on pipe-delimited annotation columns (REVEL, AM_score, MutPred2,
clinvar_sig_2025, gnomad_MAF): in the condensed file these follow the same
one-part-per-underlying-DNA-change convention as `hgvs_c`, and the parts don't
always agree (e.g. one DNA-level candidate has a REVEL score and another
doesn't). An assayed variant (or measurement row) counts as having a value, or
matching a target ClinVar significance, if *any* one of its pipe-delimited
parts does. In the expanded (DNA-variant-level) file there's only ever one
part, so this reduces to a plain presence/absence check there. "Observed in
gnomAD" follows this repo's existing convention elsewhere (see
notebooks/analysis/README_OddsPath_classifications.md) of treating any non-empty
`gnomad_MAF` as observed, regardless of the allele count.

The Dataset summary and Genes represented sections above always count
CALM1/CALM2/CALM3 as one gene target (see above). Elsewhere in the report --
the Meta-analyses gene column and the ExCALIBR calibration coverage section --
they're counted as three separate genes by default, since that's how the
underlying data labels them; pass `--merge-calm-genes` to count them as one
there too, since they encode the same calmodulin protein.

Two further sections cover ExCALIBR calibration coverage and reclassification
agreement, sourced from `--excalibr-calibrations-file` (default
Supplementary_Data_4.xlsx) and `--controls-file` (default
Supplementary_Data_5.xlsx) respectively:

- **ExCALIBR calibration coverage** (Extended Data Figure 2): how many genes
  have a row in Supplementary_Data_4's `ExCALIBR_calibrations` sheet, and how
  many of those genes have at least one dataset where ExCALIBR assigned at
  least one point of evidence in either direction (i.e. at least one of that
  row's `range_-8`..`range_8` columns is non-null). That sheet's `dataset`
  values are matched against this script's own dataset metadata (the same
  `Dataset Name` -> `Gene` mapping used elsewhere) to determine gene
  membership, after stripping a trailing `_clinvar_2018` suffix (this sheet's
  convention for marking a mixed-year duplicate calibration of the same
  dataset) and Unicode-normalizing (NFC) both sides, since accented
  characters can round-trip through Excel in different composed/decomposed
  forms across files. Any `dataset` value that still doesn't resolve raises
  `ValueError`, same as an unrecognized dataset in the main metadata file.

- **ClinGen Evidence Repository control set**: how many DNA variants have a
  determinate original ClinGen classification (`Assertion_ClinGen_repo`) and
  MAVE experimental data (after the same checkpoint-downstream exclusions as
  the filter funnel above, minus its training-variant split); how many
  retain a determinate classification after
  `src/recalculate_clingen_classification.py` recomputes it with functional
  (BS3/PS3) and predictive (BP4/PP3) evidence removed
  (`Updated_Classification_ClinGen_repo`); and how many are ultimately
  retained per predictor (the `ClinGen_Repo_*_GeneSpecific` sheets, after
  also excluding that predictor's own training variants -- these three
  totals are also each `Control concordance`'s "ClinGen" `Total` above). See
  `compute_clingen_evidence_repository_stats`.

- **Reclassification agreement** (Figure 4c): for every sheet in the controls
  file whose name starts with `controls_` (one per predictor/calibration
  combination -- REVEL, AlphaMissense, MutPred2 -- which should agree with
  each other since this doesn't depend on which predictor is active), how
  often ExCALIBR's evidence assignment (`ExC_points_2025`) and the functional
  class assignment (`OP_points`) agree with the row's ClinVar
  pathogenic-or-benign control label (`clnsig_group_18_25`). A row's points
  column is treated as assigning pathogenic evidence if positive, benign
  evidence if negative, and no evidence if zero or missing -- rows with no
  evidence assigned are excluded from the agreement percentage (reported
  separately as "no point of evidence assigned") but still counted in the
  section's total.

- **Control concordance**: for the ClinVar (`controls_*_GeneSpecific`) and
  ClinGen Evidence Repository (`ClinGen_Repo_*_GeneSpecific`) control sets
  separately, how many variants (and what percent) are concordant,
  discordant, or classified VUS -- once for OddsPath calibration evidence
  alone (`OP_points` sign, read once from the REVEL sheet since functional
  evidence is shared across predictors), and once each for the combined
  ExCALIBR/OddsPath + REVEL/AlphaMissense/MutPred2 gene-specific evidence
  (`Class_REVEL`/`Class_AM`/`Class_MP2`, each from that predictor's own
  sheet), rendered as one table per control source with a row per evidence
  source so all four can be compared at a glance. Discordant is further
  broken out by direction: control Pathogenic/Likely Pathogenic reclassified
  Benign/Likely Benign by the evidence source, vs. the reverse -- the two sum
  to the Discordant count. ClinVar's control label is `clnsig_group_18_25`;
  ClinGen's is `Updated_Classification_ClinGen_repo` (its own assertion,
  since `clnsig_group_18_25` isn't a clean ClinVar label for ClinGen-only
  controls). See `compute_control_concordance`.

- **Variant classification**: how many distinct DNA variants have a
  classification, how many of those are pathogenic or benign, and how many
  ClinVar VUS / unobserved variants are "resolved" -- reclassified/classified
  pathogenic or benign -- and what percent of that category that is, further
  split into how many (and what percent of the category) are
  Pathogenic/Likely Pathogenic vs. Benign/Likely Benign, and, of the
  Benign/Likely Benign side, how many (and what percent) have only the
  minimum possible evidence, -1 point: (re)classification of DNA variants
  using functional evidence points (ExCALIBR for most genes, OddsPath for
  F9/TP53 -- see `docs/variant_classification.md`) plus REVEL predictor
  evidence, gene-specific calibration falling back to genome-wide, from the
  controls file's `controls_REVEL_GeneSpecific`, `ClinGen_Repo_REVEL_GeneSpecific`,
  `VUS_REVEL`, `gnomAD_REVEL`, and `Unobserved_REVEL` sheets' precomputed
  `Class_REVEL` column, deduplicated to one row per distinct DNA variant
  (`VARIANT_CLASSIFICATION_COORD_COLS`) -- restricted to variants falling
  into one of those five categories, and (for `controls`/`ClinGen_Repo`)
  deduplicated by that category's own DNA-resolution-preferred rule (see
  `docs/variant_classification.md#decided-approach`).

Immediately after that section, a "Chi-squared tests" section reports, for each predictor, whether the
Pathogenic/Likely Pathogenic and Benign/Likely Benign rates cited in the
manuscript's Fig. 6d/e paragraphs actually differ between gnomAD, ClinVar
VUS, and Unobserved variants: gnomAD vs. ClinVar VUS (Pathogenic/Likely
Pathogenic rate and, separately, Benign/Likely Benign rate), and Unobserved
vs. ClinVar VUS and vs. gnomAD (Pathogenic/Likely Pathogenic rate only, the
two comparisons the manuscript draws for the Unobserved category). Each
comparison is a 2x2 Pearson's chi-squared test of independence with Yates'
continuity correction (R's `chisq.test()` default for a 2x2 table,
equivalent to `prop.test(..., correct = TRUE)`) between the two groups'
raw counts/totals from the Supplementary Data 5 version's own
`{category}_total`/`{category}_resolved_pathogenic`/`{category}_resolved_
benign` stats -- the rate over the full category, not just the "resolved"
(pathogenic-or-benign) subset. See
`compute_variant_classification_chi_squared_tests`.

A final section, "Gene-level discordance", reports the genes with the most
control variants where the ExCALIBR/OddsPath + REVEL gene-specific
classification (`Class_REVEL`, from the controls file's
`controls_REVEL_GeneSpecific` sheet, deduplicated to one row per distinct DNA
variant) disagrees with the row's ClinVar pathogenic-or-benign control label
(`clnsig_group_18_25`) -- i.e. calls a ClinVar Pathogenic/Likely pathogenic
control Benign/Likely benign, or vice versa. It shows the top 5 genes by
total discordant-variant count, each broken down by direction of discordance
(ClinVar P/LP reclassified B/LB, vs. ClinVar B/LB reclassified P/LP). See
`compute_gene_discordance_stats`.

A "Simplified consequence x SpliceAI score breakdown" section reports, for
each of the five Supplementary Data 5 variant categories (from the same
REVEL sheets as the Supplementary Data 5 variant classification section --
`controls_REVEL_GeneSpecific`, `ClinGen_Repo_REVEL_GeneSpecific`,
`VUS_REVEL`, `gnomAD_REVEL`, and `Unobserved_REVEL`, each deduplicated to
one row per distinct DNA variant the same way), a table with one row per
`simplified_consequence` value and one column per category, each cell
showing that consequence/category combination's distinct DNA variant count
split two ways: those with every SpliceAI delta score (`spliceAI_DS_AG`,
`spliceAI_DS_AL`, `spliceAI_DS_DG`, `spliceAI_DS_DL`) null or below 0.2, and
those with at least one of those four scores at or above 0.2. See
`compute_consequence_splice_breakdown`.

Both file arguments are optional and default to the paths above. Output is
written as plain text (to stdout, and optionally to `--output` as well).

A "Filtering effects on the reclassification dataset" section, immediately
after the ExCALIBR calibration coverage section, shows a sequential funnel:
starting from every DNA-level measurement row in the expanded file, it applies
-- in the pipeline's actual order -- every exclusion
`Variant_Classification_analysis.ipynb` and `src/build_variant_reclassification_dataset.py`
apply before a variant reaches the reclassification export, reporting each
step's effect two ways: distinct DNA variants (the reclassification
pipeline's own dedup key, `Gene`/`Chrom`/
`hg38_start`/`ref_allele`/`alt_allele`) and distinct assayed variants (each
row's `mavedb_variant_urn` mapped back to its parent condensed-file row's
`(hgvs_c, hgvs_p)` pair -- the same key `distinct_variants_assayed` uses, so
this matches that figure exactly over the unfiltered file -- an assayed
variant is only counted as filtered out once every one of its DNA-level
candidates has been). See `compute_reclassification_filter_funnel`'s
docstring for the exact steps and `--checkpoint-file`/`--chek2-file` for the
two extra inputs this section needs.
"""

import unicodedata
from pathlib import Path

import click
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

DEFAULT_CONDENSED_FILE = Path("data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz")
DEFAULT_EXPANDED_FILE = Path("data/output/maves/integrated_variant_effect_dataset.tsv.gz")
DEFAULT_METADATA_FILE = Path("data/output/supplementary_data/Supplementary_Data_3.xlsx")
DEFAULT_EXCALIBR_CALIBRATIONS_FILE = Path("data/output/supplementary_data/Supplementary_Data_4.xlsx")
DEFAULT_CONTROLS_FILE = Path("data/output/supplementary_data/Supplementary_Data_5.xlsx")

METADATA_SHEET = "Curation"
DATASET_COL = "Dataset"
GENE_COL = "Gene"
GENOMIC_VARIANT_COL = "hgvs_c"
PROTEIN_VARIANT_COL = "hgvs_p"
VARIANT_KEY_COLS = [GENOMIC_VARIANT_COL, PROTEIN_VARIANT_COL]
MAVEDB_VARIANT_URN_COL = "mavedb_variant_urn"

NUCLEOTIDE_OR_AA_COL = "nucleotide_or_aa"
AA_LEVEL_VALUE = "aa"
NT_LEVEL_VALUE = "nt"

RNA_SCORE_COL = "rna_score"

METADATA_DATASET_COL = "Dataset Name"
IGVF_PRODUCED_COL = "IGVF Produced?"
SCORE_SET_TYPE_COL = "Primary Score Set or Meta-analysis?"
MEASUREMENT_VALUE = "primary score set"

CALM_GENES = frozenset({"CALM1", "CALM2", "CALM3"})
CALM_MERGED_LABEL = "CALM1/2/3"

SCORE_COLUMNS = {
    "REVEL": "REVEL",
    "AlphaMissense": "AM_score",
    "MutPred2": "MutPred2",
}
CLINVAR_COL = "clinvar_sig_2025"
CLINVAR_COL_2018 = "clinvar_sig_2018"
MIXED_YEAR_GENES = frozenset({"BRCA1", "PTEN", "MSH2", "TP53"})
GNOMAD_COL = "gnomad_MAF"
TRANSCRIPT_REF_COL = "transcript_ref"
TRANSCRIPT_ALT_COL = "transcript_alt"
VUS_LABEL = "VUS (ClinVar 2025)"
PATHOGENIC_OR_BENIGN_LABEL = "Pathogenic or benign (ClinVar 2025)"
CLINVAR_CONFLICT_LABEL = "ClinVar conflict (excluded from VUS/pathogenic-or-benign)"
GNOMAD_LABEL = "Observed in gnomAD"
NO_ANNOTATION_LABEL = "No ClinVar or gnomAD annotation"
SNV_LABEL = "SNV"
SNV_ACCESSIBLE_LABEL = "SNV-accessible"
VUS_VALUES = frozenset({"Uncertain significance"})
PATHOGENIC_VALUES = frozenset({"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"})
BENIGN_VALUES = frozenset({"Benign", "Likely benign", "Benign/Likely benign"})
PATHOGENIC_OR_BENIGN_VALUES = PATHOGENIC_VALUES | BENIGN_VALUES
CONFLICTING_CLASSIFICATION_VALUE = "Conflicting classifications of pathogenicity"

EXCALIBR_CALIBRATIONS_SHEET = "ExCALIBR_calibrations"
EXCALIBR_DATASET_COL = "dataset"
EXCALIBR_RANGE_COLUMN_PREFIX = "range_"
EXCALIBR_MIXED_YEAR_DATASET_SUFFIX = "_clinvar_2018"
EXCALIBR_EXCLUDED_GENES = frozenset({"F9", "TP53", "SFPQ"})

# --- Reclassification filter funnel -----------------------------------------
# Pre-checkpoint exclusions from Variant_Classification_analysis.ipynb (cells
# 3, 4, and 42), reproduced here against the raw expanded file since the
# notebook's own checkpoint (DEFAULT_CHECKPOINT_FILE below) already has them
# baked in. See compute_reclassification_filter_funnel's docstring.
DEFAULT_CHECKPOINT_FILE = Path("data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz")
DEFAULT_CHEK2_FILE = Path("data/input/maves/CHEK2_Gebbia_2024.xlsx")

FUNNEL_GENOMIC_KEY_COLS = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele"]
LDLR_VLDL_PRIORITY_RANGES = [(66, 106), (234, 272)]  # LA module 2, LA module 6
LDLR_VLDL_PRIORITY_DATASET = "LDLR_Tabet_2025_presence_VLDL"
LDLR_VLDL_FALLBACK_DATASETS = frozenset({"LDLR_Tabet_2025_uptake", "LDLR_Tabet_2025_abundance"})
LDLR_LA1_EXCLUDE_RANGE = (25, 65)  # LA module 1
F9_TP53_RESTRICTED_DATASETS = frozenset(
    {
        "TP53_Boettcher_2019",
        "TP53_Fortuno_2021",
        "TP53_Giacomelli_2018_combined_score",
        "TP53_Giacomelli_2018_p53WT_Nutlin3",
        "TP53_Giacomelli_2018_p53null_Nutlin3",
        "TP53_Giacomelli_2018_p53null_etoposide",
        "TP53_Kato_2003_AIP1nWT",
        "TP53_Kato_2003_BAXnWT",
        "TP53_Kato_2003_GADD45nWT",
        "TP53_Kato_2003_MDM2nWT",
        "TP53_Kato_2003_NOXAnWT",
        "TP53_Kato_2003_P53R2nWT",
        "TP53_Kato_2003_WAF1nWT",
        "TP53_Kato_2003_h1433snWT",
        "F9_Popp_2025_carboxy_F9_specific",
        "F9_Popp_2025_carboxy_gla_motif",
        "F9_Popp_2025_heavy_chain",
        "F9_Popp_2025_light_chain",
        "F9_Popp_2025_strep_2",
    }
)

# Post-checkpoint exclusions -- mirrors
# src.build_variant_reclassification_dataset.apply_notebook_exclusions, split
# into named sub-steps so each one's effect can be measured separately. Kept
# as its own copy (rather than importing that function directly) so this
# report can show the funnel's intermediate steps, which that function
# doesn't expose.
# `VariantNotes` can carry more than one of these tags at once, `;`-joined
# (e.g. "splice_variant_not_measured;conflicting_fxn_data") -- see
# docs/variant_classification.md#variantnotes-values. Each individual tag
# below is its own funnel step, applied in this order, so a combined-tag row
# is removed at whichever of its tags' steps comes first.
CONFLICTING_FXN_DATA_TAG = "conflicting_fxn_data"
SPLICE_VARIANT_NOT_MEASURED_TAG = "splice_variant_not_measured"
START_LOST_VARIANT_NOT_MEASURED_TAG = "start_lost_variant_not_measured"

# The REVEL- and MutPred2-training exclusions are alternative endpoints of the
# funnel, not a chain: both apply to the same "Other flagged variants" survivor
# set independently (a variant can be a training variant for one predictor and
# not the other), so neither is applied on top of the other -- see
# compute_reclassification_filter_funnel's tail.
REVEL_TRAINING_STEP_LABEL = "- REVEL-training variant (revel_train_amino == 'Yes')"
MUTPRED2_TRAINING_STEP_LABEL = "- MutPred2-training variant (mp2_train_amino == 'Yes')"

CONTROLS_SHEET_PREFIX = "controls_"
CONTROLS_CLINVAR_GROUP_COL = "clnsig_group_18_25"
EXCALIBR_EVIDENCE_POINTS_COL = "ExC_points_2025"
FUNCTIONAL_CLASS_POINTS_COL = "OP_points"
RECLASSIFICATION_POINTS_COLUMNS = {
    "ExCALIBR evidence": EXCALIBR_EVIDENCE_POINTS_COL,
    "Functional class": FUNCTIONAL_CLASS_POINTS_COL,
}
AGREE_LABEL = "Agrees with ClinVar"
DISAGREE_LABEL = "Disagrees with ClinVar"
NO_EVIDENCE_LABEL = "No point of evidence assigned"
RECLASSIFICATION_LABELS_ORDER = [AGREE_LABEL, DISAGREE_LABEL, NO_EVIDENCE_LABEL]

VARIANT_CLASSIFICATION_COORD_COLS = ["Gene", "Chrom", "hg38_start", "ref_allele", "alt_allele"]
VARIANT_CLASSIFICATION_CATEGORY_SHEETS = {
    "controls": "controls_REVEL_GeneSpecific",
    "ClinGen_Repo": "ClinGen_Repo_REVEL_GeneSpecific",
    "VUS": "VUS_REVEL",
    "gnomAD": "gnomAD_REVEL",
    "Unobserved": "Unobserved_REVEL",
}
CLASS_REVEL_COL = "Class_REVEL"
TOTAL_POINTS_REVEL_COL = "Total_Points_REVEL"
CLASS_PATHOGENIC_VALUES = frozenset({"Pathogenic", "Likely Pathogenic"})
CLASS_BENIGN_VALUES = frozenset({"Benign", "Likely Benign"})
CLASS_PATHOGENIC_OR_BENIGN_VALUES = CLASS_PATHOGENIC_VALUES | CLASS_BENIGN_VALUES

SIMPLIFIED_CONSEQUENCE_COL = "simplified_consequence"
NO_CONSEQUENCE_LABEL = "(no consequence)"
SPLICEAI_SCORE_COLS = ["spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"]
SPLICEAI_SCORE_THRESHOLD = 0.2
SPLICEAI_LOW_LABEL = f"<{SPLICEAI_SCORE_THRESHOLD} or missing"
SPLICEAI_HIGH_LABEL = f">={SPLICEAI_SCORE_THRESHOLD}"
CONSEQUENCE_SPLICE_BREAKDOWN_TITLE = (
    "=== Simplified consequence x SpliceAI score breakdown (Supplementary Data 5 REVEL sheets) ==="
)

# Per-predictor sheet-name suffix and Class_*/Total_Points_* column names for
# `compute_variant_classification_stats`'s REVEL/AlphaMissense/MutPred2 table
# -- each predictor's own sheet excludes *that* predictor's own training
# variants (e.g. `VUS_REVEL` drops REVEL-trained variants but not AM/MP2-
# trained ones), so unlike `VARIANT_CLASSIFICATION_CATEGORY_SHEETS` above
# (REVEL only, reused by the gene-discordance/control-concordance sections),
# each predictor here reads its *own* dedicated sheet rather than reusing
# REVEL's -- the three sheets' row counts genuinely differ (e.g. `VUS_REVEL`
# 17,788 vs. `VUS_AM` 17,851 vs. `VUS_MP2` 17,364 on the real workbook).
VARIANT_CLASSIFICATION_PREDICTORS = ["REVEL", "AlphaMissense", "MutPred2"]
VARIANT_CLASSIFICATION_SHEET_SUFFIX = {"REVEL": "REVEL", "AlphaMissense": "AM", "MutPred2": "MP2"}
VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR = {
    predictor: {
        "controls": f"controls_{suffix}_GeneSpecific",
        "ClinGen_Repo": f"ClinGen_Repo_{suffix}_GeneSpecific",
        "VUS": f"VUS_{suffix}",
        "gnomAD": f"gnomAD_{suffix}",
        "Unobserved": f"Unobserved_{suffix}",
    }
    for predictor, suffix in VARIANT_CLASSIFICATION_SHEET_SUFFIX.items()
}
VARIANT_CLASSIFICATION_CLASS_COL_BY_PREDICTOR = {
    "REVEL": "Class_REVEL",
    "AlphaMissense": "Class_AM",
    "MutPred2": "Class_MP2",
}
VARIANT_CLASSIFICATION_POINTS_COL_BY_PREDICTOR = {
    "REVEL": "Total_Points_REVEL",
    "AlphaMissense": "Total_Points_AM",
    "MutPred2": "Total_Points_MP2",
}
# Fxn_points (functional evidence alone) is shared across predictors; each
# predictor has its own gene-specific/genome-wide predictor-points and
# Conflicting_* column (see docs/variant_classification.md#handling-conflicting-evidence).
# Total_Points_* == Fxn_points + Points_*_GeneSpecific_GenomeWide, so at
# LIKELY_BENIGN_POINTS_THRESHOLD (-1) exactly, a row is either "single
# source" (one side is exactly 0) or "conflicting" (opposite signs) --
# those two are mutually exclusive and exhaustive at that specific total.
FUNCTIONAL_POINTS_COL = "Fxn_points"
VARIANT_CLASSIFICATION_PREDICTOR_POINTS_COL_BY_PREDICTOR = {
    "REVEL": "Points_REVEL_GeneSpecific_GenomeWide",
    "AlphaMissense": "Points_AM_GeneSpecific_GenomeWide",
    "MutPred2": "Points_MP2_GeneSpecific_GenomeWide",
}
VARIANT_CLASSIFICATION_CONFLICTING_COL_BY_PREDICTOR = {
    "REVEL": "Conflicting_REVEL",
    "AlphaMissense": "Conflicting_AM",
    "MutPred2": "Conflicting_MP2",
}
CONFLICTING_EVIDENCE_VALUE = "Conflicting evidence"

# --- Control concordance (ClinVar vs. ClinGen, OddsPath alone vs. combined) --
# Reuses each predictor's own controls_*_GeneSpecific/ClinGen_Repo_*_GeneSpecific
# sheets (already read for the Supplementary Data 5 variant-classification
# section above) as the ClinVar/ClinGen control sets, respectively.
# "OddsPath calibration" evidence is FUNCTIONAL_CLASS_POINTS_COL (OP_points)
# alone, read once from REVEL's sheet since functional evidence is shared
# across predictors (see FUNCTIONAL_POINTS_COL above). "ExCALIBR/OddsPath +
# <predictor> gene-specific" evidence is each predictor's own sheet's
# precomputed Class_<predictor> column, which already combines functional +
# that predictor's gene-specific (falling back to genome-wide) evidence --
# see VARIANT_CLASSIFICATION_TITLE.
CLINGEN_CLASSIFICATION_COL = "Updated_Classification_ClinGen_repo"
# ClinGen's own original asserted classification, before evidence removal --
# see ASSERTION_CLINGEN_REPO_COL's use in compute_clingen_evidence_repository_stats.
ASSERTION_CLINGEN_REPO_COL = "Assertion_ClinGen_repo"
CONCORDANT_LABEL = "Concordant"
DISCORDANT_LABEL = "Discordant"
CONTROL_VUS_LABEL = "VUS"
# Directional breakdown of DISCORDANT_LABEL -- mutually exclusive and exhaustive
# over the discordant rows, since a row's control classification is either
# pathogenic- or benign-leaning (never both), and it disagrees with the
# evidence source's assignment in exactly one of these two directions.
DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL = "  ...control PLP, evidence BLB"
DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL = "  ...control BLB, evidence PLP"
CONTROL_CONCORDANCE_LABELS_ORDER = [
    CONCORDANT_LABEL,
    DISCORDANT_LABEL,
    DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL,
    DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL,
    CONTROL_VUS_LABEL,
]
CONTROL_CONCORDANCE_EVIDENCE_LABEL = "OddsPath calibration"
COMBINED_EVIDENCE_LABEL_BY_PREDICTOR = {
    predictor: f"ExCALIBR/OddsPath + {predictor} gene-specific" for predictor in VARIANT_CLASSIFICATION_PREDICTORS
}
COMBINED_REVEL_EVIDENCE_LABEL = COMBINED_EVIDENCE_LABEL_BY_PREDICTOR["REVEL"]
# {control source label: (category key into VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR,
#  control classification column, pathogenic values, benign values)}
CONTROL_CONCORDANCE_SOURCES = {
    "ClinVar": ("controls", CONTROLS_CLINVAR_GROUP_COL, PATHOGENIC_VALUES, BENIGN_VALUES),
    "ClinGen": ("ClinGen_Repo", CLINGEN_CLASSIFICATION_COL, CLASS_PATHOGENIC_VALUES, CLASS_BENIGN_VALUES),
}

VARIANT_CLASSIFICATION_TITLE = (
    "=== Variant classification (Supplementary Data 5; "
    "ExCALIBR/OddsPath + REVEL/AlphaMissense/MutPred2, gene-specific + genome-wide fallback) ==="
)

# --- Chi-squared tests on the Supplementary Data 5 variant-classification rates ---
# The four gnomAD-vs-ClinVar-VUS-vs-Unobserved PLP/BLB rate comparisons cited
# in the manuscript's Fig. 6d/e paragraphs, run per predictor against the
# same vus_total/gnomad_total/unobserved_total and _resolved_pathogenic/
# _resolved_benign counts `compute_variant_classification_stats` already
# computes -- see `compute_variant_classification_chi_squared_tests`.
VARIANT_CLASSIFICATION_CHI_SQUARED_TITLE = (
    "=== Chi-squared tests (Supplementary Data 5; PLP/BLB rate, "
    "gnomAD vs. ClinVar VUS vs. Unobserved; Fig. 6d/e) ==="
)
VARIANT_CLASSIFICATION_CHI_SQUARED_COMPARISONS = [
    {
        "label": "Pathogenic/Likely Pathogenic rate: gnomAD vs. ClinVar VUS",
        "count_key": "resolved_pathogenic",
        "group_a": ("gnomAD", "gnomad"),
        "group_b": ("ClinVar VUS", "vus"),
    },
    {
        "label": "Benign/Likely Benign rate: gnomAD vs. ClinVar VUS",
        "count_key": "resolved_benign",
        "group_a": ("gnomAD", "gnomad"),
        "group_b": ("ClinVar VUS", "vus"),
    },
    {
        "label": "Pathogenic/Likely Pathogenic rate: Unobserved vs. ClinVar VUS",
        "count_key": "resolved_pathogenic",
        "group_a": ("Unobserved", "unobserved"),
        "group_b": ("ClinVar VUS", "vus"),
    },
    {
        "label": "Pathogenic/Likely Pathogenic rate: Unobserved vs. gnomAD",
        "count_key": "resolved_pathogenic",
        "group_a": ("Unobserved", "unobserved"),
        "group_b": ("gnomAD", "gnomad"),
    },
]

LIKELY_PATHOGENIC_POINTS_THRESHOLD = 6
LIKELY_BENIGN_POINTS_THRESHOLD = -1
# Just short of LIKELY_PATHOGENIC_POINTS_THRESHOLD -- the strongest possible Uncertain
# (unresolved) calls, one or two points away from resolving Likely Pathogenic. Reported
# as an overlapping cross-cut of the unresolved evidence-source split (a near-pathogenic
# row can be concordant, discordant, experimental-only, or predictive-only), not a sixth
# mutually-exclusive bucket.
NEAR_PATHOGENIC_POINTS_VALUES = frozenset({4, 5})

# Column names shared with the reclassification export
# (`data/output/reclassification/integrated_variant_effect_reclassification.tsv.gz`,
# `src/build_variant_reclassification_dataset.py`'s output) -- reused by
# `src/ablation_variant_reclassification.py`, which reads that file directly.
RECLASSIFICATION_CLINVAR_COL = "clinvar_sig_2025"
RECLASSIFICATION_GNOMAD_COL = "gnomad_MAF"
RECLASSIFICATION_REF_COL = "ref_allele"
RECLASSIFICATION_ALT_COL = "alt_allele"


def points_are_pathogenic_or_benign(points):
    """True where combined evidence points fall outside the Uncertain range
    (0-5) -- i.e. Likely Pathogenic/Pathogenic (>=6) or Likely Benign/Benign
    (<=-1) -- per the point-to-class cutoffs documented in
    `docs/variant_classification.md` and `Variant_Classification_analysis.ipynb`.
    """
    return (points >= LIKELY_PATHOGENIC_POINTS_THRESHOLD) | (points <= LIKELY_BENIGN_POINTS_THRESHOLD)


def split_genes(gene_value):
    return [g.strip() for g in gene_value.split(",") if g.strip()]


def _pipe_parts(value):
    return [p.strip() for p in value.split("|")]


def has_any_value(series):
    """Row-wise: does any '|'-delimited part of this field hold a value?

    See the module docstring's note on pipe-delimited annotation columns.
    """
    return series.apply(lambda v: any(p != "" for p in _pipe_parts(v)))


def matches_any_value(series, targets):
    """Row-wise: does any '|'-delimited part of this field match one of `targets`?"""
    return series.apply(lambda v: any(p in targets for p in _pipe_parts(v)))


def clinvar_significance_flags(clinvar_series):
    """Row-wise pathogenic/benign/VUS/literal-conflict membership (see `matches_any_value`).

    Returns (has_pathogenic, has_benign, has_vus, has_literal_conflict) boolean Series.
    Each of these is OR-able across rows sharing a distinct-variant key (via
    `.groupby(...).any()`) before being passed to `clinvar_classification_from_flags`,
    since "does the union of parts across this group contain a pathogenic value"
    is equivalent to "does any row in the group have a pathogenic value".
    """
    return (
        matches_any_value(clinvar_series, PATHOGENIC_VALUES),
        matches_any_value(clinvar_series, BENIGN_VALUES),
        matches_any_value(clinvar_series, VUS_VALUES),
        matches_any_value(clinvar_series, {CONFLICTING_CLASSIFICATION_VALUE}),
    )


def clinvar_classification_from_flags(has_pathogenic, has_benign, has_vus, has_literal_conflict):
    """Combine pathogenic/benign/VUS/literal-conflict membership into mutually
    exclusive VUS/pathogenic-or-benign/conflict boolean flags.

    A conflict is either ClinVar's own "Conflicting classifications of
    pathogenicity" call, or disagreement between a pathogenic-leaning and
    benign-leaning call -- e.g. across several reverse-translated DNA
    candidates for one protein change, or (once the inputs are grouped via
    `.any()`) across several datasets/measurements of the same variant.
    Mirrors the conflict handling in Analysis/Curation_summary_V5_cleaned.ipynb.

    Returns (vus, pathogenic_or_benign, conflict) boolean Series.
    """
    conflict = has_literal_conflict | (has_pathogenic & has_benign)
    pathogenic_or_benign = ~conflict & (has_pathogenic | has_benign)
    vus = ~conflict & ~pathogenic_or_benign & has_vus
    return vus, pathogenic_or_benign, conflict


def is_snv_accessible(ref_series, alt_series):
    """Row-wise: does any '|'-delimited (transcript_ref, transcript_alt) pair form a SNV?

    A pair is a single-nucleotide substitution if both sides are exactly one
    base. In the expanded (DNA-variant-level) file there's only ever one
    pair, so this is simply "is this variant a SNV". In the condensed
    (assayed, protein-level) file, a row's `hgvs_c` can list several
    underlying DNA-level reverse-translation candidates in lockstep with
    `transcript_ref`/`transcript_alt`; this is True if at least one candidate
    is a SNV ("SNV-accessible").
    """

    def _any_snv(ref_value, alt_value):
        return any(len(ref) == 1 and len(alt) == 1 for ref, alt in zip(_pipe_parts(ref_value), _pipe_parts(alt_value)))

    return pd.Series([_any_snv(r, a) for r, a in zip(ref_series, alt_series)], index=ref_series.index)


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


def merge_calm_gene_names(genes):
    """Collapse CALM1/CALM2/CALM3 -- which the underlying MAVE data treats as
    three separate gene labels but which correspond to a single gene target
    (they encode the same calmodulin protein) -- into one `CALM_MERGED_LABEL`
    entry in `genes`.
    """
    if genes & CALM_GENES:
        genes = (genes - CALM_GENES) | {CALM_MERGED_LABEL}
    return genes


def genes_in(df, merge_calm_genes=False):
    """Return the set of genes represented in `df`.

    If `merge_calm_genes` is set, CALM1/CALM2/CALM3 are collapsed into one
    `CALM_MERGED_LABEL` entry -- see `merge_calm_gene_names`.
    """
    genes = set()
    for value in df[GENE_COL].unique():
        genes.update(split_genes(value))
    if merge_calm_genes:
        genes = merge_calm_gene_names(genes)
    return genes


def compute_bucket_stats(condensed, dataset_names, measurement_datasets):
    """Compute the summary stats for one bucket of datasets.

    `distinct_variants_assayed` counts distinct (hgvs_c, hgvs_p) pairs, as
    before. `distinct_protein_variants_assayed` and
    `distinct_dna_variants_assayed` break that total out by assay resolution
    (see `NUCLEOTIDE_OR_AA_COL`): distinct `hgvs_p` among rows whose variant
    was reported at protein resolution (`nucleotide_or_aa == "aa"`), and
    distinct `hgvs_c` among rows reported at DNA resolution
    (`nucleotide_or_aa == "nt"`), respectively.

    `rna_scores` counts measurement rows with a non-empty `rna_score` --
    every such row in the underlying data also carries a regular
    `auth_reported_score`, so this is a breakdown of
    `variant_effect_measurements`, not an addition to it.

    CALM1/CALM2/CALM3 are always collapsed into one gene (`CALM_MERGED_LABEL`)
    for `genes_represented`, matching the "Genes represented" list below,
    which lists the same target once rather than under three separate names.
    """
    sub = condensed[condensed[DATASET_COL].isin(dataset_names)]
    is_measurement_row = sub[DATASET_COL].isin(measurement_datasets)
    n_measurements = int(is_measurement_row.sum())
    n_composite = int((~is_measurement_row).sum())
    n_rna_scores = int((is_measurement_row & (sub[RNA_SCORE_COL] != "")).sum())
    n_variants = sub[VARIANT_KEY_COLS].drop_duplicates().shape[0]
    n_protein_variants = sub.loc[sub[NUCLEOTIDE_OR_AA_COL] == AA_LEVEL_VALUE, PROTEIN_VARIANT_COL].nunique()
    n_dna_variants = sub.loc[sub[NUCLEOTIDE_OR_AA_COL] == NT_LEVEL_VALUE, GENOMIC_VARIANT_COL].nunique()
    genes = merge_calm_gene_names(genes_in(sub))
    return {
        "datasets": sub[DATASET_COL].nunique(),
        "variant_effect_measurements": n_measurements,
        "rna_scores": n_rna_scores,
        "composite_scores": n_composite,
        "distinct_protein_variants_assayed": n_protein_variants,
        "distinct_dna_variants_assayed": n_dna_variants,
        "distinct_variants_assayed": n_variants,
        "genes_represented": len(genes),
    }, genes


def compute_all_stats_from_frame(condensed, metadata):
    """Compute all bucket stats given an already-loaded condensed frame and metadata."""
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

    combined_measurements = combined_stats["variant_effect_measurements"]
    combined_measurements_and_rna = combined_measurements + combined_stats["rna_scores"]
    for bucket_stats in (igvf_stats, non_igvf_stats, combined_stats):
        bucket_stats["pct_variant_effect_measurements"] = (
            100 * bucket_stats["variant_effect_measurements"] / combined_measurements
            if combined_measurements
            else float("nan")
        )
        bucket_stats["pct_variant_effect_measurements_with_rna_scores"] = (
            100
            * (bucket_stats["variant_effect_measurements"] + bucket_stats["rna_scores"])
            / combined_measurements_and_rna
            if combined_measurements_and_rna
            else float("nan")
        )

    stats = {
        "IGVF": igvf_stats,
        "Community (non-IGVF)": non_igvf_stats,
        "Combined (IGVF + community)": combined_stats,
    }
    gene_breakdown = compute_gene_breakdown(igvf_genes, non_igvf_genes)
    return stats, gene_breakdown


def compute_gene_breakdown(igvf_genes, non_igvf_genes):
    """Partition genes represented in the dataset by whether they're covered
    by IGVF-produced datasets, non-IGVF ("community") datasets, or both.

    CALM1/CALM2/CALM3 are always collapsed into one `CALM_MERGED_LABEL` entry
    here (see `merge_calm_gene_names`), matching `compute_bucket_stats`'s
    `genes_represented`/`genes_not_in_igvf_data` counts -- listing (or
    counting) the same calmodulin target three times under three different
    names isn't useful.

    Returns {label: sorted gene list}.
    """
    igvf_genes = merge_calm_gene_names(igvf_genes)
    non_igvf_genes = merge_calm_gene_names(non_igvf_genes)
    return {
        "IGVF only": sorted(igvf_genes - non_igvf_genes),
        "Community (non-IGVF) only": sorted(non_igvf_genes - igvf_genes),
        "Both IGVF and community (non-IGVF)": sorted(igvf_genes & non_igvf_genes),
    }


def format_gene_breakdown(gene_breakdown):
    lines = ["=== Genes represented ==="]
    for label, genes in gene_breakdown.items():
        lines.append(f"{label} ({len(genes)}): {', '.join(genes)}")
    return "\n".join(lines)


IGVF_DATASET_MEASUREMENT_COUNTS_TITLE = "=== IGVF-produced datasets: measurements, RNA scores, and composite scores ==="


def compute_igvf_dataset_measurement_counts(condensed, metadata):
    """Per-IGVF-dataset breakdown of `Dataset summary`'s
    `variant_effect_measurements`, `rna_scores`, and `composite_scores`
    columns (see `compute_bucket_stats`) -- one row per IGVF-produced dataset
    present in `condensed`.

    A dataset's rows are either all measurement rows or all composite-score
    rows (its `SCORE_SET_TYPE_COL`), so exactly one of
    `variant_effect_measurements`/`composite_scores` is nonzero per row;
    `rna_scores` is a breakdown of `variant_effect_measurements` (every row
    with a non-empty `rna_score` also has a regular `auth_reported_score`),
    not an addition to it.

    Returns a DataFrame with `Dataset`, `variant_effect_measurements`,
    `rna_scores`, and `composite_scores` columns, sorted by
    `variant_effect_measurements` descending (ties broken by dataset name).
    """
    igvf_datasets = set(metadata.index[metadata[IGVF_PRODUCED_COL].eq("Yes")]) & set(condensed[DATASET_COL].unique())
    measurement_datasets = set(metadata.index[metadata[SCORE_SET_TYPE_COL].eq(MEASUREMENT_VALUE)])

    rows = []
    for dataset in igvf_datasets:
        sub = condensed[condensed[DATASET_COL] == dataset]
        is_measurement = dataset in measurement_datasets
        n_rows = len(sub)
        rows.append(
            {
                "Dataset": dataset,
                "variant_effect_measurements": n_rows if is_measurement else 0,
                "rna_scores": int((sub[RNA_SCORE_COL] != "").sum()) if is_measurement else 0,
                "composite_scores": 0 if is_measurement else n_rows,
            }
        )
    table = pd.DataFrame(rows, columns=["Dataset", "variant_effect_measurements", "rna_scores", "composite_scores"])
    return table.sort_values(
        ["variant_effect_measurements", "Dataset"], ascending=[False, True], kind="stable"
    ).reset_index(drop=True)


def format_igvf_dataset_measurement_counts(table):
    lines = [IGVF_DATASET_MEASUREMENT_COUNTS_TITLE]
    if len(table):
        lines.append(table.to_string(index=False))
    return "\n".join(lines)


COMPOSITE_SCORE_DATASETS_TITLE = "=== Meta-analyses (composite scores) ==="


def compute_composite_score_datasets(condensed, metadata, merge_calm_genes=False):
    """Per-dataset breakdown of the composite-score rows counted in `Dataset
    summary`'s `composite_scores` column -- every dataset whose
    `SCORE_SET_TYPE_COL` isn't `MEASUREMENT_VALUE` ("primary score set"), i.e.
    every meta-analysis dataset present in `condensed`. The `Scores` column
    sums to the "Combined (IGVF + community)" bucket's `composite_scores`
    count from `compute_bucket_stats`, since a dataset's rows are either all
    measurement rows or all composite-score rows (see `compute_bucket_stats`)
    and this covers every non-measurement dataset.

    Returns a DataFrame with `Gene` (comma-joined, `merge_calm_genes` applied
    the same way as `genes_in`), `Dataset`, `IGVF / Community`, and `Scores`
    (row count in `condensed` for that dataset) columns, sorted by `Scores`
    descending (ties broken by dataset name).
    """
    is_measurement = metadata[SCORE_SET_TYPE_COL].eq(MEASUREMENT_VALUE)
    composite_datasets = set(metadata.index[~is_measurement]) & set(condensed[DATASET_COL].unique())

    rows = []
    for dataset in composite_datasets:
        sub = condensed[condensed[DATASET_COL] == dataset]
        genes = genes_in(sub, merge_calm_genes=merge_calm_genes)
        is_igvf = metadata.loc[dataset, IGVF_PRODUCED_COL] == "Yes"
        rows.append(
            {
                "Gene": ", ".join(sorted(genes)),
                "Dataset": dataset,
                "IGVF / Community": "IGVF" if is_igvf else "Community",
                "Scores": len(sub),
            }
        )
    table = pd.DataFrame(rows, columns=["Gene", "Dataset", "IGVF / Community", "Scores"])
    return table.sort_values(["Scores", "Dataset"], ascending=[False, True], kind="stable").reset_index(drop=True)


def format_composite_score_datasets(table):
    lines = [COMPOSITE_SCORE_DATASETS_TITLE, f"Total composite scores: {int(table['Scores'].sum())}"]
    if len(table):
        lines.append(table.to_string(index=False))
    return "\n".join(lines)


def compute_all_stats(condensed_path, metadata_path):
    condensed = pd.read_csv(condensed_path, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_path)
    return compute_all_stats_from_frame(condensed, metadata)


PCT_COLUMNS = ["pct_variant_effect_measurements", "pct_variant_effect_measurements_with_rna_scores"]


def stats_to_dataframe(stats):
    table = pd.DataFrame(stats).T.reindex(
        columns=[
            "datasets",
            "variant_effect_measurements",
            "rna_scores",
            "pct_variant_effect_measurements",
            "pct_variant_effect_measurements_with_rna_scores",
            "composite_scores",
            "distinct_protein_variants_assayed",
            "distinct_dna_variants_assayed",
            "distinct_variants_assayed",
            "genes_represented",
            "genes_not_in_igvf_data",
        ]
    )
    count_columns = [c for c in table.columns if c not in PCT_COLUMNS]
    table[count_columns] = table[count_columns].astype("Int64")
    table[PCT_COLUMNS] = table[PCT_COLUMNS].astype(float).round(1)
    return table


def format_genomic_variant_count(expanded_path, expanded):
    """One line reporting the row count of `expanded_path` (the DNA/genomic-
    resolution expanded file `expanded` was read from) -- every row is one
    genomic-variant measurement, so this is `len(expanded)`, i.e. `wc -l`
    minus the header line.

    Deliberately doesn't shell out to `wc -l`: this counts pandas-parsed
    rows instead, which is the only reliable way to count records in a
    field-quoted, potentially-multiline TSV -- `wc -l` counts newline
    characters, so a single quoted field containing an embedded newline
    would silently inflate its count relative to the true record count.
    """
    return f"Genomic variants (rows in {expanded_path}): {len(expanded)}"


def mixed_year_clinvar_series(df):
    """Per-row ClinVar significance: ClinVar 2018 for BRCA1/PTEN/MSH2/TP53, ClinVar 2025 otherwise.

    A row's `Gene` value is assumed never to mix a `MIXED_YEAR_GENES` gene
    with a non-member gene -- true today, since the only multi-gene `Gene`
    values are the CALM1/CALM2/CALM3 combination, which isn't one of the
    four.
    """
    use_2018 = df[GENE_COL].apply(lambda v: bool(MIXED_YEAR_GENES & set(split_genes(v))))
    return df[CLINVAR_COL_2018].where(use_2018, df[CLINVAR_COL])


def non_clinvar_row_flags(df, snv_label):
    """Per-row boolean flags that don't depend on ClinVar significance: score
    coverage, gnomAD observation, and the SNV(-accessible) flag.
    """
    flags = pd.DataFrame(index=df.index)
    for label, col in SCORE_COLUMNS.items():
        flags[label] = has_any_value(df[col])
    flags[GNOMAD_LABEL] = has_any_value(df[GNOMAD_COL])
    flags[snv_label] = is_snv_accessible(df[TRANSCRIPT_REF_COL], df[TRANSCRIPT_ALT_COL])
    return flags


def variant_flags(df, snv_label, clinvar_series=None, allow_clinvar_conflicts=False):
    """Per-row boolean flags for score coverage and clinical attributes.

    `snv_label` is the column name to store the SNV(-accessible) flag under
    -- `SNV_LABEL` for the DNA-level file, `SNV_ACCESSIBLE_LABEL` for the
    assayed (protein-level) file. See `is_snv_accessible`.

    `clinvar_series` overrides which ClinVar significance values are used for
    the VUS/pathogenic-or-benign/no-annotation flags -- defaults to
    `df[CLINVAR_COL]` (ClinVar 2025 for every gene). Pass
    `mixed_year_clinvar_series(df)` to use ClinVar 2018 for BRCA1/PTEN/MSH2/TP53
    instead.

    By default (`allow_clinvar_conflicts=False`), a row whose own
    pipe-delimited ClinVar values disagree (or literally report "Conflicting
    classifications of pathogenicity") is flagged `CLINVAR_CONFLICT_LABEL`
    instead of VUS/pathogenic-or-benign. Pass `allow_clinvar_conflicts=True`
    to fold such rows into VUS/pathogenic-or-benign via any-match instead
    (the original behavior), which also drops the conflict column entirely.
    """
    if clinvar_series is None:
        clinvar_series = df[CLINVAR_COL]
    flags = non_clinvar_row_flags(df, snv_label)
    if allow_clinvar_conflicts:
        flags[VUS_LABEL] = matches_any_value(clinvar_series, VUS_VALUES)
        flags[PATHOGENIC_OR_BENIGN_LABEL] = matches_any_value(clinvar_series, PATHOGENIC_OR_BENIGN_VALUES)
    else:
        vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(
            *clinvar_significance_flags(clinvar_series)
        )
        flags[VUS_LABEL] = vus
        flags[PATHOGENIC_OR_BENIGN_LABEL] = pathogenic_or_benign
        flags[CLINVAR_CONFLICT_LABEL] = conflict
    flags[NO_ANNOTATION_LABEL] = ~has_any_value(clinvar_series) & ~flags[GNOMAD_LABEL]
    return flags


def distinct_variant_flags(df, snv_label, clinvar_series=None, allow_clinvar_conflicts=False):
    """Collapse per-row flags to one row per distinct (hgvs_c, hgvs_p) variant.

    Score coverage, gnomAD observation, and no-annotation flags count as
    True for a distinct variant if any of its (possibly several, across
    datasets and DNA-level candidates) measurement rows does. See
    `variant_flags` for `clinvar_series` and `allow_clinvar_conflicts`;
    when conflicts aren't allowed, the pathogenic/benign/VUS/literal-conflict
    membership used to classify VUS/pathogenic-or-benign/conflict (see
    `clinvar_classification_from_flags`) is OR'd across *all* of the
    variant's rows first, so e.g. a variant called Pathogenic by one
    measurement and Benign by another is flagged a conflict rather than
    both/either.
    """
    if clinvar_series is None:
        clinvar_series = df[CLINVAR_COL]

    non_clinvar = non_clinvar_row_flags(df, snv_label)
    keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), non_clinvar.reset_index(drop=True)], axis=1)
    distinct = keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()

    has_clinvar_value = has_any_value(clinvar_series).reset_index(drop=True)
    hcv_keyed = pd.concat(
        [df[VARIANT_KEY_COLS].reset_index(drop=True), has_clinvar_value.rename("_has_clinvar_value")], axis=1
    )
    hcv_distinct = hcv_keyed.groupby(VARIANT_KEY_COLS, as_index=False)["_has_clinvar_value"].any()
    distinct = distinct.merge(hcv_distinct, on=VARIANT_KEY_COLS)

    if allow_clinvar_conflicts:
        clinvar_flags = pd.DataFrame(
            {
                VUS_LABEL: matches_any_value(clinvar_series, VUS_VALUES),
                PATHOGENIC_OR_BENIGN_LABEL: matches_any_value(clinvar_series, PATHOGENIC_OR_BENIGN_VALUES),
            }
        ).reset_index(drop=True)
        clinvar_keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), clinvar_flags], axis=1)
        clinvar_distinct = clinvar_keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()
        distinct = distinct.merge(clinvar_distinct, on=VARIANT_KEY_COLS)
    else:
        has_pathogenic, has_benign, has_vus, has_literal_conflict = clinvar_significance_flags(clinvar_series)
        significance_flags = pd.DataFrame(
            {
                "_has_pathogenic": has_pathogenic,
                "_has_benign": has_benign,
                "_has_vus": has_vus,
                "_has_literal_conflict": has_literal_conflict,
            }
        ).reset_index(drop=True)
        significance_keyed = pd.concat([df[VARIANT_KEY_COLS].reset_index(drop=True), significance_flags], axis=1)
        significance_distinct = significance_keyed.groupby(VARIANT_KEY_COLS, as_index=False).any()
        distinct = distinct.merge(significance_distinct, on=VARIANT_KEY_COLS)

        vus, pathogenic_or_benign, conflict = clinvar_classification_from_flags(
            distinct["_has_pathogenic"],
            distinct["_has_benign"],
            distinct["_has_vus"],
            distinct["_has_literal_conflict"],
        )
        distinct[VUS_LABEL] = vus
        distinct[PATHOGENIC_OR_BENIGN_LABEL] = pathogenic_or_benign
        distinct[CLINVAR_CONFLICT_LABEL] = conflict
        distinct = distinct.drop(columns=["_has_pathogenic", "_has_benign", "_has_vus", "_has_literal_conflict"])

    distinct[NO_ANNOTATION_LABEL] = ~distinct["_has_clinvar_value"] & ~distinct[GNOMAD_LABEL]
    return distinct.drop(columns=VARIANT_KEY_COLS + ["_has_clinvar_value"])


def summarize_flags(flags):
    """Return (total, DataFrame[count, pct]) for a set of boolean flag columns."""
    total = len(flags)
    counts = flags.sum().astype(int)
    pct = (100 * counts / total).round(1) if total else counts.astype(float)
    return total, pd.DataFrame({"count": counts, "pct": pct})


def format_count_table(title, total, table):
    lines = [title, f"Total: {total}"]
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


def summarize_clinical_flags(flags, snv_label):
    """Like `summarize_flags(flags[clinical_labels])`, plus two columns breaking out
    how many of the level's SNV(-accessible) variants fall into each clinical-attribute
    bucket, as a percentage of that level's total SNV(-accessible) count.

    `clinical_labels` is `CLINICAL_LABELS_ORDER` filtered to columns actually present in
    `flags` -- `CLINVAR_CONFLICT_LABEL` is only there when `variant_flags`/
    `distinct_variant_flags` were called with `allow_clinvar_conflicts=False`.

    Returns (total, snv_total, table).
    """
    clinical_labels = [label for label in CLINICAL_LABELS_ORDER if label in flags.columns]
    total, table = summarize_flags(flags[clinical_labels])
    snv_total = int(flags[snv_label].sum())
    snv_counts = flags[clinical_labels].apply(lambda col: int((col & flags[snv_label]).sum()))
    pct_snv_col = f"% of {snv_label}"
    table[snv_label] = snv_counts
    table[pct_snv_col] = (100 * snv_counts / snv_total).round(1) if snv_total else snv_counts.astype(float)
    return total, snv_total, table


def format_clinical_table(title, total, snv_total, table, snv_label):
    lines = [title, f"Total: {total} ({snv_total} {snv_label})"]
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        pct_snv_col = f"% of {snv_label}"
        body[pct_snv_col] = body[pct_snv_col].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


SCORE_LABELS = list(SCORE_COLUMNS.keys())
CLINICAL_LABELS_ORDER = [
    VUS_LABEL,
    PATHOGENIC_OR_BENIGN_LABEL,
    CLINVAR_CONFLICT_LABEL,
    GNOMAD_LABEL,
    NO_ANNOTATION_LABEL,
]


def build_variant_level_reports(condensed, expanded, condensed_path, expanded_path, allow_clinvar_conflicts=False):
    """Build the score-coverage and clinical-attribute text sections.

    Each is reported at four levels: assayed (protein-resolution, from the
    condensed file) vs. DNA-level (from the expanded file) variants, and
    distinct variants vs. variant measurements (rows). The clinical-attribute
    breakdown is additionally reported twice per level: once using ClinVar
    2025 throughout, and once substituting ClinVar 2018 for BRCA1, PTEN,
    MSH2, and TP53 (see `mixed_year_clinvar_series`). See `variant_flags` for
    `allow_clinvar_conflicts`.
    """
    levels = [
        ("assayed variants, distinct", condensed, condensed_path, SNV_ACCESSIBLE_LABEL, distinct_variant_flags),
        ("assayed variant measurements", condensed, condensed_path, SNV_ACCESSIBLE_LABEL, variant_flags),
        ("DNA variants, distinct", expanded, expanded_path, SNV_LABEL, distinct_variant_flags),
        ("DNA variant measurements", expanded, expanded_path, SNV_LABEL, variant_flags),
    ]

    score_sections = []
    clinical_sections = []
    clinical_sections_mixed_year = []
    for label, df, source, snv_label, flags_fn in levels:
        flags = flags_fn(df, snv_label, allow_clinvar_conflicts=allow_clinvar_conflicts)
        total, table = summarize_flags(flags[SCORE_LABELS])
        score_sections.append(format_count_table(f"Score coverage -- {label} (from {source})", total, table))

        clinical_total, snv_total, clinical_table = summarize_clinical_flags(flags, snv_label)
        clinical_sections.append(
            format_clinical_table(
                f"Clinical attributes -- {label} (from {source})", clinical_total, snv_total, clinical_table, snv_label
            )
        )

        mixed_flags = flags_fn(
            df, snv_label, clinvar_series=mixed_year_clinvar_series(df), allow_clinvar_conflicts=allow_clinvar_conflicts
        )
        mixed_total, mixed_snv_total, mixed_table = summarize_clinical_flags(mixed_flags, snv_label)
        clinical_sections_mixed_year.append(
            format_clinical_table(
                f"Clinical attributes -- {label} (from {source})",
                mixed_total,
                mixed_snv_total,
                mixed_table,
                snv_label,
            )
        )

    return score_sections, clinical_sections, clinical_sections_mixed_year


def excalibr_dataset_to_gene_map(calibration_datasets, metadata):
    """Map each ExCALIBR_calibrations `dataset` value to its `Gene`, via the
    same Curation metadata (`Dataset Name` -> `Gene`) used elsewhere in this
    script.

    A `dataset` value ending in `_clinvar_2018` is a mixed-year duplicate
    calibration of the same underlying dataset -- that suffix isn't part of
    Supplementary_Data_3's `Dataset Name`, so it's stripped before lookup.
    Both sides are also Unicode-normalized (NFC) before comparison, since
    accented characters (e.g. in author names) can round-trip through Excel
    in different composed/decomposed forms across files.

    Raises ValueError listing any `dataset` value that still doesn't resolve
    to a `Dataset Name` in the metadata after these adjustments.
    """
    normalized_index = {unicodedata.normalize("NFC", name): name for name in metadata.index}

    def resolve(dataset):
        stripped = dataset.removesuffix(EXCALIBR_MIXED_YEAR_DATASET_SUFFIX)
        return normalized_index.get(unicodedata.normalize("NFC", stripped))

    mapping = {}
    unmapped = []
    for dataset in calibration_datasets:
        resolved = resolve(dataset)
        if resolved is None:
            unmapped.append(dataset)
        else:
            mapping[dataset] = metadata.loc[resolved, GENE_COL]
    if unmapped:
        raise ValueError(
            f"{EXCALIBR_CALIBRATIONS_SHEET} dataset(s) not found in {METADATA_SHEET} metadata: {sorted(unmapped)}"
        )
    return mapping


def compute_excalibr_calibration_stats(calibrations, metadata, merge_calm_genes=False):
    """How many genes have an ExCALIBR calibration, and how many of those
    genes have at least one dataset where ExCALIBR assigned at least one
    point of evidence (pathogenic or benign) -- i.e. at least one non-null
    `range_-8`..`range_8` value in that dataset's row.

    A gene can have several calibration rows (one per dataset); this counts
    distinct genes, not rows. Also reports both counts a second way, excluding
    `EXCALIBR_EXCLUDED_GENES` (F9, TP53, SFPQ), since F9/TP53 use OddsPath
    rather than ExCALIBR calibration (see `docs/variant_classification.md`).
    """
    dataset_to_gene = excalibr_dataset_to_gene_map(calibrations[EXCALIBR_DATASET_COL], metadata)
    range_cols = [c for c in calibrations.columns if c.startswith(EXCALIBR_RANGE_COLUMN_PREFIX)]
    row_has_evidence = calibrations[range_cols].notna().any(axis=1)
    has_evidence_by_dataset = row_has_evidence.groupby(calibrations[EXCALIBR_DATASET_COL]).any()

    genes_with_calibration = set()
    genes_with_evidence = set()
    for dataset, gene_value in dataset_to_gene.items():
        genes = set(split_genes(gene_value))
        if merge_calm_genes and genes & CALM_GENES:
            genes -= CALM_GENES
            genes.add(CALM_MERGED_LABEL)
        genes_with_calibration |= genes
        if has_evidence_by_dataset.get(dataset, False):
            genes_with_evidence |= genes

    return {
        "genes_with_excalibr_calibrations": len(genes_with_calibration),
        "genes_with_evidence_assigned": len(genes_with_evidence),
        "genes_with_excalibr_calibrations_excl": len(genes_with_calibration - EXCALIBR_EXCLUDED_GENES),
        "genes_with_evidence_assigned_excl": len(genes_with_evidence - EXCALIBR_EXCLUDED_GENES),
    }


def format_calibration_summary(stats):
    total = stats["genes_with_excalibr_calibrations"]
    with_evidence = stats["genes_with_evidence_assigned"]
    pct = 100 * with_evidence / total if total else float("nan")
    total_excl = stats["genes_with_excalibr_calibrations_excl"]
    with_evidence_excl = stats["genes_with_evidence_assigned_excl"]
    pct_excl = 100 * with_evidence_excl / total_excl if total_excl else float("nan")
    excl_suffix = "excluding F9/TP53/SFPQ"
    return "\n".join(
        [
            "=== ExCALIBR calibration coverage (Extended Data Figure 2) ===",
            f"Genes with ExCALIBR calibrations: {total} ({total_excl} {excl_suffix})",
            (
                f"Genes with >=1 dataset assigning >=1 point of evidence (pathogenic or benign): "
                f"{with_evidence} ({pct:.1f}%) ({with_evidence_excl} ({pct_excl:.1f}%) {excl_suffix})"
            ),
        ]
    )


def funnel_distinct_dna_variants(df):
    """Count distinct DNA variants in `df` by the reclassification pipeline's
    own dedup key (`FUNNEL_GENOMIC_KEY_COLS`).

    `hg38_start` is parsed numerically and `Chrom` has any trailing `.0`
    stripped before deduping: the notebook's checkpoint file
    (`DEFAULT_CHECKPOINT_FILE`) writes both columns out from a
    float-typed pandas column, so plain string comparison against the
    expanded file's own (un-suffixed) string values would spuriously split
    what's really the same variant into two groups.
    """
    hg38_start = pd.to_numeric(df["hg38_start"], errors="coerce")
    chrom = df["Chrom"].astype(str).str.replace(r"\.0$", "", regex=True)
    key = pd.DataFrame(
        {
            "Gene": df["Gene"].values,
            "Chrom": chrom.values,
            "hg38_start": hg38_start.values,
            "ref_allele": df["ref_allele"].values,
            "alt_allele": df["alt_allele"].values,
        }
    )
    return key.drop_duplicates().shape[0]


def compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed):
    """Sequential funnel from every DNA-level measurement row in the expanded
    file down to the reclassification export, applying each exclusion in the
    pipeline's actual order and recording distinct-DNA-variant and
    distinct-assayed-variant counts after each step.

    Steps 1-3 (LDLR +VLDL assay priority, LDLR LA module 1, F9/TP53
    restricted datasets) are `Variant_Classification_analysis.ipynb` cells 3,
    4, and 42 respectively, reproduced against `expanded` since they're what
    produces `checkpoint` in the first place. The remaining steps mirror
    `src.build_variant_reclassification_dataset.apply_notebook_exclusions`,
    split into named sub-steps: Flag=='*' removal is split into the CHEK2
    QC flag specifically vs. any other pre-existing flag, since the notebook
    sets the former (`Filter_CI == 1`) immediately before removing all
    flagged rows together; and the single `VariantNotes`-tag exclusion is
    split into its three underlying tags (`CONFLICTING_FXN_DATA_TAG`/
    `SPLICE_VARIANT_NOT_MEASURED_TAG`/`START_LOST_VARIANT_NOT_MEASURED_TAG`),
    applied in that order so a combined-tag row (e.g. both conflicting and
    splice-related) is removed at whichever of its tags' steps comes first.

    A "distinct assayed variant" is *not* a distinct `mavedb_variant_urn`:
    that identifies one MaveDB score-set record (one per dataset), so the
    same variant assayed by two datasets would count twice. Instead, since
    `mavedb_variant_urn` is 1:1 with `condensed`'s rows, each expanded/
    checkpoint row's urn is mapped back to its parent condensed row's
    `VARIANT_KEY_COLS` (`hgvs_c`, `hgvs_p`) pair -- the same key
    `distinct_variants_assayed` uses -- and distinct *pairs* are counted.
    This reproduces `distinct_variants_assayed` exactly over the unfiltered
    file: splitting by `nucleotide_or_aa` instead (distinct `hgvs_p` for `aa`
    rows plus distinct `hgvs_c` for `nt` rows) overcounts by the number of
    variants assayed at *both* resolutions -- an `aa`-side reverse-translation
    candidate that happens to equal some `nt`-side row's exact `hgvs_c` for
    the same `hgvs_p` (535 such variants in the full dataset) -- since that
    split double-counts them instead of merging them the way the
    `(hgvs_c, hgvs_p)` pair key naturally does.

    The REVEL- and MutPred2-training exclusions (the last two steps) are
    siblings, not a chain: both filter the same "Other flagged variants"
    survivor set independently -- see `REVEL_TRAINING_STEP_LABEL`/
    `MUTPRED2_TRAINING_STEP_LABEL`. Each step's `*_removed` fields are the
    drop from its own logical predecessor (that survivor set, for both of
    these two), not necessarily from whichever step precedes it in the
    returned list.

    Returns a list of `{"label", "rows", "rows_removed",
    "variant_effect_measurements", "measurements_removed",
    "distinct_dna_variants", "dna_removed", "distinct_assayed_variants",
    "assayed_removed"}` dicts, one per step (including a row 0 for the
    unfiltered expanded file, whose `*_removed` fields are all 0).
    `variant_effect_measurements` is `mavedb_variant_urn`'s distinct-value
    count -- the same measurement-level count as
    `integrated_variant_effect_dataset.condensed.tsv`'s row count (`urn` is
    1:1 with `condensed`'s rows), reached here by grouping the expanded/
    checkpoint file's (possibly several per measurement, e.g. multiple
    DNA-level candidates for one protein-resolution measurement) rows by
    their shared `urn` -- unlike `rows`, which counts those DNA-level rows
    directly without collapsing by `urn`.
    """
    urn_to_assayed_key = condensed.set_index(MAVEDB_VARIANT_URN_COL)[VARIANT_KEY_COLS].apply(tuple, axis=1)

    steps = []

    def record(label, df, baseline=None):
        """Append and return a new step dict for `df`. `*_removed` fields are
        the drop from `baseline` (defaulting to the most-recently-recorded
        step) -- pass `baseline` explicitly for a step that doesn't follow
        the previous list entry (e.g. the two training-variant siblings).
        """
        if baseline is None:
            baseline = steps[-1] if steps else None
        row_count = len(df)
        measurements = int(df[MAVEDB_VARIANT_URN_COL].nunique())
        dna = funnel_distinct_dna_variants(df)
        assayed = int(df[MAVEDB_VARIANT_URN_COL].map(urn_to_assayed_key).nunique())
        step_record = {
            "label": label,
            "rows": row_count,
            "rows_removed": 0 if baseline is None else baseline["rows"] - row_count,
            "variant_effect_measurements": measurements,
            "measurements_removed": (
                0 if baseline is None else baseline["variant_effect_measurements"] - measurements
            ),
            "distinct_dna_variants": dna,
            "dna_removed": 0 if baseline is None else baseline["distinct_dna_variants"] - dna,
            "distinct_assayed_variants": assayed,
            "assayed_removed": 0 if baseline is None else baseline["distinct_assayed_variants"] - assayed,
        }
        steps.append(step_record)
        return step_record

    record("Assayed DNA-level measurements (expanded file)", expanded)

    aa_pos = pd.to_numeric(expanded["aa_pos"], errors="coerce")
    in_vldl_range = pd.Series(False, index=expanded.index)
    for lo, hi in LDLR_VLDL_PRIORITY_RANGES:
        in_vldl_range |= aa_pos.between(lo, hi)
    ldlr_vldl_range_mask = (expanded["Gene"] == "LDLR") & in_vldl_range
    vldl_rows = expanded.loc[
        ldlr_vldl_range_mask & (expanded["Dataset"] == LDLR_VLDL_PRIORITY_DATASET),
        ["aa_pos", "aa_ref", "aa_alt"],
    ]
    vldl_keys = set(map(tuple, vldl_rows.values))
    row_aa_keys = pd.Series(list(zip(expanded["aa_pos"], expanded["aa_ref"], expanded["aa_alt"])), index=expanded.index)
    superseded_mask = (
        ldlr_vldl_range_mask & expanded["Dataset"].isin(LDLR_VLDL_FALLBACK_DATASETS) & row_aa_keys.isin(vldl_keys)
    )
    step = expanded[~superseded_mask]
    record("- LDLR: +VLDL assay preferred over abundance/uptake (LA modules 2/6)", step)

    la1_mask = (step["Gene"] == "LDLR") & aa_pos.reindex(step.index).between(*LDLR_LA1_EXCLUDE_RANGE)
    step = step[~la1_mask]
    record("- LDLR: LA module 1 excluded entirely (aa 25-65)", step)

    step = step[~step["Dataset"].isin(F9_TP53_RESTRICTED_DATASETS)]
    record("- F9/TP53: restricted to the meta-analysis dataset only", step)

    step = checkpoint[checkpoint[GENE_COL] != "SFPQ"].copy()
    record("- SFPQ excluded entirely (insufficient ClinVar controls)", step)

    step = step.reset_index(drop=True)
    hgvs_p_no_transcript = step["hgvs_p"].str.replace(r"^[^:]+:", "", regex=True)
    merged = pd.merge(
        pd.DataFrame(
            {"_hgvs_p_no_transcript": hgvs_p_no_transcript, "auth_reported_score": step["auth_reported_score"]}
        ),
        chek2[["hgvs_pro", "score", "Filter_CI"]],
        left_on=["_hgvs_p_no_transcript", "auth_reported_score"],
        right_on=["hgvs_pro", "score"],
        how="left",
    )
    chek2_flag_mask = (merged["Filter_CI"] == 1).fillna(False)
    chek2_flag_mask.index = step.index
    flag = step["Flag"].where(~chek2_flag_mask, "*")
    pre_existing_flag_mask = (flag == "*") & ~chek2_flag_mask
    step = step.assign(Flag=flag)

    variant_notes = step["VariantNotes"].fillna("")
    step = step[~variant_notes.str.contains(CONFLICTING_FXN_DATA_TAG, regex=False)]
    record("- Conflicting functional data (opposite-sign Fxn_points across assays)", step)

    variant_notes = variant_notes.reindex(step.index)
    step = step[~variant_notes.str.contains(SPLICE_VARIANT_NOT_MEASURED_TAG, regex=False)]
    record("- Splice variant not measured (no functional measurement at all)", step)

    variant_notes = variant_notes.reindex(step.index)
    step = step[~variant_notes.str.contains(START_LOST_VARIANT_NOT_MEASURED_TAG, regex=False)]
    record("- Start-lost variant not measured (no functional measurement at all)", step)

    step = step[step["splice_var_amino"] != "Yes"]
    record("- Splice variant not measured at the amino-acid level (aa-level candidates sharing a splice-affecting group)", step)

    step = step[~chek2_flag_mask.reindex(step.index, fill_value=False)]
    record("- CHEK2 QC flag (Filter_CI == 1)", step)

    step = step[~pre_existing_flag_mask.reindex(step.index, fill_value=False)]
    other_flagged_step = record("- Other flagged variants (pre-existing Flag == '*')", step)

    # Siblings, not a chain -- both filter `step` (the "Other flagged variants"
    # survivors) independently, so each is baselined against that same step
    # rather than against each other.
    record(REVEL_TRAINING_STEP_LABEL, step[step["revel_train_amino"] != "Yes"], baseline=other_flagged_step)
    record(MUTPRED2_TRAINING_STEP_LABEL, step[step["mp2_train_amino"] != "Yes"], baseline=other_flagged_step)

    return steps


def _reclassification_funnel_table_row(step):
    return {
        "": step["label"],
        "rows": step["rows"],
        "rows_removed": step["rows_removed"],
        "variant_effect_measurements": step["variant_effect_measurements"],
        "measurements_removed": step["measurements_removed"],
        "distinct_dna_variants": step["distinct_dna_variants"],
        "dna_removed": step["dna_removed"],
        "distinct_assayed_variants": step["distinct_assayed_variants"],
        "assayed_removed": step["assayed_removed"],
    }


def format_reclassification_filter_funnel(steps):
    """Render `compute_reclassification_filter_funnel`'s steps as a table,
    plus a two-line "reclassified out of original" summary at each
    resolution.

    Each of `rows`/`variant_effect_measurements`/`distinct_dna_variants`/
    `distinct_assayed_variants` is paired with its own `*_removed` column
    (each step's own drop from its logical predecessor -- see
    `compute_reclassification_filter_funnel`; 0 for the first row).

    The REVEL- and MutPred2-training exclusion steps
    (`REVEL_TRAINING_STEP_LABEL`/`MUTPRED2_TRAINING_STEP_LABEL`) are rendered
    as their own table below a divider, since they're alternative endpoints
    of the funnel rather than the next step after the rest -- see
    `compute_reclassification_filter_funnel`. The summary lines below still
    describe the REVEL branch specifically.
    """
    by_label = {step["label"]: step for step in steps}
    alternative_labels = {REVEL_TRAINING_STEP_LABEL, MUTPRED2_TRAINING_STEP_LABEL}
    main_steps = [step for step in steps if step["label"] not in alternative_labels]

    table = pd.DataFrame([_reclassification_funnel_table_row(step) for step in main_steps]).set_index("")
    alternatives_table = pd.DataFrame(
        [_reclassification_funnel_table_row(by_label[label]) for label in (REVEL_TRAINING_STEP_LABEL, MUTPRED2_TRAINING_STEP_LABEL)]
    ).set_index("")

    final = by_label[REVEL_TRAINING_STEP_LABEL]
    starting_dna = steps[0]["distinct_dna_variants"]
    starting_assayed = steps[0]["distinct_assayed_variants"]
    dna_pct = 100 * final["distinct_dna_variants"] / starting_dna if starting_dna else float("nan")
    assayed_pct = 100 * final["distinct_assayed_variants"] / starting_assayed if starting_assayed else float("nan")

    lines = [
        "=== Filtering effects on the reclassification dataset ===",
        table.to_string(),
        "-" * 80,
        "Alternative endpoints (mutually exclusive, neither applied on top of the other):",
        alternatives_table.to_string(),
        "",
        (f"Distinct DNA variants reclassified: {final['distinct_dna_variants']} of {starting_dna} ({dna_pct:.1f}%)"),
        (
            f"Distinct assayed variants reclassified: {final['distinct_assayed_variants']} "
            f"of {starting_assayed} ({assayed_pct:.1f}%)"
        ),
    ]
    return "\n".join(lines)


def reclassification_flags(clinvar_group, points):
    """Per-row agree/disagree/no-evidence flags for one points column against
    the row's ClinVar pathogenic-or-benign control label.

    A row is in scope only if `clinvar_group` is a pathogenic- or
    benign-leaning control label (see `PATHOGENIC_VALUES`/`BENIGN_VALUES`);
    `points` is treated as pathogenic evidence if positive, benign evidence
    if negative, and no evidence if zero or missing. Returns (flags, in_scope),
    both boolean Series/DataFrame aligned to `clinvar_group`'s index.
    """
    is_pathogenic = clinvar_group.isin(PATHOGENIC_VALUES)
    is_benign = clinvar_group.isin(BENIGN_VALUES)
    in_scope = is_pathogenic | is_benign

    assigned_pathogenic = points > 0
    assigned_benign = points < 0
    no_evidence = ~assigned_pathogenic & ~assigned_benign

    agree = (is_pathogenic & assigned_pathogenic) | (is_benign & assigned_benign)
    disagree = (is_pathogenic & assigned_benign) | (is_benign & assigned_pathogenic)

    flags = pd.DataFrame({AGREE_LABEL: agree, DISAGREE_LABEL: disagree, NO_EVIDENCE_LABEL: no_evidence})
    return flags, in_scope


def compute_reclassification_agreement(controls_df):
    """For each points column in `RECLASSIFICATION_POINTS_COLUMNS`, the
    agree/disagree/no-evidence breakdown (see `reclassification_flags`) plus
    the agreement percentage among rows with evidence assigned.

    Returns {label: (total, determinate, agreement_pct, table)}, where `table`
    is the `summarize_flags` output over the in-scope rows.
    """
    results = {}
    for label, points_col in RECLASSIFICATION_POINTS_COLUMNS.items():
        flags, in_scope = reclassification_flags(controls_df[CONTROLS_CLINVAR_GROUP_COL], controls_df[points_col])
        flags = flags[in_scope]
        total, table = summarize_flags(flags[RECLASSIFICATION_LABELS_ORDER])
        determinate = total - int(flags[NO_EVIDENCE_LABEL].sum())
        agreement_pct = 100 * int(flags[AGREE_LABEL].sum()) / determinate if determinate else float("nan")
        results[label] = (total, determinate, agreement_pct, table)
    return results


def format_reclassification_table(title, total, determinate, agreement_pct, table):
    lines = [title, f"Total control variants: {total}", f"Determinate calls (evidence assigned): {determinate}"]
    if determinate:
        lines.append(f"Agreement with ClinVar PLP/BLB (of determinate calls): {agreement_pct:.1f}%")
    if total:
        body = table.copy()
        body["pct"] = body["pct"].map(lambda x: f"{x:.1f}%")
        lines.append(body.to_string())
    return "\n".join(lines)


def build_reclassification_report(workbook):
    """One section per (`controls_`-prefixed sheet, points column) pair in
    `workbook` (an open `pd.ExcelFile` over the controls file). See the module
    docstring's "Reclassification agreement" note for why every such sheet is
    expected to agree.
    """
    sheets = [name for name in workbook.sheet_names if name.startswith(CONTROLS_SHEET_PREFIX)]
    sections = []
    for sheet in sheets:
        controls_df = workbook.parse(sheet)
        agreement = compute_reclassification_agreement(controls_df)
        for label, (total, determinate, agreement_pct, table) in agreement.items():
            sections.append(
                format_reclassification_table(f"{label} -- {sheet}", total, determinate, agreement_pct, table)
            )
    return sections


def control_concordance_flags(control_group, pathogenic_values, benign_values, assigned_pathogenic, assigned_benign):
    """Per-row concordant/discordant/VUS flags for one evidence source against
    a control classification column.

    A row is in scope only if `control_group` is one of `pathogenic_values`/
    `benign_values`. `assigned_pathogenic`/`assigned_benign` are boolean
    Series already derived from whichever evidence source is being compared
    (e.g. `OP_points > 0`/`< 0` for OddsPath alone, or `Class_REVEL` category
    membership for the combined-with-REVEL evidence) -- a row with neither
    set counts as VUS. `DISCORDANT_LABEL` also carries two mutually-exclusive
    sub-flags (`DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL`/
    `DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL`) breaking discordant rows
    out by direction. Returns (flags, in_scope), both boolean Series/
    DataFrame aligned to `control_group`'s index.
    """
    is_pathogenic = control_group.isin(pathogenic_values)
    is_benign = control_group.isin(benign_values)
    in_scope = is_pathogenic | is_benign

    concordant = (is_pathogenic & assigned_pathogenic) | (is_benign & assigned_benign)
    discordant_plp_to_blb = is_pathogenic & assigned_benign
    discordant_blb_to_plp = is_benign & assigned_pathogenic
    discordant = discordant_plp_to_blb | discordant_blb_to_plp
    vus = ~assigned_pathogenic & ~assigned_benign

    flags = pd.DataFrame(
        {
            CONCORDANT_LABEL: concordant,
            DISCORDANT_LABEL: discordant,
            DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL: discordant_plp_to_blb,
            DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL: discordant_blb_to_plp,
            CONTROL_VUS_LABEL: vus,
        }
    )
    return flags, in_scope


def compute_control_concordance(workbook):
    """For each control source in `CONTROL_CONCORDANCE_SOURCES` (ClinVar,
    ClinGen), the concordant/discordant/VUS breakdown against OddsPath
    calibration evidence alone, plus the combined ExCALIBR/OddsPath +
    <predictor> gene-specific evidence for each of REVEL/AlphaMissense/
    MutPred2 -- see the module-level comment above `CLINGEN_CLASSIFICATION_COL`
    for exactly what each evidence source is. `Discordant` is further split
    into its two directions (control Pathogenic/Likely Pathogenic reclassified
    Benign/Likely Benign by the evidence source, or vice versa) -- see
    `control_concordance_flags`.

    Returns {(control_source_label, evidence_label): (total, table)}, where
    `table` is the `summarize_flags` output over the in-scope rows.
    """
    results = {}
    for control_label, (category, control_col, pathogenic_values, benign_values) in CONTROL_CONCORDANCE_SOURCES.items():
        revel_sheets = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR["REVEL"]
        revel_df = workbook.parse(revel_sheets[category])
        op_points = revel_df[FUNCTIONAL_CLASS_POINTS_COL]
        oddspath_flags, oddspath_in_scope = control_concordance_flags(
            revel_df[control_col], pathogenic_values, benign_values, op_points > 0, op_points < 0
        )
        results[(control_label, CONTROL_CONCORDANCE_EVIDENCE_LABEL)] = summarize_flags(
            oddspath_flags.loc[oddspath_in_scope, CONTROL_CONCORDANCE_LABELS_ORDER]
        )

        for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
            sheets = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]
            class_col = VARIANT_CLASSIFICATION_CLASS_COL_BY_PREDICTOR[predictor]
            df = workbook.parse(sheets[category])
            combined_flags, combined_in_scope = control_concordance_flags(
                df[control_col],
                pathogenic_values,
                benign_values,
                df[class_col].isin(CLASS_PATHOGENIC_VALUES),
                df[class_col].isin(CLASS_BENIGN_VALUES),
            )
            results[(control_label, COMBINED_EVIDENCE_LABEL_BY_PREDICTOR[predictor])] = summarize_flags(
                combined_flags.loc[combined_in_scope, CONTROL_CONCORDANCE_LABELS_ORDER]
            )
    return results


def format_control_concordance_report(concordance):
    """Render `compute_control_concordance`'s output as one combined table per
    control source (ClinVar, ClinGen), with one row per evidence source
    (OddsPath alone, then combined with each of REVEL/AlphaMissense/MutPred2)
    so all four evidence sources can be compared at a glance instead of
    spreading them across separate subsections. Each row's two directional
    discordance columns (`DISCORDANT_CONTROL_PLP_TO_EVIDENCE_BLB_LABEL`/
    `DISCORDANT_CONTROL_BLB_TO_EVIDENCE_PLP_LABEL`) are reported as a percent
    of that row's own total, same as `Concordant`/`Discordant`/`VUS` --
    together they sum to `Discordant`.
    """
    evidence_labels = [CONTROL_CONCORDANCE_EVIDENCE_LABEL] + [
        COMBINED_EVIDENCE_LABEL_BY_PREDICTOR[predictor] for predictor in VARIANT_CLASSIFICATION_PREDICTORS
    ]
    sections = [
        "=== Control concordance (ClinVar vs. ClinGen; OddsPath alone vs. combined with "
        "REVEL/AlphaMissense/MutPred2) ==="
    ]
    for control_label, (category, *_rest) in CONTROL_CONCORDANCE_SOURCES.items():
        sheet_pattern = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR["REVEL"][category].replace("REVEL", "*")
        rows = {}
        for evidence_label in evidence_labels:
            total, table = concordance[(control_label, evidence_label)]
            rows[evidence_label] = {"Total": total} | {
                label: _format_count_and_pct(int(table.loc[label, "count"]), total)
                for label in CONTROL_CONCORDANCE_LABELS_ORDER
            }
        report_table = pd.DataFrame(rows).T
        sections.append(f"{control_label} controls ({sheet_pattern} sheets):\n{report_table.to_string()}")
    return "\n\n".join(sections)


def distinct_dna_variants(df):
    """Collapse `df` to one row per distinct DNA variant, identified by
    genomic coordinates (`VARIANT_CLASSIFICATION_COORD_COLS`) -- the same key
    `docs/variant_classification.md` uses to identify a physical DNA variant
    (`Gene` is included alongside the genomic coordinates as a defensive
    tie-breaker, matching that doc's own dedup key).
    """
    return df.drop_duplicates(subset=VARIANT_CLASSIFICATION_COORD_COLS)


def _checkpoint_other_flagged_survivors(checkpoint, chek2):
    """Reproduce `compute_reclassification_filter_funnel`'s exclusion chain
    (SFPQ, CHEK2 QC flag, `VariantNotes` conflict/splice/start-lost tags,
    `splice_var_amino`, any other pre-existing `Flag == '*'`), collapsed to
    its final surviving population -- the same "Other flagged variants"
    step, i.e. `Variant_Classification_analysis.ipynb`'s `sankey_f` after
    cell 77, before the REVEL-/MutPred2-training split. Duplicated here
    (rather than sharing code with the funnel) because the funnel needs
    every sub-step recorded individually for its own report section, while
    this only needs the final surviving rows -- see that function's
    docstring for the step-by-step version this mirrors.
    """
    step = checkpoint[checkpoint[GENE_COL] != "SFPQ"].copy()
    step = step.reset_index(drop=True)
    hgvs_p_no_transcript = step["hgvs_p"].str.replace(r"^[^:]+:", "", regex=True)
    merged = pd.merge(
        pd.DataFrame(
            {"_hgvs_p_no_transcript": hgvs_p_no_transcript, "auth_reported_score": step["auth_reported_score"]}
        ),
        chek2[["hgvs_pro", "score", "Filter_CI"]],
        left_on=["_hgvs_p_no_transcript", "auth_reported_score"],
        right_on=["hgvs_pro", "score"],
        how="left",
    )
    chek2_flag_mask = (merged["Filter_CI"] == 1).fillna(False)
    chek2_flag_mask.index = step.index
    step = step.assign(Flag=step["Flag"].where(~chek2_flag_mask, "*"))

    variant_notes = step["VariantNotes"].fillna("")
    step = step[~variant_notes.str.contains(CONFLICTING_FXN_DATA_TAG, regex=False)]
    variant_notes = variant_notes.reindex(step.index)
    step = step[~variant_notes.str.contains(SPLICE_VARIANT_NOT_MEASURED_TAG, regex=False)]
    variant_notes = variant_notes.reindex(step.index)
    step = step[~variant_notes.str.contains(START_LOST_VARIANT_NOT_MEASURED_TAG, regex=False)]
    step = step[step["splice_var_amino"] != "Yes"]
    return step[step["Flag"] != "*"]


def compute_clingen_evidence_repository_stats(checkpoint, chek2, controls_workbook):
    """The ClinGen Evidence Repository "evidence removal" control set: how
    many DNA variants have a determinate (Pathogenic/Likely Pathogenic/
    Benign/Likely Benign) *original* ClinGen classification
    (`ASSERTION_CLINGEN_REPO_COL`) and MAVE experimental data; how many of
    those retain a determinate classification after
    `src/recalculate_clingen_classification.py` recomputes it with the
    functional (BS3/PS3) and predictive (BP4/PP3) evidence codes removed
    (`CLINGEN_CLASSIFICATION_COL`); and how many are ultimately retained for
    analysis with each of REVEL/AlphaMissense/MutPred2 (the
    `ClinGen_Repo_*_GeneSpecific` sheets, after also excluding that
    predictor's own training variants).

    Distinct DNA variants for the first two counts are by simple genomic-
    coordinate dedup (`distinct_dna_variants`), not the notebook's fuller
    nt/aa tie-break (`src/lib/dedup.py`), so those two totals may differ
    slightly from a from-scratch reproduction of the notebook's own
    per-predictor dedup -- but the three `{predictor}_total` figures below
    are read directly from the actual `ClinGen_Repo_*_GeneSpecific` sheets,
    so they match the notebook's real output exactly.
    """
    survivors = _checkpoint_other_flagged_survivors(checkpoint, chek2)
    determinate_original = survivors[survivors[ASSERTION_CLINGEN_REPO_COL].isin(CLASS_PATHOGENIC_OR_BENIGN_VALUES)]
    pre_removal = distinct_dna_variants(determinate_original)

    retained = determinate_original[
        determinate_original[CLINGEN_CLASSIFICATION_COL].isin(CLASS_PATHOGENIC_OR_BENIGN_VALUES)
    ]
    post_removal = distinct_dna_variants(retained)

    stats = {
        "pre_removal_total": len(pre_removal),
        "pre_removal_genes": int(pre_removal[GENE_COL].nunique()),
        "pre_removal_pathogenic": int((pre_removal[ASSERTION_CLINGEN_REPO_COL] == "Pathogenic").sum()),
        "pre_removal_likely_pathogenic": int((pre_removal[ASSERTION_CLINGEN_REPO_COL] == "Likely Pathogenic").sum()),
        "pre_removal_benign": int((pre_removal[ASSERTION_CLINGEN_REPO_COL] == "Benign").sum()),
        "pre_removal_likely_benign": int((pre_removal[ASSERTION_CLINGEN_REPO_COL] == "Likely Benign").sum()),
        "post_removal_total": len(post_removal),
        "post_removal_genes": int(post_removal[GENE_COL].nunique()),
        "post_removal_plp": int(post_removal[CLINGEN_CLASSIFICATION_COL].isin(CLASS_PATHOGENIC_VALUES).sum()),
        "post_removal_blb": int(post_removal[CLINGEN_CLASSIFICATION_COL].isin(CLASS_BENIGN_VALUES).sum()),
    }

    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        sheet = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]["ClinGen_Repo"]
        df = controls_workbook.parse(sheet)
        stats[f"{predictor}_total"] = len(df)
        stats[f"{predictor}_genes"] = int(df[GENE_COL].nunique())
        stats[f"{predictor}_plp"] = int(df[CLINGEN_CLASSIFICATION_COL].isin(CLASS_PATHOGENIC_VALUES).sum())
        stats[f"{predictor}_blb"] = int(df[CLINGEN_CLASSIFICATION_COL].isin(CLASS_BENIGN_VALUES).sum())

    return stats


CLINGEN_EVIDENCE_REPOSITORY_TITLE = "=== ClinGen Evidence Repository control set (evidence removed & reclassified) ==="


def format_clingen_evidence_repository_summary(stats, title=CLINGEN_EVIDENCE_REPOSITORY_TITLE):
    lines = [
        title,
        "Determinate (Pathogenic/Likely Pathogenic/Benign/Likely Benign) original ClinGen Evidence Repository "
        "classification and MAVE experimental data: "
        f"{stats['pre_removal_total']} distinct DNA variants across {stats['pre_removal_genes']} genes",
        f"  Pathogenic: {stats['pre_removal_pathogenic']}",
        f"  Likely Pathogenic: {stats['pre_removal_likely_pathogenic']}",
        f"  Benign: {stats['pre_removal_benign']}",
        f"  Likely Benign: {stats['pre_removal_likely_benign']}",
        "",
        "After removing the functional (BS3/PS3) and predictive (BP4/PP3) evidence used in the original "
        "classification and recalculating from the remaining evidence codes:",
        "  Retained a determinate classification: "
        + _format_count_and_pct(stats["post_removal_total"], stats["pre_removal_total"])
        + f" distinct DNA variants across {stats['post_removal_genes']} genes",
        f"    Pathogenic or Likely Pathogenic: {stats['post_removal_plp']}",
        f"    Benign or Likely Benign: {stats['post_removal_blb']}",
        "",
        "Post filtering (excluding each predictor's own training variants), retained for analysis:",
    ]
    for predictor in VARIANT_CLASSIFICATION_PREDICTORS:
        lines.append(
            f"  {predictor}: {stats[f'{predictor}_total']} variants across {stats[f'{predictor}_genes']} genes "
            f"(Pathogenic or Likely Pathogenic: {stats[f'{predictor}_plp']}, "
            f"Benign or Likely Benign: {stats[f'{predictor}_blb']})"
        )
    return "\n".join(lines)


def _compute_variant_classification_stats_for_predictor(workbook, predictor):
    """`compute_variant_classification_stats` for one predictor, reading that
    predictor's own dedicated category sheets (`VARIANT_CLASSIFICATION_
    CATEGORY_SHEETS_BY_PREDICTOR`) and its own `Class_*`/`Total_Points_*`/
    `Conflicting_*` columns (plus the shared `Fxn_points`).
    """
    category_sheets = VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR[predictor]
    class_col = VARIANT_CLASSIFICATION_CLASS_COL_BY_PREDICTOR[predictor]
    points_col = VARIANT_CLASSIFICATION_POINTS_COL_BY_PREDICTOR[predictor]
    predictor_points_col = VARIANT_CLASSIFICATION_PREDICTOR_POINTS_COL_BY_PREDICTOR[predictor]
    conflicting_col = VARIANT_CLASSIFICATION_CONFLICTING_COL_BY_PREDICTOR[predictor]

    sheets = {label: distinct_dna_variants(workbook.parse(sheet_name)) for label, sheet_name in category_sheets.items()}

    combined = distinct_dna_variants(pd.concat(sheets.values(), ignore_index=True))
    classified = combined[class_col].notna()
    pathogenic_or_benign = combined[class_col].isin(CLASS_PATHOGENIC_OR_BENIGN_VALUES)

    def category_counts(df):
        class_values = df[class_col]
        is_pathogenic = class_values.isin(CLASS_PATHOGENIC_VALUES)
        is_benign = class_values.isin(CLASS_BENIGN_VALUES)
        is_unresolved = ~(is_pathogenic | is_benign)
        at_threshold = is_benign & (df[points_col] == LIKELY_BENIGN_POINTS_THRESHOLD)
        single_source = (df[FUNCTIONAL_POINTS_COL] == 0) | (df[predictor_points_col] == 0)
        conflicting = df[conflicting_col] == CONFLICTING_EVIDENCE_VALUE
        has_experimental = df[FUNCTIONAL_POINTS_COL] != 0
        has_predictive = df[predictor_points_col] != 0
        only_experimental = has_experimental & ~has_predictive
        only_predictive = ~has_experimental & has_predictive
        both_evidence = has_experimental & has_predictive
        return {
            "resolved": int(class_values.isin(CLASS_PATHOGENIC_OR_BENIGN_VALUES).sum()),
            "resolved_pathogenic": int(is_pathogenic.sum()),
            "resolved_pathogenic_only_experimental": int((is_pathogenic & only_experimental).sum()),
            "resolved_pathogenic_only_predictive": int((is_pathogenic & only_predictive).sum()),
            "resolved_pathogenic_both_evidence": int((is_pathogenic & both_evidence).sum()),
            "resolved_benign": int(is_benign.sum()),
            "resolved_benign_only_experimental": int((is_benign & only_experimental).sum()),
            "resolved_benign_only_predictive": int((is_benign & only_predictive).sum()),
            "resolved_benign_both_evidence": int((is_benign & both_evidence).sum()),
            "resolved_benign_at_threshold": int(at_threshold.sum()),
            "resolved_benign_at_threshold_single_source": int((at_threshold & single_source).sum()),
            "resolved_benign_at_threshold_conflicting": int((at_threshold & conflicting).sum()),
            "unresolved": int(is_unresolved.sum()),
            "unresolved_concordant": int((is_unresolved & both_evidence & ~conflicting).sum()),
            "unresolved_discordant": int((is_unresolved & both_evidence & conflicting).sum()),
            "unresolved_only_experimental": int((is_unresolved & only_experimental).sum()),
            "unresolved_only_predictive": int((is_unresolved & only_predictive).sum()),
            "unresolved_neither": int((is_unresolved & ~has_experimental & ~has_predictive).sum()),
            "unresolved_zero_or_one_source": int((is_unresolved & ~both_evidence).sum()),
            "unresolved_near_pathogenic": int((is_unresolved & df[points_col].isin(NEAR_PATHOGENIC_POINTS_VALUES)).sum()),
        }

    vus_counts = category_counts(sheets["VUS"])
    gnomad_counts = category_counts(sheets["gnomAD"])
    unobserved_counts = category_counts(sheets["Unobserved"])

    def _prefixed(counts, prefix):
        return {f"{prefix}_{key}": value for key, value in counts.items()}

    return {
        "total_classified": int(classified.sum()),
        "total_pathogenic_or_benign": int(pathogenic_or_benign.sum()),
        "vus_total": len(sheets["VUS"]),
        **_prefixed(vus_counts, "vus"),
        "gnomad_total": len(sheets["gnomAD"]),
        **_prefixed(gnomad_counts, "gnomad"),
        "unobserved_total": len(sheets["Unobserved"]),
        **_prefixed(unobserved_counts, "unobserved"),
    }


def compute_variant_classification_stats(workbook):
    """Reclassification summary (functional evidence from ExCALIBR -- or
    OddsPath for F9/TP53, see `docs/variant_classification.md` -- plus
    predictor evidence from REVEL, AlphaMissense, and MutPred2, each gene-
    specific calibration falling back to genome-wide) across the five
    Supplementary_Data_5 variant categories, for each predictor.

    For each predictor, reads that predictor's own category sheets
    (`VARIANT_CLASSIFICATION_CATEGORY_SHEETS_BY_PREDICTOR`) from `workbook`
    and reports:

    - Across every category combined (a variant appearing in more than one
      category's sheet -- e.g. a ClinVar control also observed in gnomAD --
      is counted once): how many distinct DNA variants have a `Class_*`
      classification, and how many of those are Pathogenic/Likely
      Pathogenic/Benign/Likely Benign rather than Uncertain.
    - Within `VUS` alone: how many ClinVar Uncertain-significance variants
      are "resolved" -- reclassified Pathogenic/Likely Pathogenic/Benign/
      Likely Benign by this pipeline -- and what percent of VUS that is, plus
      how many (and what percent) of VUS land on each side (Pathogenic/Likely
      Pathogenic vs. Benign/Likely Benign). Each side is further split by
      which evidence source(s) contributed a nonzero value (`Fxn_points` for
      experimental, the predictor's own gene-specific/genome-wide points for
      predictive): only experimental, only predictive, or both -- mutually
      exclusive and exhaustive for any resolved (nonzero-total) row, since a
      zero total (both sides zero, or opposite-sign values that cancel to
      zero) always falls in the Uncertain range instead. Separately, the
      Benign/Likely Benign side reports how many have only the minimum
      possible evidence (`Total_Points_* == LIKELY_BENIGN_POINTS_THRESHOLD`,
      i.e. -1 point), further split into how many of *those* got there from
      only one evidence source (`Fxn_points == 0` or the predictor's own
      points are 0 -- the other side contributed nothing) vs. from
      conflicting functional-vs-predictor evidence (`Conflicting_* ==
      CONFLICTING_EVIDENCE_VALUE`, i.e. the two sides disagreed in sign but
      still summed to -1). At exactly -1 point these two are mutually
      exclusive and exhaustive: `Total_Points_* == Fxn_points +
      Points_*_GeneSpecific_GenomeWide`, and two same-signed nonzero terms
      can't sum to a magnitude-1 total, so every -1-point row is either
      single-source or conflicting, never both or neither.
      Separately, `VUS` also reports how many variants are "unresolved" --
      still `Uncertain` after this pipeline -- and what percent of VUS that
      is, further split by the same evidence-source combination used above,
      with the "both" case itself split by sign agreement (`Conflicting_* ==
      CONFLICTING_EVIDENCE_VALUE`): concordant (both present, same sign),
      discordant (both present, opposite signs -- i.e. `Conflicting_*`
      flagged this row, but the evidence still summed to something in the
      Uncertain range), experimental only, predictive only, or neither --
      mutually exclusive and exhaustive over the unresolved rows -- plus a
      combined "zero or one source" total (experimental-only + predictive-
      only + neither, i.e. everything except concordant/discordant). Each of
      these is reported as a percent of the *unresolved* count, not of
      `vus_total`. Also reports how many unresolved rows have
      `Total_Points_* in NEAR_PATHOGENIC_POINTS_VALUES` (4 or 5) -- the
      strongest possible Uncertain calls, one or two points short of
      resolving Likely Pathogenic. This overlaps the categories above (a
      near-pathogenic row can be concordant, discordant, experimental-only,
      or predictive-only) rather than adding a sixth mutually-exclusive
      bucket.
    - Within `gnomAD` alone (population variants with a `gnomad_MAF` value,
      per `notebooks/analysis/README_Variant_Classification_analysis.md` --
      not necessarily lacking a ClinVar call, so this can overlap `VUS`) and
      within `Unobserved` alone (variants absent from both ClinVar and
      gnomAD): the same resolved count/percent, P/LP vs. B/LB split (with
      its experimental/predictive/both evidence-source split), minimum-
      evidence B/LB count (with its single-source vs. conflicting split),
      and unresolved count/percent (with its own concordant/discordant/
      experimental-only/predictive-only/neither split, plus the near-
      pathogenic overlap), computed exactly the same way as `VUS` above.

    Both category-specific counts are computed over each category's own
    distinct DNA variants, independent of the cross-category dedup used for
    the combined total above. Each predictor's population differs slightly
    from the others' (each sheet excludes only *that* predictor's own
    training variants), so totals are not expected to match across
    predictors.

    Returns `{predictor: stats}` for `predictor` in
    `VARIANT_CLASSIFICATION_PREDICTORS`; see
    `format_variant_classification_table`.
    """
    return {
        predictor: _compute_variant_classification_stats_for_predictor(workbook, predictor)
        for predictor in VARIANT_CLASSIFICATION_PREDICTORS
    }


def _format_count_and_pct(count, total):
    pct = 100 * count / total if total else float("nan")
    return f"{count} of {total} ({pct:.1f}%)"


def format_variant_classification_table(stats_by_predictor, title=VARIANT_CLASSIFICATION_TITLE):
    """One row per predictor in `stats_by_predictor` (see
    `compute_variant_classification_stats`) per section -- follows this
    file's established `<DataFrame>.to_string()` convention (see
    `format_count_table`/`format_gene_discordance_summary`).
    """
    predictors = list(stats_by_predictor)

    overall = pd.DataFrame(
        {
            "Distinct DNA variants classified": [stats_by_predictor[p]["total_classified"] for p in predictors],
            "Pathogenic or benign": [
                _format_count_and_pct(
                    stats_by_predictor[p]["total_pathogenic_or_benign"], stats_by_predictor[p]["total_classified"]
                )
                for p in predictors
            ],
        },
        index=predictors,
    )

    def _pct_col(count_key, denom_key):
        return [_format_count_and_pct(stats_by_predictor[p][count_key], stats_by_predictor[p][denom_key]) for p in predictors]

    def _category_table(prefix):
        total_key = f"{prefix}_total"
        resolved_key = f"{prefix}_resolved"
        pathogenic_key = f"{prefix}_resolved_pathogenic"
        benign_key = f"{prefix}_resolved_benign"
        at_threshold_key = f"{prefix}_resolved_benign_at_threshold"
        return pd.DataFrame(
            {
                "Total": [stats_by_predictor[p][total_key] for p in predictors],
                "Resolved": _pct_col(resolved_key, total_key),
                "Pathogenic/Likely Pathogenic": _pct_col(pathogenic_key, total_key),
                "  P/LP: experimental only": _pct_col(f"{pathogenic_key}_only_experimental", pathogenic_key),
                "  P/LP: predictive only": _pct_col(f"{pathogenic_key}_only_predictive", pathogenic_key),
                "  P/LP: both": _pct_col(f"{pathogenic_key}_both_evidence", pathogenic_key),
                "Benign/Likely Benign": _pct_col(benign_key, total_key),
                "  B/LB: experimental only": _pct_col(f"{benign_key}_only_experimental", benign_key),
                "  B/LB: predictive only": _pct_col(f"{benign_key}_only_predictive", benign_key),
                "  B/LB: both": _pct_col(f"{benign_key}_both_evidence", benign_key),
                "Benign/Likely Benign (-1 pt only)": _pct_col(at_threshold_key, total_key),
                "  ...from only one source": _pct_col(f"{at_threshold_key}_single_source", at_threshold_key),
                "  ...conflicting functional/predictive data": _pct_col(
                    f"{at_threshold_key}_conflicting", at_threshold_key
                ),
            },
            index=predictors,
        )

    def _unresolved_table(prefix):
        total_key = f"{prefix}_total"
        unresolved_key = f"{prefix}_unresolved"
        return pd.DataFrame(
            {
                "Total": [stats_by_predictor[p][total_key] for p in predictors],
                "Unresolved": _pct_col(unresolved_key, total_key),
                "Concordant": _pct_col(f"{unresolved_key}_concordant", unresolved_key),
                "Discordant": _pct_col(f"{unresolved_key}_discordant", unresolved_key),
                "Experimental only": _pct_col(f"{unresolved_key}_only_experimental", unresolved_key),
                "Predictive only": _pct_col(f"{unresolved_key}_only_predictive", unresolved_key),
                "Neither": _pct_col(f"{unresolved_key}_neither", unresolved_key),
                "0 or 1 source (total)": _pct_col(f"{unresolved_key}_zero_or_one_source", unresolved_key),
                "+4/+5 points (overlaps above)": _pct_col(f"{unresolved_key}_near_pathogenic", unresolved_key),
            },
            index=predictors,
        )

    vus_table = _category_table("vus")
    gnomad_table = _category_table("gnomad")
    unobserved_table = _category_table("unobserved")
    vus_unresolved_table = _unresolved_table("vus")
    gnomad_unresolved_table = _unresolved_table("gnomad")
    unobserved_unresolved_table = _unresolved_table("unobserved")

    lines = [
        title,
        overall.to_string(),
        "",
        "ClinVar VUS resolved (reclassified pathogenic or benign):",
        vus_table.to_string(),
        "",
        "ClinVar VUS unresolved:",
        vus_unresolved_table.to_string(),
        "",
        "gnomAD variants resolved (classified pathogenic or benign):",
        gnomad_table.to_string(),
        "",
        "gnomAD variants unresolved:",
        gnomad_unresolved_table.to_string(),
        "",
        "Unobserved variants resolved (classified pathogenic or benign):",
        unobserved_table.to_string(),
        "",
        "Unobserved variants unresolved:",
        unobserved_unresolved_table.to_string(),
    ]
    return "\n".join(lines)


def _chi_squared_2x2(count_a, total_a, count_b, total_b):
    """Pearson's chi-squared test of independence on a 2x2 contingency table
    comparing two independent proportions, `count_a/total_a` vs.
    `count_b/total_b`, with Yates' continuity correction -- R's
    `chisq.test()` default for a 2x2 table (equivalent to
    `prop.test(..., correct = TRUE)`), the standard test for "is this rate
    different between two groups" cited as "chi-squared test" throughout the
    manuscript.

    Returns `(chi2, dof, p, table)`, where `table` is the 2x2
    `[[count_a, total_a - count_a], [count_b, total_b - count_b]]` array
    `chi2_contingency` was run on.
    """
    table = np.array([[count_a, total_a - count_a], [count_b, total_b - count_b]])
    chi2, p, dof, _ = chi2_contingency(table, correction=True)
    return chi2, dof, p, table


def compute_variant_classification_chi_squared_tests(stats_by_predictor):
    """For each predictor in `stats_by_predictor` (see
    `compute_variant_classification_stats`), run `_chi_squared_2x2` on each
    comparison in `VARIANT_CLASSIFICATION_CHI_SQUARED_COMPARISONS`: the
    Pathogenic/Likely Pathogenic and Benign/Likely Benign rate comparisons
    (gnomAD vs. ClinVar VUS, and Unobserved vs. ClinVar VUS / vs. gnomAD for
    Pathogenic/Likely Pathogenic) cited in the manuscript's Fig. 6d/e
    paragraphs.

    Each comparison's two groups' counts/totals are read from that
    predictor's own `{prefix}_total` and `{prefix}_{count_key}` stats (e.g.
    `gnomad_resolved_pathogenic` out of `gnomad_total`) -- the rate over the
    full category, not just the "resolved" (pathogenic-or-benign) subset --
    where `prefix` is each comparison's `group_a`/`group_b` second element.

    Returns `{predictor: [{label, group_a: {label, count, total}, group_b:
    {label, count, total}, table, chi2, dof, p}, ...]}`, one dict per
    comparison in `VARIANT_CLASSIFICATION_CHI_SQUARED_COMPARISONS`, same
    order.
    """
    results = {}
    for predictor, stats in stats_by_predictor.items():
        predictor_results = []
        for comparison in VARIANT_CLASSIFICATION_CHI_SQUARED_COMPARISONS:
            count_key = comparison["count_key"]
            a_label, a_prefix = comparison["group_a"]
            b_label, b_prefix = comparison["group_b"]
            count_a, total_a = stats[f"{a_prefix}_{count_key}"], stats[f"{a_prefix}_total"]
            count_b, total_b = stats[f"{b_prefix}_{count_key}"], stats[f"{b_prefix}_total"]
            chi2, dof, p, table = _chi_squared_2x2(count_a, total_a, count_b, total_b)
            predictor_results.append(
                {
                    "label": comparison["label"],
                    "group_a": {"label": a_label, "count": count_a, "total": total_a},
                    "group_b": {"label": b_label, "count": count_b, "total": total_b},
                    "table": table,
                    "chi2": chi2,
                    "dof": dof,
                    "p": p,
                }
            )
        results[predictor] = predictor_results
    return results


def format_variant_classification_chi_squared_tests(
    results_by_predictor, title=VARIANT_CLASSIFICATION_CHI_SQUARED_TITLE
):
    """Text report for `compute_variant_classification_chi_squared_tests`'s
    output: one block per predictor, one comparison per block, each showing
    both groups' count/total/pct, the 2x2 contingency table the test was run
    on, and the resulting chi2 statistic, degrees of freedom, and p-value.
    """
    lines = [
        title,
        "Pearson's chi-squared test of independence, 2x2 contingency table "
        "([[group A count, group A total - count], [group B count, group B total - count]]), "
        "with Yates' continuity correction (R's chisq.test() default for a 2x2 table, "
        "equivalent to prop.test(..., correct = TRUE)).",
    ]
    for predictor, comparisons in results_by_predictor.items():
        lines.append("")
        lines.append(f"-- {predictor} --")
        for comparison in comparisons:
            group_a, group_b = comparison["group_a"], comparison["group_b"]
            table = comparison["table"]
            lines.append(comparison["label"] + ":")
            lines.append(f"  {group_a['label']}: " + _format_count_and_pct(group_a["count"], group_a["total"]))
            lines.append(f"  {group_b['label']}: " + _format_count_and_pct(group_b["count"], group_b["total"]))
            lines.append(f"  2x2 table: [[{table[0, 0]}, {table[0, 1]}], [{table[1, 0]}, {table[1, 1]}]]")
            lines.append(f"  chi2 = {comparison['chi2']:.4f}, df = {comparison['dof']}, p = {comparison['p']:.4g}")
    return "\n".join(lines)


GENE_DISCORDANCE_SHEET = VARIANT_CLASSIFICATION_CATEGORY_SHEETS["controls"]
GENE_DISCORDANCE_TOP_N = 5
DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL = "ClinVar P/LP -> Class_REVEL B/LB"
DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL = "ClinVar B/LB -> Class_REVEL P/LP"
GENE_DISCORDANCE_TITLE = (
    "=== Gene-level discordance: ExCALIBR/OddsPath + REVEL gene-specific classification "
    f"vs. ClinVar control label ({GENE_DISCORDANCE_SHEET}) ==="
)


def compute_gene_discordance_stats(workbook, top_n=GENE_DISCORDANCE_TOP_N):
    """Per-gene counts of control variants where the ExCALIBR/OddsPath + REVEL
    gene-specific classification (`Class_REVEL`) disagrees with the row's
    ClinVar pathogenic-or-benign control label (`clnsig_group_18_25`) --
    calling a ClinVar Pathogenic/Likely pathogenic control Benign/Likely
    benign, or vice versa.

    Reads `GENE_DISCORDANCE_SHEET` (`controls_REVEL_GeneSpecific`) from
    `workbook` (an open `pd.ExcelFile` over the controls file), deduplicated
    to one row per distinct DNA variant the same way
    `compute_variant_classification_stats` does (`distinct_dna_variants`).

    Returns `(by_gene, total_discordant, total_controls)`: `by_gene` is a
    DataFrame indexed by gene, with `total`,
    `DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL`, and
    `DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL` columns, covering every gene with
    at least one discordant variant, sorted by `total` descending (ties
    broken alphabetically by gene). `format_gene_discordance_summary` takes
    the top `top_n` rows of this.
    """
    controls_df = distinct_dna_variants(workbook.parse(GENE_DISCORDANCE_SHEET))
    clinvar_group = controls_df[CONTROLS_CLINVAR_GROUP_COL]
    class_revel = controls_df[CLASS_REVEL_COL]

    pathogenic_to_benign = clinvar_group.isin(PATHOGENIC_VALUES) & class_revel.isin(CLASS_BENIGN_VALUES)
    benign_to_pathogenic = clinvar_group.isin(BENIGN_VALUES) & class_revel.isin(CLASS_PATHOGENIC_VALUES)

    flags = pd.DataFrame(
        {
            DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL: pathogenic_to_benign,
            DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL: benign_to_pathogenic,
        }
    )
    discordant = flags.any(axis=1)
    by_gene = flags[discordant].groupby(controls_df.loc[discordant, GENE_COL]).sum().astype(int)
    by_gene["total"] = by_gene.sum(axis=1)
    by_gene = by_gene.sort_index().sort_values("total", ascending=False, kind="stable")
    by_gene = by_gene[["total", DISCORDANT_PATHOGENIC_TO_BENIGN_LABEL, DISCORDANT_BENIGN_TO_PATHOGENIC_LABEL]]

    return by_gene, int(discordant.sum()), len(controls_df)


def format_gene_discordance_summary(by_gene, total_discordant, total_controls, top_n=GENE_DISCORDANCE_TOP_N):
    lines = [
        GENE_DISCORDANCE_TITLE,
        f"Total discordant control variants: {total_discordant} of {total_controls}",
        f"Top {top_n} genes by discordant-variant count:",
    ]
    top = by_gene.head(top_n)
    if len(top):
        lines.append(top.to_string())
    return "\n".join(lines)


def has_high_spliceai_score(df, splice_score_cols=SPLICEAI_SCORE_COLS, threshold=SPLICEAI_SCORE_THRESHOLD):
    """True where any of `splice_score_cols` is >= `threshold`; a row with all
    four scores null (a complex delins near a splice junction can score `NaN`
    on all four, per `docs/splice_variant_filtering_pipeline.md`) is False,
    not unknown -- matches that doc's own `splice_variant` derivation.
    """
    return (df[splice_score_cols] >= threshold).any(axis=1)


def compute_consequence_splice_breakdown(workbook, category_sheets=VARIANT_CLASSIFICATION_CATEGORY_SHEETS):
    """For each of the five Supplementary Data 5 variant categories
    (`category_sheets`, default the REVEL sheets in
    `VARIANT_CLASSIFICATION_CATEGORY_SHEETS`), read that category's sheet from
    `workbook` (an open `pd.ExcelFile` over the controls file), deduplicate to
    one row per distinct DNA variant (`distinct_dna_variants`), and count
    those variants by `simplified_consequence` (null values grouped under
    `NO_CONSEQUENCE_LABEL`), split into those with every SpliceAI delta score
    null or below `SPLICEAI_SCORE_THRESHOLD` ("low") and those with at least
    one at or above it ("high") -- see `has_high_spliceai_score`.

    Returns `{category: DataFrame}`, each DataFrame indexed by
    `simplified_consequence` value with `low`/`high` integer columns. See
    `format_consequence_splice_breakdown_table`.
    """
    breakdown = {}
    for category, sheet_name in category_sheets.items():
        df = distinct_dna_variants(workbook.parse(sheet_name))
        consequence = df[SIMPLIFIED_CONSEQUENCE_COL].fillna(NO_CONSEQUENCE_LABEL)
        high = has_high_spliceai_score(df)
        counts = high.groupby(consequence).agg(["sum", "count"])
        counts["low"] = counts["count"] - counts["sum"]
        breakdown[category] = counts.rename(columns={"sum": "high"})[["low", "high"]].astype(int)
    return breakdown


def format_consequence_splice_breakdown_table(breakdown, title=CONSEQUENCE_SPLICE_BREAKDOWN_TITLE):
    """Table form of `compute_consequence_splice_breakdown`'s output: one row
    per `simplified_consequence` value (in `NO_CONSEQUENCE_LABEL`-last order),
    one column per category, each cell "{low} / {high}" -- distinct DNA
    variants with every SpliceAI score null or below `SPLICEAI_SCORE_
    THRESHOLD`, vs. with at least one at or above it. A trailing "Total" row
    sums each category's column.
    """
    categories = list(breakdown)
    consequences = sorted(
        set().union(*(set(table.index) for table in breakdown.values())),
        key=lambda consequence: (consequence == NO_CONSEQUENCE_LABEL, consequence),
    )

    def _cell(table, consequence):
        if consequence not in table.index:
            return "0 / 0"
        row = table.loc[consequence]
        return f"{row['low']} / {row['high']}"

    result = pd.DataFrame(
        {category: [_cell(breakdown[category], consequence) for consequence in consequences] for category in categories},
        index=consequences,
    )
    result.loc["Total"] = {
        category: f"{breakdown[category]['low'].sum()} / {breakdown[category]['high'].sum()}"
        for category in categories
    }
    lines = [
        title,
        f"Each cell: distinct DNA variants with every SpliceAI score ({', '.join(SPLICEAI_SCORE_COLS)}) "
        f"{SPLICEAI_LOW_LABEL} / with at least one SpliceAI score {SPLICEAI_HIGH_LABEL}.",
        result.to_string(),
    ]
    return "\n".join(lines)


def build_report_text(
    table,
    gene_breakdown,
    genomic_variant_summary,
    igvf_dataset_measurement_counts_summary,
    composite_score_datasets_summary,
    score_sections,
    clinical_sections,
    clinical_sections_mixed_year,
    calibration_summary,
    filter_funnel_summary,
    clingen_evidence_repository_summary,
    reclassification_sections,
    control_concordance_summary,
    variant_classification_summary,
    variant_classification_chi_squared_summary,
    gene_discordance_summary,
    consequence_splice_breakdown_summary,
    allow_clinvar_conflicts=False,
):
    conflict_note = (
        "conflicting/ambiguous ClinVar calls folded in via any-match"
        if allow_clinvar_conflicts
        else "conflicting/ambiguous ClinVar calls excluded"
    )
    parts = [
        "=== Dataset summary ===",
        table.to_string(),
        genomic_variant_summary,
        format_gene_breakdown(gene_breakdown),
        igvf_dataset_measurement_counts_summary,
        composite_score_datasets_summary,
        "=== Score coverage (REVEL, AlphaMissense, MutPred2) ===",
        *score_sections,
        f"=== Clinical attributes (ClinVar 2025, gnomAD; {conflict_note}) ===",
        *clinical_sections,
        f"=== Clinical attributes (ClinVar 2025, except ClinVar 2018 for BRCA1/PTEN/MSH2/TP53; gnomAD; {conflict_note}) ===",
        *clinical_sections_mixed_year,
        calibration_summary,
        filter_funnel_summary,
        clingen_evidence_repository_summary,
        "=== Reclassification agreement (Figure 4c) ===",
        *reclassification_sections,
        control_concordance_summary,
        variant_classification_summary,
        variant_classification_chi_squared_summary,
        gene_discordance_summary,
        consequence_splice_breakdown_summary,
    ]
    return "\n\n".join(parts)


@click.command(help=__doc__)
@click.argument(
    "condensed_file",
    required=False,
    default=DEFAULT_CONDENSED_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "metadata_file",
    required=False,
    default=DEFAULT_METADATA_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.argument(
    "expanded_file",
    required=False,
    default=DEFAULT_EXPANDED_FILE,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
)
@click.option(
    "--excalibr-calibrations-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_EXCALIBR_CALIBRATIONS_FILE,
    help=(
        f"Path to the ExCALIBR calibrations workbook (default {DEFAULT_EXCALIBR_CALIBRATIONS_FILE}), "
        f"read from its '{EXCALIBR_CALIBRATIONS_SHEET}' sheet, for the calibration-coverage section."
    ),
)
@click.option(
    "--controls-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CONTROLS_FILE,
    help=(
        f"Path to the controls workbook (default {DEFAULT_CONTROLS_FILE}), read from every sheet "
        f"whose name starts with '{CONTROLS_SHEET_PREFIX}', for the reclassification-agreement section."
    ),
)
@click.option(
    "--checkpoint-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CHECKPOINT_FILE,
    help=(
        f"Path to Variant_Classification_analysis.ipynb's pre-category-split checkpoint "
        f"(default {DEFAULT_CHECKPOINT_FILE}), for the filter-funnel section."
    ),
)
@click.option(
    "--chek2-file",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_CHEK2_FILE,
    help=f"Path to the CHEK2 QC workbook (default {DEFAULT_CHEK2_FILE}), for the filter-funnel section.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Optional path to also write the full report as a text file",
)
@click.option(
    "--merge-calm-genes",
    is_flag=True,
    default=False,
    help=(
        "Count CALM1, CALM2, and CALM3 as a single gene target ('CALM1/2/3') "
        "instead of three separate genes in the Meta-analyses gene column and "
        "the ExCALIBR calibration coverage section (the Dataset summary and "
        "Genes represented sections always merge them)."
    ),
)
@click.option(
    "--allow-clinvar-conflicts",
    is_flag=True,
    default=False,
    help=(
        "Fold variants with a conflicting or ambiguous ClinVar call -- "
        "disagreement between a pathogenic-leaning and benign-leaning call "
        "across a variant's measurements/DNA candidates, or ClinVar's own "
        "'Conflicting classifications of pathogenicity' call -- into the "
        "VUS/pathogenic-or-benign clinical-attribute buckets via any-match. "
        "By default such variants are excluded from both buckets into their "
        "own 'ClinVar conflict' bucket instead, mirroring the conflict "
        "handling in Analysis/Curation_summary_V5_cleaned.ipynb."
    ),
)
def main(
    condensed_file,
    metadata_file,
    expanded_file,
    excalibr_calibrations_file,
    controls_file,
    checkpoint_file,
    chek2_file,
    output,
    merge_calm_genes,
    allow_clinvar_conflicts,
):
    condensed = pd.read_csv(condensed_file, sep="\t", dtype=str, keep_default_na=False)
    metadata = load_dataset_metadata(metadata_file)

    try:
        stats, gene_breakdown = compute_all_stats_from_frame(condensed, metadata)
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    table = stats_to_dataframe(stats)
    igvf_dataset_measurement_counts_summary = format_igvf_dataset_measurement_counts(
        compute_igvf_dataset_measurement_counts(condensed, metadata)
    )
    composite_score_datasets_summary = format_composite_score_datasets(
        compute_composite_score_datasets(condensed, metadata, merge_calm_genes=merge_calm_genes)
    )

    expanded = pd.read_csv(expanded_file, sep="\t", dtype=str, keep_default_na=False)
    genomic_variant_summary = format_genomic_variant_count(expanded_file, expanded)
    score_sections, clinical_sections, clinical_sections_mixed_year = build_variant_level_reports(
        condensed, expanded, condensed_file, expanded_file, allow_clinvar_conflicts=allow_clinvar_conflicts
    )

    calibrations = pd.read_excel(excalibr_calibrations_file, sheet_name=EXCALIBR_CALIBRATIONS_SHEET)
    try:
        calibration_stats = compute_excalibr_calibration_stats(
            calibrations, metadata, merge_calm_genes=merge_calm_genes
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc
    calibration_summary = format_calibration_summary(calibration_stats)

    controls_workbook = pd.ExcelFile(controls_file)
    reclassification_sections = build_reclassification_report(controls_workbook)
    control_concordance_summary = format_control_concordance_report(compute_control_concordance(controls_workbook))
    variant_classification_stats_by_predictor = compute_variant_classification_stats(controls_workbook)
    variant_classification_summary = format_variant_classification_table(variant_classification_stats_by_predictor)
    variant_classification_chi_squared_summary = format_variant_classification_chi_squared_tests(
        compute_variant_classification_chi_squared_tests(variant_classification_stats_by_predictor)
    )

    gene_discordance_summary = format_gene_discordance_summary(*compute_gene_discordance_stats(controls_workbook))
    consequence_splice_breakdown_summary = format_consequence_splice_breakdown_table(
        compute_consequence_splice_breakdown(controls_workbook)
    )

    # Deliberately not dtype=str: the CHEK2 merge below matches
    # auth_reported_score/score by exact numeric equality, which only lines
    # up if both are parsed as their natural (float) type -- see
    # compute_reclassification_filter_funnel's docstring.
    checkpoint = pd.read_csv(checkpoint_file, low_memory=False)
    chek2 = pd.read_excel(chek2_file, header=0)
    funnel_steps = compute_reclassification_filter_funnel(expanded, checkpoint, chek2, condensed)
    filter_funnel_summary = format_reclassification_filter_funnel(funnel_steps)

    clingen_evidence_repository_summary = format_clingen_evidence_repository_summary(
        compute_clingen_evidence_repository_stats(checkpoint, chek2, controls_workbook)
    )

    report = build_report_text(
        table,
        gene_breakdown,
        genomic_variant_summary,
        igvf_dataset_measurement_counts_summary,
        composite_score_datasets_summary,
        score_sections,
        clinical_sections,
        clinical_sections_mixed_year,
        calibration_summary,
        filter_funnel_summary,
        clingen_evidence_repository_summary,
        reclassification_sections,
        control_concordance_summary,
        variant_classification_summary,
        variant_classification_chi_squared_summary,
        gene_discordance_summary,
        consequence_splice_breakdown_summary,
        allow_clinvar_conflicts=allow_clinvar_conflicts,
    )
    click.echo(report)

    if output:
        output.write_text(report + "\n")
        click.echo(f"\nWrote report to {output}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Ablation analysis for the biobank-analysis reclassification method
(`src/build_variant_reclassification_dataset.py`): what does functional
evidence (ExCALIBR/OddsPath) or predictor evidence (REVEL/AlphaMissense/
MutPred2, gene-specific calibration falling back to genome-wide) accomplish
on its own, compared to the pipeline's actual combined-evidence method?

Starts from the same checkpoint file and the same downstream exclusions as
`build_variant_reclassification_dataset.py` (`apply_notebook_exclusions`),
then builds three "arms" per predictor:

- **Functional-only**: `Functional_points` alone.
- **Predictor-only**: that predictor's own `Points_<predictor>_GeneSpecific_
  GenomeWide` alone (gene-specific calibration, falling back to genome-wide,
  falling back to 0).
- **Combined**: `Functional_points + Points_<predictor>_GeneSpecific_
  GenomeWide` -- for REVEL this is exactly `build_variant_reclassification_
  dataset.py`'s own `Combined_points` column, byte-for-byte.

Since `Functional_points` doesn't depend on which predictor is active, the
functional-only arm is computed once and reused across all three predictor
sections.

**Each arm gets its own one-row-per-DNA-variant deduplication**, using that
arm's own points column via `dedup_by_max_abs_points`, rather than reusing
`build_variant_reclassification_dataset.py`'s Combined_points-based winner.
This matters: the representative assay row for a given variant can differ
between arms (e.g. one assay has the strongest functional evidence, a
different assay has the strongest REVEL evidence for the same variant), and
reusing the combined-evidence winner for a single-evidence-source arm would
silently understate what that source alone could show. This also means the
functional-only arm here is *not* a novel dedup rule: it's the same
single-pass, uniform, no-category-restriction `abs`-max-plus-`Dataset`-name-
tiebreak dedup `build_variant_reclassification_dataset.py` already uses for
its combined arm, just keyed on `Functional_points` instead of
`Combined_points` -- which also happens to match the *notebook's own*
dedup key (`Fxn_points`) for its category-specific `VUS`/`gnomAD`/
`Unobserved`/`controls`/`ClinGen_Repo` outputs (see
`docs/variant_classification.md`), even though the notebook's dedup is
otherwise more elaborate (nt/aa resolution preference, category splits) than
this script's single uniform pass.

Each arm is classified pathogenic-or-benign vs. uncertain using the same
fixed point cutoffs as the rest of this pipeline
(`points_are_pathogenic_or_benign`, from `src/mave_dataset_stats.py`: >=6 or
<=-1 total points), and summarized the same way as `mave_dataset_stats.py`'s
reclassification-export section: how many distinct DNA variants are
pathogenic-or-benign, how many ClinVar VUS are resolved, how many unobserved
(SNV) variants are resolved.

Finally, for each predictor, a three-way agree/disagree contingency table
(functional-only vs. predictor-only vs. combined, each resolved in
`--direction` or not) is reported twice: once over every variant, and once
restricted to a configurable scope of variant categories (`--scope`,
repeatable; default `vus`+`unobserved`, matching this pipeline's headline
claim):

- `vus`/`gnomad`/`unobserved` -- ClinVar VUS, gnomAD population variants, and
  unobserved variants, the previously-uncertain populations the pipeline's
  headline claim is about.
- `clinvar_control`/`clingen_control` -- ClinVar or ClinGen Evidence
  Repository *control* variants (already-known-truth calibration variants,
  reproducing `Variant_Classification_analysis.ipynb`'s own `controls`/
  `clingen` category membership filters -- see `is_clinvar_control`/
  `is_clingen_control`), for checking each arm's agreement with a known
  answer rather than counting newly-resolved uncertain variants. Not part of
  the default scope -- pass explicitly.
- `clinvar_control_missense_only`/`clingen_control_missense_only` -- the
  same two control populations, restricted to strictly missense variants
  (`is_missense_only`, i.e. `condensed_consequence == "missense_variant"`,
  excluding start-loss) -- one of the two consequence types REVEL/
  AlphaMissense/MutPred2 actually score, so the predictor-only and combined
  arms are a meaningful comparison rather than diluted by consequence types
  no predictor ever scores. Also not part of the default scope. (Not to be
  confused with `--consequence missense_only`, which restricts the *entire*
  analysis this way -- see below.)

The second table isolates how many of that scope's variants get resolved (in
`--direction`; default `any`, either direction -- pass `pathogenic` or
`benign` to restrict the contingency table/chart to just Pathogenic/Likely
Pathogenic or just Benign/Likely Benign outcomes), and in particular how many
are resolved by *combining* evidence but not by either source alone (true
synergy, not just one source doing all the work). Restricting `--direction`
also reshapes what counts as the "lost to combining" conflict bucket/chart
bar: with `--direction pathogenic`, it becomes variants a single arm alone
already classified Pathogenic/Likely Pathogenic *before* combining evidence
pulled the combined score back out of that range (into Uncertain, or even
across to Benign-leaning) -- not just any-direction sign conflicts.

Usage: `python -m src.ablation_variant_reclassification [checkpoint_file]
[--chek2-file path] [--consequence all|missense|missense_only]
[--predictor REVEL] [--predictor AlphaMissense] [--scope vus]
[--scope gnomad] [--scope unobserved] [--scope clinvar_control]
[--scope clingen_control] [--scope clinvar_control_missense_only]
[--scope clingen_control_missense_only]
[--direction any|pathogenic|benign] [--output path]
[--plot path] [--plot-comparison path] [--plot-redundant-gain path]
[--document path] [--document-chart ablation] [--document-split-dir dir]
[--document-split-format pdf|png|svg]`.
`--consequence missense` restricts the entire analysis -- every arm, scope,
direction, chart, and document section, not just the `*_missense_only`
scopes -- to missense/start-loss variants (`is_missense_or_start_loss`), the
only two consequence types REVEL/AlphaMissense/MutPred2 actually score;
`--consequence missense_only` restricts it further, to strictly
`condensed_consequence == "missense_variant"` (excluding start-loss, see
`is_missense_only` -- the same function the `*_missense_only` `--scope`
values use, but applied to the whole analysis rather than just those two
control scopes); omitting it (`all`, the default) applies no restriction.
`--predictor` is repeatable; omitting it runs all three.
`--scope` is repeatable; omitting it defaults to `vus`+`unobserved`.
`--direction` is single-select; omitting it defaults to `any`. `--plot`
renders the "added value of combining evidence" chart (one stacked bar per
predictor, restricted to the same scope/direction -- see
`build_predictor_section`/`save_synergy_chart`) to an image file (format
inferred from the extension, e.g. `.png`/`.svg`/`.pdf`). `--plot-comparison`
renders that same chart for Pathogenic/Likely Pathogenic and Benign/Likely
Benign side by side in one image, independently of `--direction` (see
`build_direction_comparison_chart_data`/`save_synergy_chart_comparison`).
`--plot-redundant-gain` drills into `--plot`'s "Each alone sufficient
(redundant)" band: a histogram, per predictor, of how many points combining
added on top of each single already-sufficient source (see
`build_redundant_gain_chart_data`/`save_redundant_gain_chart`).
`--document-split-dir` writes `--document`'s per-section chart blocks out as
separate files (one per chart-type/scope pair, each already carrying its own
legend) instead of one combined multi-section image -- see
`save_document_charts_as_files` and Extended Data Figure 10 in
`docs/figures.md`.
"""

from pathlib import Path

import click
import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from src.build_variant_reclassification_dataset import (
    DEFAULT_CHECKPOINT_FILE,
    DEFAULT_CHEK2_FILE,
    add_points_columns,
    apply_notebook_exclusions,
)
from src.lib.dedup import (
    GENOMIC_KEY_COLS,
    aa_dedup_or_mark,
    catch_mis_2,
    clingen_aa_sort_key,
    controls_aa_sort_key,
    dedup_by_max_abs_points,
)
from src.mave_dataset_stats import (
    LIKELY_BENIGN_POINTS_THRESHOLD,
    LIKELY_PATHOGENIC_POINTS_THRESHOLD,
    MIXED_YEAR_GENES,
    RECLASSIFICATION_ALT_COL,
    RECLASSIFICATION_CLINVAR_COL,
    RECLASSIFICATION_GNOMAD_COL,
    RECLASSIFICATION_REF_COL,
    VUS_VALUES,
    points_are_pathogenic_or_benign,
    split_genes,
)

DEFAULT_OUTPUT_FILE = None

FUNCTIONAL_POINTS_COL = "Functional_points"
FUNCTIONAL_ARM_LABEL = "Functional evidence only"

# Predictor label -> its own gene-specific-falling-back-to-genome-wide points
# column in the checkpoint file. `REVEL`'s combined arm reuses
# `build_variant_reclassification_dataset.py`'s own `Combined_points` column
# directly rather than recomputing it, so the two scripts can never drift
# apart on that one column.
PREDICTOR_POINTS_COLUMNS = {
    "REVEL": "Points_REVEL_GeneSpecific_GenomeWide",
    "AlphaMissense": "Points_AM_GeneSpecific_GenomeWide",
    "MutPred2": "Points_MP2_GeneSpecific_GenomeWide",
}
REVEL_COMBINED_POINTS_COL = "Combined_points"

# Which resolution direction counts as "resolved" for the contingency
# table/chart -- see `resolved_mask`. `"any"` (both directions) is this
# pipeline's usual pathogenic-or-benign framing and this script's default;
# `"pathogenic"`/`"benign"` restrict to one direction only, which also
# reshapes what counts as a "lost to combining" conflict (see module
# docstring).
DIRECTION_ANY = "any"
DIRECTION_PATHOGENIC = "pathogenic"
DIRECTION_BENIGN = "benign"
DIRECTION_CHOICES = [DIRECTION_ANY, DIRECTION_PATHOGENIC, DIRECTION_BENIGN]
DIRECTION_LABELS = {
    DIRECTION_ANY: "pathogenic-or-benign",
    DIRECTION_PATHOGENIC: "Pathogenic/Likely Pathogenic",
    DIRECTION_BENIGN: "Benign/Likely Benign",
}

# --- "Added value of combining evidence" chart (see save_synergy_chart) ---------------------

FUNCTIONAL_ALONE_LABEL = "Functional alone sufficient"
PREDICTOR_ALONE_LABEL = "Predictive alone sufficient"
BOTH_ALONE_LABEL = "Either alone sufficient"
SYNERGY_LABEL = "Neither alone sufficient"
CONFLICT_LABEL = "Lost to combining (evidence conflict)"
SYNERGY_CATEGORY_ORDER = [FUNCTIONAL_ALONE_LABEL, PREDICTOR_ALONE_LABEL, BOTH_ALONE_LABEL, SYNERGY_LABEL]
CONFLICT_ORDER = [CONFLICT_LABEL]

# Chart chrome, and the first four slots of the dataviz skill's validated
# categorical palette (`references/palette.md`) -- fixed order, never cycled.
# `CONFLICT_LABEL` deliberately uses the skill's reserved "critical" status
# color instead of a fifth categorical slot: it isn't a resolution category
# (combining evidence did *not* resolve those variants), so it's drawn as a
# separate below-the-baseline bar rather than stacked in with the four that
# are.
CHART_SURFACE = "#fcfcfb"
CHART_INK_PRIMARY = "#0b0b0b"
CHART_INK_SECONDARY = "#52514e"
CHART_INK_MUTED = "#898781"
CHART_GRIDLINE = "#e1e0d9"
CHART_BASELINE = "#c3c2b7"
CONFLICT_COLOR = "#d03b3b"
SYNERGY_CATEGORY_COLORS = {
    FUNCTIONAL_ALONE_LABEL: "#2a78d6",
    PREDICTOR_ALONE_LABEL: "#eb6834",
    BOTH_ALONE_LABEL: "#1baf7a",
    SYNERGY_LABEL: "#eda100",
    CONFLICT_LABEL: CONFLICT_COLOR,
}

# The three "alone sufficient" categories eligible for the `--show-class-
# upgrade` light/dark/hatched split (see `class_upgrade_status`) --
# `SYNERGY_LABEL` and `CONFLICT_LABEL` are unaffected, since "upgraded to a
# stronger class" only makes sense for a variant a single source (or, for
# `BOTH_ALONE_LABEL`, the stronger of the two) already resolved on its own.
CLASS_UPGRADE_CATEGORIES = (FUNCTIONAL_ALONE_LABEL, PREDICTOR_ALONE_LABEL, BOTH_ALONE_LABEL)

# The three possible outcomes `class_upgrade_status` assigns a variant,
# relative to the ACMG tier (`class_tier`) its already-sufficient single
# source(s) established alone: combining moved it up a tier, left it in the
# same tier, or -- despite the variant remaining `--direction`-resolved
# overall, so it's *not* `CONFLICT_LABEL` -- eroded it down a tier. See
# `--show-class-upgrade`'s help text for the visual encoding of each.
CLASS_UPGRADED = "upgraded"
CLASS_UNCHANGED = "unchanged"
CLASS_DOWNGRADED = "downgraded"
CLASS_UPGRADE_STATUS_ORDER = [CLASS_DOWNGRADED, CLASS_UNCHANGED, CLASS_UPGRADED]
CLASS_DOWNGRADE_HATCH = "///"


def _lighten(hex_color, amount=0.55):
    """Mix `hex_color` toward white by `amount` (0-1). Used by `--show-class-
    upgrade` for the "unchanged"/"downgraded" sub-shades; the unlightened
    base color (`SYNERGY_CATEGORY_COLORS`) represents "upgraded to a stronger
    class" -- so the chart looks identical to the default when the flag is
    off, and only gains the lighter tint (plus, for "downgraded", a hatch --
    `CLASS_DOWNGRADE_HATCH`) as an extra distinction.
    """
    r, g, b = (int(hex_color.lstrip("#")[i : i + 2], 16) for i in (0, 2, 4))
    r = round(r + (255 - r) * amount)
    g = round(g + (255 - g) * amount)
    b = round(b + (255 - b) * amount)
    return f"#{r:02x}{g:02x}{b:02x}"


CLASS_UPGRADE_LIGHT_COLORS = {
    category: _lighten(SYNERGY_CATEGORY_COLORS[category]) for category in CLASS_UPGRADE_CATEGORIES
}
CLASS_UPGRADE_CAPTION = 'Hatched = downgraded a class (still resolved) -- same light shade as "(unchanged)"'


def add_ablation_points_columns(df):
    """Add `Functional_points` (via `add_points_columns`, coerced numeric
    with missing treated as 0) plus one `_combined_points_<predictor>` column
    per `PREDICTOR_POINTS_COLUMNS` entry (`Functional_points` plus that
    predictor's own points, missing treated as 0). REVEL's combined column is
    aliased to the existing `Combined_points` rather than recomputed.
    """
    df = add_points_columns(df)
    df[FUNCTIONAL_POINTS_COL] = pd.to_numeric(df[FUNCTIONAL_POINTS_COL], errors="coerce").fillna(0)
    for predictor, points_col in PREDICTOR_POINTS_COLUMNS.items():
        df[points_col] = pd.to_numeric(df[points_col], errors="coerce").fillna(0)
        if predictor == "REVEL":
            continue
        df[combined_points_col(predictor)] = df[FUNCTIONAL_POINTS_COL] + df[points_col]
    return df


def combined_points_col(predictor):
    return REVEL_COMBINED_POINTS_COL if predictor == "REVEL" else f"_combined_points_{predictor}"


def build_arm(df, points_col):
    """Deduplicate `df` to one row per DNA variant by `abs(points_col)` (see
    `dedup_by_max_abs_points`) -- this arm's own winner, independent of any
    other arm's.
    """
    return dedup_by_max_abs_points(df, points_col=points_col, genomic_key_cols=GENOMIC_KEY_COLS)


def is_unobserved(df):
    """Same definition `mave_dataset_stats.py`'s reclassification-export
    section uses: an SNV with no ClinVar significance and no gnomAD MAF.
    """
    is_snv = (df[RECLASSIFICATION_REF_COL].str.len() == 1) & (df[RECLASSIFICATION_ALT_COL].str.len() == 1)
    return df[RECLASSIFICATION_CLINVAR_COL].isna() & df[RECLASSIFICATION_GNOMAD_COL].isna() & is_snv


def summarize_arm(deduped, points_col):
    """Distinct-DNA-variant classification summary for one arm, matching the
    shape of `mave_dataset_stats.compute_variant_classification_stats_from_
    reclassification_file`: total classified, how many pathogenic-or-benign,
    and ClinVar-VUS/unobserved resolved counts.
    """
    pathogenic_or_benign = points_are_pathogenic_or_benign(deduped[points_col])
    vus = deduped[RECLASSIFICATION_CLINVAR_COL].isin(VUS_VALUES)
    unobserved = is_unobserved(deduped)
    return {
        "total": len(deduped),
        "pathogenic_or_benign": int(pathogenic_or_benign.sum()),
        "vus_total": int(vus.sum()),
        "vus_resolved": int((vus & pathogenic_or_benign).sum()),
        "unobserved_total": int(unobserved.sum()),
        "unobserved_resolved": int((unobserved & pathogenic_or_benign).sum()),
    }


def genomic_key_index(deduped):
    """A string-typed MultiIndex over `GENOMIC_KEY_COLS`, safe to align/merge
    across arms -- `Chrom` mixes `int`/`str` values for the same chromosome
    across rows in practice (see `src/lib/dedup.py`), which silently defeats
    an exact-type comparison.
    """
    return pd.MultiIndex.from_frame(deduped[GENOMIC_KEY_COLS].astype(str))


def resolved_mask(points, direction=DIRECTION_ANY):
    """Boolean Series: does `points` cross this pipeline's fixed
    classification threshold in `direction` -- `DIRECTION_ANY` (either
    direction, `points_are_pathogenic_or_benign`'s >=6-or-<=-1 test),
    `DIRECTION_PATHOGENIC` (`>= LIKELY_PATHOGENIC_POINTS_THRESHOLD` only), or
    `DIRECTION_BENIGN` (`<= LIKELY_BENIGN_POINTS_THRESHOLD` only).
    """
    if direction == DIRECTION_PATHOGENIC:
        return points >= LIKELY_PATHOGENIC_POINTS_THRESHOLD
    if direction == DIRECTION_BENIGN:
        return points <= LIKELY_BENIGN_POINTS_THRESHOLD
    return points_are_pathogenic_or_benign(points)


def resolved_flags(deduped, points_col, direction=DIRECTION_ANY):
    """Boolean `direction`-classified flag (see `resolved_mask`) for each
    variant in `deduped`, indexed by `genomic_key_index` for cross-arm
    alignment.
    """
    return pd.Series(resolved_mask(deduped[points_col], direction).values, index=genomic_key_index(deduped))


def points_series(deduped, points_col):
    """Numeric points for each variant in `deduped`, indexed by
    `genomic_key_index` for cross-arm alignment -- the raw-value analog of
    `resolved_flags`, used by `redundant_gain_counts` to measure how many
    points combining added on top of an already-sufficient single source.
    """
    return pd.Series(pd.to_numeric(deduped[points_col], errors="coerce").values, index=genomic_key_index(deduped))


def control_truth_direction_series(arm, scopes):
    """Combined known pathogenic/benign truth direction across `scopes` (a
    subset of `CONTROL_TRUTH_FUNCTIONS`'s keys), indexed by
    `genomic_key_index` for cross-arm alignment -- the truth-direction
    analog of `points_series`. Where more than one requested scope assigns a
    variant a truth direction, `SCOPE_ORDER`'s earlier scope wins (ClinVar
    before ClinGen) -- a defensive, deterministic tiebreak for the rare case
    of a variant with conflicting ClinVar/ClinGen assertions; the two
    sources are independent, so this isn't expected to matter in practice.
    `scopes` containing no control scope (or `arm` having no matching
    variants) leaves every value `pandas.NA`, degrading gracefully to an
    all-empty concordance chart rather than an error.
    """
    direction = pd.Series(pd.NA, index=arm.index, dtype="object")
    for scope in SCOPE_ORDER:
        if scope not in scopes or scope not in CONTROL_TRUTH_FUNCTIONS:
            continue
        scope_direction = CONTROL_TRUTH_FUNCTIONS[scope](arm)
        direction = direction.where(direction.notna(), scope_direction)
    return pd.Series(direction.values, index=genomic_key_index(arm))


def class_tier(points):
    """Ordinal ACMG evidence-strength tier for `points` (see `points_series`),
    sign-agnostic: `2` = Pathogenic/Benign (points >= 10 or <= -7), `1` =
    Likely Pathogenic/Likely Benign (6-9 or -6..-1), `0` = Uncertain (0-5) or
    otherwise unresolved. Mirrors `Variant_Classification_analysis.ipynb`
    cell 38's `points_class` cutoffs (see `docs/variant_classification.md`'s
    evidence-scoring table), collapsed to "how many tiers deep" rather than
    which pathogenic/benign side -- used by `class_upgrade_flags` to detect
    when combining pushes a variant into a *stronger* tier than a single
    already-sufficient source established alone.
    """
    tier = pd.Series(0, index=points.index, dtype="Int64")
    tier = tier.mask(((points >= 6) & (points <= 9)) | ((points >= -6) & (points <= -1)), 1)
    tier = tier.mask((points >= 10) | (points <= -7), 2)
    return tier


def is_vus(df):
    return df[RECLASSIFICATION_CLINVAR_COL].isin(VUS_VALUES)


def is_gnomad(df):
    """gnomAD population variants: observed in gnomAD, regardless of ClinVar
    status -- reproduces `Variant_Classification_analysis.ipynb` cell 104's
    `gnomad = sankey_f[sankey_f['gnomad_MAF'].notna()]` exactly (a bare
    `gnomad_MAF.notna()` test, with no ClinVar-null requirement).

    This is **not** mutually exclusive with `is_vus`: a variant with both a
    ClinVar "Uncertain significance" call and a gnomAD MAF satisfies both
    (confirmed on the real checkpoint: 10,666 such rows) and legitimately
    appears in both `VUS_REVEL` and `gnomAD_REVEL` in the real Supplementary
    Data 5. An earlier version of this function additionally required
    `RECLASSIFICATION_CLINVAR_COL.isna()`, which doesn't match the notebook
    and silently undercounted the gnomAD population by tens of thousands of
    variants -- any ClinVar-annotated-but-MAF-observed row was wrongly
    excluded.
    """
    return df[RECLASSIFICATION_GNOMAD_COL].notna()


CLNSIG_GROUP_COL = "clnsig_group_18_25"
CLINVAR_CONFLICT_COL = "clinvar_conflict_flag_18_25"
CLINVAR_CONFLICT_VALUE = "has clinvar conflict"
CLINVAR_STAR_2018_COL = "clinvar_star_2018"
CLINVAR_STAR_2025_COL = "clinvar_star_2025"
CLINGEN_CLASSIFICATION_COL = "Updated_Classification_ClinGen_repo"
CLINGEN_VUS_VALUE = "VUS"

# Already on the checkpoint (`src/annotate_simplified_consequence.py`'s
# single-value, pipe-collapsed consequence term per row -- `"conflict_
# consequence"` where the underlying `simplified_consequence` had more than
# one disagreeing term). Used by `is_missense_or_start_loss` (for
# `--consequence missense`) and `is_missense_only` (for `--consequence
# missense_only` and the `*_missense_only` control scopes) to restrict to
# the consequence type(s) REVEL/AlphaMissense/MutPred2 actually score
# (predictor scoring is missense-SNV-only by construction -- see
# `vendor/variant-annotation/src/annotate_predictors.py` -- so every other
# consequence type is guaranteed unscored, making a predictor-only or
# combined arm meaningless for it).
CONDENSED_CONSEQUENCE_COL = "condensed_consequence"
MISSENSE_STARTLOSS_CONSEQUENCE_VALUES = frozenset({"missense_variant", "start_lost"})

CONTROLS_PATHOGENIC_VALUES = frozenset({"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"})
CONTROLS_BENIGN_VALUES = frozenset({"Benign", "Likely benign", "Benign/Likely benign"})
CONTROLS_CLINVAR_VALUES = CONTROLS_PATHOGENIC_VALUES | CONTROLS_BENIGN_VALUES

# Cell 79/87 of Variant_Classification_analysis.ipynb's `one_plus_stars` list
# -- doesn't include "practice guideline" (unobserved in this dataset), see
# src/lib/dedup.py's CLINVAR_REVIEW_STATUS_RANK for the equivalent full rank.
ONE_PLUS_STAR_VALUES = frozenset(
    {
        "criteria provided, single submitter",
        "criteria provided, multiple submitters, no conflicts",
        "reviewed by expert panel",
        "criteria provided, conflicting classifications",
    }
)


def mixed_year_star_series(df):
    """Mixed-year ClinVar review status: `clinvar_star_2018` for BRCA1/PTEN/
    MSH2/TP53, `clinvar_star_2025` otherwise -- the same per-gene switch
    `mave_dataset_stats.mixed_year_clinvar_series` applies to `clinvar_sig_*`,
    reproducing `Variant_Classification_analysis.ipynb` cell 79's
    `clinvar_star_18_25`.
    """
    use_2018 = df["Gene"].apply(lambda v: bool(MIXED_YEAR_GENES & set(split_genes(v))))
    return df[CLINVAR_STAR_2018_COL].where(use_2018, df[CLINVAR_STAR_2025_COL])


def is_clinvar_control(df):
    """ClinVar control variant, reproducing `Variant_Classification_
    analysis.ipynb` cells 78/80/85/88's row-level `controls` membership
    test: a mixed-year, conflict-resolved ClinVar Pathogenic/Likely
    Pathogenic/Benign/Likely Benign call (`clnsig_group_18_25`, already on
    the checkpoint -- plain per-row for nt-resolution candidates,
    per-amino-acid-group for aa-resolution ones), no ClinVar conflict
    (`clinvar_conflict_flag_18_25`), and 1+ star review status
    (`mixed_year_star_series`).

    This reproduces the notebook's row-level membership test only -- not its
    additional aa-group-level `has_clinvar_star_conflict` exclusion (which
    drops an entire amino-acid-position group if some of its rows are 1+
    star and others 0 star) or its controls-specific dedup. Like this
    script's other scopes, membership is checked per candidate row ahead of
    this script's own independent per-arm dedup, not reconstructed from the
    notebook's controls-specific pipeline.
    """
    star = mixed_year_star_series(df)
    return (
        df[CLNSIG_GROUP_COL].isin(CONTROLS_CLINVAR_VALUES)
        & (df[CLINVAR_CONFLICT_COL] != CLINVAR_CONFLICT_VALUE)
        & star.isin(ONE_PLUS_STAR_VALUES)
    )


def is_clingen_control(df):
    """ClinGen Evidence Repository control variant, reproducing
    `Variant_Classification_analysis.ipynb` cell 107's `clingen` filter: a
    non-null, non-VUS `Updated_Classification_ClinGen_repo` (already on the
    checkpoint).
    """
    classification = df[CLINGEN_CLASSIFICATION_COL]
    return classification.notna() & (classification != CLINGEN_VUS_VALUE)


def is_missense_or_start_loss(df):
    """Missense or start-loss variant (`condensed_consequence` in
    `MISSENSE_STARTLOSS_CONSEQUENCE_VALUES`) -- the only two consequence
    types REVEL/AlphaMissense/MutPred2 actually score. Used to restrict the
    `clinvar_control`/`clingen_control` scopes to a population where the
    predictor-only and combined arms are actually meaningful, rather than
    diluted by consequence types no predictor ever scores (which would make
    every such variant either `FUNCTIONAL_ALONE_LABEL` or unresolved by
    construction, not a genuine test of predictor or combined performance).
    """
    return df[CONDENSED_CONSEQUENCE_COL].isin(MISSENSE_STARTLOSS_CONSEQUENCE_VALUES)


MISSENSE_ONLY_CONSEQUENCE_VALUE = "missense_variant"


def is_missense_only(df):
    """Strictly `condensed_consequence == "missense_variant"` -- unlike
    `is_missense_or_start_loss`, this excludes start-loss variants entirely
    rather than grouping them in with missense. Used by `--consequence
    missense_only`, for a run that wants a strictly missense-only candidate
    universe (e.g. Extended Data Figure 10's missense-only panel set) rather
    than the broader missense-or-start-loss population `--consequence
    missense` restricts to.
    """
    return df[CONDENSED_CONSEQUENCE_COL] == MISSENSE_ONLY_CONSEQUENCE_VALUE


CONSEQUENCE_ALL = "all"
CONSEQUENCE_MISSENSE = "missense"
CONSEQUENCE_MISSENSE_ONLY = "missense_only"
CONSEQUENCE_CHOICES = [CONSEQUENCE_ALL, CONSEQUENCE_MISSENSE, CONSEQUENCE_MISSENSE_ONLY]
CONSEQUENCE_LABELS = {
    CONSEQUENCE_ALL: "all consequences",
    CONSEQUENCE_MISSENSE: "missense/start-loss only",
    CONSEQUENCE_MISSENSE_ONLY: "missense only (excluding start-loss)",
}


def restrict_to_consequence(df, consequence):
    """Restrict `df` to `is_missense_or_start_loss` rows when `consequence`
    is `CONSEQUENCE_MISSENSE`, to `is_missense_only` rows (excluding
    start-loss) when it's `CONSEQUENCE_MISSENSE_ONLY`, `df` unchanged when
    it's `CONSEQUENCE_ALL` (the default). Applied once, upfront in `main`,
    before any arm is built -- so unlike the `clinvar_control_missense_only`/
    `clingen_control_missense_only` `--scope` values, which restrict only
    their own control population, this restricts every scope, direction,
    chart, and document section to the same candidate universe.

    Only under `CONSEQUENCE_MISSENSE_ONLY` do the `clinvar_control_missense_
    only`/`clingen_control_missense_only` scopes become redundant with their
    unrestricted `clinvar_control`/`clingen_control` counterparts (the whole
    population is already restricted to strictly `is_missense_only`, the
    same test those two scopes apply on their own) -- expected, not a bug.
    Under `CONSEQUENCE_MISSENSE`, they stay meaningfully different: the whole
    population still includes start-loss variants (`is_missense_or_start_
    loss` is broader than `is_missense_only`), so `clinvar_control_missense_
    only`/`clingen_control_missense_only` still further exclude those from
    `clinvar_control`/`clingen_control`.
    """
    if consequence == CONSEQUENCE_MISSENSE:
        return df[is_missense_or_start_loss(df)]
    if consequence == CONSEQUENCE_MISSENSE_ONLY:
        return df[is_missense_only(df)]
    return df


def is_clinvar_control_missense_only(df):
    """`is_clinvar_control`, restricted to `is_missense_only` -- strictly
    `condensed_consequence == "missense_variant"`, excluding start-loss
    (unlike `--consequence missense`'s broader `is_missense_or_start_loss`).
    """
    return is_clinvar_control(df) & is_missense_only(df)


def is_clingen_control_missense_only(df):
    """`is_clingen_control`, restricted to `is_missense_only` -- see
    `is_clinvar_control_missense_only`.
    """
    return is_clingen_control(df) & is_missense_only(df)


CLINGEN_PATHOGENIC_VALUES = frozenset({"Pathogenic", "Likely Pathogenic"})
CLINGEN_BENIGN_VALUES = frozenset({"Benign", "Likely Benign"})


def clinvar_control_truth_direction(df):
    """Known pathogenic/benign truth direction for ClinVar control variants
    (see `is_clinvar_control`): `DIRECTION_PATHOGENIC`/`DIRECTION_BENIGN`
    (reused here as plain "which side" sentinels, not `--direction`'s
    resolved-threshold meaning) per row, `pandas.NA` where `is_clinvar_
    control` doesn't hold. Used by `concordance_status` to check whether an
    arm's own resolved direction agrees with this variant's known
    ClinVar-control truth.
    """
    in_control = is_clinvar_control(df)
    sig = df[CLNSIG_GROUP_COL]
    direction = pd.Series(pd.NA, index=df.index, dtype="object")
    direction[in_control & sig.isin(CONTROLS_PATHOGENIC_VALUES)] = DIRECTION_PATHOGENIC
    direction[in_control & sig.isin(CONTROLS_BENIGN_VALUES)] = DIRECTION_BENIGN
    return direction


def clingen_control_truth_direction(df):
    """Known pathogenic/benign truth direction for ClinGen Evidence
    Repository control variants (see `is_clingen_control`) --
    `clinvar_control_truth_direction`'s counterpart, from
    `Updated_Classification_ClinGen_repo`.
    """
    in_control = is_clingen_control(df)
    classification = df[CLINGEN_CLASSIFICATION_COL]
    direction = pd.Series(pd.NA, index=df.index, dtype="object")
    direction[in_control & classification.isin(CLINGEN_PATHOGENIC_VALUES)] = DIRECTION_PATHOGENIC
    direction[in_control & classification.isin(CLINGEN_BENIGN_VALUES)] = DIRECTION_BENIGN
    return direction


def clinvar_control_missense_only_truth_direction(df):
    """`clinvar_control_truth_direction`, restricted to `is_missense_only`
    (see `is_clinvar_control_missense_only`)."""
    return clinvar_control_truth_direction(df).where(is_missense_only(df))


def clingen_control_missense_only_truth_direction(df):
    """`clingen_control_truth_direction`, restricted to `is_missense_only`
    (see `is_clingen_control_missense_only`)."""
    return clingen_control_truth_direction(df).where(is_missense_only(df))


# `--scope` key -> known-truth-direction function, for the scopes that carry
# a ground truth to check concordance against (see `control_truth_direction_
# series`/`concordance_status`) -- a subset of `SCOPE_MASKS`'s keys, since
# `vus`/`gnomad`/`unobserved` have no known answer to check against by
# definition (that's the whole reason they need resolving in the first
# place).
CONTROL_TRUTH_FUNCTIONS = {
    "clinvar_control": clinvar_control_truth_direction,
    "clinvar_control_missense_only": clinvar_control_missense_only_truth_direction,
    "clingen_control": clingen_control_truth_direction,
    "clingen_control_missense_only": clingen_control_missense_only_truth_direction,
}


# --- `--special-control-dedup`: Supplementary Data 5's own ClinVar/ClinGen -----------------
# control dedup, reproduced here instead of this script's own generic,
# whole-population `dedup_by_max_abs_points` (see `special_control_dedup`).

CLINVAR_DATE_2018_COL = "clinvar_date_last_reviewed_2018"
CLINVAR_DATE_2025_COL = "clinvar_date_last_reviewed_2025"


def mixed_year_review_date_series(df):
    """Mixed-year ClinVar review date: `clinvar_date_last_reviewed_2018` for
    BRCA1/PTEN/MSH2/TP53, `clinvar_date_last_reviewed_2025` otherwise --
    `mixed_year_star_series`'s counterpart for review date, reproducing
    `Variant_Classification_analysis.ipynb` cell 79's `clinvar_date_last_
    reviewed_18_25`. Used only by `special_control_dedup`'s aa-stage
    recency tie-break (`src.lib.dedup.controls_aa_sort_key`).
    """
    use_2018 = df["Gene"].apply(lambda v: bool(MIXED_YEAR_GENES & set(split_genes(v))))
    return df[CLINVAR_DATE_2018_COL].where(use_2018, df[CLINVAR_DATE_2025_COL])


ZERO_STAR_VALUES = frozenset(
    {"no classification for the single variant", "no classification provided", "no assertion criteria provided"}
)


def summarize_clinvar_star_group(star_values):
    """Per-amino-acid-group ClinVar review-status classification,
    reproducing `Variant_Classification_analysis.ipynb` cell 87's
    `summarize_clnstar`: `"Unseen"` (every member's review status missing),
    `"has_clinvar_star_conflict"` (the group mixes 1+-star and 0-star
    members), `"one_plus_star"` (every rated member is 1+ star),
    `"zero_star"` (every rated member is 0 star). Only `"one_plus_star"`
    groups survive `special_control_dedup`'s ClinVar aa-branch (cell 88) --
    a group with even one 0-star member alongside a 1+-star one is dropped
    *entirely*, not just its 0-star rows: the group's amino-acid-level
    evidence can't be trusted to represent a reliable control when its own
    members disagree this sharply on review confidence.
    """
    values = {v for v in star_values if pd.notna(v)}
    if not values:
        return "Unseen"
    has_one_plus = any(v in ONE_PLUS_STAR_VALUES for v in values)
    has_zero = any(v in ZERO_STAR_VALUES for v in values)
    if has_one_plus and has_zero:
        return "has_clinvar_star_conflict"
    if has_one_plus:
        return "one_plus_star"
    if has_zero:
        return "zero_star"
    raise ValueError(f"Unexpected ClinVar review-status values encountered: {values}")


def controls_base(df):
    """Base ClinVar-control population *before* the nt/aa-branch-specific
    review-status filtering `special_control_dedup` applies next --
    reproduces `Variant_Classification_analysis.ipynb` cells 78+80's
    `controls` exactly (mixed-year `clnsig_group_18_25` Pathogenic/Benign
    call, no ClinVar conflict). Deliberately *not* `is_clinvar_control`,
    which additionally requires 1+ star at the row level: the notebook only
    applies that requirement per-branch below (a plain filter on the nt
    branch, a group-level exclusion on the aa branch), not uniformly on the
    base population the way `is_clinvar_control` does for this script's
    generic, non-special dedup path.
    """
    return df[df[CLNSIG_GROUP_COL].isin(CONTROLS_CLINVAR_VALUES) & (df[CLINVAR_CONFLICT_COL] != CLINVAR_CONFLICT_VALUE)]


NUCLEOTIDE_OR_AA_COL = "nucleotide_or_aa"
RESEQ_TRANSCRIPT_COL = "RefSeq Transcript ID"
RESEQ_TRANSCRIPT_STRIPPED_COL = "Ref_seq_transcript_ID_stripped"
AA_GROUP_COLS = ["Gene", "aa_pos", "aa_ref", "aa_alt", RESEQ_TRANSCRIPT_STRIPPED_COL]
MAX_FXN_PTS_VALUE = "max_fxn_pts"
MAX_PRED_PTS_VALUE = "max_pred_pts"
FXN_POINTS_COL = "Fxn_points"  # verbatim source of `Functional_points` (see `add_points_columns`)
CONTROL_DEDUP_STRATEGY = "nt_then_abs_max"  # Variant_Classification_analysis.ipynb's actual production setting
VARIANT_ROLE_COL = "Variant_Role"
VARIANT_ROLE_PRIMARY = "primary"

PREDICTOR_MAX_FLAG_COLUMNS = {
    "REVEL": "GeneSpecific_REVEL_max",
    "AlphaMissense": "GeneSpecific_AM_max",
    "MutPred2": "GeneSpecific_MP2_max",
}
# No AlphaMissense entry: the notebook excludes REVEL-/MutPred2-training
# variants from their own controls/ClinGen_Repo sheets (cells 91/113) but
# never does the equivalent for AlphaMissense.
PREDICTOR_TRAIN_AMINO_COLUMNS = {"REVEL": "revel_train_amino", "MutPred2": "mp2_train_amino"}
TRAIN_AMINO_YES_VALUE = "Yes"


def special_control_dedup(df, predictor, control_kind):
    """One row per DNA variant for `control_kind` (`"clinvar"` or
    `"clingen"`), restricted to `predictor`'s own gene-specific-maximum
    population -- reproduces `Variant_Classification_analysis.ipynb`'s
    actual production `controls_<predictor>_GeneSpecific`/`ClinGen_Repo_
    <predictor>_GeneSpecific` Supplementary Data 5 sheets (`CONTROLS_
    CLINGEN_DEDUP_STRATEGY = "nt_then_abs_max"`, cells 78-92/107-114)
    byte-for-byte, via the same shared `src.lib.dedup` mechanics the
    notebook itself uses.

    This script's own `dedup_by_max_abs_points` (used everywhere else, via
    `build_arm`) has no equivalent to this pipeline's amino-acid-level
    double-counting guard: several distinct NT variants can share one
    amino-acid change and therefore one shared amino-acid-level functional/
    predictor measurement, which a coordinate-only dedup would otherwise
    let each of those NT variants independently claim as if it were unique
    supporting evidence, inflating control-population counts.

    Two branches, split by `nucleotide_or_aa`:

    - **nt branch**: rows with a direct, NT-resolution measurement.
      Collapse exact-coordinate duplicates by `VariantNotes` tag order
      (same tie-break as `dedup_vus_gnomad_unobserved`'s `"v1"` strategy),
      then -- ClinVar only -- keep only 1+-star review-status rows
      (`ONE_PLUS_STAR_VALUES`).
    - **aa branch**: rows whose evidence only resolves to a shared
      amino-acid change (`AA_GROUP_COLS`). ClinVar only: drop any
      amino-acid group whose ClinVar review-status membership conflicts
      across its rows (`summarize_clinvar_star_group`) -- a group can't be
      trusted as a control when its own members disagree this sharply on
      review confidence. Then, restricted to rows already flagged as both
      the shared functional maximum (`VariantNotes == "max_fxn_pts"`) *and*
      `predictor`'s own per-variant maximum (`GeneSpecific_<predictor>_max
      == "max_pred_pts"`), resolve which single NT variant represents each
      remaining amino-acid group (`aa_dedup_or_mark`, `controls_aa_sort_
      key`/`clingen_aa_sort_key`) -- ClinVar prefers the highest-confidence,
      most-recently-reviewed record; ClinGen prefers the non-retracted,
      most-recently-approved/published record. Under `"nt_then_abs_max"`,
      `aa_dedup_or_mark` itself keeps every candidate, tagged `Variant_Role`
      "primary" (the chosen representative) or "secondary" -- this function
      then keeps only "primary" rows, matching the canonical, distributed
      `Supplementary_Data_5.xlsx` exactly (confirmed directly: 12,433 rows
      on `controls_REVEL_GeneSpecific` on the real dataset, vs. 23,322 in
      the separate `Supplementary_Data_5.with_secondary_variants.xlsx`
      export, which keeps every "secondary" row too) -- this is the actual
      double-counting guard: each amino-acid group contributes exactly one
      representative row, not one per constituent NT variant.

    The two branches are concatenated, `predictor`'s own training variants
    excluded (`PREDICTOR_TRAIN_AMINO_COLUMNS` -- no exclusion for
    AlphaMissense, matching the notebook), and finally `catch_mis_2`
    resolves the remaining case where the nt branch's own survivor and the
    aa branch's chosen representative land on the *exact same* genomic
    coordinate (an nt-resolution record always wins over an aa-resolution
    one for the same variant, under `"nt_then_abs_max"`).

    Returns a DataFrame with one row per surviving DNA variant, carrying
    `Functional_points`/that predictor's own points column/its combined
    points column (`add_ablation_points_columns` already added these
    upstream, unaffected by which dedup path a row survives through).
    """
    if control_kind == "clinvar":
        base = controls_base(df)
    elif control_kind == "clingen":
        base = df[is_clingen_control(df)]
    else:
        raise ValueError(f"Unknown control_kind: {control_kind!r} (expected 'clinvar' or 'clingen')")

    nt_rows = base[base[NUCLEOTIDE_OR_AA_COL] == "nt"]
    aa_rows = base[base[NUCLEOTIDE_OR_AA_COL] == "aa"].copy()

    nt_survivors = nt_rows.sort_values(by="VariantNotes", na_position="last").drop_duplicates(
        subset=GENOMIC_KEY_COLS, keep="first"
    )
    if control_kind == "clinvar":
        nt_survivors = nt_survivors[mixed_year_star_series(nt_survivors).isin(ONE_PLUS_STAR_VALUES)]

    aa_rows[RESEQ_TRANSCRIPT_STRIPPED_COL] = aa_rows[RESEQ_TRANSCRIPT_COL].str.replace(r"\.\d+$", "", regex=True)
    if control_kind == "clinvar":
        aa_rows["clinvar_star_18_25"] = mixed_year_star_series(aa_rows)
        group_status = aa_rows.groupby(AA_GROUP_COLS)["clinvar_star_18_25"].transform(summarize_clinvar_star_group)
        aa_rows = aa_rows[group_status == "one_plus_star"].copy()
        aa_rows["clinvar_date_last_reviewed_18_25"] = mixed_year_review_date_series(aa_rows)
        aa_sort_by, aa_ascending = controls_aa_sort_key(aa_rows, CONTROL_DEDUP_STRATEGY)
    else:
        aa_sort_by, aa_ascending = clingen_aa_sort_key(aa_rows, CONTROL_DEDUP_STRATEGY)

    max_flag_col = PREDICTOR_MAX_FLAG_COLUMNS[predictor]
    predictor_aa_candidates = aa_rows[
        (aa_rows["VariantNotes"] == MAX_FXN_PTS_VALUE) & (aa_rows[max_flag_col] == MAX_PRED_PTS_VALUE)
    ]
    aa_survivors = aa_dedup_or_mark(
        predictor_aa_candidates, AA_GROUP_COLS, GENOMIC_KEY_COLS, aa_sort_by, aa_ascending, CONTROL_DEDUP_STRATEGY
    )
    aa_survivors = aa_survivors[aa_survivors[VARIANT_ROLE_COL] == VARIANT_ROLE_PRIMARY]

    combined = pd.concat([nt_survivors, aa_survivors])
    train_amino_col = PREDICTOR_TRAIN_AMINO_COLUMNS.get(predictor)
    if train_amino_col is not None:
        combined = combined[combined[train_amino_col] != TRAIN_AMINO_YES_VALUE]

    return catch_mis_2(combined, GENOMIC_KEY_COLS, points_col=FXN_POINTS_COL, strategy=CONTROL_DEDUP_STRATEGY)


def special_control_population(df, predictor, scopes):
    """One-row-per-DNA-variant population covering whichever ClinVar/
    ClinGen control scopes are present in `scopes`, built from `special_
    control_dedup` instead of this script's generic, whole-population
    dedup -- the shared building block behind `--special-control-dedup`'s
    effect on both the ablation contingency table/chart and the control-
    concordance chart. Scopes outside the four ClinVar/ClinGen control keys
    (`vus`/`gnomad`/`unobserved`) are ignored here entirely -- callers
    combine this with the generic path's own result for those, so this
    function never needs to know about them. Returns `None` if `scopes`
    contains none of the four control keys.

    Requesting a control category's full scope (e.g. `clinvar_control`)
    together with its own `_missense_only` restriction is redundant,
    not additive -- exactly like `scope_mask`'s union semantics elsewhere in
    this script -- so this only ever runs `special_control_dedup` once per
    control kind actually needed: unrestricted if the full scope is
    present, missense-only-restricted (`is_missense_only`, excluding
    start-loss) only if just the restricted scope is requested without its
    parent.

    If both ClinVar and ClinGen control scopes are requested together, the
    two dedup passes are concatenated and, in the vanishingly rare case a
    variant satisfies both, deduplicated by genomic coordinate with ClinVar
    taking priority -- the same `SCOPE_ORDER` ClinVar-before-ClinGen
    tie-break `control_truth_direction_series` already uses elsewhere in
    this script, applied here so every downstream `.reindex()` call sees a
    unique index.
    """
    frames = []
    for control_kind, full_scope, missense_scope in (
        ("clinvar", "clinvar_control", "clinvar_control_missense_only"),
        ("clingen", "clingen_control", "clingen_control_missense_only"),
    ):
        if full_scope in scopes:
            deduped = special_control_dedup(df, predictor, control_kind)
        elif missense_scope in scopes:
            deduped = special_control_dedup(df, predictor, control_kind)
            deduped = deduped[is_missense_only(deduped)]
        else:
            continue
        frames.append(deduped)

    if not frames:
        return None
    combined = pd.concat(frames)
    genomic_key = combined[GENOMIC_KEY_COLS].astype(str)
    return combined[~genomic_key.duplicated(keep="first")]


def special_control_flags_for_scopes(df, predictor, scopes, direction):
    """The `flags` frame (see `build_contingency`) for `predictor`, built
    from `special_control_population`'s specially-deduped population
    instead of this script's generic per-arm dedup. Returns `None` if
    `scopes` contains no ClinVar/ClinGen control scope.
    """
    population = special_control_population(df, predictor, scopes)
    if population is None:
        return None
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)
    return pd.concat(
        {
            FUNCTIONAL_ARM_LABEL: resolved_flags(population, FUNCTIONAL_POINTS_COL, direction),
            f"{predictor} only": resolved_flags(population, points_col, direction),
            f"Functional + {predictor}": resolved_flags(population, combined_col, direction),
        },
        axis=1,
    )


SPECIAL_DEDUP_CONTROL_SCOPES = frozenset(
    {"clinvar_control", "clinvar_control_missense_only", "clingen_control", "clingen_control_missense_only"}
)


def restricted_flags_for_scope(
    df,
    predictor,
    predictor_arm,
    combined_arm,
    functional_flags,
    scopes,
    scope_flags,
    direction,
    special_control_dedup=False,
):
    """The scope-restricted `flags` frame (see `build_contingency`) for one
    predictor. By default (`special_control_dedup=False`), this is exactly
    the previous behavior: `predictor_flags_for_direction`'s flags,
    restricted to the precomputed `scope_flags` mask.

    When `special_control_dedup` is set, whichever of `scopes` is a
    ClinVar/ClinGen control scope (`SPECIAL_DEDUP_CONTROL_SCOPES`) instead
    uses `special_control_flags_for_scopes`'s own specially-deduped
    population for those scopes (see `special_control_dedup` for why the
    two differ) -- `vus`/`gnomad`/`unobserved` scopes are entirely
    unaffected either way, still computed via the generic path. The two
    results are safe to union (grouping by every `genomic_key_index` level
    and taking `.any()`): a variant's genomic-coordinate index means the
    same thing regardless of which population resolved it.
    """
    flags = predictor_flags_for_direction(predictor, predictor_arm, combined_arm, functional_flags, direction)
    if not special_control_dedup:
        return flags[scope_flags.reindex(flags.index, fill_value=False)]

    control_scopes = [s for s in scopes if s in SPECIAL_DEDUP_CONTROL_SCOPES]
    non_control_scopes = [s for s in scopes if s not in SPECIAL_DEDUP_CONTROL_SCOPES]

    frames = []
    if non_control_scopes:
        non_control_mask = scope_mask(predictor_arm, non_control_scopes)
        frames.append(flags[non_control_mask.reindex(flags.index, fill_value=False)])
    if control_scopes:
        special_flags = special_control_flags_for_scopes(df, predictor, control_scopes, direction)
        if special_flags is not None:
            frames.append(special_flags)
    if not frames:
        return flags.iloc[0:0]
    unioned = pd.concat(frames)
    return unioned.groupby(level=list(range(unioned.index.nlevels))).any()


# Category key -> membership predicate, and the fixed order/prose used to
# describe a `--scope` selection (see `describe_scope`). `clinvar_control`/
# `clingen_control` use the notebook's own mixed-year `clnsig_group_18_25`
# (not the plain `clinvar_sig_2025` the other three scopes use), so a
# variant can in principle satisfy both a control scope and `vus` at once
# for the four mixed-year genes (BRCA1/PTEN/MSH2/TP53) -- harmless under
# `scope_mask`'s union semantics, but worth knowing if combining a control
# scope with `vus`. Separately, `vus` and `gnomad` are *not* mutually
# exclusive with each other either -- see `is_gnomad`'s docstring.
#
# `clinvar_control_missense_only`/`clingen_control_missense_only`
# are strict subsets of `clinvar_control`/`clingen_control` (see `is_
# missense_only`): the same control population, restricted to strictly
# missense variants (excluding start-loss) -- one of the two consequence
# types REVEL/AlphaMissense/MutPred2 actually score, so the predictor-only
# and combined arms are a meaningful comparison rather than diluted by
# consequence types no predictor ever scores.
SCOPE_MASKS = {
    "vus": is_vus,
    "gnomad": is_gnomad,
    "unobserved": is_unobserved,
    "clinvar_control": is_clinvar_control,
    "clinvar_control_missense_only": is_clinvar_control_missense_only,
    "clingen_control": is_clingen_control,
    "clingen_control_missense_only": is_clingen_control_missense_only,
}
SCOPE_ORDER = [
    "vus",
    "gnomad",
    "unobserved",
    "clinvar_control",
    "clinvar_control_missense_only",
    "clingen_control",
    "clingen_control_missense_only",
]
SCOPE_LABELS = {
    "vus": "ClinVar VUS",
    "gnomad": "gnomAD population",
    "unobserved": "unobserved",
    "clinvar_control": "ClinVar control",
    "clinvar_control_missense_only": "ClinVar control (missense only)",
    "clingen_control": "ClinGen control",
    "clingen_control_missense_only": "ClinGen control (missense only)",
}
DEFAULT_SCOPES = ("vus", "unobserved")


def describe_scope(scopes):
    """English description of a `--scope` selection, e.g. `("unobserved",
    "vus")` -> `"ClinVar VUS or unobserved"`, `("vus", "gnomad",
    "unobserved")` -> `"ClinVar VUS, gnomAD population, or unobserved"`.
    Always rendered in `SCOPE_ORDER`, regardless of input order.
    """
    labels = [SCOPE_LABELS[s] for s in SCOPE_ORDER if s in scopes]
    if len(labels) == 1:
        return labels[0]
    if len(labels) == 2:
        return f"{labels[0]} or {labels[1]}"
    return ", ".join(labels[:-1]) + f", or {labels[-1]}"


def scope_mask(deduped, scopes):
    """Boolean flag for each variant in `deduped`, True if it falls into any
    of the requested `scopes` (subset of `SCOPE_MASKS`'s keys), indexed by
    `genomic_key_index`. ClinVar/gnomAD annotations are properties of the DNA
    variant's coordinates, not of whichever assay row an arm's dedup happened
    to pick, so this is the same regardless of which arm's deduped frame it's
    computed from.
    """
    mask = pd.Series(False, index=deduped.index)
    for scope in scopes:
        mask |= SCOPE_MASKS[scope](deduped)
    return pd.Series(mask.values, index=genomic_key_index(deduped))


def format_arm_summary(label, stats):
    lines = [
        f"--- {label} ---",
        f"Distinct DNA variants classified: {stats['total']}",
        "Pathogenic or benign: " + _format_count_and_pct(stats["pathogenic_or_benign"], stats["total"]),
        "ClinVar VUS resolved (reclassified pathogenic or benign): "
        + _format_count_and_pct(stats["vus_resolved"], stats["vus_total"]),
        "Unobserved variants resolved (classified pathogenic or benign): "
        + _format_count_and_pct(stats["unobserved_resolved"], stats["unobserved_total"]),
    ]
    return "\n".join(lines)


def _format_count_and_pct(count, total):
    pct = 100 * count / total if total else float("nan")
    return f"{count} of {total} ({pct:.1f}%)"


def build_contingency(title, flags):
    """Three-way agree/disagree table over `flags` (a boolean DataFrame with
    one column per arm, one row per variant) -- how many variants fall into
    each True/False combination across the arms, dropping the all-False
    (unresolved-by-anything) combination since it isn't informative here.
    """
    total = len(flags)
    lines = [title, f"Variants in scope: {total}"]
    if not total:
        return "\n".join(lines)
    counts = flags.groupby(list(flags.columns), dropna=False).size()
    counts = counts[counts.index.to_frame().any(axis=1)] if len(counts) else counts
    counts = counts.sort_values(ascending=False)
    table = counts.rename("count").reset_index()
    table["pct"] = (100 * table["count"] / total).round(1)
    lines.append(table.to_string(index=False))
    return "\n".join(lines)


def categorize_synergy(flags, functional_col, predictor_col, combined_col):
    """Map each row of a 3-boolean-column flags frame (see `build_contingency`)
    to one of `SYNERGY_CATEGORY_ORDER`/`CONFLICT_LABEL`'s five mutually
    exclusive labels, for `save_synergy_chart`. A variant resolved by nothing
    (all three flags `False`) gets `pandas.NA`, dropped from the chart --
    same population `build_contingency` already drops.
    """
    f, p, c = flags[functional_col], flags[predictor_col], flags[combined_col]
    category = pd.Series(pd.NA, index=flags.index, dtype="object")
    category[c & f & ~p] = FUNCTIONAL_ALONE_LABEL
    category[c & p & ~f] = PREDICTOR_ALONE_LABEL
    category[c & f & p] = BOTH_ALONE_LABEL
    category[c & ~f & ~p] = SYNERGY_LABEL
    category[~c & (f | p)] = CONFLICT_LABEL
    return category


def synergy_counts(flags, functional_col, predictor_col, combined_col):
    """Per-category variant counts (see `categorize_synergy`), reindexed to
    `SYNERGY_CATEGORY_ORDER + CONFLICT_ORDER` with missing categories filled
    with 0.
    """
    category = categorize_synergy(flags, functional_col, predictor_col, combined_col)
    return category.value_counts().reindex(SYNERGY_CATEGORY_ORDER + CONFLICT_ORDER, fill_value=0)


def class_upgrade_status(functional_points, predictor_points, combined_points, category):
    """For variants in one of `CLASS_UPGRADE_CATEGORIES` (a `categorize_
    synergy` label), did the *combined* score reach a stronger, the same, or
    a weaker ACMG evidence tier (`class_tier`) than the tier its already-
    sufficient single source(s) established alone? `BOTH_ALONE_LABEL`
    compares against the stronger (higher-tier) of the two individual
    sources -- "whichever was higher." A "weaker" (`CLASS_DOWNGRADED`)
    result is still `--direction`-resolved overall -- that's what makes it
    `CLASS_UPGRADE_CATEGORIES` rather than `CONFLICT_LABEL` -- it's only
    weaker than the single source's own tier, e.g. an opposing-but-
    insufficient-to-declassify second source eroding Pathogenic (tier 2)
    down to Likely Pathogenic (tier 1) while the variant stays resolved
    either way.

    Returns a Series of `CLASS_UPGRADED`/`CLASS_UNCHANGED`/`CLASS_DOWNGRADED`
    strings, `pandas.NA` outside those three categories.

    `functional_points`/`predictor_points`/`combined_points` must share the
    same index as `category` (see `points_series`/`categorize_synergy` --
    every arm dedups the same candidate universe, so they always do).
    """
    points = pd.concat(
        {"functional": functional_points, "predictor": predictor_points, "combined": combined_points}, axis=1
    )
    functional_tier = class_tier(points["functional"])
    predictor_tier = class_tier(points["predictor"])
    combined_tier = class_tier(points["combined"])
    higher_single_tier = pd.concat([functional_tier, predictor_tier], axis=1).max(axis=1)

    baseline = pd.Series(pd.NA, index=category.index, dtype="Int64")
    baseline = baseline.mask(category == FUNCTIONAL_ALONE_LABEL, functional_tier)
    baseline = baseline.mask(category == PREDICTOR_ALONE_LABEL, predictor_tier)
    baseline = baseline.mask(category == BOTH_ALONE_LABEL, higher_single_tier)

    status = pd.Series(pd.NA, index=category.index, dtype="object")
    in_scope = baseline.notna()
    status[in_scope & (combined_tier > baseline)] = CLASS_UPGRADED
    status[in_scope & (combined_tier == baseline)] = CLASS_UNCHANGED
    status[in_scope & (combined_tier < baseline)] = CLASS_DOWNGRADED
    return status


def class_upgrade_counts(functional_points, predictor_points, combined_points, category):
    """Per-`CLASS_UPGRADE_CATEGORIES`, per-`CLASS_UPGRADE_STATUS_ORDER`
    variant counts (see `class_upgrade_status`). Returns `{category_label:
    {"upgraded": int, "unchanged": int, "downgraded": int}}` for the three
    "alone sufficient" categories.
    """
    status = class_upgrade_status(functional_points, predictor_points, combined_points, category)
    counts = {}
    for label in CLASS_UPGRADE_CATEGORIES:
        in_category = category == label
        counts[label] = {
            outcome: int((in_category & (status == outcome)).sum()) for outcome in CLASS_UPGRADE_STATUS_ORDER
        }
    return counts


def predictor_class_upgrade_counts(predictor, functional_arm, predictor_arm, combined_arm, restricted_flags):
    """`class_upgrade_counts` for one predictor, from already-built arms and
    that predictor's already scope/direction-restricted `flags` frame (see
    `predictor_flags_for_direction`) -- the category labels
    (`categorize_synergy`) and the raw points are computed fresh here since
    `restricted_flags` only carries the boolean resolved flags.
    """
    category = categorize_synergy(
        restricted_flags, FUNCTIONAL_ARM_LABEL, f"{predictor} only", f"Functional + {predictor}"
    )
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)
    functional_points = points_series(functional_arm, FUNCTIONAL_POINTS_COL).reindex(category.index)
    predictor_points = points_series(predictor_arm, points_col).reindex(category.index)
    combined_points = points_series(combined_arm, combined_col).reindex(category.index)

    return class_upgrade_counts(functional_points, predictor_points, combined_points, category)


CONCORDANT = "concordant"
DISCORDANT = "discordant"
UNRESOLVED = "unresolved"
CONCORDANCE_STATUS_ORDER = [CONCORDANT, DISCORDANT, UNRESOLVED]

CONCORDANT_COLOR = "#0ca30c"  # dataviz skill's reserved "good" status color
DISCORDANT_COLOR = "#d03b3b"  # dataviz skill's reserved "critical" status color
UNRESOLVED_CONCORDANCE_COLOR = CHART_INK_MUTED
CONCORDANCE_STATUS_COLORS = {
    CONCORDANT: CONCORDANT_COLOR,
    DISCORDANT: DISCORDANT_COLOR,
    UNRESOLVED: UNRESOLVED_CONCORDANCE_COLOR,
}

CONCORDANCE_ARM_ORDER = ["functional", "predictor", "combined"]
CONCORDANCE_ARM_LABELS = {"functional": "Functional", "predictor": "Predictor", "combined": "Combined"}


def concordance_status(points, truth_direction):
    """Per-variant concordance status for one arm's points against a known
    control truth direction (see `control_truth_direction_series`):
    `CONCORDANT` (points resolved pathogenic/benign per `LIKELY_PATHOGENIC_
    POINTS_THRESHOLD`/`LIKELY_BENIGN_POINTS_THRESHOLD` and matching `truth_
    direction`), `DISCORDANT` (resolved but the opposite direction from
    truth), `UNRESOLVED` (didn't cross either threshold). `pandas.NA` where
    `truth_direction` itself is NA (not a control variant, or not in the
    requested control scope(s)).
    """
    resolved_pathogenic = points >= LIKELY_PATHOGENIC_POINTS_THRESHOLD
    resolved_benign = points <= LIKELY_BENIGN_POINTS_THRESHOLD

    status = pd.Series(pd.NA, index=points.index, dtype="object")
    has_truth = truth_direction.notna()
    truth_pathogenic = truth_direction == DIRECTION_PATHOGENIC
    truth_benign = truth_direction == DIRECTION_BENIGN

    status[has_truth & ((resolved_pathogenic & truth_pathogenic) | (resolved_benign & truth_benign))] = CONCORDANT
    status[has_truth & ((resolved_pathogenic & truth_benign) | (resolved_benign & truth_pathogenic))] = DISCORDANT
    status[has_truth & ~resolved_pathogenic & ~resolved_benign] = UNRESOLVED
    return status


def concordance_counts(functional_points, predictor_points, combined_points, truth_direction):
    """Per-arm (functional/predictor/combined) concordant/discordant/
    unresolved variant counts against `truth_direction` (see `concordance_
    status`). Returns `{"functional": {...}, "predictor": {...}, "combined":
    {...}}`, each a `{"concordant": int, "discordant": int, "unresolved":
    int}` dict.
    """
    arms = {"functional": functional_points, "predictor": predictor_points, "combined": combined_points}
    counts = {}
    for arm_label, points in arms.items():
        status = concordance_status(points, truth_direction)
        counts[arm_label] = {outcome: int((status == outcome).sum()) for outcome in CONCORDANCE_STATUS_ORDER}
    return counts


def predictor_control_concordance_counts(predictor, functional_arm, predictor_arm, combined_arm, truth_direction):
    """`concordance_counts` for one predictor, from already-built arms and a
    precomputed `truth_direction` (see `control_truth_direction_series`) --
    `predictor_class_upgrade_counts`'s counterpart for control concordance.
    """
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)
    functional_points = points_series(functional_arm, FUNCTIONAL_POINTS_COL).reindex(truth_direction.index)
    predictor_points = points_series(predictor_arm, points_col).reindex(truth_direction.index)
    combined_points = points_series(combined_arm, combined_col).reindex(truth_direction.index)
    return concordance_counts(functional_points, predictor_points, combined_points, truth_direction)


def build_predictor_arms(predictor, df):
    """Dedup the predictor-only and combined arms for `predictor` (see
    `build_arm`) once -- shared across every `--direction`, since arm dedup
    keys on `abs(points)`, which doesn't depend on direction; only the
    resolved-flag computation (`predictor_flags_for_direction`) does. Reused
    by `build_direction_comparison_chart_data` to avoid re-deduping per
    direction when rendering a multi-direction comparison chart.
    """
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)
    return build_arm(df, points_col), build_arm(df, combined_col)


def predictor_flags_for_direction(predictor, predictor_arm, combined_arm, functional_flags, direction):
    """The `flags` frame (see `build_contingency`) for one predictor, from
    already-built arms (`build_predictor_arms`) and `direction`.
    """
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)
    predictor_flags = resolved_flags(predictor_arm, points_col, direction)
    combined_flags = resolved_flags(combined_arm, combined_col, direction)
    return pd.concat(
        {
            FUNCTIONAL_ARM_LABEL: functional_flags,
            f"{predictor} only": predictor_flags,
            f"Functional + {predictor}": combined_flags,
        },
        axis=1,
    )


def build_predictor_section(
    df,
    predictor,
    predictor_arm,
    combined_arm,
    functional_flags,
    scopes,
    scope_flags,
    scope_label,
    direction=DIRECTION_ANY,
    special_control_dedup=False,
):
    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
    combined_col = combined_points_col(predictor)

    predictor_stats = summarize_arm(predictor_arm, points_col)
    combined_stats = summarize_arm(combined_arm, combined_col)

    flags = predictor_flags_for_direction(predictor, predictor_arm, combined_arm, functional_flags, direction)
    restricted_flags = restricted_flags_for_scope(
        df,
        predictor,
        predictor_arm,
        combined_arm,
        functional_flags,
        scopes,
        scope_flags,
        direction,
        special_control_dedup,
    )

    direction_label = DIRECTION_LABELS[direction]
    text = "\n\n".join(
        [
            f"=== Ablation: {predictor} ===",
            format_arm_summary(f"{predictor} only", predictor_stats),
            format_arm_summary(f"Functional + {predictor} (combined; the pipeline's method)", combined_stats),
            build_contingency(f"{predictor}: {direction_label} agreement, all variants", flags),
            build_contingency(
                f"{predictor}: {direction_label} agreement, {scope_label} variants only "
                "(added value of combining evidence)",
                restricted_flags,
            ),
        ]
    )
    return text, restricted_flags


def _chart_data_for_scope_direction(
    df,
    predictors,
    arms_by_predictor,
    functional_arm,
    functional_flags,
    scopes,
    scope_flags,
    scope_label,
    scope_total,
    direction,
    show_class_upgrade=False,
    special_control_dedup=False,
):
    """Build one `save_synergy_chart`-shaped `chart_data` dict (`{"scope_total",
    "scope_label", "direction_label", "predictor_flags", ["class_upgrade"]}`)
    for one scope/direction pair, from already-built arms. Shared by
    `build_direction_comparison_chart_data` and `build_document_chart_data`
    so arms are only ever deduped once regardless of how many scopes/
    directions/chart types end up needing a `chart_data` built from them.
    """
    predictor_flags = {}
    class_upgrade = {} if show_class_upgrade else None
    for predictor in predictors:
        predictor_arm, combined_arm = arms_by_predictor[predictor]
        restricted_flags = restricted_flags_for_scope(
            df,
            predictor,
            predictor_arm,
            combined_arm,
            functional_flags,
            scopes,
            scope_flags,
            direction,
            special_control_dedup,
        )
        predictor_flags[predictor] = restricted_flags
        if show_class_upgrade:
            class_upgrade[predictor] = predictor_class_upgrade_counts(
                predictor, functional_arm, predictor_arm, combined_arm, restricted_flags
            )

    chart_data = {
        "scope_total": scope_total,
        "scope_label": scope_label,
        "direction_label": DIRECTION_LABELS[direction],
        "predictor_flags": predictor_flags,
    }
    if show_class_upgrade:
        chart_data["class_upgrade"] = class_upgrade
    return chart_data


def _gain_data_for_predictors(
    df,
    predictors,
    arms_by_predictor,
    functional_points,
    functional_flags,
    scopes,
    scope_flags,
    direction,
    special_control_dedup=False,
):
    """Build `{predictor: redundant_gain_counts(...)}` (see
    `build_redundant_gain_chart_data`) for one scope/direction pair, from
    already-built arms. Shared with `build_document_chart_data` for the same
    reason as `_chart_data_for_scope_direction`.
    """
    gain_data = {}
    for predictor in predictors:
        points_col = PREDICTOR_POINTS_COLUMNS[predictor]
        combined_col = combined_points_col(predictor)
        predictor_arm, combined_arm = arms_by_predictor[predictor]

        restricted_flags = restricted_flags_for_scope(
            df,
            predictor,
            predictor_arm,
            combined_arm,
            functional_flags,
            scopes,
            scope_flags,
            direction,
            special_control_dedup,
        )
        category = categorize_synergy(
            restricted_flags, FUNCTIONAL_ARM_LABEL, f"{predictor} only", f"Functional + {predictor}"
        )
        both_alone_index = category[category == BOTH_ALONE_LABEL].index

        predictor_points = points_series(predictor_arm, points_col)
        combined_points = points_series(combined_arm, combined_col)
        gain_data[predictor] = redundant_gain_counts(
            functional_points, predictor_points, combined_points, both_alone_index
        )
    return gain_data


def build_ablation_report(
    df,
    predictors,
    scopes=DEFAULT_SCOPES,
    direction=DIRECTION_ANY,
    show_class_upgrade=False,
    special_control_dedup=False,
):
    """Returns `(report_text, chart_data)`. `chart_data` (see
    `save_synergy_chart`) is `{"scope_total": int, "scope_label": str,
    "direction_label": str, "predictor_flags": {predictor: restricted_flags}}`,
    reusing each predictor's already-built scope-restricted contingency flags
    rather than rebuilding the (expensive, per-arm-dedup) arms a second time.
    `scopes` is a subset of `SCOPE_MASKS`'s keys restricting the "added value
    of combining evidence" contingency table/chart to variants in the union
    of those categories. `direction` (`DIRECTION_ANY`/`DIRECTION_PATHOGENIC`/
    `DIRECTION_BENIGN`, see `resolved_mask`) restricts what counts as
    "resolved" for that same contingency table/chart -- the per-arm summary
    blocks (`format_arm_summary`) are unaffected, always reporting combined
    pathogenic-or-benign totals regardless of `direction`. `show_class_upgrade`
    additionally computes `chart_data["class_upgrade"]` (see
    `predictor_class_upgrade_counts`), for `--show-class-upgrade`'s light/dark
    chart split -- skipped by default since it costs an extra pass per
    predictor. `special_control_dedup` -- see `restricted_flags_for_scope`.
    """
    functional_arm = build_arm(df, FUNCTIONAL_POINTS_COL)
    functional_stats = summarize_arm(functional_arm, FUNCTIONAL_POINTS_COL)
    functional_flags = resolved_flags(functional_arm, FUNCTIONAL_POINTS_COL, direction)
    scope_flags = scope_mask(functional_arm, scopes)
    scope_label = describe_scope(scopes)

    sections = [format_arm_summary(FUNCTIONAL_ARM_LABEL, functional_stats)]
    predictor_flags = {}
    class_upgrade = {}
    for predictor in predictors:
        predictor_arm, combined_arm = build_predictor_arms(predictor, df)
        section_text, restricted_flags = build_predictor_section(
            df,
            predictor,
            predictor_arm,
            combined_arm,
            functional_flags,
            scopes,
            scope_flags,
            scope_label,
            direction,
            special_control_dedup,
        )
        sections.append(section_text)
        predictor_flags[predictor] = restricted_flags
        if show_class_upgrade:
            class_upgrade[predictor] = predictor_class_upgrade_counts(
                predictor, functional_arm, predictor_arm, combined_arm, restricted_flags
            )

    chart_data = {
        "scope_total": int(scope_flags.sum()),
        "scope_label": scope_label,
        "direction_label": DIRECTION_LABELS[direction],
        "predictor_flags": predictor_flags,
    }
    if show_class_upgrade:
        chart_data["class_upgrade"] = class_upgrade
    return "\n\n".join(sections), chart_data


def build_direction_comparison_chart_data(
    df,
    predictors,
    scopes=DEFAULT_SCOPES,
    directions=(DIRECTION_PATHOGENIC, DIRECTION_BENIGN),
    show_class_upgrade=False,
    special_control_dedup=False,
):
    """Like `build_ablation_report`'s `chart_data`, but for more than one
    `--direction` at once, for `save_synergy_chart_comparison`'s side-by-side
    chart. Reuses each predictor's arms (`build_predictor_arms`) across every
    requested direction rather than rebuilding (re-deduping) them per
    direction, since arm dedup doesn't depend on direction -- only the
    resolved-flag computation does. Returns `{direction: chart_data}`, one
    `chart_data` dict per `save_synergy_chart`'s shape. `show_class_upgrade`
    -- see `build_ablation_report`. `special_control_dedup` -- see
    `restricted_flags_for_scope`.
    """
    functional_arm = build_arm(df, FUNCTIONAL_POINTS_COL)
    scope_flags = scope_mask(functional_arm, scopes)
    scope_label = describe_scope(scopes)
    scope_total = int(scope_flags.sum())

    arms_by_predictor = {predictor: build_predictor_arms(predictor, df) for predictor in predictors}

    chart_data_by_direction = {}
    for direction in directions:
        functional_flags = resolved_flags(functional_arm, FUNCTIONAL_POINTS_COL, direction)
        chart_data_by_direction[direction] = _chart_data_for_scope_direction(
            df,
            predictors,
            arms_by_predictor,
            functional_arm,
            functional_flags,
            scopes,
            scope_flags,
            scope_label,
            scope_total,
            direction,
            show_class_upgrade,
            special_control_dedup,
        )
    return chart_data_by_direction


def redundant_gain_counts(functional_points, predictor_points, combined_points, both_alone_index):
    """For the variants at `both_alone_index` (the `BOTH_ALONE_LABEL` subset
    of a `categorize_synergy` call -- functional-only and predictor-only
    already independently resolve these), how many points combining added on
    top of each single source, as two integer-keyed histograms.

    "Points gained" is `abs(combined_points) - abs(single_source_points)`:
    the actual net change in evidence strength from adding the other source,
    not just that other source's own points restated. The two agree exactly
    when the two sources are same-signed (the usual case for a `--direction
    pathogenic`/`benign` run, where both are constrained to one sign by
    construction) -- but under `--direction any`, a variant can land in
    `BOTH_ALONE_LABEL` with the two sources pointing in *opposite*
    directions (each independently resolved, in different directions, while
    the combined score still happens to resolve too), in which case this
    formula correctly returns a small or negative "gain" (a partial
    erosion), where naively restating the other source's own points would
    not.

    Returns `{"gained_via_predictor": Series(int -> variant count),
    "gained_via_functional": Series(int -> variant count), "total": int}`.
    """
    functional_points = functional_points.reindex(both_alone_index)
    predictor_points = predictor_points.reindex(both_alone_index)
    combined_points = combined_points.reindex(both_alone_index)

    gained_via_predictor = (combined_points.abs() - functional_points.abs()).round().astype("Int64")
    gained_via_functional = (combined_points.abs() - predictor_points.abs()).round().astype("Int64")

    return {
        "gained_via_predictor": gained_via_predictor.value_counts().sort_index(),
        "gained_via_functional": gained_via_functional.value_counts().sort_index(),
        "total": len(both_alone_index),
    }


def build_redundant_gain_chart_data(
    df, predictors, scopes=DEFAULT_SCOPES, direction=DIRECTION_ANY, special_control_dedup=False
):
    """For each predictor, the points-gained histograms (`redundant_gain_
    counts`) for the `BOTH_ALONE_LABEL` band of the "added value of
    combining evidence" contingency table -- restricted to the same
    `scopes`/`direction`-restricted population `build_ablation_report`'s
    chart uses. Returns `{predictor: redundant_gain_counts(...)}`, for
    `save_redundant_gain_chart`. `special_control_dedup` -- see
    `restricted_flags_for_scope`.
    """
    functional_arm = build_arm(df, FUNCTIONAL_POINTS_COL)
    functional_flags = resolved_flags(functional_arm, FUNCTIONAL_POINTS_COL, direction)
    functional_points = points_series(functional_arm, FUNCTIONAL_POINTS_COL)
    scope_flags = scope_mask(functional_arm, scopes)
    arms_by_predictor = {predictor: build_predictor_arms(predictor, df) for predictor in predictors}

    return _gain_data_for_predictors(
        df,
        predictors,
        arms_by_predictor,
        functional_points,
        functional_flags,
        scopes,
        scope_flags,
        direction,
        special_control_dedup,
    )


CONTROL_SCOPES = tuple(CONTROL_TRUTH_FUNCTIONS)


def build_control_concordance_chart_data(df, predictors, scopes=CONTROL_SCOPES, special_control_dedup=False):
    """For each predictor, per-arm (functional/predictor/combined)
    concordant/discordant/unresolved counts (`predictor_control_concordance_
    counts`) against the known ClinVar/ClinGen control truth direction
    (`control_truth_direction_series`), restricted to `scopes` -- which
    control population(s) to check against; defaults to both. `scopes`
    outside `CONTROL_SCOPES` (e.g. `vus`) contribute no truth direction and
    so no concordance counts, by construction of `control_truth_direction_
    series`.

    `special_control_dedup` swaps the generic, whole-population dedup
    (`build_arm`) for `special_control_population`'s specially-deduped
    population (see `special_control_dedup` the function, for why the two
    differ) -- this chart is always control-scope-only already (`scopes`
    defaults to `CONTROL_SCOPES`), so unlike `restricted_flags_for_scope`
    there's no generic-path population to merge this against.

    Returns `{predictor: concordance_counts(...)}`, for `save_control_
    concordance_chart`.
    """
    if not special_control_dedup:
        functional_arm = build_arm(df, FUNCTIONAL_POINTS_COL)
        truth_direction = control_truth_direction_series(functional_arm, scopes)
        arms_by_predictor = {predictor: build_predictor_arms(predictor, df) for predictor in predictors}

        concordance_data = {}
        for predictor in predictors:
            predictor_arm, combined_arm = arms_by_predictor[predictor]
            concordance_data[predictor] = predictor_control_concordance_counts(
                predictor, functional_arm, predictor_arm, combined_arm, truth_direction
            )
        return concordance_data

    empty_counts = {arm: dict.fromkeys(CONCORDANCE_STATUS_ORDER, 0) for arm in CONCORDANCE_ARM_ORDER}
    concordance_data = {}
    for predictor in predictors:
        population = special_control_population(df, predictor, scopes)
        if population is None:
            concordance_data[predictor] = empty_counts
            continue
        truth_direction = control_truth_direction_series(population, scopes)
        points_col = PREDICTOR_POINTS_COLUMNS[predictor]
        combined_col = combined_points_col(predictor)
        concordance_data[predictor] = concordance_counts(
            points_series(population, FUNCTIONAL_POINTS_COL),
            points_series(population, points_col),
            points_series(population, combined_col),
            truth_direction,
        )
    return concordance_data


# Chart types `--document`/`build_document_chart_data`/`save_ablation_document`
# can include in each of a document's seven sections -- one per `SCOPE_ORDER`
# category. "ablation" is `--plot`'s chart (restricted to `direction`),
# "comparison" is `--plot-comparison`'s P/LP-vs-B/LB chart (independent of
# `direction`, always both), "concordance" is `--plot-control-concordance`'s
# per-arm concordant/discordant/unresolved chart (independent of `direction`;
# only nonempty for the `clinvar_control`/`clinvar_control_missense_only`/
# `clingen_control`/`clingen_control_missense_only` sections -- see
# `control_truth_direction_series`), "gain" is `--plot-redundant-
# gain`'s histogram (restricted to `direction`).
DOCUMENT_CHART_TYPES = ("ablation", "comparison", "concordance", "gain")


def build_document_chart_data(
    df,
    predictors,
    chart_types=DOCUMENT_CHART_TYPES,
    direction=DIRECTION_ANY,
    show_class_upgrade=False,
    special_control_dedup=False,
):
    """Build the chart data for `save_ablation_document`'s multi-section
    document: one section per `SCOPE_ORDER` category -- `vus`, `gnomad`,
    `unobserved`, `clinvar_control`, `clinvar_control_missense_only`,
    `clingen_control`, `clingen_control_missense_only` -- used
    *individually* here, unlike `--scope` elsewhere in this script, which
    unions whichever categories are requested into a single population. Each
    section gets whichever of `chart_types` were requested.

    Reuses each predictor's arms (`build_predictor_arms`) across every scope
    and chart type -- dedup doesn't depend on scope or direction, only the
    resolved-flag computation does -- and each direction's resolved-flags
    Series is itself computed at most once and reused across every section
    that needs it, so a full seven-section, three-chart-type document costs
    the same one-time arm-dedup pass as any single chart, not twenty-one.

    Returns `{scope: {"scope_label": str, "scope_total": int, "ablation":
    chart_data | None, "comparison": {direction: chart_data} | None,
    "concordance": {predictor: concordance_counts(...)} | None, "gain":
    gain_data | None}}` -- the value for a chart type not in `chart_types` is
    always `None`, so `save_ablation_document` can tell "not requested" apart
    from "requested but empty." "concordance" is only ever nonempty for the
    `clinvar_control`/`clingen_control` sections (see `control_truth_
    direction_series`) -- every other scope's section still gets a
    `"concordance"` dict when requested, just with every count zero.
    """
    functional_arm = build_arm(df, FUNCTIONAL_POINTS_COL)
    functional_points = points_series(functional_arm, FUNCTIONAL_POINTS_COL)
    arms_by_predictor = {predictor: build_predictor_arms(predictor, df) for predictor in predictors}

    needs_ablation = "ablation" in chart_types
    needs_comparison = "comparison" in chart_types
    needs_concordance = "concordance" in chart_types
    needs_gain = "gain" in chart_types

    directions_needed = set()
    if needs_ablation or needs_gain:
        directions_needed.add(direction)
    if needs_comparison:
        directions_needed.update((DIRECTION_PATHOGENIC, DIRECTION_BENIGN))
    functional_flags_by_direction = {
        d: resolved_flags(functional_arm, FUNCTIONAL_POINTS_COL, d) for d in directions_needed
    }

    document_data = {}
    for scope in SCOPE_ORDER:
        scope_flags = scope_mask(functional_arm, (scope,))
        scope_label = SCOPE_LABELS[scope]
        scope_total = int(scope_flags.sum())
        section = {
            "scope_label": scope_label,
            "scope_total": scope_total,
            "ablation": None,
            "comparison": None,
            "concordance": None,
            "gain": None,
        }

        if needs_ablation:
            section["ablation"] = _chart_data_for_scope_direction(
                df,
                predictors,
                arms_by_predictor,
                functional_arm,
                functional_flags_by_direction[direction],
                (scope,),
                scope_flags,
                scope_label,
                scope_total,
                direction,
                show_class_upgrade,
                special_control_dedup,
            )
        if needs_comparison:
            section["comparison"] = {
                d: _chart_data_for_scope_direction(
                    df,
                    predictors,
                    arms_by_predictor,
                    functional_arm,
                    functional_flags_by_direction[d],
                    (scope,),
                    scope_flags,
                    scope_label,
                    scope_total,
                    d,
                    show_class_upgrade,
                    special_control_dedup,
                )
                for d in (DIRECTION_PATHOGENIC, DIRECTION_BENIGN)
            }
        if needs_concordance:
            if special_control_dedup:
                empty_counts = {arm: dict.fromkeys(CONCORDANCE_STATUS_ORDER, 0) for arm in CONCORDANCE_ARM_ORDER}
                section["concordance"] = {}
                for predictor in predictors:
                    population = special_control_population(df, predictor, (scope,))
                    if population is None:
                        section["concordance"][predictor] = empty_counts
                        continue
                    truth_direction = control_truth_direction_series(population, (scope,))
                    points_col = PREDICTOR_POINTS_COLUMNS[predictor]
                    combined_col = combined_points_col(predictor)
                    section["concordance"][predictor] = concordance_counts(
                        points_series(population, FUNCTIONAL_POINTS_COL),
                        points_series(population, points_col),
                        points_series(population, combined_col),
                        truth_direction,
                    )
            else:
                truth_direction = control_truth_direction_series(functional_arm, (scope,))
                section["concordance"] = {
                    predictor: predictor_control_concordance_counts(
                        predictor,
                        functional_arm,
                        arms_by_predictor[predictor][0],
                        arms_by_predictor[predictor][1],
                        truth_direction,
                    )
                    for predictor in predictors
                }
        if needs_gain:
            section["gain"] = _gain_data_for_predictors(
                df,
                predictors,
                arms_by_predictor,
                functional_points,
                functional_flags_by_direction[direction],
                (scope,),
                scope_flags,
                direction,
                special_control_dedup,
            )

        document_data[scope] = section
    return document_data


def _relative_luminance(hex_color):
    """WCAG relative luminance of a `#rrggbb` color, 0 (black) to 1 (white)."""
    r, g, b = (int(hex_color.lstrip("#")[i : i + 2], 16) / 255 for i in (0, 2, 4))

    def linearize(channel):
        return channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4

    r, g, b = linearize(r), linearize(g), linearize(b)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _label_ink_for(hex_color):
    """White or `CHART_INK_PRIMARY` text, whichever contrasts more against `hex_color`."""
    luminance = _relative_luminance(hex_color)
    white_contrast = (1.0 + 0.05) / (luminance + 0.05)
    black_contrast = (luminance + 0.05) / (0.0 + 0.05)
    return "#ffffff" if white_contrast >= black_contrast else CHART_INK_PRIMARY


def synergy_chart_series(predictor, restricted_flags):
    return synergy_counts(restricted_flags, FUNCTIONAL_ARM_LABEL, f"{predictor} only", f"Functional + {predictor}")


def _draw_stacked_segment(ax, x_positions, heights, bottoms, color, label, totals, bar_width, hatch=None):
    """Draw one stacked-bar segment (one `ax.bar` call plus its >=3%-of-
    `totals` value labels, see `_draw_synergy_panel`) and return the updated
    `bottoms` for the next segment. `label=None` omits this segment from the
    legend -- used only for the "downgraded" sub-shade under `--show-class-
    upgrade`, which reuses the "unchanged" shade's own light color (plus a
    hatch) rather than adding a third legend entry; the "upgraded"/
    "unchanged" sub-shades each get their own explicit dark/light legend
    entry (see `_draw_synergy_panel`). `hatch` (e.g. `CLASS_DOWNGRADE_HATCH`)
    overlays a texture on the segment without needing a third color.
    """
    ax.bar(
        x_positions,
        heights,
        bottom=bottoms,
        width=bar_width,
        color=color,
        label=label,
        edgecolor=CHART_SURFACE,
        linewidth=1.5,
        hatch=hatch,
        zorder=3,
    )
    for xi, height, bottom, total in zip(x_positions, heights, bottoms, totals):
        if total and height / total >= 0.03:
            ax.text(
                xi,
                bottom + height / 2,
                f"{height:,}",
                ha="center",
                va="center",
                fontsize=9,
                color=_label_ink_for(color),
                zorder=4,
            )
    return [b + h for b, h in zip(bottoms, heights)]


def _draw_synergy_panel(ax, chart_data, show_ylabel=True, title=None):
    """Draw one `save_synergy_chart` panel onto `ax` -- see that function's
    docstring for the chart's content. Doesn't create a figure, set a
    figure-level title, add a legend, or save; callers handle those (see
    `save_synergy_chart`/`save_synergy_chart_comparison`) so several panels
    can share one title/legend in a multi-panel figure. `show_ylabel=False`
    omits the y-axis label -- used for every panel but the first in a
    side-by-side comparison, since it would otherwise repeat. `title`
    overrides the default full "Of N ... in scope, reclassified ..." subtitle
    -- used by `save_synergy_chart_comparison` to keep each panel's own title
    short (the "Of N ... in scope" part is identical across panels there, so
    it's hoisted into one shared figure-level line instead; the full text
    doesn't fit legibly within a single panel's width in a multi-panel
    figure).
    """
    predictor_counts = {
        predictor: synergy_chart_series(predictor, flags) for predictor, flags in chart_data["predictor_flags"].items()
    }
    scope_total = chart_data["scope_total"]
    scope_label = chart_data["scope_label"]
    direction_label = chart_data["direction_label"]
    class_upgrade = chart_data.get("class_upgrade")
    predictors = list(predictor_counts)

    ax.set_facecolor(CHART_SURFACE)

    x_positions = list(range(len(predictors)))
    bar_width = 0.55
    bottoms = [0.0] * len(predictors)
    resolved_totals = [sum(int(predictor_counts[p][c]) for c in SYNERGY_CATEGORY_ORDER) for p in predictors]

    for category in SYNERGY_CATEGORY_ORDER:
        if class_upgrade is not None and category in CLASS_UPGRADE_CATEGORIES:
            for outcome in CLASS_UPGRADE_STATUS_ORDER:
                outcome_heights = [int(class_upgrade[p][category][outcome]) for p in predictors]
                bottoms = _draw_stacked_segment(
                    ax,
                    x_positions,
                    outcome_heights,
                    bottoms,
                    CLASS_UPGRADE_LIGHT_COLORS[category]
                    if outcome != CLASS_UPGRADED
                    else SYNERGY_CATEGORY_COLORS[category],
                    f"{category} ({outcome})" if outcome != CLASS_DOWNGRADED else None,
                    resolved_totals,
                    bar_width,
                    hatch=CLASS_DOWNGRADE_HATCH if outcome == CLASS_DOWNGRADED else None,
                )
        else:
            heights = [int(predictor_counts[p][category]) for p in predictors]
            bottoms = _draw_stacked_segment(
                ax,
                x_positions,
                heights,
                bottoms,
                SYNERGY_CATEGORY_COLORS[category],
                category,
                resolved_totals,
                bar_width,
            )

    for xi, total in zip(x_positions, resolved_totals):
        pct = 100 * total / scope_total if scope_total else float("nan")
        ax.text(
            xi,
            total + max(scope_total, 1) * 0.015,
            f"{total:,} resolved\n({pct:.1f}%)",
            ha="center",
            va="bottom",
            fontsize=9,
            color=CHART_INK_SECONDARY,
            zorder=4,
        )

    conflict_counts = [int(predictor_counts[p][CONFLICT_LABEL]) for p in predictors]
    max_conflict = max(conflict_counts, default=0)
    ax.bar(
        x_positions,
        [-count for count in conflict_counts],
        width=bar_width * 0.45,
        color=CONFLICT_COLOR,
        label=CONFLICT_LABEL,
        zorder=3,
    )
    for xi, count in zip(x_positions, conflict_counts):
        if count:
            # Offset scales with max_conflict (the below-axis range this
            # panel's ylim bottom is itself set from), not scope_total (the
            # resolved-side range) -- under --direction pathogenic/benign,
            # conflicts are rare relative to the scope population, and a
            # scope_total-scaled offset pushed this label far below the tiny
            # conflict bar, past the axis and down near the x-tick labels.
            ax.text(
                xi,
                -count - max(max_conflict, 1) * 0.06,
                f"{count:,}",
                ha="center",
                va="top",
                fontsize=8,
                color=CONFLICT_COLOR,
                zorder=4,
            )

    max_resolved = max(resolved_totals, default=0)
    ax.set_ylim(bottom=-(max_conflict * 1.8 or 1), top=max_resolved * 1.3 or 1)

    ax.axhline(0, color=CHART_BASELINE, linewidth=1, zorder=2)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(predictors, fontsize=11, color=CHART_INK_PRIMARY)
    if show_ylabel:
        ax.set_ylabel(f"Distinct DNA variants ({scope_label})", fontsize=10, color=CHART_INK_SECONDARY)
    ax.yaxis.grid(True, color=CHART_GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", which="both", length=0, colors=CHART_INK_MUTED, labelsize=9)
    if title is None:
        title = f"Of {scope_total:,} {scope_label} DNA variants in scope, reclassified {direction_label}"
    ax.set_title(title, fontsize=10, color=CHART_INK_SECONDARY, pad=14)


def save_synergy_chart(chart_data, output_path):
    """Render the "added value of combining evidence" chart to `output_path`
    (format inferred from its extension) -- see the module docstring and
    `docs/ablation_variant_reclassification.md`. One stacked bar per
    predictor, broken into `SYNERGY_CATEGORY_ORDER`'s four resolution
    categories; a separate below-the-baseline bar per predictor for
    `CONFLICT_LABEL` (variants combining evidence did *not* resolve, even
    though a single source alone would have).

    Segment labels are shown only for categories reaching >=3% of that
    predictor's resolved total, to avoid clutter from near-zero slivers; the
    full per-category counts remain available in the text report regardless.
    """
    predictors = list(chart_data["predictor_flags"])
    fig, ax = plt.subplots(figsize=(2.2 + 1.8 * len(predictors), 6.5), dpi=200)
    fig.patch.set_facecolor(CHART_SURFACE)
    _draw_synergy_panel(ax, chart_data)

    title = "Added value of combining functional + predictor evidence"
    title_y = 0.98
    if "class_upgrade" in chart_data:
        title += f"\n{CLASS_UPGRADE_CAPTION}"
        title_y = 1.08
    fig.suptitle(title, fontsize=13, color=CHART_INK_PRIMARY, y=title_y)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.14),
        ncol=3,
        frameon=False,
        fontsize=9,
        labelcolor=CHART_INK_SECONDARY,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor=CHART_SURFACE)
    plt.close(fig)


def save_synergy_chart_comparison(
    chart_data_by_direction, output_path, directions=(DIRECTION_PATHOGENIC, DIRECTION_BENIGN)
):
    """Render `save_synergy_chart`'s chart for two or more directions as
    side-by-side panels in one image (format inferred from `output_path`'s
    extension), with one shared title/legend below -- see
    `build_direction_comparison_chart_data`. Each panel keeps its own
    independent y-axis scale: Pathogenic/Likely Pathogenic and Benign/Likely
    Benign resolved counts can differ by an order of magnitude in practice
    (Benign/Likely Benign calls are far more common in this dataset), and a
    shared scale would crush the smaller panel flat.
    """
    panels = [chart_data_by_direction[direction] for direction in directions]
    predictors = list(panels[0]["predictor_flags"])
    # scope_total/scope_label are identical across every panel here (same
    # scope population, different direction), so they're shown once as a
    # shared figure-level line rather than repeated per panel -- see
    # `_draw_synergy_panel`'s `title` param.
    scope_total = panels[0]["scope_total"]
    scope_label = panels[0]["scope_label"]

    fig, axes = plt.subplots(1, len(panels), figsize=((2.2 + 1.8 * len(predictors)) * len(panels), 6.5), dpi=200)
    fig.patch.set_facecolor(CHART_SURFACE)
    axes = list(axes) if hasattr(axes, "__len__") else [axes]

    for i, (ax, chart_data) in enumerate(zip(axes, panels)):
        _draw_synergy_panel(ax, chart_data, show_ylabel=(i == 0), title=f"Reclassified {chart_data['direction_label']}")

    title_lines = [
        "Added value of combining functional + predictor evidence, by resolution direction",
        f"Of {scope_total:,} {scope_label} DNA variants in scope",
    ]
    if "class_upgrade" in panels[0]:
        title_lines.append(CLASS_UPGRADE_CAPTION)
    fig.suptitle("\n".join(title_lines), fontsize=13, color=CHART_INK_PRIMARY, y=1.06)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=3,
        frameon=False,
        fontsize=9,
        labelcolor=CHART_INK_SECONDARY,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor=CHART_SURFACE)
    plt.close(fig)


def _gain_chart_x_positions(gain_data):
    """The sorted, deduplicated set of integer "points gained" values across
    every predictor/series in `gain_data`, plus their x-axis slot index --
    shared between `save_redundant_gain_chart` and `save_ablation_document`
    so a document section's rows share one consistent x-axis.
    """
    all_values = sorted(
        {
            int(value)
            for predictor in gain_data
            for series_name in ("gained_via_predictor", "gained_via_functional")
            for value in gain_data[predictor][series_name].index
        }
    )
    return all_values, {value: i for i, value in enumerate(all_values)}


def _gain_legend_handles():
    """Explicit legend handles for the points-gained histogram's two series.

    Built from `Patch` objects rather than `ax.get_legend_handles_labels()`:
    when a section/scope has zero "both alone sufficient" variants,
    `_draw_gain_panel`'s `ax.bar([], [], ...)` calls produce an empty
    `BarContainer` with no patches to read a color from, and matplotlib's
    auto-generated legend proxy for an empty container silently falls back
    to a default color instead of the one actually passed to `color=` --
    both series would render as the same (wrong) color in the legend. Fixed
    colors sidestep that regardless of whether any given panel's bars are
    empty.
    """
    return [
        Patch(color=SYNERGY_CATEGORY_COLORS[FUNCTIONAL_ALONE_LABEL], label="Gained via functional evidence"),
        Patch(color=SYNERGY_CATEGORY_COLORS[PREDICTOR_ALONE_LABEL], label="Gained via predictor evidence"),
    ]


def _draw_gain_panel(ax, predictor, data, all_values, x_positions, bar_width=0.38):
    """Draw one predictor's points-gained histogram panel onto `ax` -- see
    `save_redundant_gain_chart`'s docstring for the chart's content. Doesn't
    set a figure-level title, x-axis label, or legend; callers handle those
    (see `save_redundant_gain_chart`/`save_ablation_document`) so several
    panels can share one legend/title in a multi-panel figure.
    """
    ax.set_facecolor(CHART_SURFACE)
    functional_heights = [int(data["gained_via_functional"].get(value, 0)) for value in all_values]
    predictor_heights = [int(data["gained_via_predictor"].get(value, 0)) for value in all_values]
    xs = [x_positions[value] for value in all_values]

    ax.bar(
        [x - bar_width / 2 for x in xs],
        functional_heights,
        width=bar_width,
        color=SYNERGY_CATEGORY_COLORS[FUNCTIONAL_ALONE_LABEL],
        label="Gained via functional evidence",
        zorder=3,
    )
    ax.bar(
        [x + bar_width / 2 for x in xs],
        predictor_heights,
        width=bar_width,
        color=SYNERGY_CATEGORY_COLORS[PREDICTOR_ALONE_LABEL],
        label="Gained via predictor evidence",
        zorder=3,
    )

    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels([str(value) for value in all_values], fontsize=9, color=CHART_INK_MUTED)
    ax.set_ylabel("Variants", fontsize=9, color=CHART_INK_SECONDARY)
    ax.yaxis.grid(True, color=CHART_GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", which="both", length=0, colors=CHART_INK_MUTED, labelsize=9)
    ax.set_title(
        f"{predictor} -- {data['total']:,} variants where either source alone already sufficed",
        fontsize=10,
        color=CHART_INK_PRIMARY,
        loc="left",
        pad=10,
    )


def save_redundant_gain_chart(gain_data, output_path):
    """Render the points-gained histogram panel (see
    `build_redundant_gain_chart_data`) to `output_path` (format inferred
    from its extension): one row per predictor, x-axis = exact integer
    points gained, y-axis = variant count, two dodged bars per x-value --
    how many points combining added on top of the functional-only arm
    (blue, `FUNCTIONAL_ALONE_LABEL`'s color) or the predictor-only arm
    (orange, `PREDICTOR_ALONE_LABEL`'s color) -- for the "Each alone
    sufficient (redundant)" band of the main synergy chart.
    """
    predictors = list(gain_data)
    all_values, x_positions = _gain_chart_x_positions(gain_data)

    fig, axes = plt.subplots(len(predictors), 1, figsize=(9, 3.4 * len(predictors)), dpi=200, squeeze=False)
    fig.patch.set_facecolor(CHART_SURFACE)
    axes = [row[0] for row in axes]

    for ax, predictor in zip(axes, predictors):
        _draw_gain_panel(ax, predictor, gain_data[predictor], all_values, x_positions)

    axes[-1].set_xlabel("Points gained by adding the other evidence source", fontsize=9, color=CHART_INK_SECONDARY)

    fig.suptitle(
        'Points gained by combining, within the "each alone sufficient" band',
        fontsize=13,
        color=CHART_INK_PRIMARY,
        y=0.995,
    )
    handles = _gain_legend_handles()
    fig.legend(
        handles,
        [handle.get_label() for handle in handles],
        loc="lower center",
        bbox_to_anchor=(0.5, -0.02),
        ncol=2,
        frameon=False,
        fontsize=9,
        labelcolor=CHART_INK_SECONDARY,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor=CHART_SURFACE)
    plt.close(fig)


def _draw_concordance_panel(ax, predictor, arm_counts, bar_width=0.6):
    """Draw one predictor's control-concordance panel onto `ax` -- see
    `save_control_concordance_chart`'s docstring for the chart's content.
    Doesn't set a figure-level title or legend; callers handle those so
    several panels can share one legend/title in a multi-panel figure.
    """
    ax.set_facecolor(CHART_SURFACE)
    x_positions = list(range(len(CONCORDANCE_ARM_ORDER)))
    totals = [sum(arm_counts[arm][status] for status in CONCORDANCE_STATUS_ORDER) for arm in CONCORDANCE_ARM_ORDER]
    bottoms = [0.0] * len(CONCORDANCE_ARM_ORDER)
    for status in CONCORDANCE_STATUS_ORDER:
        heights = [arm_counts[arm][status] for arm in CONCORDANCE_ARM_ORDER]
        bottoms = _draw_stacked_segment(
            ax,
            x_positions,
            heights,
            bottoms,
            CONCORDANCE_STATUS_COLORS[status],
            status.capitalize(),
            totals,
            bar_width,
        )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(
        [CONCORDANCE_ARM_LABELS[arm] for arm in CONCORDANCE_ARM_ORDER], fontsize=10, color=CHART_INK_PRIMARY
    )
    ax.set_ylabel("Control variants", fontsize=9, color=CHART_INK_SECONDARY)
    ax.yaxis.grid(True, color=CHART_GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", which="both", length=0, colors=CHART_INK_MUTED, labelsize=9)
    ax.set_title(predictor, fontsize=11, color=CHART_INK_PRIMARY, loc="left", pad=10)


def save_control_concordance_chart(concordance_data, output_path):
    """Render the control-concordance chart to `output_path` (format
    inferred from its extension) -- see `build_control_concordance_chart_
    data`: one panel per predictor, each a three-bar cluster (functional-
    only / predictor-only / combined arm), stacked into concordant (green,
    the arm's resolved direction matches the known ClinVar/ClinGen control
    truth), discordant (red, resolved but the opposite direction from
    truth), and unresolved (gray, didn't cross either resolution threshold)
    segments.

    Uses the dataviz skill's reserved status palette (good/critical), not
    the categorical palette the ablation/comparison/--show-class-upgrade
    charts use -- this chart encodes a correctness dimension (did the call
    match a known answer), not a category identity, so reusing the
    categorical colors here would collide with their existing meaning
    elsewhere in this script's charts.
    """
    predictors = list(concordance_data)
    first_arm_counts = concordance_data[predictors[0]]["functional"]
    total = sum(first_arm_counts[status] for status in CONCORDANCE_STATUS_ORDER)

    fig, axes = plt.subplots(1, len(predictors), figsize=(3.4 * len(predictors), 5.5), dpi=200, squeeze=False)
    fig.patch.set_facecolor(CHART_SURFACE)
    axes = list(axes[0])

    for i, (ax, predictor) in enumerate(zip(axes, predictors)):
        _draw_concordance_panel(ax, predictor, concordance_data[predictor])
        if i > 0:
            ax.set_ylabel("")

    fig.suptitle(
        "Control concordance: does each arm's call agree with the known ClinVar/ClinGen truth?\n"
        f"Of {total:,} ClinVar/ClinGen control variants in scope",
        fontsize=13,
        color=CHART_INK_PRIMARY,
        y=1.06,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.08),
        ncol=3,
        frameon=False,
        fontsize=9,
        labelcolor=CHART_INK_SECONDARY,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor=CHART_SURFACE)
    plt.close(fig)


# Approximate inches each document block's chart itself needs, matching
# that chart type's own standalone figure proportions (see save_synergy_
# chart/save_synergy_chart_comparison/save_control_concordance_chart/
# save_redundant_gain_chart) -- "gain" scales with the number of predictor
# rows, the other three don't. Each block's total height additionally
# reserves `_document_block_legend_height`'s room for that block's own
# legend (see `_draw_document_block`).
DOCUMENT_BLOCK_LEGEND_HEIGHT = 0.9
# "ablation"/"comparison" under --show-class-upgrade split each class-
# upgrade-eligible category into its own light/dark legend entry (see
# `_document_legend_handles`), roughly doubling that block's legend row
# count -- wrapping to two rows at the same `ncol` needs more room than one.
DOCUMENT_BLOCK_LEGEND_HEIGHT_CLASS_UPGRADE = 1.5
DOCUMENT_ABLATION_CHART_HEIGHT = 6.0
DOCUMENT_COMPARISON_CHART_HEIGHT = 6.0
DOCUMENT_CONCORDANCE_CHART_HEIGHT = 5.5
DOCUMENT_GAIN_ROW_HEIGHT = 2.6
DOCUMENT_GAIN_BASE_WIDTH = 9.0


def _document_block_legend_height(block_type, show_class_upgrade):
    """Inches reserved for one block's own legend band (see
    `_draw_document_block`) -- `DOCUMENT_BLOCK_LEGEND_HEIGHT_CLASS_UPGRADE`
    for "ablation"/"comparison" when `show_class_upgrade` is set (their
    legend gains extra light/dark entries then), `DOCUMENT_BLOCK_LEGEND_
    HEIGHT` otherwise.
    """
    if show_class_upgrade and block_type in ("ablation", "comparison"):
        return DOCUMENT_BLOCK_LEGEND_HEIGHT_CLASS_UPGRADE
    return DOCUMENT_BLOCK_LEGEND_HEIGHT


def _document_block_heights(chart_types, num_predictors, show_class_upgrade=False):
    """`(block_types, block_heights)` for one document section -- the subset
    of `DOCUMENT_CHART_TYPES` present in `chart_types`, in that fixed order,
    with each one's approximate total height in inches: its chart proportion
    (`DOCUMENT_*_CHART_HEIGHT`/`DOCUMENT_GAIN_ROW_HEIGHT`) plus its own
    legend band (`_document_block_legend_height`).
    """
    block_types = [t for t in DOCUMENT_CHART_TYPES if t in chart_types]
    chart_heights = {
        "ablation": DOCUMENT_ABLATION_CHART_HEIGHT,
        "comparison": DOCUMENT_COMPARISON_CHART_HEIGHT,
        "concordance": DOCUMENT_CONCORDANCE_CHART_HEIGHT,
        "gain": DOCUMENT_GAIN_ROW_HEIGHT * max(num_predictors, 1),
    }
    return block_types, [chart_heights[t] + _document_block_legend_height(t, show_class_upgrade) for t in block_types]


def _document_width(chart_types, num_predictors):
    """Figure width in inches -- wide enough for whichever requested chart
    type needs the most horizontal room. "comparison" needs two side-by-side
    ablation-chart-width panels; "ablation" needs one; "concordance" needs
    one predictor-width panel per predictor (see `save_control_concordance_
    chart`); "gain" uses its own fixed width (see `save_redundant_gain_
    chart`).
    """
    ablation_width = 2.2 + 1.8 * num_predictors
    concordance_width = 3.2 * num_predictors
    candidates = [DOCUMENT_GAIN_BASE_WIDTH]
    if "ablation" in chart_types:
        candidates.append(ablation_width)
    if "comparison" in chart_types:
        candidates.append(2 * ablation_width)
    if "concordance" in chart_types:
        candidates.append(concordance_width)
    return max(candidates)


def _document_legend_handles(chart_types, show_class_upgrade=False):
    """The fixed legend handles needed for `chart_types` -- the five
    `SYNERGY_CATEGORY_ORDER`/`CONFLICT_LABEL` swatches if "ablation" or
    "comparison" is requested (both use the same five), plus the three
    `CONCORDANCE_STATUS_COLORS` swatches if "concordance" is, plus
    `_gain_legend_handles`'s two if "gain" is. Called once per block with a
    single-element `chart_types` (see `_draw_document_block`) so each block
    only shows the handles relevant to its own chart type, not the union
    every other block type would also need.

    `show_class_upgrade` splits each `CLASS_UPGRADE_CATEGORIES` swatch into
    its own light ("(unchanged)") and dark ("(upgraded)") entries, mirroring
    exactly what `_draw_synergy_panel` draws in that case -- these handles
    are built independently of the axes (unlike a standalone chart's
    `ax.get_legend_handles_labels()`), so they must stay in sync by hand.
    """
    handles = []
    if "ablation" in chart_types or "comparison" in chart_types:
        for category in SYNERGY_CATEGORY_ORDER:
            if show_class_upgrade and category in CLASS_UPGRADE_CATEGORIES:
                handles.append(
                    Patch(color=CLASS_UPGRADE_LIGHT_COLORS[category], label=f"{category} ({CLASS_UNCHANGED})")
                )
                handles.append(Patch(color=SYNERGY_CATEGORY_COLORS[category], label=f"{category} ({CLASS_UPGRADED})"))
            else:
                handles.append(Patch(color=SYNERGY_CATEGORY_COLORS[category], label=category))
        handles.append(Patch(color=CONFLICT_COLOR, label=CONFLICT_LABEL))
    if "concordance" in chart_types:
        handles += [
            Patch(color=CONCORDANCE_STATUS_COLORS[status], label=status.capitalize())
            for status in CONCORDANCE_STATUS_ORDER
        ]
    if "gain" in chart_types:
        handles += _gain_legend_handles()
    return handles


def _draw_document_block(block_fig, block_type, block_height, section_data, predictors, show_class_upgrade=False):
    """Draw one chart-type block (see `DOCUMENT_CHART_TYPES`) onto its own
    subfigure within a document section -- see `save_ablation_document` --
    plus that block's own legend (just the handles its own chart type needs,
    see `_document_legend_handles`), rather than one legend shared across
    the whole document.

    The legend is reserved in a bottom margin carved out of `block_fig`
    itself via `subplots_adjust` (`block_height` includes `_document_block_
    legend_height`'s share for this), then anchored at `block_fig`'s own
    bottom edge (`bbox_to_anchor=(0.5, 0.0)`, figure-fraction coordinates) --
    so it never needs to extend past `block_fig`'s own bounds. An axes-
    anchored legend positioned outside a nested subfigure's own [0, 1]
    coordinate range doesn't have this guarantee: `bbox_inches="tight"` has
    been confirmed to silently clip such a legend at the subfigure's edge
    instead of expanding the saved figure to fit it.
    """
    block_fig.patch.set_facecolor(CHART_SURFACE)

    if block_type == "ablation":
        ax = block_fig.subplots(1, 1)
        _draw_synergy_panel(ax, section_data["ablation"])
    elif block_type == "comparison":
        axes = block_fig.subplots(1, 2)
        comparison_data = section_data["comparison"]
        for i, direction in enumerate((DIRECTION_PATHOGENIC, DIRECTION_BENIGN)):
            chart_data = comparison_data[direction]
            _draw_synergy_panel(
                axes[i], chart_data, show_ylabel=(i == 0), title=f"Reclassified {chart_data['direction_label']}"
            )
    elif block_type == "concordance":
        concordance_data = section_data["concordance"]
        axes = block_fig.subplots(1, len(predictors), squeeze=False)
        axes = list(axes[0])
        for i, (ax, predictor) in enumerate(zip(axes, predictors)):
            _draw_concordance_panel(ax, predictor, concordance_data[predictor])
            if i > 0:
                ax.set_ylabel("")
    else:
        gain_data = section_data["gain"]
        all_values, x_positions = _gain_chart_x_positions(gain_data)
        axes = block_fig.subplots(len(predictors), 1, squeeze=False)
        axes = [row[0] for row in axes]
        for ax, predictor in zip(axes, predictors):
            _draw_gain_panel(ax, predictor, gain_data[predictor], all_values, x_positions)
        axes[-1].set_xlabel("Points gained by adding the other evidence source", fontsize=9, color=CHART_INK_SECONDARY)

    legend_fraction = _document_block_legend_height(block_type, show_class_upgrade) / block_height
    block_fig.subplots_adjust(bottom=legend_fraction + 0.03)
    legend_handles = _document_legend_handles((block_type,), show_class_upgrade)
    block_fig.legend(
        legend_handles,
        [handle.get_label() for handle in legend_handles],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.0),
        ncol=min(len(legend_handles), 5),
        frameon=False,
        fontsize=8,
        labelcolor=CHART_INK_SECONDARY,
    )


DOCUMENT_HEADER_HEIGHT = 0.6


def save_ablation_document(
    document_data, output_path, chart_types=DOCUMENT_CHART_TYPES, predictors=None, show_class_upgrade=False
):
    """Render `build_document_chart_data`'s multi-section document to
    `output_path` (format inferred from its extension, e.g.
    `.png`/`.svg`/`.pdf`): a header (document title only, spanning the full
    document width) followed by one section per `SCOPE_ORDER` category --
    all seven, always, regardless of population size -- each containing
    whichever of `chart_types` were requested, stacked top to bottom in
    `DOCUMENT_CHART_TYPES` order (ablation, comparison, concordance, gain).
    Each block carries its own legend (see `_draw_document_block`) rather
    than the document showing one shared legend in its header.

    `predictors` orders the predictors within each section -- required if
    `document_data` might have empty (all-`None`) sections for some scope
    (nothing to infer predictor order from otherwise); safe to omit
    otherwise, since it's then inferred from the first section with content.

    `show_class_upgrade` must match whatever was passed to `build_document_
    chart_data` when building `document_data` -- it only controls whether
    each "ablation"/"comparison" block's own legend additionally splits into
    light/dark entries per class-upgrade-eligible category (see `_document_
    legend_handles`) and reserves the taller legend band that split needs
    (see `_document_block_legend_height`); it has no effect on the bars
    themselves; those are already split (or not) based on whether `document_
    data`'s own `chart_data` entries carry a `"class_upgrade"` key.
    """
    if predictors is None:
        for section_data in document_data.values():
            if section_data.get("ablation") is not None:
                predictors = list(section_data["ablation"]["predictor_flags"])
            elif section_data.get("comparison") is not None:
                first_direction_data = next(iter(section_data["comparison"].values()))
                predictors = list(first_direction_data["predictor_flags"])
            elif section_data.get("concordance") is not None:
                predictors = list(section_data["concordance"])
            elif section_data.get("gain") is not None:
                predictors = list(section_data["gain"])
            if predictors:
                break
    num_predictors = len(predictors)

    block_types, block_heights = _document_block_heights(chart_types, num_predictors, show_class_upgrade)
    section_title_height = 0.5
    section_height = section_title_height + sum(block_heights)
    width = _document_width(chart_types, num_predictors)
    body_height = section_height * len(SCOPE_ORDER)

    fig = plt.figure(figsize=(width, DOCUMENT_HEADER_HEIGHT + body_height), dpi=150)
    fig.patch.set_facecolor(CHART_SURFACE)
    header_fig, body_fig = fig.subfigures(nrows=2, ncols=1, height_ratios=[DOCUMENT_HEADER_HEIGHT, body_height])
    header_fig.patch.set_facecolor(CHART_SURFACE)
    body_fig.patch.set_facecolor(CHART_SURFACE)

    header_fig.suptitle(
        "Ablation analysis: added value of combining functional + predictor evidence",
        fontsize=16,
        color=CHART_INK_PRIMARY,
        fontweight="bold",
        y=0.55,
    )

    section_figs = body_fig.subfigures(nrows=len(SCOPE_ORDER), ncols=1)
    for section_fig, scope in zip(section_figs, SCOPE_ORDER):
        section_fig.patch.set_facecolor(CHART_SURFACE)
        section_data = document_data[scope]
        section_fig.suptitle(
            f"{section_data['scope_label']} -- {section_data['scope_total']:,} variants",
            fontsize=15,
            color=CHART_INK_PRIMARY,
            fontweight="bold",
            y=0.99,
        )

        block_figs = section_fig.subfigures(nrows=len(block_types), ncols=1, height_ratios=block_heights)
        if len(block_types) == 1:
            block_figs = [block_figs]
        for block_fig, block_type, block_height in zip(block_figs, block_types, block_heights):
            _draw_document_block(block_fig, block_type, block_height, section_data, predictors, show_class_upgrade)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, bbox_inches="tight", facecolor=CHART_SURFACE)
    plt.close(fig)


def document_chart_filename(chart_type, scope, fmt):
    return f"{chart_type}_{scope}.{fmt}"


def save_document_charts_as_files(document_data, output_dir, chart_types=DOCUMENT_CHART_TYPES, fmt="pdf"):
    """Render each (chart type, scope) block of `build_document_chart_data`'s
    per-scope sections to its own file under `output_dir`, one file per
    block (e.g. `ablation_vus.pdf`, `concordance_clinvar_control.pdf`) --
    an alternative to `save_ablation_document`'s single combined
    multi-section image, for a figure assembled by hand from individual
    panels (see Extended Data Figure 10, `docs/figures.md`).

    Each file is produced by that chart type's own standalone `save_*`
    function (`save_synergy_chart`/`save_synergy_chart_comparison`/
    `save_control_concordance_chart`/`save_redundant_gain_chart`), so it
    already carries its own legend exactly as it would as a standalone
    `--plot`/`--plot-comparison`/`--plot-control-concordance`/`--plot-
    redundant-gain` chart -- there's no document-specific legend layout to
    reproduce here.

    A "concordance" block is skipped (no file written) for a scope outside
    `CONTROL_TRUTH_FUNCTIONS`'s keys (`vus`/`gnomad`/`unobserved`): those
    sections carry no known truth to check against, so `--document`'s own
    equivalent block there is uninformative (every count zero, see
    `control_truth_direction_series`) -- this just omits the empty file
    rather than writing one nobody would use.

    Returns the list of paths written, in `SCOPE_ORDER` x `DOCUMENT_CHART_
    TYPES` order.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for scope in SCOPE_ORDER:
        section_data = document_data[scope]
        if "ablation" in chart_types and section_data["ablation"] is not None:
            path = output_dir / document_chart_filename("ablation", scope, fmt)
            save_synergy_chart(section_data["ablation"], path)
            written.append(path)
        if "comparison" in chart_types and section_data["comparison"] is not None:
            path = output_dir / document_chart_filename("comparison", scope, fmt)
            save_synergy_chart_comparison(section_data["comparison"], path)
            written.append(path)
        if (
            "concordance" in chart_types
            and section_data["concordance"] is not None
            and scope in CONTROL_TRUTH_FUNCTIONS
        ):
            path = output_dir / document_chart_filename("concordance", scope, fmt)
            save_control_concordance_chart(section_data["concordance"], path)
            written.append(path)
        if "gain" in chart_types and section_data["gain"] is not None:
            path = output_dir / document_chart_filename("gain", scope, fmt)
            save_redundant_gain_chart(section_data["gain"], path)
            written.append(path)
    return written


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
    "--consequence",
    type=click.Choice(CONSEQUENCE_CHOICES),
    default=CONSEQUENCE_ALL,
    help=(
        "Restrict the whole analysis to variants of a given simplified consequence, applied once "
        "upfront before any arm is built: 'missense' keeps only missense/start-loss variants "
        "(condensed_consequence in missense_variant/start_lost -- the only consequence types "
        "REVEL/AlphaMissense/MutPred2 actually score), 'missense_only' keeps strictly "
        "condensed_consequence == missense_variant (excluding start-loss), 'all' (default) applies no "
        "restriction. Unlike the 'clinvar_control_missense_only'/'clingen_control_missense_only' "
        "--scope values, which restrict only those two control populations, this restricts every scope, "
        "direction, chart, and document section."
    ),
)
@click.option(
    "--predictor",
    "predictors",
    multiple=True,
    type=click.Choice(list(PREDICTOR_POINTS_COLUMNS)),
    help="Restrict to one predictor (repeatable). Default: all three (REVEL, AlphaMissense, MutPred2).",
)
@click.option(
    "--scope",
    "scopes",
    multiple=True,
    type=click.Choice(SCOPE_ORDER),
    help=(
        "Restrict the 'added value of combining evidence' contingency table/chart to one or more "
        "variant categories (repeatable): 'vus' (ClinVar VUS), 'gnomad' (gnomAD population variants, "
        "no ClinVar call), 'unobserved' (no ClinVar call, not observed in gnomAD), 'clinvar_control' "
        "(ClinVar Pathogenic/Likely Pathogenic/Benign/Likely Benign control variants, 1+ star, no "
        "ClinVar conflict), 'clingen_control' (ClinGen Evidence Repository control variants), "
        "'clinvar_control_missense_only'/'clingen_control_missense_only' (the same two "
        "control populations, restricted to strictly missense variants, excluding start-loss -- one "
        "of the two consequence types REVEL/AlphaMissense/MutPred2 actually score). Default: vus, "
        "unobserved."
    ),
)
@click.option(
    "--direction",
    type=click.Choice(DIRECTION_CHOICES),
    default=DIRECTION_ANY,
    help=(
        "Restrict the 'added value of combining evidence' contingency table/chart to one resolution "
        "direction: 'pathogenic' (Pathogenic/Likely Pathogenic only), 'benign' (Benign/Likely Benign "
        "only), or 'any' (either direction, the default). Restricting this also reshapes the 'lost to "
        "combining' conflict bucket/chart bar to variants a single arm alone classified in that "
        "direction before combining evidence pulled the result back out of it. The per-arm summary "
        "blocks always report combined pathogenic-or-benign totals regardless of this option."
    ),
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_FILE,
    help="Optional path to also write the full report as a text file.",
)
@click.option(
    "--plot",
    "plot_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to also render the 'added value of combining evidence' chart, restricted to "
        "--direction (format inferred from the extension, e.g. .png/.svg/.pdf)."
    ),
)
@click.option(
    "--plot-comparison",
    "plot_comparison_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to render a side-by-side Pathogenic/Likely Pathogenic vs. Benign/Likely Benign "
        "comparison chart -- two independently-scaled panels in one image, one shared legend (format "
        "inferred from the extension). Independent of --direction, which only affects --plot/the text "
        "report; this always compares both directions regardless of --direction."
    ),
)
@click.option(
    "--plot-control-concordance",
    "plot_control_concordance_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to render a control-concordance chart: one panel per predictor, each a "
        "functional-only/predictor-only/combined bar cluster showing whether that arm's resolved call "
        "agrees (concordant), disagrees (discordant), or doesn't resolve (unresolved) against the known "
        "ClinVar/ClinGen control truth direction (format inferred from the extension). Restricted to "
        "whichever of 'clinvar_control'/'clinvar_control_missense_only'/'clingen_control'/"
        "'clingen_control_missense_only' are present in --scope, defaulting to all four if none "
        "is given via --scope."
    ),
)
@click.option(
    "--plot-redundant-gain",
    "plot_redundant_gain_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to render a supplementary chart drilling into the 'Each alone sufficient "
        "(redundant)' band of --plot's chart: for the variants there, a histogram of how many points "
        "combining added on top of each single already-sufficient source (format inferred from the "
        "extension). Restricted to the same --scope/--direction as --plot."
    ),
)
@click.option(
    "--show-class-upgrade",
    is_flag=True,
    default=False,
    help=(
        "In --plot/--plot-comparison, split the 'Functional/Predictive alone sufficient' and "
        "'Either alone sufficient' bars by what combining did to the variant's ACMG evidence class: a "
        "lighter shade, labeled '(unchanged)' in the legend (stayed in the same class, e.g. Likely "
        "Pathogenic), a darker shade, labeled '(upgraded)' (combining moved it to a stronger class, e.g. "
        "Pathogenic), or a hatched light shade (combining eroded it to a weaker class -- e.g. Pathogenic "
        "down to Likely Pathogenic -- while it stayed --direction-resolved overall, so it's not counted "
        "as 'lost to combining'). For 'Either alone sufficient', the comparison baseline is whichever "
        "single source's own class was higher."
    ),
)
@click.option(
    "--special-control-dedup",
    is_flag=True,
    default=False,
    help=(
        "For the 'clinvar_control'/'clinvar_control_missense_only'/'clingen_control'/"
        "'clingen_control_missense_only' scopes only, use the same amino-acid-level "
        "double-counting-aware dedup Supplementary Data 5 itself uses (nt-resolution evidence "
        "preferred outright over amino-acid-resolution evidence, and one representative NT variant "
        "chosen per shared amino-acid change) instead of this script's own generic, whole-population "
        "abs-max dedup -- see special_control_dedup's docstring. Applies to the text report's "
        "restricted contingency table, --plot/--plot-comparison/--plot-redundant-gain, "
        "--plot-control-concordance, and --document's per-section charts. Does not affect the "
        "'vus'/'gnomad'/'unobserved' scopes, the all-variants contingency table, or the per-arm "
        "summary stats blocks, which are unaffected either way."
    ),
)
@click.option(
    "--document",
    "document_path",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional path to render a single multi-section document (format inferred from the extension, "
        "e.g. .png/.svg/.pdf): one section per variant category -- vus, gnomad, unobserved, "
        "clinvar_control, clinvar_control_missense_only, clingen_control, "
        "clingen_control_missense_only -- all seven, individually, regardless of --scope (which "
        "only affects --plot/--plot-comparison/--plot-redundant-gain/the text report). Each section "
        "shows whichever --document-chart types were requested."
    ),
)
@click.option(
    "--document-chart",
    "document_chart_types",
    multiple=True,
    type=click.Choice(DOCUMENT_CHART_TYPES),
    help=(
        "Which chart type(s) to include in each --document section (repeatable): 'ablation' (--plot's "
        "chart, restricted to --direction), 'comparison' (--plot-comparison's P/LP-vs-B/LB chart, "
        "always both directions), 'concordance' (--plot-control-concordance's chart, always all "
        "control scopes; nonempty only in the clinvar_control/clinvar_control_missense_only/"
        "clingen_control/clingen_control_missense_only sections), 'gain' "
        "(--plot-redundant-gain's points-gained histogram, restricted to --direction). Default: all four. "
        "Also controls which chart types --document-split-dir writes."
    ),
)
@click.option(
    "--document-split-dir",
    "document_split_dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help=(
        "Optional directory to render --document's per-section chart blocks as separate files instead "
        "of one combined multi-section image: one file per (chart type, scope) pair, e.g. "
        "ablation_vus.pdf, concordance_clinvar_control.pdf -- each already carrying its own legend, "
        "identical to a standalone --plot/--plot-comparison/--plot-control-concordance/--plot-redundant-"
        "gain chart. Uses --document-chart to pick which chart types (default: all four) and always "
        "covers all seven SCOPE_ORDER categories, same as --document; --scope has no effect here. A "
        "'concordance' file is skipped for the vus/gnomad/unobserved scopes, which carry no known truth "
        "to check concordance against. See docs/figures.md's Extended Data Figure 10 section."
    ),
)
@click.option(
    "--document-split-format",
    type=click.Choice(["pdf", "png", "svg"]),
    default="pdf",
    help="File format for --document-split-dir's output files (default: pdf).",
)
def main(
    checkpoint_file,
    chek2_file,
    consequence,
    predictors,
    scopes,
    direction,
    output,
    plot_path,
    plot_comparison_path,
    plot_control_concordance_path,
    plot_redundant_gain_path,
    show_class_upgrade,
    special_control_dedup,
    document_path,
    document_chart_types,
    document_split_dir,
    document_split_format,
):
    df = pd.read_csv(checkpoint_file)
    df = apply_notebook_exclusions(df, chek2_file)
    df = add_ablation_points_columns(df)
    df = restrict_to_consequence(df, consequence)

    predictors = list(predictors) or list(PREDICTOR_POINTS_COLUMNS)
    control_scopes = tuple(s for s in scopes if s in CONTROL_SCOPES) or CONTROL_SCOPES
    scopes = tuple(s for s in SCOPE_ORDER if s in scopes) or DEFAULT_SCOPES
    report, chart_data = build_ablation_report(
        df, predictors, scopes, direction, show_class_upgrade, special_control_dedup
    )
    if consequence != CONSEQUENCE_ALL:
        report = f"Consequence restriction: {CONSEQUENCE_LABELS[consequence]}\n\n" + report
    click.echo(report)

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report + "\n")
        click.echo(f"\nWrote report to {output}")

    if plot_path:
        save_synergy_chart(chart_data, plot_path)
        click.echo(f"Wrote chart to {plot_path}")

    if plot_comparison_path:
        chart_data_by_direction = build_direction_comparison_chart_data(
            df, predictors, scopes, show_class_upgrade=show_class_upgrade, special_control_dedup=special_control_dedup
        )
        save_synergy_chart_comparison(chart_data_by_direction, plot_comparison_path)
        click.echo(f"Wrote comparison chart to {plot_comparison_path}")

    if plot_control_concordance_path:
        concordance_data = build_control_concordance_chart_data(df, predictors, control_scopes, special_control_dedup)
        save_control_concordance_chart(concordance_data, plot_control_concordance_path)
        click.echo(f"Wrote control-concordance chart to {plot_control_concordance_path}")

    if plot_redundant_gain_path:
        gain_data = build_redundant_gain_chart_data(df, predictors, scopes, direction, special_control_dedup)
        save_redundant_gain_chart(gain_data, plot_redundant_gain_path)
        click.echo(f"Wrote redundant-gain chart to {plot_redundant_gain_path}")

    if document_path or document_split_dir:
        chart_types = tuple(t for t in DOCUMENT_CHART_TYPES if t in document_chart_types) or DOCUMENT_CHART_TYPES
        document_data = build_document_chart_data(
            df, predictors, chart_types, direction, show_class_upgrade, special_control_dedup
        )
        if document_path:
            save_ablation_document(document_data, document_path, chart_types, predictors, show_class_upgrade)
            click.echo(f"Wrote document to {document_path}")
        if document_split_dir:
            written = save_document_charts_as_files(
                document_data, document_split_dir, chart_types, document_split_format
            )
            click.echo(f"Wrote {len(written)} chart files to {document_split_dir}")


if __name__ == "__main__":
    main()

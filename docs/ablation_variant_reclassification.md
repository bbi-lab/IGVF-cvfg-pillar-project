# Ablation: variant reclassification by evidence source

`src/ablation_variant_reclassification.py` answers a question
`src/build_variant_reclassification_dataset.py` doesn't: how much of the
reclassification pipeline's effect comes from functional evidence
(ExCALIBR/OddsPath) alone, how much from predictor evidence (REVEL/
AlphaMissense/MutPred2, gene-specific calibration falling back to
genome-wide) alone, and how much is genuine added value from combining the
two.

## Inputs

Same as `build_variant_reclassification_dataset.py`:

- `checkpoint_file` (positional, optional) -- `Variant_Classification_
  analysis.ipynb`'s cell-69 checkpoint, default
  `data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`.
- `--chek2-file` -- default `data/input/maves/CHEK2_Gebbia_2024.xlsx`.

The same downstream exclusions (`SFPQ`, the CHEK2 QC flag,
conflicting/unmeasured-splice `VariantNotes` tags -- see
`src/build_variant_reclassification_dataset.py`'s own module docstring) are
re-applied before any ablation arm is built, so every arm starts from the
same candidate row set.

## `--consequence`: restricting the whole analysis to one simplified consequence

`--consequence missense` (`restrict_to_consequence`) restricts the entire
candidate row set to missense/start-loss variants (`condensed_consequence`
in `{"missense_variant", "start_lost"}`, via `is_missense_or_start_loss`) --
the only two consequence types REVEL/AlphaMissense/MutPred2 actually score --
applied once, upfront, before any arm is built. `--consequence missense_only`
restricts it further still, to strictly `condensed_consequence ==
"missense_variant"` (`is_missense_only`) -- excluding start-loss variants
entirely, unlike `missense`'s broader missense-or-start-loss population; this
is the same test the `*_missense_only` `--scope` values (see below) use, but
applied to the whole analysis rather than just the two control scopes.
`--consequence all` (the default) applies no restriction.

This is different from, and broader than, the existing
`clinvar_control_missense_only`/`clingen_control_missense_only`
`--scope` values (see below): those restrict *only* their own control
population to strictly missense (excluding start-loss), leaving
`vus`/`gnomad`/`unobserved`/`clinvar_control`/`clingen_control` unrestricted.
`--consequence missense`/`--consequence missense_only` instead restrict every
scope, every direction, every chart, and every `--document` section to the
same candidate universe -- e.g. to see how the "added value of combining
evidence" chart (`--plot`) looks for missense variants only, rather than
diluted by consequence types no predictor ever scores.

Only under `--consequence missense_only` do the `*_missense_only` scopes
become redundant with their unrestricted `clinvar_control`/`clingen_control`
counterparts (the whole population is already restricted to strictly
`is_missense_only`, the same test those two scopes apply on their own) --
expected, not a bug; there's no reason to combine the two. Under
`--consequence missense`, they stay meaningfully different: the whole
population still includes start-loss variants (`is_missense_or_start_loss`
is broader than `is_missense_only`), so `clinvar_control_missense_only`/
`clingen_control_missense_only` still further exclude those from
`clinvar_control`/`clingen_control`.

## Method

For each predictor (REVEL, AlphaMissense, MutPred2 -- restrict with
repeatable `--predictor`), three "arms":

| Arm | Points column | Dedup key |
|---|---|---|
| Functional-only | `Functional_points` (ExCALIBR or OddsPath, by gene) | `abs(Functional_points)` |
| Predictor-only | `Points_<predictor>_GeneSpecific_GenomeWide` | `abs(predictor points)` |
| Combined | `Functional_points + predictor points` | `abs(combined points)` |

**Each arm redoes the one-row-per-DNA-variant deduplication using its own
points column**, rather than reusing the combined-evidence winner for the
single-source arms. This is the crux of the ablation and worth being
explicit about: `build_variant_reclassification_dataset.py`'s own dedup
picks the candidate row with the greatest `abs(Combined_points)` for a given
variant. If a functional-only ablation reused that same winning row and just
looked at its `Functional_points`, it would understate what functional
evidence alone can show whenever a *different* candidate row has stronger
functional evidence but weaker predictor evidence (and therefore lost the
combined-points contest). Redoing dedup per arm avoids that -- each arm
reports what its own evidence source, at its strongest available reading,
actually achieves.

Since `Functional_points` doesn't depend on which predictor is active, the
functional-only arm is computed once and shared across all three predictor
sections (same numbers each time).

For REVEL specifically, the combined arm's points column is the *same*
`Combined_points` column `build_variant_reclassification_dataset.py`
produces -- not a recomputed value -- so the two scripts can never drift
apart on that one arm.

## Classification and reporting

Each arm is classified pathogenic-or-benign vs. uncertain using this
pipeline's existing fixed cutoffs (`points_are_pathogenic_or_benign` in
`src/mave_dataset_stats.py`: total points >= 6 or <= -1), and summarized the
same way as that script's reclassification-export section:

- Distinct DNA variants classified (always the same across arms for a given
  predictor -- every arm dedups the same candidate universe, just picking
  different winners)
- How many are pathogenic-or-benign
- How many ClinVar VUS are resolved (reclassified pathogenic-or-benign)
- How many unobserved (SNV, no ClinVar call, no gnomAD MAF) variants are
  resolved

Then, per predictor, a three-way contingency table over each variant's
resolved-or-not flag (see `--direction` below) in the functional-only,
predictor-only, and combined arms -- reported twice: once over every
variant, and once restricted to a configurable scope of variant categories
(`--scope`, repeatable; default `vus`+`unobserved`):

| `--scope` value | Category |
|---|---|
| `vus` | ClinVar VUS (`clinvar_sig_2025` is "Uncertain significance") |
| `gnomad` | gnomAD population variant (`gnomad_MAF` non-null, regardless of ClinVar status) |
| `unobserved` | No ClinVar call, not observed in gnomAD, restricted to SNVs |
| `clinvar_control` | ClinVar control variant (see below) |
| `clinvar_control_missense_only` | `clinvar_control`, restricted to strictly missense variants, excluding start-loss (see below) |
| `clingen_control` | ClinGen Evidence Repository control variant (see below) |
| `clingen_control_missense_only` | `clingen_control`, restricted to strictly missense variants, excluding start-loss (see below) |

`vus`/`gnomad`/`unobserved` are the previously-uncertain populations this
pipeline's headline claim is about (default scope, "how many previously-
uncertain variants get resolved"). `clinvar_control`/`clingen_control` (and
their `_missense_only` restrictions) are a different kind of population
entirely -- already-known-truth calibration variants, useful for checking
each arm's agreement with a known answer rather than counting newly-resolved
uncertain variants -- so they're opt-in only, not part of the default scope.
See `--plot-control-concordance` below for a chart built specifically around
that check.

A variant is in scope if it falls into *any* of the requested categories
(union, not intersection). `unobserved` is mutually exclusive with the other
two (it requires *no* ClinVar call and *no* gnomAD MAF, the complement of
both). **`vus` and `gnomad` are not mutually exclusive**, though: `gnomad`
(`is_gnomad`) is a bare `gnomad_MAF.notna()` test with no ClinVar-null
requirement, reproducing `Variant_Classification_analysis.ipynb` cell 104
exactly, so a variant with both a ClinVar "Uncertain significance" call and a
gnomAD MAF satisfies both `vus` and `gnomad` -- and legitimately appears in
both `VUS_REVEL` and `gnomAD_REVEL` in the real Supplementary Data 5 (10,666
such rows on the real checkpoint). Selecting more than one of `vus`/`gnomad`/
`unobserved` still just widens the scope (union semantics dedup any overlap
automatically), it just doesn't mean the three partition the population into
disjoint thirds the way `clinvar_control`/`clingen_control` do relative to
each other.

**`clinvar_control`** reproduces `Variant_Classification_analysis.ipynb`
cells 78/80/85/88's row-level `controls` membership test (`is_clinvar_
control`): a mixed-year, conflict-resolved ClinVar Pathogenic/Likely
Pathogenic/Benign/Likely Benign call (`clnsig_group_18_25`, already on the
checkpoint), no ClinVar conflict (`clinvar_conflict_flag_18_25`), and 1+
star review status (`mixed_year_star_series`, the same BRCA1/PTEN/MSH2/TP53
2018-vs-2025 gene switch `mave_dataset_stats.mixed_year_clinvar_series` uses
for `clinvar_sig_*`, applied to `clinvar_star_*` instead). This reproduces
the notebook's row-level membership test only -- not its additional
aa-group-level `has_clinvar_star_conflict` exclusion or its controls-specific
dedup; like every other scope here, membership is checked per candidate row
ahead of this script's own independent per-arm dedup. Because it keys off
the *mixed-year* `clnsig_group_18_25` rather than the plain `clinvar_sig_2025`
the other three scopes use, a variant can in principle satisfy both
`clinvar_control` and `vus` at once for the four mixed-year genes -- harmless
under the union semantics above, just worth knowing if combining the two.

**`clingen_control`** reproduces cell 107's `clingen` filter (`is_clingen_
control`): a non-null, non-VUS `Updated_Classification_ClinGen_repo`
(already on the checkpoint) -- no further filtering, matching the notebook
exactly.

**`clinvar_control_missense_only`/`clingen_control_missense_only`**
are the same two control populations, additionally restricted to strictly
missense variants (`condensed_consequence == "missense_variant"`, via
`is_missense_only` -- excluding start-loss, unlike `--consequence
missense`'s broader `is_missense_or_start_loss`) -- one of the two
consequence types REVEL/AlphaMissense/MutPred2 actually score (predictor
scoring is missense-SNV-only by construction upstream, in
`vendor/variant-annotation`; start-loss is the other scored type, but these
two scopes deliberately exclude it). Every other consequence type is
guaranteed to have no predictor score, so a predictor-only or combined arm
over the full `clinvar_control`/`clingen_control` population is diluted by
variants where those arms can never contribute anything; these two
restricted scopes isolate the population where the predictor-only and
combined arms are actually a meaningful comparison against a known answer.

The all-`False` row (unresolved by anything) is dropped since it isn't
informative here. Reading the restricted table (`Resolved` = crosses the
`--direction` threshold, see below):

- `Functional=True, Predictor=False, Combined=True` -- functional evidence
  was already enough on its own.
- `Functional=False, Predictor=True, Combined=True` -- predictor evidence
  was already enough on its own.
- `Functional=False, Predictor=False, Combined=True` -- **true synergy**:
  neither source alone reaches the threshold, but their sum does.
- `Functional=True, Predictor=True, Combined=False` -- **evidence conflict**:
  both sources individually cross the threshold, but in opposite directions
  (under `--direction any`, the default), and their sum falls back into the
  uncertain range. Small in practice, but real -- confirms why redoing dedup
  per arm (rather than reusing one representative row) matters: this row
  only shows up because the two arms' own winning candidates genuinely
  disagree.

### `--direction`: restricting to one resolution direction

By default (`--direction any`), "resolved" means crossing this pipeline's
combined threshold in *either* direction (`points_are_pathogenic_or_benign`:
total points >=6 or <=-1). `--direction pathogenic` restricts every flag in
the contingency table/chart to just the Pathogenic/Likely Pathogenic test
(points >= `LIKELY_PATHOGENIC_POINTS_THRESHOLD`); `--direction benign` to
just Benign/Likely Benign (points <= `LIKELY_BENIGN_POINTS_THRESHOLD`).

This reshapes the "evidence conflict" row/bucket without any special-casing
-- it falls out of the same `Functional=True, Predictor=True, Combined=False`
-style logic (`~combined & (functional | predictor)`) applied to the
direction-restricted flags. Under `--direction pathogenic`, that row becomes
*variants a single arm alone already classified Pathogenic/Likely Pathogenic,
before combining evidence pulled the combined score back out of that range*
-- whether that's because the combined score dropped into Uncertain, or even
crossed over toward Benign-leaning. The other source doesn't need to have
independently reached the *opposite* direction's own threshold for this to
happen; it only needs enough opposing magnitude to drag the sum back down.

The per-arm summary blocks (`Distinct DNA variants classified`, `Pathogenic
or benign`, `ClinVar VUS resolved`, `Unobserved variants resolved`) are
**not** affected by `--direction` -- they always report combined
pathogenic-or-benign totals, regardless of which direction the contingency
table/chart below them is restricted to.

Variant identity is matched across arms by genomic coordinates (`Gene`,
`Chrom`, `hg38_start`, `ref_allele`, `alt_allele`, compared as strings --
see `src/lib/dedup.py`'s note on mixed `int`/`str` `Chrom` values), which is
safe because every arm dedups the identical candidate row set and ClinVar/
gnomAD annotations are properties of the variant's coordinates, not of
whichever assay row a given arm's dedup happened to pick.

#### Reconciling "lost to combining" counts across directions

**The `pathogenic` and `benign` conflict counts don't sum to the `any`
count, and that's expected -- "lost to combining" isn't a fixed property of
a variant.** It's recomputed per `--direction`, applying that direction's
own threshold to the functional, predictor, *and* combined points
independently (`categorize_synergy`'s `~c & (f | p)`, where `f`/`p`/`c` are
each `resolved_mask(..., direction)`). The same variant's raw
`Functional_points`/predictor points/`Combined_points` never change across
directions -- only which of them count as "resolved" does -- so the same
variant can be flagged conflict under one direction and cleanly resolved
under another.

Worked example, confirmed on the real integrated dataset (REVEL,
`--scope vus`): `--direction pathogenic` reports 18 conflicts, `--direction
benign` reports 1,363, `--direction any` reports 1,352 -- nowhere close to
`18 + 1363 = 1381`. Reconciled exactly:

- The 18 `pathogenic`-conflict variants are a strict subset of both the
  `benign`-conflict variants and the `any`-conflict variants. These are
  genuine "tug-of-war" cases, e.g. `Functional_points=6`,
  `REVEL_points=-1`, `Combined_points=5` (CHEK2 chr22:28694063): functional
  weakly pathogenic, predictor weakly benign, nearly cancelling to leave
  the combined score stuck in the Uncertain range (0-5) -- genuinely
  unresolved regardless of which direction you're asking about, so all
  three views agree it's a conflict.
- `any`'s 1,352 = those 18 dual-triggering variants + 1,334 more visible
  only through the benign lens (neither source alone reached the
  *pathogenic* threshold, so `pathogenic`'s report never sees them).
- The remaining 11 of the 1,363 `benign`-conflict variants are **not**
  real conflicts at all: e.g. `Functional_points=8`, `REVEL_points=-1`,
  `Combined_points=7` (BRCA1 chr17:43063914) -- functional strongly
  pathogenic, and combining produces a clean, confident Likely Pathogenic
  call (`Combined_points=7`). The *only* reason this counts as a "benign
  conflict" is that predictor's weak opposing score (`-1`) happened to just
  barely cross its own benign threshold in isolation, before being
  completely dominated by functional evidence once combined. Under `any`
  (or `pathogenic`), this variant is correctly not a conflict at all --
  it's `FUNCTIONAL_ALONE_LABEL`.

So `1334 + 18 = 1352` (`any`) and `1334 + 18 + 11 = 1363` (`benign`) --
consistent. **`--direction any`'s conflict count is the trustworthy,
direction-agnostic one.** Each single-direction chart's conflict bar
answers a narrower question ("did combining undo *this specific
direction's* call"), which can flag a handful of variants the `any` view
correctly recognizes as clean resolutions in the *other* direction --
summing or comparing conflict counts across separate `--direction` runs
isn't meaningful without accounting for this.

## Usage

```bash
poetry run python -m src.ablation_variant_reclassification
poetry run python -m src.ablation_variant_reclassification --predictor REVEL
poetry run python -m src.ablation_variant_reclassification --scope vus --scope gnomad --scope unobserved
poetry run python -m src.ablation_variant_reclassification --scope clinvar_control
poetry run python -m src.ablation_variant_reclassification --consequence missense --document ablation_missense.png
poetry run python -m src.ablation_variant_reclassification --direction pathogenic --plot ablation_plp.png
poetry run python -m src.ablation_variant_reclassification --direction benign --plot ablation_blb.png
poetry run python -m src.ablation_variant_reclassification --plot-comparison ablation_plp_vs_blb.png
poetry run python -m src.ablation_variant_reclassification --direction benign --plot-redundant-gain ablation_gain_blb.png
poetry run python -m src.ablation_variant_reclassification --show-class-upgrade --plot ablation_class_upgrade.png
poetry run python -m src.ablation_variant_reclassification --plot-control-concordance ablation_concordance.png
poetry run python -m src.ablation_variant_reclassification --document ablation_document.png
poetry run python -m src.ablation_variant_reclassification --document ablation_document_gain.png --document-chart gain
poetry run python -m src.ablation_variant_reclassification --output data/output/reclassification/ablation_variant_reclassification.txt
poetry run python -m src.ablation_variant_reclassification --plot data/output/reclassification/ablation_variant_reclassification.png
poetry run python -m src.ablation_variant_reclassification --consequence missense_only --output ablation_missense_only.txt
poetry run python -m src.ablation_variant_reclassification --document-split-dir data/output/figures/extended_data_figure_10
```

Output is plain text to stdout, and optionally also written to `--output`.
No files are written under `data/output/` by default.

## Chart (`--plot`)

`--plot <path>` (format inferred from the extension, e.g. `.png`/`.svg`/`.pdf`)
renders one stacked bar per predictor -- restricted to the same `--scope`
population and `--direction` as the text report's restricted contingency
table -- broken into four resolution categories (in a fixed color order,
following this project's data-viz conventions):

- **Functional alone sufficient** / **Predictive alone sufficient** -- one
  source alone already crosses the threshold (in `--direction`, default
  either direction).
- **Either alone sufficient (redundant)** -- each source independently crosses
  it; combining just confirms.
- **Synergy (neither alone sufficient)** -- the actual added value of
  combining: resolved only because the two sources' points were summed.

A separate bar below the zero line, **Lost to combining (evidence
conflict)**, shows variants where a single source alone crossed the
threshold but combining evidence reversed that -- these are *not* part of
the "resolved" total above the line, since the pipeline's actual combined
method does not classify them that way. Under `--direction pathogenic` (or
`benign`), this bar is specifically variants a single arm alone already
classified Pathogenic/Likely Pathogenic (or Benign/Likely Benign) *before*
combining evidence pulled the combined score back out of that range -- see
`--direction` above. Each predictor's total resolved count/percent (the sum
of the four stacked categories) is annotated above its bar.
`src.ablation_variant_reclassification.categorize_synergy` is the single
source of truth mapping the text report's contingency rows to these five
categories, so the chart and text report can never disagree with each other.

Generate separate P/LP-only and B/LB-only charts with two runs
(`--direction pathogenic --plot ...` and `--direction benign --plot ...`,
each to a different output path) -- a single `--plot` chart only ever shows
one `--direction` at a time.

### Light/dark/hatched evidence-class split (`--show-class-upgrade`)

`--show-class-upgrade` (a flag, works with both `--plot` and
`--plot-comparison`) splits the "Functional alone sufficient", "Predictive
alone sufficient", and "Either alone sufficient" bars into three sub-shades
each, without changing their total height:

- **Darker shade, labeled "(upgraded)" in the legend** (the category's
  normal base color -- so the chart looks identical to the default when the
  flag is off) -- combining pushed the variant into the *next, stronger*
  class (Likely Pathogenic -> Pathogenic, or Likely Benign -> Benign).
- **Lighter shade, labeled "(unchanged)" in the legend** -- the variant
  stayed in the same ACMG evidence-strength class (e.g. Likely Pathogenic,
  or Likely Benign) after combining as it was from the single sufficient
  source alone.
- **Lighter shade, hatched, no separate legend entry** (same color as
  "(unchanged)", explained by the caption above the chart instead) --
  combining *eroded* the variant down a class (e.g. Pathogenic -> Likely
  Pathogenic) without dropping it out of `--direction`-resolved range
  entirely -- an opposing-but-insufficient-to-declassify second source. This
  is why it's in this band at all rather than `CONFLICT_LABEL`/"Lost to
  combining": that label only fires when the *combined* score stops being
  `--direction`-resolved altogether (see `categorize_synergy`'s `~c & (f |
  p)` condition); an erosion that stays resolved, just in a weaker class,
  was previously invisible, silently folded into the plain lighter shade.

For "Either alone sufficient", the comparison baseline is whichever of the
two individual sources established the higher (stronger) class on its own --
"whichever was higher," not just functional's or just predictor's.
`class_tier` implements the ACMG cutoffs (`>=10`/`<=-7` = the stronger class,
`6-9`/`-6..-1` = the weaker one, mirroring `Variant_Classification_
analysis.ipynb` cell 38 and `docs/variant_classification.md`'s evidence-
scoring table); `class_upgrade_status` compares `abs(combined_points)`'s tier
against that baseline tier per variant and returns `CLASS_UPGRADED`/
`CLASS_UNCHANGED`/`CLASS_DOWNGRADED`. This is a genuine net-change
comparison, not just "did the other source have nonzero points."

**The legend explicitly shows both the light and dark shade of each
class-upgrade-eligible category** -- e.g. "Functional alone sufficient
(upgraded)" (dark) and "Functional alone sufficient (unchanged)" (light) as
two separate entries -- rather than one base-colored swatch per category. A
caption above the chart still explains the hatched sub-shade, which reuses
the "(unchanged)" light color rather than adding a third legend entry.

**Real findings worth knowing before reading too much into an all-light(-and-
unhatched) band**: whether upgrades can happen at all in a given band
depends on the asymmetry between this pipeline's pathogenic (`>=6`) and
benign (`<=-1`) resolution thresholds -- but only as a rough guide, not a
guarantee, for the "Functional/Predictive alone sufficient" bands
specifically. The naive argument is that the *other*, not-independently-
sufficient source's own points are bounded on one side by definition --
under `--direction benign`, "not independently sufficient" means that
source's points are `>= 0` -- so summing it with the already-sufficient
source's points could only weaken or hold the class, never strengthen it.
That argument implicitly assumes "combined" means "this same row's
functional points plus this same row's predictor points," which is exactly
what each arm's *independent* per-DNA-variant dedup (see `build_arm` above)
does not guarantee: the combined arm picks whichever underlying assay/
annotation row maximizes `abs(Combined_points)` for that variant, which can
be a *different* row than the one the functional-only (or predictor-only)
arm independently picked as its own best. When a variant has more than one
functional measurement, the row that wins "Functional alone sufficient"
(largest `abs(Functional_points)` on its own) need not be the row that
combines best with a fixed predictor score -- so a small number of
variants can still show a real upgrade in these bands. Confirmed on the
real integrated dataset (default `vus`+`unobserved` scope): one concrete
example is CHEK2 22:28689191 G>C -- the functional-only arm's winning row
has `Functional_points = +5` (pathogenic-leaning, not benign-resolving on
its own), REVEL's own points are `-3` (Likely Benign, tier 1), and the
combined arm's own winning row (a *different* functional measurement for
the same variant) combines to `Combined_points = -7` (Benign, tier 2) -- an
upgrade, even though naive addition of the two "winning" rows shown above
(`+5` and `-3`) would suggest otherwise.

Counts on the current checkpoint, `--direction benign --show-class-upgrade`
(vus+unobserved scope, summed): "Either alone sufficient" shows the bulk of
upgrades, as expected (REVEL 8,548 of 26,173; AlphaMissense 7,085 of 25,518;
MutPred2 7,804 of 26,516 -- roughly a third of each). "Functional alone
sufficient" shows zero upgrades for all three predictors in this dataset
(REVEL 0 of 37,473; AlphaMissense 0 of 39,186; MutPred2 0 of 37,316) --
consistent with the naive bound, though not guaranteed to hold in general by
the mechanism above. "Predictive alone sufficient" does show a handful of
real upgrades from exactly that per-arm-dedup mechanism -- small, but
nonzero: REVEL 2 of 4,888, AlphaMissense 6 of 6,377, MutPred2 6 of 5,812.
Under `--direction pathogenic`, "Predictive alone sufficient" and "Either
alone sufficient" are entirely empty for all three predictors -- no variant
is ever categorized into either band, since no predictor alone reaches the
Pathogenic/Likely Pathogenic threshold in this dataset (see the note on this
asymmetry elsewhere in this doc). All of `--direction pathogenic`'s upgrades
(and there are many -- REVEL 1,718 of 4,554 functional-alone variants
combining `vus`+`unobserved`) show up in "Functional alone sufficient"
itself. All three predictors also show real (if small -- well under 1% of
each band) downgraded-but-still-
resolved counts in both "Functional alone sufficient" and "Predictive alone
sufficient" under `--direction benign` -- e.g. REVEL: 178 of 37,473
functional-alone variants, 22 of 4,888 predictor-alone variants --
confirming the erosion case isn't just theoretical, and was invisible
before this flag existed.

## Side-by-side comparison chart (`--plot-comparison`)

`--plot-comparison <path>` renders the Pathogenic/Likely Pathogenic and
Benign/Likely Benign charts as two panels in one image instead of two
separate files -- independent of `--direction`, which only affects `--plot`
and the text report; `--plot-comparison` always shows both directions
regardless of what `--direction` is set to. Each panel keeps its own
independent y-axis scale (Benign/Likely Benign calls are typically far more
common than Pathogenic/Likely Pathogenic ones, often by an order of
magnitude, so a shared scale would crush the smaller panel flat), and the
two panels share one legend and one "Of N ... variants in scope" line (the
scope population is identical between panels -- only the direction differs)
rather than repeating it per panel.

Internally (`build_direction_comparison_chart_data`), each predictor's arms
are deduplicated once and reused for both directions, rather than rebuilding
them per direction -- arm dedup keys on `abs(points)`, which doesn't depend
on direction, only the resolved-flag computation does.

`--show-class-upgrade` (see above) works here too -- pass both flags together
to get the light/dark/hatched split on both panels at once, one shared
caption.

## Points-gained histogram (`--plot-redundant-gain`)

`--plot-redundant-gain <path>` drills into `--plot`'s "Each alone sufficient
(redundant)" band -- the variants where *both* the functional-only and
predictor-only arms already independently resolve the variant (in
`--direction`), so combining them wasn't strictly *necessary*. It doesn't
answer whether combining was necessary (the main chart already does); it
answers *how much extra evidence combining added anyway*, for that specific
population.

For each predictor, one panel; x-axis is the exact integer number of points
gained, y-axis is variant count, two dodged bars per value:

- **Gained via functional evidence** -- `abs(combined_points) -
  abs(predictor_points)`: how much the combined score exceeds what predictor
  evidence alone already established.
- **Gained via predictor evidence** -- `abs(combined_points) -
  abs(functional_points)`: the same, the other way around.

This is a net-change calculation, not just a restatement of the other
source's own points -- the two agree exactly when both sources are
same-signed (guaranteed under `--direction pathogenic`/`benign`, since both
are constrained to one sign by construction), but under `--direction any` a
variant can land in this band with the two sources resolved in *opposite*
directions (each independently crosses its own threshold, while the summed
score still happens to clear a threshold too); there, this calculation
correctly reports a small or negative "gain" -- a partial erosion -- rather
than overstating it. `redundant_gain_counts` is the function to read for the
exact formula.

On the real integrated dataset, this reveals a structural difference between
evidence sources: REVEL/AlphaMissense/MutPred2 points only take a handful of
discrete ACMG evidence-strength values (1/2/3/4/8 -- Supporting/Moderate/
Moderate+/Strong/Very Strong, skipping 5-7 entirely, since those aren't valid
PP3/BP4 codes), while ExCALIBR's functional points come from a continuous
score-interval calibration and can land on any integer in range -- so
"gained via predictor" clusters tightly at 1/2/3/4/8 while "gained via
functional" spreads smoothly across the full range.

Restricted to the same `--scope`/`--direction` as `--plot`. Since the "each
alone sufficient" band can be empty for some predictor/direction/scope
combinations (e.g. no predictor in this dataset independently reaches the
Pathogenic/Likely Pathogenic threshold on its own, so `--direction
pathogenic` produces an all-zero, uninformative chart for every predictor;
`--direction benign` or `any` are the more informative choices in practice),
check the "Each alone sufficient (redundant)" count on `--plot`'s chart
first to confirm there's something to drill into.

## Control-concordance chart (`--plot-control-concordance`)

Every chart above answers "how many previously-uncertain variants get
resolved." `--plot-control-concordance <path>` answers a different question:
for variants whose true classification is already known (`clinvar_control`/
`clinvar_control_missense_only`/`clingen_control`/`clingen_control_
missense_only`), does each arm's resolved call actually *agree* with
that known answer? The ablation contingency table has no way to express
this -- it only tracks whether an arm resolved a variant, not whether the
resolution was correct.

One panel per predictor, each a three-bar cluster -- functional-only,
predictor-only, combined -- stacked into:

- **Concordant** (blue) -- the arm resolved the variant, and in the same
  direction as the known ClinVar/ClinGen truth.
- **Discordant** (red) -- the arm resolved the variant, but in the *opposite*
  direction from the known truth.
- **Unresolved** (gray) -- the arm didn't cross either resolution threshold
  for this variant at all.

`concordance_status` computes this per variant/arm; `concordance_counts`/
`predictor_control_concordance_counts` roll it up per arm/predictor.

Restricted to whichever of `clinvar_control`/`clinvar_control_missense_only`/
`clingen_control`/`clingen_control_missense_only` are passed
via `--scope` (repeatable); if none is given, all four are used by default --
independent of `--scope`'s own default (`vus`+`unobserved`), since the other
three carry no known truth to check against (`control_truth_direction_series`
returns `pandas.NA` for every row outside the four control scopes). Where a
variant satisfies more than one control scope with disagreeing truth values
(rare -- the two sources are independent in practice), `SCOPE_ORDER`'s
earlier scope wins (ClinVar before ClinGen), a deterministic tiebreak with no
expected real-world impact.

**Color choice is deliberate**: this chart uses a blue/red status pair
(`CONCORDANT_COLOR`/`DISCORDANT_COLOR`) rather than the categorical palette
`--show-class-upgrade`/the ablation/comparison charts use elsewhere in this
script. Concordance is a *correctness* dimension (did the call match a known
answer), not a category identity, so reusing the categorical colors here
would collide with their existing meaning in this script's other charts.
Blue/red rather than the dataviz skill's default green/red "good"/"critical"
status colors, too: green paired with red is a classic colorblind confusion,
so this reuses the project's own recurring blue-vs-red convention for
opposite states instead (e.g. `figure_4/plot_utils.py`'s
`BENIGN_THRESHOLD_COLOR`/`PATHOGENIC_THRESHOLD_COLOR`, `Figure5_6.Rmd`'s
"Benign"/"Pathogenic"). The same reasoning applies to `BOTH_ALONE_LABEL`'s
categorical color (see its definition): that teal-green shared every
ablation/comparison chart with this same red, so it's now a blue-violet
instead (a warm brown was tried first and rejected -- too close in hue to
`PREDICTOR_ALONE_LABEL`'s orange).

On the real integrated dataset, the functional-only and combined arms are
both highly concordant (few percent discordant at most) and mostly resolved,
while each predictor-only arm leaves the large majority of control variants
*unresolved* -- predictor evidence alone rarely reaches this pipeline's
resolution threshold on its own, so functional evidence, not predictor
evidence, is what does most of the resolving work among controls. This
mirrors -- and helps explain -- the "Functional alone sufficient" band
dominating `--plot`'s chart.

## `--special-control-dedup`: Supplementary Data 5's own control dedup

By default, every scope in this script -- including `clinvar_control`/
`clingen_control` -- is deduplicated the same generic way: one row per DNA
variant, keyed on `abs(points)` (`dedup_by_max_abs_points`, via `build_arm`).
Supplementary Data 5 itself doesn't dedup its `controls_<predictor>_
GeneSpecific`/`ClinGen_Repo_<predictor>_GeneSpecific` sheets this way --  it
uses a more elaborate, amino-acid-level-aware pipeline
(`Variant_Classification_analysis.ipynb` cells 78-92/107-114,
`CONTROLS_CLINGEN_DEDUP_STRATEGY = "nt_then_abs_max"`) that this script's
generic dedup has no equivalent to: several distinct NT variants can share
one amino-acid change and therefore one shared amino-acid-level functional/
predictor measurement, which a coordinate-only dedup would otherwise let
each of those NT variants independently claim as if it were unique
supporting evidence.

`--special-control-dedup` (`special_control_dedup`) switches `clinvar_
control`/`clinvar_control_missense_only`/`clingen_control`/`clingen_
control_missense_only` over to that same pipeline, reusing the exact
shared `src.lib.dedup` mechanics (`aa_dedup_or_mark`, `catch_mis_2`,
`controls_aa_sort_key`/`clingen_aa_sort_key`) the notebook itself calls, so
the two can never drift apart independently. `vus`/`gnomad`/`unobserved`,
the all-variants contingency table, and the per-arm summary stats blocks are
completely unaffected either way.

Two branches, split by `nucleotide_or_aa`: an **nt branch** (direct,
NT-resolution measurements, collapsed by `VariantNotes` tag order, ClinVar
additionally restricted to 1+-star review status) and an **aa branch**
(evidence that only resolves to a shared amino-acid change -- ClinVar
additionally drops any amino-acid group whose review-status membership
conflicts across its rows). Restricted to rows already flagged as both the
shared functional maximum and that predictor's own per-variant maximum, one
NT variant is chosen to represent each amino-acid group (ClinVar prefers the
highest-confidence, most-recently-reviewed record; ClinGen prefers the
non-retracted, most-recently-approved/published one) -- this is the actual
double-counting guard: each amino-acid group contributes exactly one
representative row, not one per constituent NT variant. The two branches are
concatenated, that predictor's own training variants excluded (no exclusion
for AlphaMissense, matching the notebook), and `catch_mis_2` resolves the
remaining case where the nt branch's own survivor and the aa branch's chosen
representative land on the exact same genomic coordinate.

**Confirmed against the real, distributed `Supplementary_Data_5.xlsx`**: on
the current integrated dataset, `special_control_dedup(df, "REVEL",
"clinvar")` produces exactly 12,433 rows, matching `controls_REVEL_
GeneSpecific`'s own row count exactly. AlphaMissense/MutPred2 match within
~0.3% (a handful of amino-acid groups have fully tied sort keys -- identical
review rank, review date, and functional points -- whose winning "primary"
row then depends on incidental DataFrame row order, which this script
doesn't guarantee matches the notebook's own historical row order for such
ties). A separate `Supplementary_Data_5.with_secondary_variants.xlsx` export
keeps every candidate (tagged `Variant_Role` "primary"/"secondary") rather
than filtering to primary rows only -- `special_control_dedup` always
filters to primary, matching the canonical, distributed file.

Since a variant's genomic coordinate means the same thing regardless of
which population resolved it, this option composes with everything else in
this script (`--predictor`, `--direction`, `--consequence`, `--document`,
`--plot-control-concordance`) -- including a `--scope` that mixes control
and non-control categories, where only the control portion switches dedup
strategy.

Usage: `--special-control-dedup --scope clinvar_control --plot
ablation_special_dedup.png`, or `--special-control-dedup
--plot-control-concordance ablation_concordance_special.png`.

## Multi-section document (`--document`)

`--document <path>` (format inferred from the extension, e.g.
`.png`/`.svg`/`.pdf`) renders a single combined document instead of separate
chart files: a header (document title, spanning the full width) followed by
**one section per variant category, all seven by default** -- `vus`, `gnomad`,
`unobserved`, `clinvar_control`,
`clinvar_control_missense_only`, `clingen_control`,
`clingen_control_missense_only` -- regardless of `--scope`, which only
affects `--plot`/`--plot-comparison`/`--plot-redundant-gain`/the text
report. `--scope` unions whichever categories you pass into one population;
`--document` instead shows each one side by side, individually, so you can
compare them directly -- pass `--document-scope` (repeatable) to restrict
which of the seven get a section at all (see below). Each section shows a
header line ("_scope label_ -- _N_ variants") followed by whichever
`--document-chart` types were requested, stacked top to bottom in a fixed
order:

1. **`ablation`** -- `--plot`'s chart, restricted to `--direction`.
2. **`comparison`** -- `--plot-comparison`'s side-by-side P/LP-vs-B/LB chart
   (always both directions, independent of `--direction`, same as
   `--plot-comparison` itself).
3. **`concordance`** -- `--plot-control-concordance`'s chart (always all four
   control scopes, independent of `--direction` and `--scope`). Included only
   for the `clinvar_control`/`clinvar_control_missense_only`/
   `clingen_control`/`clingen_control_missense_only` sections -- the other
   three sections carry no known truth to check against
   (`control_truth_direction_series`), so their section simply omits this
   block entirely (`_document_chart_types_for_scope`) rather than rendering
   an always-empty one.
4. **`gain`** -- `--plot-redundant-gain`'s points-gained histogram,
   restricted to `--direction`.

`--document-chart` is repeatable; pass it once per chart type you want
(`--document-chart ablation --document-chart gain` for just those two), or
omit it for all four. `--show-class-upgrade` applies to the document's
`ablation`/`comparison` blocks the same way it applies to `--plot`/
`--plot-comparison` (it has no effect on `concordance`, which uses the
status palette, not the categorical one).

### Restricting which sections appear (`--document-scope`)

`--document-scope` (repeatable, choices are `SCOPE_ORDER`'s seven category
keys) restricts `--document`/`--document-split-dir` to just the named
sections instead of all seven; omitting it keeps the default (all seven).
Combined with `--document-chart`, this lets a single document mix chart
types by scope without any chart type ever rendering an empty block: e.g.

```bash
poetry run python -m src.ablation_variant_reclassification \
  --document-scope vus --document-scope gnomad --document-scope unobserved \
  --document-scope clinvar_control --document-scope clingen_control \
  --document-chart comparison --document-chart concordance \
  --document ablation_5comparison_2concordance.pdf
```

renders exactly five sections (`vus`, `gnomad`, `unobserved`,
`clinvar_control`, `clingen_control` -- the two `_missense_only` control
scopes and the plain `--scope` default population omitted), each with a
`comparison` block, and -- because `concordance` is only ever included for
`CONTROL_SCOPES` sections -- a `concordance` block appears under just the
two control ones. `--special-control-dedup` (see above) is the natural
companion here, since these two sections' numbers are exactly where it
matters. `document_scopes` must be a subset of whatever was passed to
`build_document_chart_data` when `--document-split-dir` is used
programmatically rather than via this CLI (naming an unbuilt scope raises a
`KeyError`); the CLI always builds and renders the same `--document-scope`
selection, so this doesn't come up from the command line.

**Each block carries its own legend** (just the handles relevant to that
block's own chart type -- `ablation`/`comparison` show the five
`SYNERGY_CATEGORY_ORDER`/`CONFLICT_LABEL` swatches, `concordance` shows its
three concordant/discordant/unresolved swatches, `gain` shows its two --
rather than one legend shared across the whole document. Each block reserves
a bottom margin sized to fit its own legend (`DOCUMENT_BLOCK_LEGEND_HEIGHT`)
via `subplots_adjust`, then anchors the legend at that block's own
subfigure bottom edge (figure-fraction coordinates, so it never needs to
extend past the block's own bounds) -- a legend positioned outside a nested
subfigure's own `[0, 1]` range doesn't have this guarantee:
`bbox_inches="tight"` has been confirmed to silently clip such a legend at
the subfigure's edge instead of expanding the saved figure to fit it.
`build_document_chart_data`/`save_ablation_document`/`_draw_document_block`
are the functions to read for the full mechanics.

**Runtime**: a document reuses each predictor's arms (dedup) across every
scope, chart type, and direction it needs, the same one-time-dedup
principle as `--plot-comparison`/`--show-class-upgrade` -- but rendering
seven sections' worth of charts is still substantially more work than a
single chart. All four chart types, all three predictors, produces an image
roughly 27,100px tall (proportionally up from five sections' ~19,400px) and
takes several minutes end to end on the real integrated dataset (dominated
by rendering, not the dedup pass) -- expected, not a hang. Pass
`--predictor`/`--document-chart` to narrow it down if you don't need
everything.

Restricting `--document-chart` to a single type still produces this
function's seven-section layout (just with fewer blocks per section) --
functionally a "seven populations side by side" view of whichever one chart
type you asked for, e.g. `--document-chart gain` alone shows the
points-gained histogram for `vus`/`gnomad`/`unobserved`/`clinvar_control`/
`clinvar_control_missense_only`/`clingen_control`/`clingen_control_
missense_only` stacked in one file for direct comparison, something
none of `--plot`/`--plot-comparison`/`--plot-redundant-gain`/`--plot-
control-concordance` can do on their own since each of those is restricted
to a single (possibly unioned) `--scope`.

## Splitting the document into individual chart files (`--document-split-dir`)

`--document` renders every section/chart-type block into one combined
multi-section image -- convenient for browsing the whole analysis at once,
but not how a print figure gets assembled: a manuscript figure is built from
individual panels, each its own file. `--document-split-dir <dir>`
(`save_document_charts_as_files`) renders the exact same per-section chart
data as `--document`, but writes each `(chart type, scope)` block to its own
file under `<dir>` instead -- one file per block, named
`<chart type>_<scope>.<ext>`, e.g. `ablation_vus.pdf`,
`concordance_clinvar_control.pdf`, `gain_unobserved.pdf`. This is what
[`docs/figures.md`](figures.md)'s Extended Data Figure 10 section uses.

Each file is rendered by that chart type's own standalone `save_synergy_
chart`/`save_synergy_chart_comparison`/`save_control_concordance_chart`/
`save_redundant_gain_chart` function -- the same functions `--plot`/
`--plot-comparison`/`--plot-control-concordance`/`--plot-redundant-gain`
call directly -- so **every file already carries its own legend**; there's
no shared or document-level legend to reconstruct.

A `concordance` file is skipped (not written) for the `vus`/`gnomad`/
`unobserved` scopes: those sections carry no known ClinVar/ClinGen truth to
check concordance against (see `control_truth_direction_series`), so
`--document`'s own equivalent block there is an all-"unresolved", zero-value
chart nobody would use as a figure panel -- the split output omits the file
entirely rather than writing an empty one.

`--document-chart` (repeatable; default all four) controls which chart types
are split out, and `--document-scope` (repeatable; default all seven)
controls which scopes get any files at all, exactly as both do for
`--document`. `--document-split-dir` composes with `--document` -- passing
both writes the combined image *and* the split files from one
`build_document_chart_data` pass, rather than rebuilding it twice. `--document-split-format` (`pdf`/`png`/`svg`, default
`pdf`) sets the file format for every split file, since (unlike every other
`--plot*`/`--document` path) there's no single output path to infer an
extension from.

Usage:

```bash
poetry run python -m src.ablation_variant_reclassification \
  --document-split-dir data/output/figures/extended_data_figure_10
poetry run python -m src.ablation_variant_reclassification \
  --document-chart ablation --document-chart gain \
  --document-split-dir data/output/figures/extended_data_figure_10/ablation_and_gain_only
```

See [`docs/figures.md`](figures.md) for the full Extended Data Figure 10
invocation this project actually uses.

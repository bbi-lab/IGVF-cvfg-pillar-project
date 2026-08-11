# Variant classification methodology

`Variant_Classification_analysis.ipynb` combines MAVE functional-assay
evidence with computational-predictor evidence (REVEL/AlphaMissense/
MutPred2) to classify variants into ACMG/AMP-style Pathogenic/Benign
buckets, for five variant categories: ClinVar controls, ClinGen Evidence
Repository controls, VUS, gnomAD population variants, and variants
unobserved in either. This doc records the methodology behind that
process and the specific decisions/gaps found while documenting it,
organized by topic rather than by notebook cell order. See
`notebooks/analysis/README_Variant_Classification_analysis.md` for
setup/usage, and [`docs/assay_priority_questions.md`](assay_priority_questions.md)
for the `ASSAY_PRIORITY_LIST`-specific background this doc builds on.

## Overview

- **[Evidence scoring](#evidence-scoring)** -- how functional
  (`Fxn_points`) and predictor evidence become `Total_Points` and a
  `Class`.
- **[Handling conflicting evidence](#handling-conflicting-evidence)** --
  what happens when assays disagree with each other, or functional
  evidence disagrees with predictor evidence.
- **[Gene-specific special-casing](#gene-specific-special-casing)** --
  `F9`/`TP53`'s alternate calibration, `SFPQ`'s exclusion, `CHEK2`'s QC
  flag, `LDLR`'s assay-priority-by-position override, and the
  2018-vs-2025 ClinVar vintage switch.
- **[Deduplication strategy](#deduplication-strategy)** -- how the
  pipeline picks one representative record when the same variant is
  scored by more than one assay. By far the largest section, and the
  reason this doc originally existed -- it's also where the most
  implementation/intent gaps were found.
- **[Outputs](#outputs)** -- what each category/predictor combination
  produces.
- **[Open questions](#open-questions)** -- gaps found while writing this
  doc that haven't been resolved yet.

## Evidence scoring

For each predictor (REVEL / AlphaMissense / MutPred2):

```
Total_Points_GeneSpecific_<predictor> = Fxn_points + Points_<predictor>_GeneSpecific_GenomeWide
```

where the predictor points prefer a gene-specific calibration, falling
back to the genome-wide calibration when no gene-specific one exists,
falling back to `0` when neither does:

```python
Points_REVEL_GeneSpecific_GenomeWide = (
    Points_REVEL_GeneSpecific.fillna(Points_REVEL_GenomeWide).fillna(0)
)
```

`Total_Points_*` is then mapped to a `Class_*` bucket by fixed cutoffs:

| `Total_Points_*` | `Class_*` |
|---|---|
| ≥ 10 | Pathogenic |
| 6 to 9 | Likely Pathogenic |
| 0 to 5 | Uncertain |
| −6 to −1 | Likely Benign |
| ≤ −7 | Benign |

This isn't a full account of the evidence-point methodology -- see
`notebooks/analysis/README_Variant_Classification_analysis.md` and
Tavtigian et al. 2020 / Bergquist et al. for that. `Fxn_points` itself
comes from ExCALIBR score-interval calibrations for most genes, with
several gene-specific exceptions -- see
[Gene-specific special-casing](#gene-specific-special-casing) below.

A given evidence-strength point value can correspond to more than one
disjoint score interval -- currently only `DDX3X_Radford_2023`, whose `+1`
and `-1` points each cover two separate ranges. `Supplementary_Data_4`'s
`ExCALIBR_calibrations` sheet still has one `range_<point>` column per point
value, but a cell can hold a comma-separated list of `"<low> <high>"` ranges
instead of just one (`src/load_excalibr_calibrations.py`'s `format_range`);
`get_points_from_intervals` in `Variant_Classification_analysis.ipynb`
(the cell computing `ExC_points_2025`/`ExC_points_2018`) checks a score
against every comma-separated range in a cell and returns that point value
if any of them contains it.

## Handling conflicting evidence

The deduplication strategies described below resolve *which record to
keep* when multiple assays score the same variant. Separately, the
notebook has two distinct mechanisms for evidence that actively
*disagrees* -- one upstream of dedup that excludes the variant outright,
and one downstream of it that flags the result without excluding
anything. They address different questions and neither one is affected by
the choice of dedup strategy.

### Inter-assay functional conflict (excluded, not resolved)

If different assays scoring the *same variant* disagree in direction --
one reports pathogenic-direction (`Fxn_points` > 0) evidence, another
reports benign-direction (`Fxn_points` < 0) evidence -- the variant is
dropped from the dataset entirely, before any category split or dedup
strategy runs. This is not something `"v1"`/`"abs_max"`/`"nt_then_abs_max"`
resolve by preference; it happens earlier and identically regardless of
which strategy is configured.

Mechanism (`has_opposite_signs`, applied three times):

1. Within `sankey_nuc` (nt-type rows only), grouped by genomic coordinates
   (`Gene`/`Chrom`/`hg38_start`/`ref_allele`/`alt_allele`).
2. Within `sankey_aa` (aa-type rows only), grouped by
   `Gene`/`aa_pos`/`aa_ref`/`aa_alt`/transcript.
3. On the combined `sankey_f` (nt and aa together), grouped by genomic
   coordinates again -- this catches *cross-resolution* conflicts an
   nt-only or aa-only check would miss (e.g. an nt assay reports +5, an aa
   assay reports −3, for the same physical variant).

Each check ignores `NaN` and exact-zero values, and flags a group only if
it contains both a strictly-positive and a strictly-negative value. Any
row belonging to a flagged group has `VariantNotes` set to
`'conflicting_fxn_data'` (appended if it already carried a different tag).
Every row carrying that tag -- alone or combined with the splice-variant
tags -- is then dropped from `sankey_f` in one place, before `controls`,
`ClinGen_Repo`, `VUS`, `gnomAD`, and `Unobserved` are split off, so all
five categories share the same exclusion.

Practical consequence for the dedup strategies: by the time any of them
run, every surviving candidate group is guaranteed to have same-signed (or
zero) `Fxn_points` values, so a signed-value comparison and an
absolute-value comparison can only disagree on *magnitude* among
same-signed values -- never on a genuine sign conflict. That's consistent
with the documented [stage-3 sign bias](assay_priority_questions.md#known-problems-by-category)
examples, which are all same-sign (e.g. −5 vs. −2), not opposite-sign.

### Functional-vs-predictor conflict (flagged, not excluded)

Separately, each individual assay row is checked for whether its own
*functional* evidence (`Fxn_points`, or `OP_points` for the OddsPath-based
genes) disagrees in direction with its own *computational predictor*
evidence (REVEL/AlphaMissense/MutPred2 points, genome-wide or
gene-specific) -- e.g. the functional assay says pathogenic-direction
while REVEL says benign-direction for that same row. This is a **per-row,
per-predictor** check (`split_zero`), run once, early, on the full
un-deduplicated per-assay dataset (`pp_ex_OP`) -- well before the
opposite-sign exclusion above, the category split, or either dedup
parameter.

For each predictor/scope combination (e.g. gene-specific REVEL), a new
`Conflicting_<predictor>_<scope>` column is created alongside the existing
`Total_Points_<predictor>_<scope>` column:

- If either side is missing, or both sides are `NaN`, or the two sides
  agree in sign (or either is zero): `Conflicting_*` just copies the
  numeric `Total_Points_*` value (or reports `"No evidence"` if both sides
  are missing).
- If the two sides are strictly opposite in sign: `Conflicting_*` is set
  to the string `"Conflicting evidence"` instead of a number.

Unlike the inter-assay conflict above, **this never removes the variant
and never changes `Total_Points_*`/`Class_*`** -- classification is
computed from the numeric total regardless. `Conflicting_*` is a parallel,
informational column meant to flag the disagreement for manual review
without discarding the underlying evidence or classification. Whichever
row a dedup strategy selects as the representative for a variant carries
forward *that row's own*, already-computed `Conflicting_*` value -- the
dedup strategy has no say in it.

Only the gene-specific `Conflicting_REVEL`/`Conflicting_AM`/`Conflicting_MP2`
columns (renamed from the `_GeneSpecific` suffix) make it into the final
Supplementary Data 5 output; the genome-wide and OddsPath-suffixed
variants of this column are computed but not exported.

## Gene-specific special-casing

A handful of genes get bespoke handling beyond the general pipeline
described above. All of it is applied before category split, conflict
checks, or dedup, so none of it is affected by dedup strategy -- but all
of it affects what a strategy ever gets a chance to see.

### `F9` and `TP53`: alternate calibration and restricted datasets

`F9` and `TP53` are special-cased twice over, unrelated to dedup strategy
but relevant to understanding why they rarely (`TP53`) or never (`F9`)
show up needing a dedup decision at all.

**First, their `Fxn_points` don't come from ExCALIBR at all.** Every other
gene's `Fxn_points` is `ExC_points_2025` (the ExCALIBR score-interval
calibration, with the 2018-vintage override described
[below](#brca1ptenmsh2tp53-2018-vs-2025-clinvar-vintage) for a few genes).
`F9` and `TP53` instead get `Fxn_points = OP_points` -- points derived from
`OddsPath` likelihood ratios on a fixed evidence-strength scale (BS3/PS3),
computed from `OddsNormal`/`OddsAbnormal` via fixed thresholds (e.g.
`OddsAbnormal > 350` scores `8`, the PS3 maximum; `OddsNormal < 0.053`
scores `-4`, the BS3 maximum) rather than a per-gene score-interval
calibration. This is a completely different points methodology for these
two genes, not just a dataset restriction.

**Second, a single cell drops every `F9` and `TP53` dataset except one
designated survivor per gene:**

```python
#remove F9 and TP53 datasets that are not the meta analysis

pp_ex_OP = pp_ex_OP[~pp_ex_OP['Dataset'].isin([
    'TP53_Boettcher_2019', 'TP53_Fortuno_2021',
    'TP53_Giacomelli_2018_combined_score', 'TP53_Giacomelli_2018_p53WT_Nutlin3',
    'TP53_Giacomelli_2018_p53null_Nutlin3', 'TP53_Giacomelli_2018_p53null_etoposide',
    'TP53_Kato_2003_AIP1nWT', 'TP53_Kato_2003_BAXnWT', 'TP53_Kato_2003_GADD45nWT',
    'TP53_Kato_2003_MDM2nWT', 'TP53_Kato_2003_NOXAnWT', 'TP53_Kato_2003_P53R2nWT',
    'TP53_Kato_2003_WAF1nWT', 'TP53_Kato_2003_h1433snWT',
    'F9_Popp_2025_carboxy_F9_specific', 'F9_Popp_2025_carboxy_gla_motif',
    'F9_Popp_2025_heavy_chain', 'F9_Popp_2025_light_chain', 'F9_Popp_2025_strep_2',
])]
```

| Gene | Kept (the "meta-analysis" dataset) | Removed | Removed datasets |
|---|---|---|---|
| `TP53` | `TP53_Fayer_2021_meta` | 14 datasets | `TP53_Boettcher_2019`; `TP53_Giacomelli_2018_combined_score`/`_p53WT_Nutlin3`/`_p53null_Nutlin3`/`_p53null_etoposide`; `TP53_Kato_2003_AIP1nWT`/`_BAXnWT`/`_GADD45nWT`/`_MDM2nWT`/`_NOXAnWT`/`_P53R2nWT`/`_WAF1nWT`/`_h1433snWT`; `TP53_Fortuno_2021` |
| `F9` | `F9_Popp_2025_model` | 5 datasets | `F9_Popp_2025_carboxy_F9_specific`/`_carboxy_gla_motif`/`_heavy_chain`/`_light_chain`/`_strep_2` |

All 19 removed datasets are confirmed present with real data upstream of
this filter (8,313-31,193 rows each, all aa-resolution) -- this is a
deliberate exclusion of real, substantial datasets, not a filter for
already-empty or malformed ones. The pipeline README documents the
practical rule this implements: *"TP53 and F9 datasets are limited to
meta-analysis versions (individual assay datasets excluded)."*

Judging from the dataset names, most of what's removed looks like several
correlated readouts from a *single* underlying study rather than
independent assays -- the eight removed `TP53_Kato_2003_*` datasets are
different transcriptional-reporter targets (`AIP1`, `BAX`, `GADD45`,
`MDM2`, `NOXA`, `P53R2`, `WAF1`, `h1433s`) measured in what looks like the
same experiment; the four removed `TP53_Giacomelli_2018_*` datasets are
different experimental conditions from one study; the five removed
`F9_Popp_2025_*` datasets look like different structural
domains/epitopes from one DMS study. Treating each of those as an
independently-competing "assay" in the general dedup/`ASSAY_PRIORITY_LIST`
framework would effectively let one underlying study's evidence compete
against itself many times over -- consistent with this pipeline's general
avoid-double-counting theme (see
[The two categories of variant, and why they're handled differently](#the-two-categories-of-variant-and-why-theyre-handled-differently))
even though the exact rationale isn't spelled out in the notebook itself.
`TP53_Fortuno_2021` doesn't fit that pattern as neatly (it isn't obviously
a sub-readout of another kept/removed `TP53` dataset), but it's excluded
by the same rule regardless.

**Consequence for dedup:** by the time any dedup strategy runs, `TP53` has
at most one surviving dataset (`TP53_Fayer_2021_meta`) and `F9` has
exactly one (`F9_Popp_2025_model`) -- neither gene can ever present a
multi-assay dedup decision to `"v1"`/`"abs_max"`/`"nt_then_abs_max"`, or to
`ASSAY_PRIORITY_LIST`, regardless of strategy. This is why `TP53` appears
in the [known problems](#known-problems-with-v1) and
[empirical comparison](#comparing-v1-to-the-decided-approach) sections at
all only via its *other* aa dataset interactions (none, currently -- `TP53`
has exactly one surviving aa assay, so it was never actually a candidate
for `assay_priority` in the [table in the deduplication section below](#where-assay_priority_list-actually-decides-an-outcome)),
and why `F9` doesn't appear anywhere in this doc's tables at all.

### `SFPQ`: excluded entirely

`SFPQ` is dropped from the dataset in a single line, before any other
category-level filtering: `sankey_f = sankey_f[sankey_f['Gene'] != 'SFPQ']`.
Per the pipeline README, this is because of insufficient controls -- not
enough ClinVar-classified `SFPQ` variants to support the `controls`-based
calibration this pipeline otherwise relies on. `SFPQ` doesn't appear in any
of this doc's output categories at all.

### `CHEK2`: an external QC flag

Before the general conflict/flag cleanup, `sankey_f` is left-joined
against `CHEK2_Gebbia_2024.xlsx` on (`hgvs_p`, `auth_reported_score`) /
(`hgvs_pro`, `score`), and any row where that file's `Filter_CI == 1` has
its `Flag` set to `'*'` -- which then gets removed by the same
"remove flagged variants" step (cell 75) that drops everything else
flagged for any other reason, regardless of gene. A second column from
that file, `Filter_Hypercomplement`, is merged in alongside `Filter_CI`
but never used for anything.

### `LDLR`: `+VLDL` uptake assay prioritized in two LA modules, LA module 1 excluded entirely

Within LA module 2 (`aa_pos` 66-106) and LA module 6 (`aa_pos` 234-272), a
row filter runs immediately after `pp` is loaded, before any other
processing: for any amino-acid substitution in those two ranges where
`LDLR_Tabet_2025_presence_VLDL` (LDL uptake measured in the *presence* of
excess VLDL, i.e. "+VLDL") has a measurement, the corresponding rows from
`LDLR_Tabet_2025_uptake` ("−VLDL" uptake, without VLDL) and
`LDLR_Tabet_2025_abundance` (a cell-surface abundance assay) are dropped
for that same substitution. Substitutions in those ranges that `+VLDL`
doesn't cover keep their `uptake`/`abundance` data unchanged -- this is a
priority with fallback, not a blanket exclusion of the other two assays.
Per the original investigator, `+VLDL` is less subject to a blind spot the
other two assays have specifically in these two modules; the same
prioritization was deliberately *not* extended to LA module 1, which
instead gets the blanket exclusion described below.

Unlike the other special-casing in this section, this rule runs on raw,
per-assay rows before dedup, category split, or `ASSAY_PRIORITY_LIST` ever
run, so it applies identically regardless of dedup strategy and to all
five variant categories alike -- `controls`/`ClinGen_Repo` and
`VUS`/`gnomAD`/`Unobserved` equally. Confirmed directly against the
current data: 10,780 rows (5,392 `abundance`, 5,388 `uptake`) are dropped
by this rule, and every remaining LA2/LA6 substitution that had `+VLDL`
coverage has `LDLR_Tabet_2025_presence_VLDL` as its only surviving assay
in that position range -- `uptake`/`abundance` are retained only for the
substitutions `+VLDL` didn't cover (21 substitutions in each module, per
the reclassification-stage data).

Immediately after that priority rule, in the same pre-dedup step, LA
module 1 (`aa_pos` 25-65) is dropped entirely: every `LDLR` row in that
range is removed outright, regardless of assay, rather than having one
assay prioritized over the others. This is a blanket exclusion, not a
priority-with-fallback -- no LA module 1 substitution survives into
Supplementary Data 5.

### `BRCA1`/`PTEN`/`MSH2`/`TP53`: 2018 vs. 2025 ClinVar vintage

Two different gene lists switch specific genes from the 2025 ClinVar
snapshot this pipeline otherwise uses to the 2018 snapshot, for two
different purposes:

- **`Fxn_points`** (cell 17): `BRCA1`, `PTEN`, `MSH2` use `ExC_points_2018`
  (falling back to `ExC_points_2025` if missing) instead of
  `ExC_points_2025`. `TP53` is *not* in this list -- its `Fxn_points`
  already comes from `OP_points` (see
  [above](#f9-and-tp53-alternate-calibration-and-restricted-datasets)),
  which the code comments note "already considers 2018 clinvar controls."
- **ClinVar star count/significance used for `controls` selection**
  (`clinvar_star_18_25`/`clinvar_18_25`, cells 56/62/77): `BRCA1`, `PTEN`,
  `MSH2`, *and* `TP53` all use the 2018 values instead of 2025.

Per the pipeline README, this reflects "legacy calibrations" for these
genes -- their calibrations were originally derived against 2018 ClinVar
controls, so classification stays internally consistent with that
calibration rather than mixing a 2018-derived calibration with
2025-vintage control selection.

### Predictor training-variant exclusion is group-propagated on the aa side

`revel_train_amino`/`mp2_train_amino` (used to exclude REVEL/MutPred2
training variants from their respective per-predictor outputs) aren't
computed per row the same way on both sides. On the nt side (cell 62)
it's a plain per-row map: `np.where(sankey_nuc['REVEL_train'] == True,
"Yes", "No")`. On the aa side (cells 60-61), it's group-propagated:
`sankey_aa.groupby(group_cols_aa)['REVEL_train'].transform(...)` marks an
*entire* aa-coordinate group `"Yes"` if *any* row within it was a REVEL
training variant -- so if one assay's row for a given amino-acid
substitution happens to be a REVEL training variant, every *other*
assay's row for that same substitution is also excluded from the
REVEL-specific outputs, even if that other assay's own row was never used
to train REVEL. Another nt/aa asymmetry in the same family as the ones
elsewhere in this doc, though in the opposite direction (the aa side is
*more* exclusionary here, not less).

## Deduplication strategy

When the same physical variant is scored by more than one assay --
possibly a mix of nt- and aa-resolution ones -- the pipeline must collapse
that to a single representative record before computing points and a
classification. This section records the decided deduplication
methodology, why it differs between variant categories, and the `"v1"`
legacy behavior it replaced -- the largest and most-revised part of this
notebook's methodology, and the reason this doc originally existed.

### The two categories of variant, and why they're handled differently

**ClinVar controls and ClinGen Evidence Repository controls** are used as
calibration/validation sets -- the same physical variant contributing
evidence twice (once via an nt-resolution assay's score, once via an
aa-resolution assay's score for the same substitution) would double-count
that variant's evidence. These categories must pick exactly one
representative record per variant, and the *resolution* of the assay it
came from matters.

**VUS, gnomAD, and Unobserved variants** are evaluated at DNA resolution
only -- however many assays (nt- or aa-type) scored a given variant, only
one DNA-level classification is ever produced for it, so there's no
double-counting risk. The only question is which of possibly several
scores to trust for that variant.

### What double-counting avoidance means in practice

"Avoid double-counting" is easy to state abstractly but has real,
sometimes surprising consequences for which records actually survive into
the output. Two questions pin it down concretely. Both were confirmed by
running the actual `controls_nuc`/`controls_aa`/`catch_mis_2`/
`dedup_vus_gnomad_unobserved` code (not just read from source) against
small constructed inputs representing a hypothetical amino-acid change
`Gene:p.A100V` reachable by two different nucleotide substitutions --
call them **SNV1** (say chr1:1000 G>A) and **SNV2** (chr1:1003 C>T). These
are two distinct ClinVar records / two distinct rows of input data. (That
more than one nucleotide-level event can map to the same amino-acid
outcome is independently confirmed real in this dataset -- see the
`BRCA1_Adamovich_2022_HDR`/`G6PD_IGVF`/`LDLR_Tabet_2025` examples in
[Comparing `v1` to the decided approach](#comparing-v1-to-the-decided-approach)
and this doc's [Open questions](#open-questions).)

**A single assay never gives SNV1 and SNV2 different aa evidence in this
dataset.** Checked directly: of every (aa-group, `Dataset`) combination
where one assay's rows span more than one genomic coordinate for the same
amino-acid outcome (229,083 such combinations), **none** report more than
one distinct `Fxn_points` value across those coordinates -- makes sense,
since an aa-resolution assay is measuring the protein consequence, and
different codons reaching the same one get scored identically. So for
SNV1's and SNV2's aa evidence to genuinely differ, as the worked example
below needs, that evidence has to come from **two different aa-resolution
assays** scoring the same amino-acid change -- not the same assay
reporting inconsistent results for it.

#### Question 1: two distinct SNVs producing the same protein change

| Assay coverage for this gene | `controls`/`ClinGen_Repo` result | `VUS`/`gnomAD`/`Unobserved` result |
|---|---|---|
| (a) nt-resolution only | Both SNV1 and SNV2 survive, each via its own nt evidence. No interaction between them -- nt-level dedup groups by exact genomic coordinates, so two different SNVs are never in the same group. | Same: both survive independently, for the same reason. |
| (b) nt- *and* aa-resolution | **Both SNV1 and SNV2 still survive -- but only because each has its own nt-resolution record to fall back on.** Under `"nt_then_abs_max"` (the default), each one's *own* aa evidence may never even be compared to its own nt evidence -- see the worked trace below. | Both survive, **each independently choosing its own strongest evidence** (nt or aa, whichever has the greater absolute value under the default `"abs_max"`). SNV1 and SNV2 never interact, because dedup here always groups by exact genomic coordinates, never by amino-acid coordinates. |
| (c) aa-resolution only | **Only one of SNV1/SNV2 survives.** The aa-level dedup stage groups by (`Gene`, `aa_pos`, `aa_ref`, `aa_alt`, transcript) -- not by genomic coordinates -- so SNV1's and SNV2's aa-assay rows collide in the *same* group and get reduced to one row before stage 3 ever runs. **The other SNV is completely absent from the output**, not merely outscored -- there is no row for it at all. | **Both SNV1 and SNV2 survive.** Dedup groups by genomic coordinates only, so an aa-only variant still gets its own group; two different SNVs sharing an aa outcome never collide. |

**Case (b), worked through precisely.** Construct SNV1 with nt evidence
`Fxn_points=3` and aa evidence `5` from `Dataset P`; SNV2 with nt evidence
`6` and aa evidence `8` from `Dataset Q` (a *different* aa-resolution
assay -- see above for why it has to be a different one) -- i.e. for
*each* SNV individually, its own aa evidence is the stronger one:

- **`controls`/`ClinGen_Repo`, any strategy:** stage 2 (aa dedup) first
  picks a single winner *across both SNVs' aa rows*, using the configured
  strategy (assay priority under `"v1"`; magnitude under `"abs_max"`/
  `"nt_then_abs_max"`) -- say SNV2's aa row (`8`) wins that contest. SNV1's
  aa row (`5`) is discarded *right there*, before stage 3 -- it never gets
  a chance to compete against SNV1's own nt row (`3`). Stage 3 then merges:
  SNV2's surviving aa row (`8`) collides with SNV2's own nt row (`6`) at
  the same coordinates -- under `"nt_then_abs_max"` the nt row wins
  regardless of magnitude, so the final SNV2 record uses `6`, not `8`.
  SNV1's nt row (`3`) has no competitor left (its aa row is already gone)
  and survives untouched. **Final: SNV1 → `3`, SNV2 → `6` -- both from nt
  evidence, and neither variant ever got the benefit of comparing its own
  nt and aa evidence to each other**, because the choice was made one
  level up, across SNVs, first.
- **`VUS`/`gnomAD`/`Unobserved`, `"abs_max"` (the default):** each SNV's
  own nt and aa rows are in the *same* group (same coordinates), and
  nothing from the other SNV is in that group. SNV1 compares `3` vs. `5`
  and keeps `5`; SNV2 compares `6` vs. `8` and keeps `8`. **Final: SNV1 →
  `5`, SNV2 → `8`** -- see below for why "correctly compared" needs a
  caveat here too.

**A sharper way to see the cost.** Since SNV1 and SNV2 cause the exact
same protein change, `Dataset P`'s `5` and `Dataset Q`'s `8` aren't really
"SNV1's evidence" and "SNV2's evidence" as separate, competing
quantities -- they're two independent measurements *of the same thing*
(whether the shared amino-acid substitution disrupts function), and the
strongest available evidence for that shared consequence is `8`, full
stop, regardless of which SNV a given ClinVar record happens to describe.
`controls`/`ClinGen_Repo`'s aa-stage dedup effectively recognizes this --
it compares every row sharing the amino-acid identity, across both
datasets, and correctly identifies `8` as the strongest one available.
What it does *not* do is recognize that this shared, strongest evidence is
equally relevant to *any* SNV producing that protein change. Instead, the
winning row carries forward whatever specific genomic coordinates its
assay happened to use (`Dataset Q`'s, which match SNV2) -- so only SNV2
ever gets to weigh `8` against its own nt evidence at the merge stage.
**SNV1 doesn't lose a fair fight against stronger evidence; it never gets
to make its case with the strongest available evidence for its own
protein consequence at all**, because that evidence happened to be filed
under someone else's genomic coordinates. Stated more precisely than "aa
loses to nt when they collide": a variant's aa evidence can be taken out
of contention by a *different* variant's measurement of the identical
protein change, before its own nt-vs-aa comparison ever happens. This
holds under **all three strategies** (`"v1"`, `"abs_max"`,
`"nt_then_abs_max"`) -- switching strategy only changes *which* SNV
inherits the shared evidence and *which* type wins the merge, never
whether one SNV gets frozen out of it entirely.

**This may well be the intended trade-off, not a cost to eliminate.**
`controls`/`ClinGen_Repo` feed calibration/validation exercises where each
control is meant to contribute one independent data point. If SNV1's and
SNV2's aa evidence were *both* separately counted as "the assay correctly
scored a pathogenic control" whenever they share a protein consequence,
a single underlying measurement (`Dataset Q`'s `8`) would inflate the
assay's apparent ClinVar concordance by being credited twice for what is,
functionally, one observation. Collapsing to one row per amino-acid group
before computing any calibration-relevant count is a defensible way to
prevent exactly that. Framed this way, the aa-stage dedup is arguably
doing the *right* thing at the amino-acid level -- the open question is
what it does downstream of that decision (below).

**Which SNV wins, exactly, and does ClinVar/ClinGen evidence quality ever
factor in?** Checked directly: **no.** The aa-stage sort key is
`assay_priority` (`"v1"`) or `abs(Fxn_points)` (`"abs_max"`/
`"nt_then_abs_max"`) -- see [Three-stage pipeline](#three-stage-pipeline-for-controlsclingen_repo-under-v1)
below -- neither one references ClinVar star count, review status, or
ClinGen curation quality anywhere. The *only* ClinVar-quality check that
happens before the aa-stage dedup is a coarse, group-level, binary gate
(`summarize_clnstar`, cell 84): it requires *at least one* row in the
aa-group to have a "1+ star" review status (`criteria provided, single
submitter` counts exactly the same as `reviewed by expert panel` for this
purpose), and it excludes the *entire* group if some rows are "1+ star"
and others are "0 star" (tagged `has_clinvar_star_conflict`). Once a group
clears that bar, **which specific SNV's row goes on to win the aa-stage
dedup, and therefore which SNV's own ClinVar star count/review
status/significance ends up representing the group in the final output,
is decided entirely by the *functional* assay data** -- with no
preference for the SNV backed by the more authoritative ClinVar/ClinGen
record. If SNV1 is `reviewed by expert panel` with weaker functional
evidence and SNV2 is `criteria provided, single submitter` with stronger
functional evidence, SNV2 wins and its (weaker) ClinVar review status is
what the output reports -- SNV1's stronger clinical evidence doesn't
factor into the choice at all. **Whether it should** -- e.g. preferring
the more authoritatively-reviewed ClinVar record when functional evidence
doesn't clearly discriminate, or weighting it into the tie-break directly
-- is an open question; see [Open questions](#open-questions).

**`VUS`/`gnomAD`/`Unobserved` has the same underlying limitation, just in
milder form.** Its single-pass dedup avoids the cross-SNV *collision*
problem above (SNV1 and SNV2 never compete for a single output row), but
it doesn't pool aa evidence across genomic representations of the same
amino-acid change either: grouping strictly by each SNV's own coordinates
means SNV1 only ever sees `Dataset P`'s `5`, never `Dataset Q`'s stronger
`8` -- even though `Q`'s measurement is just as much about the protein
consequence SNV1 causes as it is about SNV2's. No part of this pipeline
pools aa-level evidence across the genomic representations of a shared
amino-acid change; `controls`/`ClinGen_Repo`'s aa-stage dedup is the
*closest* thing to that, and even it only benefits whichever one SNV
happens to inherit the winning row's coordinates.

#### Question 2: one ClinVar record, scored by both an nt- and an aa-resolution assay

This is the simple case that [Decided approach](#decided-approach) below
is about. Confirmed directly:

| | `controls`/`ClinGen_Repo` | `VUS`/`gnomAD`/`Unobserved` |
|---|---|---|
| Result | Exactly one record survives: the **nt**-resolution one, under `"nt_then_abs_max"` (the default) -- *even if the aa-resolution evidence has the larger magnitude*. | Exactly one record survives: **whichever has the greater absolute value**, nt or aa -- under `"abs_max"` (the default), resolution doesn't matter, only magnitude. |

Worked example (nt evidence `3`, aa evidence `5`, one physical variant):
`controls`/`ClinGen_Repo` under `"nt_then_abs_max"` keeps the nt row (`3`),
discarding the larger aa value (`5`) -- the deliberate, by-design cost of
avoiding double-counting for this category (see
[Why "greatest absolute value" rather than signed value](#why-greatest-absolute-value-rather-than-signed-value)
for why the tie-break itself uses magnitude, and
[Decided approach](#decided-approach) for why nt is preferred outright
here in the first place). `VUS`/`gnomAD`/`Unobserved` under `"abs_max"`
keeps the aa row (`5`) instead, since there's no double-counting risk to
guard against for these categories and the larger value wins outright.

### Decided approach

#### ClinVar controls / ClinGen Evidence Repository controls

1. **Always prefer a DNA-resolution (nt) assay's record over an
   amino-acid-resolution (aa) assay's record for the same variant**,
   regardless of either one's point value.
2. Among remaining candidates of the same resolution (multiple nt assays,
   or multiple aa assays when no nt assay covers the variant), **take the
   record with the greatest absolute point value**.

This is implemented as the `"nt_then_abs_max"` value of the
`CONTROLS_CLINGEN_DEDUP_STRATEGY` notebook parameter (now the default).

#### VUS / gnomAD / Unobserved

Simply **take the record with the greatest absolute point value**, with no
preference between nt- and aa-resolution assays.

This is implemented as the `"abs_max"` value of the
`VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY` notebook parameter (now the
default).

### Why "greatest absolute value" rather than signed value

Within either category's tie-break step, the point value with the larger
*magnitude* is preferred over the one closer to zero, regardless of sign.
The pipeline's original (pre-parameter) tie-break for the controls/ClinGen
nt+aa merge instead compared *signed* points, which has a real, documented
failure mode: a positive value always beat a negative one no matter the
magnitude, and between two negative values the one *closer to zero* (the
weaker benign-direction evidence) won -- see
[Stage 3's sign bias](assay_priority_questions.md#known-problems-by-category)
for the mechanism and confirmed examples (e.g. `PALB2` chr16:23623123,
where the original logic kept `Fxn_points` -1 over -5 because -1 > -5,
discarding the stronger evidence). Comparing absolute values instead
avoids that failure mode by construction.

### The `"v1"` strategy in detail

`"v1"` isn't a single, uniform rule -- it's shorthand for "whatever the
pipeline did before either dedup parameter existed," and what that was
differs by category. This section documents it in full, including where it
doesn't match either its own documentation or a coherent design intent.
The evidence and cell-level detail behind all of this lives in
[`docs/assay_priority_questions.md`](assay_priority_questions.md); this
section carries over the tables most relevant to understanding `"v1"` as a
whole.

#### Three-stage pipeline for `controls`/`ClinGen_Repo` under `"v1"`

`controls` and `ClinGen_Repo` are split by assay resolution and go through
three dedup passes, not one -- nt and aa candidates are deduped
separately, then the two survivors are merged and deduped again:

| Stage | Input | Sort key | Uses `assay_priority`? |
|---|---|---|---|
| 1. nt subset dedup | `controls_nuc` / `clingen_nuc` | `VariantNotes` tag (same weak sort as `VUS`) | no |
| 2. aa subset dedup | `controls_aa` / `clingen_aa` | `assay_priority` | **yes** -- the one place `ASSAY_PRIORITY_LIST` is used |
| 3. merge nt+aa survivors, dedup again | `pd.concat([stage-1 output, stage-2 output])` | `catch_mis_2`, *signed* `Fxn_points` | no |

Stage 3 exists because stages 1 and 2 each only guarantee one row *per
assay type* -- if the same physical variant was scored by both an nt assay
and an aa assay, both survive stages 1-2 and collide again on genomic
coordinates, so stage 3 dedups the combined result a second time.

#### All five categories compared

| Category | Split by assay type? | Dedup stages | Uses `assay_priority`? | Supp Data 5 tabs |
|---|---|---|---|---|
| `controls` | yes | 3 (nt / aa / merge, table above) | yes (stage 2 only) | `controls_REVEL/AM/MP2_GeneSpecific` |
| `ClinGen_Repo` | yes | 3 (nt / aa / merge, table above) | yes (stage 2 only) | `ClinGen_Repo_REVEL/AM/MP2_GeneSpecific` |
| `VUS` | no | 1, `VariantNotes` over all rows | no | `VUS_REVEL/AM/MP2` |
| `gnomAD` | no | 1, `VariantNotes` over all rows | no | `gnomAD_REVEL/AM/MP2` |
| `Unobserved` | no | 1, `VariantNotes` over all rows | no | `Unobserved_REVEL/AM/MP2` |

#### Where `ASSAY_PRIORITY_LIST` actually decides an outcome

`ASSAY_PRIORITY_LIST` is consulted only inside the `controls_aa`/`clingen_aa`
dedup (stage 2 above), and only *after* a pre-filter has already run:
`sankey_aa` rows are tagged `VariantNotes = 'max_fxn_pts'` only if their
`abs(Fxn_points)` equals the *maximum* `abs(Fxn_points)` within their own
(`Gene`, `aa_pos`, `aa_ref`, `aa_alt`, transcript) group, and only rows
carrying that tag (plus a matching predictor-max tag) reach the
`assay_priority` sort at all.

**This means `assay_priority` does not simply pick "whichever listed assay
scored the variant" or "the first listed assay that contains a variant" --
magnitude filters first.** If two assays score the same aa-level variant
and one has a strictly larger `abs(Fxn_points)`, that one wins outright;
the smaller one never reaches the `assay_priority` sort, regardless of
either assay's rank (or absence from the list). `assay_priority` only ever
breaks a genuine **tie**: two or more candidate rows that already share
the group's maximum `abs(Fxn_points)` value. That's also why an *unlisted*
assay (falling back to `9999`) can still beat a *listed* one outright, if
its magnitude is strictly larger -- being on the list is not a trump card
over evidence strength, only over other evidence of equal strength.

Checked directly against every aa-type variant for the 9 genes with 2+
relevant aa assays (the 8 genes with 2+ *listed* aa assays, per the table
below, plus `CHEK2`, which has one listed assay competing against one
unlisted one). These counts span the full annotated dataset -- i.e. every
variant these assays score, before the category split and before
`controls`/`ClinGen_Repo`-specific filters like ClinVar star count --
not just the subset that ultimately lands in `controls`/`ClinGen_Repo`,
so they illustrate how often the situation arises for these assays in
general, not a live per-category count:

| Gene | Listed aa assays (rank) | Variants scored by 2+ of them | Decided by a genuine tie (`assay_priority` matters) | Decided by magnitude alone (`assay_priority` irrelevant) |
|---|---|---|---|---|
| `ASPA` | `_abundance` (18), `_toxicity` (19) | 6,151 | 3,914 | 2,237 |
| `BRCA1` | `Adamovich_..._Cisplatin_Resistance` (20), `Adamovich_..._HDR` (21) | 1,086 | 106 | 980 |
| `GCK` | `Gersing_2023_complementation` (26), `Gersing_2024_abundance` (27) | 9,019 | 4,472 | 4,547 |
| `KCNE1` | `_trafficking` (28), `_potassium_flux` (29), `_trafficking_WT_background_DN` (30) | 2,553 | 340 | 2,213 |
| `KCNH2` | `Jiang_2022` (31), `Kozek_Glazer_2020` (32), `O_Neill_2024_surface_expression` (33) | 302 | 50 | 252 |
| `PTEN` | `Matreyek_2018` (37), `Mighell_2018` (38) | 4,252 | 2,590 | 1,662 |
| `SCN5A` | `Glazer_2020` (39), `Ma_2024` (40) | 0 | 0 | 0 |
| `KCNQ4` | `_current_homozygous` (45), `_v12_homozygous` (46) | 3,348 | 75 | 3,273 |
| `CHEK2` | `Gebbia_2024` (22) vs. unlisted `CHEK2_McCarthy-Leo_2024` | 3,485 | 318 | 3,167 |

In every tie, the winner is exactly the tied candidate with the best
(lowest-numbered) rank, confirmed directly -- e.g. `BRCA1`'s 106 ties are
all won by `Adamovich_..._Cisplatin_Resistance` (rank 20) over
`Adamovich_..._HDR` (rank 21); `CHEK2`'s 318 ties are all won by the listed
`Gebbia_2024` over the unlisted (`9999`-fallback) `McCarthy-Leo_2024`.

Two things worth noting:

- **`SCN5A` never actually exercises `assay_priority`, despite having 2
  listed aa assays.** `SCN5A_Glazer_2020` and `SCN5A_Ma_2024` simply never
  score the same variant in this dataset, so the rank between them (39 vs.
  40) has decided zero real outcomes so far -- structurally eligible, but
  empirically inert.
- **For `CHEK2`, the unlisted assay has the larger magnitude more often
  than the listed one does**, among the cases magnitude alone decides
  (2,298 vs. 869 -- not shown separately above, folded into the 3,167
  total). Being on `ASSAY_PRIORITY_LIST` doesn't correlate with an assay
  actually producing the stronger evidence for `CHEK2`; it only wins the
  cases where the two exactly tie.

#### Known problems with `"v1"`

Each category's `"v1"` dedup has a different, specific way of not keeping
the strongest evidence -- none of them is a uniform coin-flip; each is
wrong in one particular, describable direction:

| Category | Problem | Confirmed impact |
|---|---|---|
| `controls` / `ClinGen_Repo` | (1) `assay_priority` only covers the aa half | 170/1,147 nt-type multi-assay `controls` groups picked a different `Dataset` than `assay_priority` would; 0 REVEL/MP2, 1 AM classification flip |
| `controls` / `ClinGen_Repo` | (2) stage 3's sign bias | discards stronger evidence whenever the two survivors disagree in sign, or agree on a negative sign |
| `VUS` | systematic nt-over-aa bias | 921 nt/aa collisions, 46/84/68 REVEL/AM/MP2 flips |
| `gnomAD`, `Unobserved` | same nt-over-aa bias (same code as `VUS`) | not quantified |
| all five | string-sort fragility (tag rename/addition/tie) | none confirmed live -- a latent risk, not an active bug |

(1) and (2)'s counts predate today's `ASSAY_PRIORITY_LIST` trim to
section-1-only (see [`docs/assay_priority_questions.md`](assay_priority_questions.md));
they illustrate the mechanism, not a currently-live count.

#### `VariantNotes` values

Not an input column -- built up inside `Variant_Classification_analysis.ipynb`
on the full per-assay-row dataset, before it's split into
`controls`/`ClinGen_Repo`/`VUS`/`gnomAD`/`Unobserved`:

| Value | Meaning |
|---|---|
| `splice_variant_not_measured` | Splice variant with no functional measurement |
| `start_lost_variant_not_measured` | Start-lost aa variant with no functional measurement |
| `conflicting_fxn_data` | This coordinate group has rows with opposite-sign `Fxn_points` across assays -- see [Handling conflicting evidence](#handling-conflicting-evidence) |
| `First_max_fxn_pts` | Row has the max absolute `Fxn_points` in its nt-level coordinate group |
| `max_fxn_pts` | Same, aa-level coordinate group |
| `''` / `NaN` | None of the above -- not flagged, not its group's max |
| e.g. `splice_variant_not_measured;conflicting_fxn_data` | Multiple conditions, joined with `;` |

`"v1"`'s `VariantNotes`-tag sort for `VUS`/`gnomAD`/`Unobserved` relies on
`First_max_fxn_pts` and `max_fxn_pts` being the only two non-blank values
left by the time dedup runs (rows with the other tags were already dropped
-- see [Handling conflicting evidence](#handling-conflicting-evidence)).
That's also the root of the nt-over-aa bias below: `F` (`First_max_fxn_pts`)
sorts ahead of `m` (`max_fxn_pts`) alphabetically, which has nothing to do
with evidence strength.

#### Gaps between what `"v1"` does and what was intended

- **The pipeline README describes one conceptual strategy that matches
  neither actual implementation.** `README_Variant_Classification_analysis.md`'s
  "Deduplication Logic" section says the pipeline "select[s] variants with
  the highest absolute functional points." That matches none of `"v1"`'s
  three sort keys for `controls`/`ClinGen_Repo` (only one even looks at
  points, and that one uses the *signed* value) nor
  `VUS`/`gnomAD`/`Unobserved`'s alphabetical tag sort.
- **`catch_mis_2`'s own name and original docstring said "highest
  functional points," but it sorted the signed value.** `-1` was treated
  as higher than `-5` because `-1 > -5`, not because `-1` is stronger
  evidence -- the opposite of what "highest absolute...points" would mean,
  and the opposite of what the function's own name (`catch_mis_2` --
  "catch missed [duplicates]") suggests it was meant to do: clean up
  whatever the earlier stages missed, not introduce a new, different bias.
- **`ASSAY_PRIORITY_LIST` ranks nt-type assays too, but was only ever
  consulted for the aa half.** Its own docstring flagged `BRCA2_Hu_2024`
  vs. `BRCA2_IGVF` as an open question implying the nt side mattered, yet
  stage 1 (nt) of the controls/ClinGen_Repo pipeline never looked at it --
  only stage 2 (aa) did.
- **The tag-based `VariantNotes` sort was never an evidence comparison at
  all**, despite standing in for one across three of the five categories.
  It happens to agree with "greatest absolute value" *within* a single
  assay type (exactly one row per group gets tagged, and that tag always
  sorts ahead of blank) -- purely because of how the tagging was set up,
  not because the sort itself compares evidence. The moment two different
  *types* of tag both need comparing (nt vs. aa), it silently falls back
  to alphabetical order.

### Implementation

Both parameters live in `Variant_Classification_analysis.ipynb`, defined
together in the cell just before the `controls`/`ClinGen_Repo` nt/aa split
(cell 79 as of this writing -- cell numbers shift as the notebook is
edited, so search for `CONTROLS_CLINGEN_DEDUP_STRATEGY` if the number is
stale). Each takes the same three values, though what each value means is
implemented separately per category since their original (pre-parameter)
behaviors genuinely differed:

| Value | `controls`/`ClinGen_Repo` (`catch_mis_2`) | `VUS`/`gnomAD`/`Unobserved` (`dedup_vus_gnomad_unobserved`) |
|---|---|---|
| `"v1"` | Original behavior: aa-level ties broken by `ASSAY_PRIORITY_LIST` rank, nt+aa merge broken by *signed* `Fxn_points`. | Original behavior: sorted by `VariantNotes` tag, which has an accidental nt-over-aa bias (`First_max_fxn_pts` sorts ahead of `max_fxn_pts` alphabetically) -- see [VUS: systematic nt-over-aa bias](assay_priority_questions.md#known-problems-by-category). |
| `"abs_max"` | Greatest absolute `Fxn_points` wins, nt and aa candidates treated identically. | Same. |
| `"nt_then_abs_max"` | An nt-type record always wins over an aa-type record; ties within a type broken by greatest absolute `Fxn_points`. **Current default.** | Same rule, offered for completeness/comparison but not the default for this category -- see [above](#the-two-categories-of-variant-and-why-theyre-handled-differently). |

`"v1"` is retained (not removed) so the original behavior can still be
reproduced for audit/comparison purposes; it is not the recommended
setting for either parameter going forward.

### Comparing `v1` to the decided approach

Ran the full pipeline twice, both times against the same input data as
production (`data/output/maves/integrated_variant_effect_dataset.tsv.gz`)
but with the eight datasets added after `v1` was written excluded from the
very start (`BRCA2_Huang_2025_SGE`, `CHEK2_McCarthy-Leo_2024`,
`LDLR_Tabet_2025_uptake`, `LDLR_Tabet_2025_abundance`,
`LDLR_Tabet_2025_presence_VLDL`, `PALB2_Boonen_2026`,
`PALB2_Boonen_2026_SGE`, `TP53_Funk_2025`) -- once with both
`CONTROLS_CLINGEN_DEDUP_STRATEGY` and `VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY`
set to `"v1"` (results at
`data/output/reclassification/integrated_variant_effect_dataset_analysis.v1.csv.gz`
and `data/output/supplementary_data/Supplementary_Data_5.v1.xlsx`), and once
with both left at current defaults (`nt_then_abs_max` / `abs_max`; results at
`data/output/reclassification/integrated_variant_effect_dataset_analysis.current_no_new_data.csv.gz`
and `data/output/supplementary_data/Supplementary_Data_5.current_no_new_data.xlsx`).
Both runs share the identical restricted input and category-assignment code
path -- the eight-dataset exclusion happens in a cell immediately after `pp`
is loaded, before any other processing -- so the only thing that differs
between them is the two dedup-strategy parameters.

**`VUS`, `gnomAD`, and `Unobserved` have zero unmatched rows in either
direction, in every predictor.** This is exactly what the theoretical
argument in [Question 1](#question-1-two-distinct-snvs-producing-the-same-protein-change)
predicts: those three categories' dedup (`dedup_vus_gnomad_unobserved`)
groups by genomic coordinates only, never amino-acid coordinates, so it
structurally can't produce an unmatched row from a strategy difference --
confirmed here by identical row counts and a 1:1 key match between `v1`
and current defaults across the board.

| Category | `v1` rows | current rows (excl. new datasets) | matched | unmatched (`v1`-only / current-only) | dataset pick differs | `Class_*` flips |
|---|---|---|---|---|---|---|
| `controls` × REVEL | 11,026 | 11,026 | 10,119 | 907 / 907 | 115 | 16 |
| `controls` × MP2 | 9,969 | 9,975 | 8,688 | 1,281 / 1,287 | 69 | 8 |
| `controls` × AM | 11,059 | 11,056 | 9,916 | 1,143 / 1,140 | 62 | 15 |
| `ClinGen_Repo` × REVEL | 293 | 293 | 290 | 3 / 3 | 5 | 1 |
| `ClinGen_Repo` × MP2 | 97 | 97 | 97 | 0 / 0 | 8 | 0 |
| `ClinGen_Repo` × AM | 299 | 299 | 295 | 4 / 4 | 22 | 1 |
| `VUS` × REVEL | 16,711 | 16,711 | 16,711 | 0 / 0 | 303 | 9 |
| `VUS` × MP2 | 16,367 | 16,367 | 16,367 | 0 / 0 | 288 | 14 |
| `VUS` × AM | 16,760 | 16,760 | 16,760 | 0 / 0 | 306 | 18 |
| `gnomAD` × REVEL | 32,737 | 32,737 | 32,737 | 0 / 0 | 721 | 10 |
| `gnomAD` × MP2 | 31,759 | 31,759 | 31,759 | 0 / 0 | 678 | 11 |
| `gnomAD` × AM | 32,849 | 32,849 | 32,849 | 0 / 0 | 723 | 14 |
| `Unobserved` × REVEL | 88,270 | 88,270 | 88,270 | 0 / 0 | 3,160 | 0 |
| `Unobserved` × MP2 | 87,354 | 87,354 | 87,354 | 0 / 0 | 3,136 | 0 |
| `Unobserved` × AM | 88,308 | 88,308 | 88,308 | 0 / 0 | 3,161 | 0 |

**`controls`/`ClinGen_Repo`'s unmatched rows are genuinely about the two
dedup strategies differing**, not a category-assignment side effect: both
runs share the identical category-assignment code path, differing only in
`CONTROLS_CLINGEN_DEDUP_STRATEGY`/`VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY`.
Rows are matched between the two runs by (`Gene`, `Chrom`, `hg38_start`, `ref_allele`,
`alt_allele`) -- i.e. by genomic position, which is also each row's final
dedup key. When `v1`'s `assay_priority`-ranked aa-stage winner and the
`abs(Fxn_points)`-ranked aa-stage winner for the same amino-acid
substitution come from two *different* SNVs (see
[Question 1](#question-1-two-distinct-snvs-producing-the-same-protein-change)),
each run reports a different SNV's genomic coordinates as that
substitution's representative row, and a coordinate-keyed join sees one row
vanish from each side's key set rather than "the dataset pick changed."
Confirmed on real data: `HMBS` chr11:119093175 G>A (`HMBS_van_Loggerenberg_2023_combined`,
aa, `Fxn_points`=1, `NP_000181.2:p.Gly326=`) is `v1`'s pick for that
synonymous substitution; current defaults instead pick a *different* SNV at
the same position, chr11:119093175 G>T (`HMBS_van_Loggerenberg_2023_erythroid`,
aa, `Fxn_points`=1, same `hgvs_p`) -- so the coordinate-keyed join sees
neither row as a match for the other, even though both strategies agree on
the winning evidence magnitude. This also explains why `controls`/`ClinGen_Repo`
row *counts* differ only slightly between the two runs: a handful of these
SNV-representative swaps change which row's ClinVar star count/summary
makes a group eligible for `controls` at all, so the two strategies can each
pick up a few rows the other doesn't -- `controls` × MP2 has 9,975 current
rows vs. 9,969 `v1` rows, `controls` × AM has 11,056 vs. 11,059, roughly
symmetric either way rather than one strategy consistently yielding more
rows than the other.

**Nearly all of `controls`/`ClinGen_Repo`'s unmatched rows are this same
re-keying, not a genuinely different set of variants.** Matching each
category's unmatched rows across runs by amino-acid identity (`Gene`,
`hgvs_p`) instead of genomic position finds a same-substitution counterpart
for almost every one of them -- confirming they're the HMBS-style
mechanism above, not unrelated variants that happen to appear on only one
side:

| Category | unmatched (`v1` / current) | re-keyed to the same substitution | true orphans (`v1` / current) | of re-keyed: dataset differs | of re-keyed: `Class_*` flips |
|---|---|---|---|---|---|
| `controls` × REVEL | 907 / 907 | 905 | 2 / 2 | 30 | 0 |
| `controls` × MP2 | 1,281 / 1,287 | 1,276 | 5 / 11 | 54 | 0 |
| `controls` × AM | 1,143 / 1,140 | 1,138 | 5 / 2 | 24 | 0 |
| `ClinGen_Repo` × REVEL | 3 / 3 | 3 | 0 / 0 | 0 | 0 |
| `ClinGen_Repo` × MP2 | 0 / 0 | 0 | 0 / 0 | 0 | 0 |
| `ClinGen_Repo` × AM | 4 / 4 | 4 | 0 / 0 | 0 | 0 |

Folding the re-keyed rows' own dataset/`Class_*` comparison into the main
table's "dataset pick differs" / "`Class_*` flips" columns gives a fuller
picture of how often the two strategies actually disagree on evidence for
the *same* substitution, independent of which SNV ends up representing it:

| Category | dataset pick differs (genomic-key + re-keyed) | `Class_*` flips (genomic-key + re-keyed) |
|---|---|---|
| `controls` × REVEL | 115 + 30 = 145 | 16 + 0 = 16 |
| `controls` × MP2 | 69 + 54 = 123 | 8 + 0 = 8 |
| `controls` × AM | 62 + 24 = 86 | 15 + 0 = 15 |
| `ClinGen_Repo` × REVEL | 5 + 0 = 5 | 1 + 0 = 1 |
| `ClinGen_Repo` × MP2 | 8 + 0 = 8 | 0 + 0 = 0 |
| `ClinGen_Repo` × AM | 22 + 0 = 22 | 1 + 0 = 1 |

Even after folding in every re-keyed row, **zero additional classification
flips appear** -- the re-keyed rows change which dataset/SNV is on record
for a substitution far more often (24-54 of them) than they change its
evidence tier. What's left after re-keying -- 0-11 rows per
category/predictor -- are true orphans with no same-substitution
counterpart on the other side at all: cases where the winning SNV's *own*
ClinVar star count/summary is what made its aa-group eligible for
`controls`/`ClinGen_Repo` in the first place, so switching which SNV wins
changes category membership itself, not just which row represents an
already-shared substitution. `VUS`/`gnomAD`/`Unobserved` need no such
reconciliation -- they already have zero unmatched rows to begin with.

**As found previously, tie artifacts account for most "dataset pick
differs" cases, and every classification flip is genuine.** Splitting
each category's differing picks into tie artifacts (`abs(Fxn_points)`
unchanged, only *which* tied row got reported changed) versus genuine
differences (magnitude actually changed):

| Category | tie-artifact / genuine picks | flips within tie-artifact / genuine |
|---|---|---|
| `controls` × REVEL | 93 / 22 | 0 / 16 |
| `controls` × MP2 | 58 / 11 | 0 / 8 |
| `controls` × AM | 41 / 21 | 0 / 15 |
| `ClinGen_Repo` × REVEL | 4 / 1 | 0 / 1 |
| `ClinGen_Repo` × MP2 | 8 / 0 | 0 / 0 |
| `ClinGen_Repo` × AM | 21 / 1 | 0 / 1 |
| `VUS` × REVEL | 279 / 24 | 0 / 9 |
| `VUS` × MP2 | 264 / 24 | 0 / 14 |
| `VUS` × AM | 282 / 24 | 0 / 18 |
| `gnomAD` × REVEL | 703 / 18 | 0 / 10 |
| `gnomAD` × MP2 | 661 / 17 | 0 / 11 |
| `gnomAD` × AM | 705 / 18 | 0 / 14 |
| `Unobserved` × REVEL | 3,160 / 0 | 0 / 0 |
| `Unobserved` × MP2 | 3,136 / 0 | 0 / 0 |
| `Unobserved` × AM | 3,161 / 0 | 0 / 0 |

Zero flips came from a tie-artifact pick in any of the fifteen
category/predictor combinations. (`Unobserved` has no ClinVar/ClinGen
classification, so `Class_*` doesn't apply there -- all its dataset-pick
differences are tie artifacts and none can flip anything.)

**Almost every genuine difference and every flip is `BRCA1`, with one
`SCN5A` exception in `gnomAD`.** Grouping genuine (non-tie-artifact) flips
by gene:

| Category | Genuine flips | Gene breakdown |
|---|---|---|
| `controls` × REVEL / MP2 / AM | 16 / 8 / 15 | 100% `BRCA1` |
| `ClinGen_Repo` × REVEL / AM | 1 / 1 | 100% `BRCA1` |
| `VUS` × REVEL / MP2 / AM | 9 / 14 / 18 | 100% `BRCA1` |
| `gnomAD` × REVEL / MP2 / AM | 10 / 11 / 14 | `BRCA1` (9/10/13) + `SCN5A` (1/1/1) |

Two mechanisms account for essentially the entire effect:

1. **`catch_mis_2`'s nt/aa merge, `BRCA1`.** `BRCA1_Findlay_2018` (nt) and
   `BRCA1_Adamovich_2022_HDR`/`_Cisplatin_Resistance` (aa) frequently
   cover the same genomic position with comparable-magnitude
   `Fxn_points`. `v1`'s signed-descending sort lets an aa-type row beat
   an nt-type row whenever its signed value is greater -- including a
   less-negative aa value beating a more-negative nt value. E.g.
   chr17:43070959 A>G: `v1` keeps `BRCA1_Adamovich_2022_HDR` (aa,
   `Fxn_points`=0) over `BRCA1_Findlay_2018` (nt, `Fxn_points`=−5) since
   0 > −5, classifying it *Uncertain*; `nt_then_abs_max` keeps Findlay's
   nt row regardless of sign, classifying it *Likely Benign*. This is the
   same sign-bias mechanism described under [Why "greatest absolute
   value" rather than signed value](#why-greatest-absolute-value-rather-than-signed-value),
   now confirmed to be almost entirely responsible for the real-world
   `BRCA1` classification differences in every category that includes
   this gene's variants.
2. **`dedup_vus_gnomad_unobserved`'s magnitude comparison, `SCN5A`.**
   Unlike `controls`/`ClinGen_Repo`, `VUS`/`gnomAD`/`Unobserved` dedup all
   candidates for a genomic position in one pass regardless of
   resolution. One `gnomAD` position has two aa-type candidates from
   different assays: `SCN5A_Glazer_2020` (`Fxn_points`=0) and
   `SCN5A_Ma_2024` (`Fxn_points`=3). `v1`'s `VariantNotes`-alphabetical
   sort picks Glazer (*Uncertain*); `abs_max` correctly picks the larger
   magnitude, Ma_2024 (*Likely Pathogenic*).

`PALB2` and `TP53` -- flagged as differing genes in the earlier
`assay_priority_questions.md` comparison -- don't appear in either the
raw dataset-pick-differs or the genuine-flip breakdown here: their new
post-`v1` datasets (`PALB2_Boonen_2026`, `PALB2_Boonen_2026_SGE`,
`TP53_Funk_2025`) are excluded from both sides in this comparison, and
whatever remaining evidence they have from `v1`-era datasets doesn't
produce a dataset-pick difference under either strategy. `LDLR` has zero
rows in this comparison at all -- every `LDLR` dataset postdates `v1` and
is excluded from both sides. `GCK`, `HMBS`, `BRCA1`, `PAX6`, and `ASPA`
dominate the *raw* `controls`/`ClinGen_Repo` pick-differs counts (1-54
changed picks each, across predictors -- `GCK` alone accounts for 54/34/7
of `controls`'s REVEL/MP2/AM picks-differ counts) but, aside from `BRCA1`'s
already-counted genuine differences above, contribute zero genuine
differences and zero flips -- every other gene's pick change is a tie
artifact.

## Outputs

Each of the five categories × three predictors (REVEL, AlphaMissense,
MutPred2) produces one output table, assembled into a `dfs` dict (cell
113) with 15 entries: `controls_REVEL/AM/MP2_GeneSpecific`,
`ClinGen_Repo_REVEL/AM/MP2_GeneSpecific`, `VUS_REVEL/AM/MP2`,
`gnomAD_REVEL/AM/MP2`, `Unobserved_REVEL/AM/MP2`. Every one of them has its
intermediate working columns (`assay_priority`, `VariantNotes`, the
training-variant helper columns, and others) dropped via a shared
`COLUMNS_TO_DROP` list, then gets a shared `rename_dict` applied
(`Total_Points_GeneSpecific_REVEL` → `Total_Points_REVEL`, etc.) and is
reordered to a fixed `column_order`. Each table is written to its own CSV
under `data/output/predictor_calibration/gene_specific/`, then all 15 are
bundled into `Supplementary_Data_5.xlsx`/`.xlsx.gz` (one sheet per table).
See `notebooks/analysis/README_Variant_Classification_analysis.md` for the
full output file list and per-column descriptions.

## See also

- [`docs/assay_priority_questions.md`](assay_priority_questions.md) --
  background on `nucleotide_or_aa`, the `controls_aa`/`clingen_aa` split,
  the `ASSAY_PRIORITY_LIST` open questions, and the empirical comparison
  referenced above.
- `notebooks/analysis/README_Variant_Classification_analysis.md` --
  overall pipeline inputs/outputs and usage.

## Open questions

### Does the ClinVar submitter-conflict exclusion need an nt-side check too?

`controls` excludes variants where ClinVar submitters disagree on
pathogenicity/benignity (cell 78: `controls[controls['clinvar_conflict_flag_18_25']
!= 'has clinvar conflict']`). The flag it checks, `clinvar_conflict_flag_18_25`,
is computed **only for `sankey_aa`** (cell 57):

```python
sankey_aa['clinvar_conflict_flag_18_25'] = (
    sankey_aa.groupby(group_cols_aa)['clinvar_18_25']
    .transform(lambda x: 'has clinvar conflict' if has_conflict(x) else np.nan)
)
```

`sankey_nuc` never computes this column at all (compare to cell 62, its
aa-side counterpart for the analogous `clinvar_18_25`/`clnsig_group_18_25`
setup, which has no equivalent conflict check). After
`sankey_f = pd.concat([sankey_nuc, sankey_aa])`, nt-type rows get `NaN`
for `clinvar_conflict_flag_18_25` -- and `NaN != 'has clinvar conflict'`
evaluates to `True` in pandas, so **nt-type `controls` variants can never
be excluded for a ClinVar submitter conflict, regardless of whether they
actually have one.** This is the same shape as the already-documented
[`assay_priority` only covers the aa half](assay_priority_questions.md#known-problems-by-category)
gap -- an aa-only mechanism silently not applying to nt-type rows -- just
for a different exclusion rule.

**Confirmed directly against the current data:**

- `0` of `168,237` nt-type rows in the full annotated dataset have a
  non-null `clinvar_conflict_flag_18_25` -- the column is never populated
  for any nt-type row, confirming the gap is real and total, not partial.
- Applying the aa-side `has_conflict` logic to nt-type `controls`-eligible
  rows directly (grouped by genomic coordinates instead of aa coordinates)
  finds `0` variant groups that would actually be flagged. So, as with the
  `assay_priority` gap, **this is confirmed live but not currently
  consequential** -- no nt-type `controls` variant in the current dataset
  is known to have a genuine ClinVar submitter conflict that's slipping
  through.

**Open questions:**

1. Should `sankey_nuc` get the same `has_conflict` check, grouped by
   genomic coordinates, so nt-type `controls` variants are protected the
   same way aa-type ones are -- even though it isn't changing any output
   today, `BRCA2`'s nt-assay count (and any other gene's) keeps growing,
   the same way it does for the `assay_priority` gap?
2. If so, should it reuse `ClinGen_Repo`'s parallel construction (cells
   104-111 largely mirror `controls`'s aa/nt split) or be factored into one
   shared helper, given both categories would then need the identical
   nt-side conflict check that already exists for their aa side?
3. Is there a reason this check was aa-only by design (rather than by
   omission) that isn't apparent from the code -- e.g. something about how
   nt-level ClinVar conflicts are already handled upstream, elsewhere in
   the pipeline, that would make a redundant check here unnecessary?

### Should ClinVar/ClinGen evidence quality factor into which SNV represents a shared protein consequence?

When two distinct ClinVar-classified SNVs produce the same amino-acid
change and both have aa-resolution evidence, `controls`/`ClinGen_Repo`'s
aa-stage dedup collapses them to one representative row -- plausibly the
right thing to do, so a single functional measurement doesn't get counted
as two independent calibration hits (see
[What double-counting avoidance means in practice](#what-double-counting-avoidance-means-in-practice)).
That representative is chosen **purely by the functional assay data**
(`assay_priority` under `"v1"`, `abs(Fxn_points)` under `"abs_max"`/
`"nt_then_abs_max"`) -- confirmed directly, neither sort key references
ClinVar/ClinGen evidence at all. The only clinical-quality check that runs
first is coarse and binary: does *at least one* row in the amino-acid
group have a "1+ star" ClinVar review status, with no mix of "1+ star" and
"0 star" rows in the same group. Within the "1+ star" bucket, a
single-submitter record and an expert-panel-reviewed record are
indistinguishable to this gate.

So today, if SNV1 carries the stronger ClinVar/ClinGen evidence (say,
expert-panel-reviewed) but SNV2's assay row has the larger functional
score, SNV2 wins -- and SNV2's own (weaker) ClinVar star count/review
status/significance is what ends up representing the group in the final
`controls`/`ClinGen_Repo` output, not SNV1's.

**Open questions:**

1. Should the choice instead prefer the SNV with the more authoritative
   ClinVar/ClinGen record -- e.g. as a tie-break when functional evidence
   is close, or unconditionally, on the theory that calibration integrity
   depends on the *clinical* truth label being as reliable as possible,
   independent of which functional measurement happens to be strongest?
2. Or is functional-evidence-only selection actually preferable -- e.g.
   because a stronger functional signal is itself a proxy for a
   cleaner/more reliable measurement, and ClinVar review status doesn't
   necessarily track how *functionally* informative a given variant's
   assay data is?
3. Does this matter in practice today, or -- like the `assay_priority` and
   `clinvar_conflict_flag_18_25` gaps above -- is it a real mechanism with
   negligible current impact? Nobody's checked how often SNV1/SNV2-style
   collisions with meaningfully different ClinVar review status actually
   occur in the current `controls`/`ClinGen_Repo` data.

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
  flag, and the 2018-vs-2025 ClinVar vintage switch.
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

> **These numbers are carried over as-is from the empirical comparison in**
> [`docs/assay_priority_questions.md`, section 3](assay_priority_questions.md#3-controls_clingen_dedup_strategy-implemented-options-and-empirical-comparison)
> **and are stale -- pending a rerun.** They were computed before
> `ASSAY_PRIORITY_LIST` was trimmed to section-1-only, before
> `VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY` existed, and before `"current"`
> was renamed to `"v1"` (the rename doesn't itself change any number here
> -- same code, same data -- but the `ASSAY_PRIORITY_LIST` trim does, since
> it changes what `"v1"` actually does for `controls`/`ClinGen_Repo`'s
> aa-stage tie-break). We'll refresh these counts and add the
> VUS/gnomAD/Unobserved comparison later.

Ran the full pipeline once per `CONTROLS_CLINGEN_DEDUP_STRATEGY` value
against the real dataset (1,354,282 input rows) and compared the resulting
`controls`/`ClinGen_Repo` `*_GeneSpecific` tables pairwise, matching rows
across runs by (`Gene`, `Chrom`, `hg38_start`, `ref_allele`, `alt_allele`).
Reproducibility was checked directly -- `v1` was run twice and the six
output tables were byte-identical both times.

| Category | Rows matched across all 3 runs | Dataset pick differs: v1→abs_max / v1→nt_then_abs_max / abs_max→nt_then_abs_max | `Class_GeneSpecific_*` flips: v1→abs_max / v1→nt_then_abs_max / abs_max→nt_then_abs_max |
|---|---|---|---|
| `controls` × REVEL | 11,359 | 350 / 355 / 17 | 29 / 31 / 12 |
| `controls` × MP2 | 9,387 | 175 / 165 / 12 | 19 / 11 / 8 |
| `controls` × AM | 11,190 | 238 / 244 / 14 | 16 / 23 / 13 |
| `ClinGen_Repo` × REVEL | 435 | 99 / 100 / 1 | 3 / 4 / 1 |
| `ClinGen_Repo` × MP2 | 128 | 26 / 26 / 0 | 0 / 0 / 0 |
| `ClinGen_Repo` × AM | 442 | 102 / 103 / 1 | 2 / 3 / 1 |

**Most "dataset pick differs" cases turned out to be tie artifacts, not
policy differences.** Splitting each `v1`-vs-other comparison into rows
where the two strategies picked the same `abs(Fxn_points)` value (a tie
artifact -- switching sort key only changed *which* tied row got reported)
versus rows where the picked value's magnitude genuinely differed:

| Category | v1→abs_max: tie-artifact / genuine | flips within tie-artifact / genuine | v1→nt_then_abs_max: tie-artifact / genuine | flips within tie-artifact / genuine |
|---|---|---|---|---|
| `controls` × REVEL | 306 / 44 | 0 / 29 | 305 / 50 | 0 / 31 |
| `controls` × MP2 | 141 / 34 | 0 / 19 | 140 / 25 | 0 / 11 |
| `controls` × AM | 202 / 36 | 0 / 16 | 201 / 43 | 0 / 23 |
| `ClinGen_Repo` × REVEL | 96 / 3 | 0 / 3 | 96 / 4 | 0 / 4 |
| `ClinGen_Repo` × MP2 | 25 / 1 | 0 / 0 | 25 / 1 | 0 / 0 |
| `ClinGen_Repo` × AM | 99 / 3 | 0 / 2 | 99 / 4 | 0 / 3 |

Every classification flip, in every category, fell in the "genuine"
column -- zero flips came from a tie-artifact row. Raw "dataset pick
differs" counts substantially overstate how much the strategy choice
actually matters.

**The genuine differences and all flips concentrated in `PALB2`, `BRCA1`,
and `TP53` -- nowhere else.** `LDLR` and `GCK` dominated the *raw*
pick-differs counts (73-180 and 23-71 changed picks respectively, across
predictors) but contributed **zero** genuine differences and zero flips in
every category checked -- every one of their pick changes was a tie
artifact, confirming the [open question](assay_priority_questions.md#2-open-questions-on-assay-priority-order-for-oddspath-calibration)
that `LDLR`'s assays have no reviewed priority order and frequently score
identically for the same variant.

| v1→abs_max, genuine diffs / flips | REVEL | MP2 | AM |
|---|---|---|---|
| `PALB2` | 27 / 16 | 22 / 13 | 21 / 11 |
| `BRCA1` | 12 / 11 | 6 / 6 | 9 / 5 |
| `TP53` | 5 / 2 | 6 / 0 | 6 / 0 |

Two `PALB2` differences were verified directly against the raw
per-assay-row data (clean cases -- exactly one nt- and one aa-type
candidate, no degenerate-codon duplicates), confirming the sign-bias
mechanism above reproduces on real `controls` data: chr16:23623123 A>G
(`F948L`) -- `v1` keeps `PALB2_Boonen_2026` (aa, `Fxn_points` −1) over
`PALB2_IGVF` (nt, `Fxn_points` −5) since −1 > −5 → *Likely Benign*;
`abs_max` correctly keeps the stronger −5 → *Benign*. `BRCA1`'s
genuine-diff rows were harder to attribute cleanly (see the full writeup
in `docs/assay_priority_questions.md` for why), and `TP53`'s weren't
individually audited.

No VUS/gnomAD/Unobserved comparison exists yet -- `VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY`
didn't exist when this comparison was run. That, and a rerun against the
trimmed `ASSAY_PRIORITY_LIST`, are still pending.

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

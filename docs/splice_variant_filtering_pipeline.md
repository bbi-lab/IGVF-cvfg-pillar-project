# Splice variant filtering pipeline

How a splice-related variant is annotated, flagged, and -- for datasets not
curated as able to detect splicing effects -- excluded, from
`integrated_variant_effect_dataset.tsv.gz` through Supplementary Data 5 and
6, the reclassification export, and Figures 4, 5, and 6.

**Status: this document covers both the preprint's original splice-variant
handling and the final (current) handling after the `splice_measure`
exception was implemented.** Before the fix, every consumer but Figure 4
excluded splice-flagged rows unconditionally, regardless of whether the
assay could detect splicing. See [`docs/splice_variant_assay_consistency.md`](splice_variant_assay_consistency.md)
for the motivation behind the change and its observed impact on each
output; Section 7 here describes what changed mechanically, and Section 8's
typology table shows preprint vs. final side by side, per category and per
output.

## Field glossary

Five fields recur throughout this trace:

| Field | Level | What it is |
|---|---|---|
| `splice_variant` | per DNA-candidate row | `Yes` if any of `spliceAI_DS_AG/DS_AL/DS_DG/DS_DL` &ge; 0.2, or `simplified_consequence == "splice_site_variant"`; else `No`. |
| `splice_var_amino` | one value per amino-acid position | The field every downstream exclusion filter actually checks (see below). |
| `splice_measure` | dataset/assay | The "Detects Splicing Variants?" curation flag from `Supplementary_Data_3.xlsx` ("Curation" sheet). Now checked by every consumer in this trace. |
| `nucleotide_or_aa` | row resolution | How MaveDB reports the score (`nt` or `aa`). Independent of `splice_measure` -- e.g. `PALB2_Boonen_2026_SGE` is `aa`-resolution yet `splice_measure == "Yes"`. |
| `VariantNotes` | row, audit/exclude tag | A second, parallel exclusion path alongside `splice_var_amino`. |

## 1. Source annotation -- `vendor/variant-annotation/` pipeline

Every DNA-candidate row gets:

- **`consequence`** -- the raw VEP term, e.g. `splice_donor_region_variant^intron_variant`. No filtering happens at this stage.
- **`simplified_consequence`** -- mapped via `data/input/reference/extended_ensembl_consequence.csv.gz` into one of two distinct SO summary terms:

  | `simplified_consequence` | Impact | Source VEP terms |
  |---|---|---|
  | `splice_site_variant` | HIGH | `splice_acceptor_variant`, `splice_donor_variant`, `"Canonical splice"` |
  | `splicing_variant` | LOW | `splice_region_variant`, `splice_donor_region_variant`, `splice_donor_5th_base_variant`, `splice_polypyrimidine_tract_variant`, `"splice region"` / `"Splice region"` |

  By deliberate, documented choice, only the exact string `splice_site_variant`
  is ever treated as splice-affecting downstream (Section 2) -- the whole
  `splicing_variant` category is not, regardless of dataset or context. This
  choice is unaffected by the `splice_measure` fix described below.
- **`spliceAI_DS_AG` / `DS_AL` / `DS_DG` / `DS_DL`** (+ `DP_*`) -- raw SpliceAI delta scores. Can be `NaN` (a complex delins near a splice junction may score `NaN` on all four). No filtering happens here either, and by the same deliberate choice, a missing score is never treated as "unknown" -- only as "not splice-affecting."
- **`splice_measure`** -- merged from `Supplementary_Data_3.xlsx` ("Curation" sheet), column "Detects Splicing Variants?". Dataset-level, not row-level.
- **`nucleotide_or_aa`** -- a resolution fact, not a splicing-detection fact.

Output: `integrated_variant_effect_dataset.tsv.gz` (+ `.condensed.tsv.gz`).

## 2. Flag computation -- `src/annotate_simplified_consequence.py` (Step 17)

`compute_splice_variant` (per DNA-candidate, pipe-delimited across positions within a row):

```
Yes  if  any(DS_AG, DS_AL, DS_DG, DS_DL) >= 0.2
       OR simplified_consequence == "splice_site_variant"
else No
```

`compute_splice_var_amino` collapses the candidates at one amino-acid
position to a single value: blank if none, the shared value if all agree,
`Yes` if any is `Yes` or if they disagree. **This is the field every
exclusion filter downstream actually checks** -- `splice_variant` itself is
only reused (recomputed identically) for `VariantNotes` tagging in the
classification and OddsPath notebooks (Section 3). Unchanged by the fix --
this step still doesn't know about `splice_measure` at all; the exception is
applied everywhere this field is later used to exclude a row.

Neither [`docs/annotate_simplified_consequence.md`](annotate_simplified_consequence.md) nor
[`docs/variant_annotation_pipeline.md`](variant_annotation_pipeline.md) (the two docs that describe this step)
document this threshold/filter logic -- it's documented only in the script's
own docstring and here.

## 3. Classification notebook -- `notebooks/analysis/Variant_Classification_analysis.ipynb`

Input: `integrated_variant_effect_dataset.tsv.gz` + `Supplementary_Data_4.xlsx` calibrations.

1. **Build `sankey` (`pp_ex_OP`)** -- merge in ExCALIBR + gene-specific + OddsPath calibrations; recompute `splice_variant` with the identical rule from Section 2. `splice_var_amino` is then rebuilt, *differently per resolution*:
   - `nt` &rarr; `= splice_variant` (a literal passthrough -- an nt-resolution row can never disagree with itself).
   - `aa` &rarr; `any()` across nucleotide realizations sharing one amino-acid substitution (grouped by `Gene, aa_ref, aa_pos, aa_alt, Ref_seq_transcript_ID_stripped`) -- a sibling realization's `Yes` can flag a row that carries no splicing evidence of its own (see typology category 7, Section 8). **This broadcast is unchanged by the fix** -- once a `splice_measure == 'Yes'` dataset's flagged row is spared exclusion (steps 2 and 6 below), every sibling nucleotide realization in its amino-acid group is spared along with it, including ones with no individual evidence of their own. This was a known, explicit decision at implementation time, not an oversight.
2. **Tag `VariantNotes`** (fixed):
   ```
   splice_variant_not_measured  <-  splice_variant == 'Yes'  AND  splice_measure != 'Yes'
   ```
   Previously gated on `nucleotide_or_aa == 'aa'` instead of `splice_measure` -- an aa-resolution, splice-aware dataset (`PALB2_Boonen_2026_SGE`) got tagged "not measured" regardless of what the assay could actually detect, while nt-resolution non-splice-aware datasets (7 datasets, 281 rows: `LARGE1_Ma_2024`, `JAG1_Gilbert_2024`, `CBS_Sun_2020_high_B6`, `CBS_Sun_2020_low_B6`, `TPK1_Weile_2017`, `BRCA2_Hu_2024`, `FKRP_Ma_2024`, `RHO_Wan_2019`) never got tagged at all. Both are corrected by keying on `splice_measure` directly. `start_lost_variant_not_measured` is untouched (`nucleotide_or_aa == 'aa' AND consequence == 'start_lost'`) -- start-lost handling was out of scope for this fix.
3. **Split `nt` / `aa` &rarr; per-position dedup across assays**: `conflicting_fxn_data` (opposite-sign `Fxn_points` at the same `{Gene, Chrom, pos, ref, alt}` across datasets); `First_max_fxn_pts` (nt) / `max_fxn_pts` (aa) picks a representative row per position by max `|Fxn_points|`. A priority tag only lands on a row whose `VariantNotes` is still blank -- a row already carrying a "not_measured" or conflict tag can block its whole group from getting a representative. Splice-aware rows no longer carry the "not_measured" tag (step 2), so they no longer block their group here either.
4. **Recombine** -- `sankey_f = concat(sankey_nuc, sankey_aa)`, plus one more `conflicting_fxn_data` pass across the full set. `Gene == 'SFPQ'` is dropped entirely here, unrelated to splicing.
5. **Checkpoint** -- write `sankey_g` (dropping ExCALIBR-only columns) to `integrated_variant_effect_dataset_analysis.csv.gz`, then immediately reload it as `sankey_f`. This file is written *before* the exclusion filter below runs, so every category in the typology survives here -- it's also read independently by the reclassification script (Section 4).
6. **The filter** (fixed) -- exclude rows where
   ```
   VariantNotes in {conflicting_fxn_data, splice_variant_not_measured,
                     splice_variant_not_measured;conflicting_fxn_data,
                     start_lost_variant_not_measured}
   OR (splice_var_amino == 'Yes' AND splice_measure != 'Yes')
   ```
   then, separately, drop `Flag == '*'`. Previously `splice_var_amino == 'Yes'` excluded unconditionally; now a row from a `splice_measure == 'Yes'` dataset survives this line regardless of its own `splice_var_amino` value.
7. **Category split**: `controls` (ClinVar Benign/Pathogenic groups; re-checks `clinvar_conflict` + the same `splice_var_amino`/`splice_measure` condition a second time, redundantly, since `sankey_f` is already filtered), `ClinGen_Repo` (`Updated_Classification_ClinGen_repo` present, not VUS), `VUS`, `gnomAD`, `Unobserved`. Deduplication when the same variant is scored by more than one assay is governed by two independently configurable strategies (`v1` / `abs_max` / `nt_then_abs_max`; see [`docs/variant_classification.md`](variant_classification.md)): `CONTROLS_CLINGEN_DEDUP_STRATEGY` for `controls`/`ClinGen_Repo`, `VUS_GNOMAD_UNOBSERVED_DEDUP_STRATEGY` for the other three. In `v1`/`nt_then_abs_max`, nt-resolution evidence is preferred outright over aa-resolution evidence for the same variant.
8. **Output**: `Supplementary_Data_5.xlsx` -- 15 sheets = 5 categories &times; {REVEL, AM, MP2} (+ a `.with_secondary_variants` sibling file, same rule).

## 4. Reclassification export -- `src/build_variant_reclassification_dataset.py`

Reads the *same* checkpoint (`integrated_variant_effect_dataset_analysis.csv.gz`)
independently of Section 3's own filter, and applies its own (fixed):

```
exclude rows where
  VariantNotes in {splice_variant_not_measured, splice_variant_not_measured;conflicting_fxn_data}
  OR (splice_var_amino == 'Yes' AND splice_measure != 'Yes')
```

Narrower disallowed-set than Section 3 -- doesn't independently exclude bare
`conflicting_fxn_data` or `start_lost_variant_not_measured` rows here (though
in practice `splice_var_amino`/`splice_measure` alone still catches every
non-splice-aware splice-flagged row the same way Section 3 does).

Output: `integrated_variant_effect_reclassification.tsv.gz` (the biobank-facing export).

Covered by `tests/test_build_variant_reclassification_dataset.py::test_splice_var_amino_kept_when_splice_measure_yes`.

## 5. OddsPath calibration -- `OddsPath_calculations.ipynb` &rarr; `OddsPath_classifications.ipynb`

Rebuilds `pp`/`pp_ex_OP` independently of Section 3 (same source tables, its
own merge), then:

- `OddsPath_calculations.ipynb`: `Flag != '*' AND (splice_var_amino != 'Yes' OR splice_measure == 'Yes')` --
  fixed, and this feeds the OddsPath/ACMG point calibration itself, not just
  Data 6.
- `OddsPath_classifications.ipynb`: two fixes --
  - `VariantNotes_OP` tagging (this notebook's own equivalent of Section 3
    step 2) now reads `(splice_var_amino == 'Yes') AND (splice_measure != 'Yes')`,
    previously just `splice_var_amino == 'Yes'` with no `splice_measure`
    check at all.
  - the same disallowed-`VariantNotes_OP` + `splice_var_amino`/`splice_measure`
    filter as Section 3's step 6, applied in two places (`OP_sankey_full` and
    `controls_OP_18_x2`).

Output: `Supplementary_Data_6.xlsx` (`controls` / `ClinGen_Repo` / etc. &times;
{REVEL, AM, MP2}, `_OP` tabs).

## 6. Figures

| Consumer | Reads | Filter |
|---|---|---|
| **Figure 4** (`src/build_figure4_data.py` &rarr; `figure_4/data_utils.py: Scoreset.splicing_filter`) | `integrated_variant_effect_dataset.tsv.gz` directly | `if splice_measure == 'Yes': keep splice variants; else: drop as usual`. The original assay-aware filter this fix generalized to the rest of the pipeline. |
| **Figures 5 & 6** (`Figure5_6.Rmd`) | `Supplementary_Data_5.xlsx` sheets only (`controls`/`ClinGen_Repo`/`VUS`/`gnomAD`/`Unobserved` &times; REVEL/AM/MP2) | Inherits Data 5's fixed filter as-is -- no changes needed in this file itself. |
| **Extended Data figures 4, 6, 7, 8, 9** (`Extended_data_figures.Rmd`) | Both `Supplementary_Data_5.xlsx` and `Supplementary_Data_6.xlsx` sheets | Inherits both fixed filters -- no changes needed in this file itself. |

## 7. What was fixed

1. **Data 5 & 6 now check `splice_measure`.** Both the direct
   `splice_var_amino`-based filter and the parallel `VariantNotes`
   disallowed-tag filter spare a row when its dataset is curated
   `splice_measure == 'Yes'` -- the 22 SGE datasets able to detect splicing
   effects (per `Supplementary_Data_3.xlsx`) now keep their splice-flagged
   variants in Data 5, Data 6, the reclassification export, and every figure
   built from them.
2. **`VariantNotes` tagging now keys on the right axis.** It reads
   `splice_measure` (assay capability) instead of `nucleotide_or_aa` (score
   resolution) -- fixing both directions of the old mismatch: an
   aa-resolution splice-aware dataset is no longer wrongly excluded, and
   nt-resolution non-splice-aware datasets are now correctly tagged too
   (though the direct `splice_var_amino` filter already caught those either
   way).
3. **Figure 4's pattern is now applied everywhere.** `Scoreset.splicing_filter`
   was, before this fix, the only place in the codebase that skipped the
   exclusion when `splice_measure == "Yes"`. Every other consumer (Data 5,
   Data 6, the reclassification export, and everything built from them) now
   applies the same exception.

**An explicit, un-reconsidered decision carried over from planning:** the
aa-resolution "guilt by association" broadcast (Section 3, step 1) is
untouched. A sibling nucleotide realization's `Yes` still spares the whole
amino-acid-substitution group once `splice_measure == 'Yes'`, even for
realizations with no individual splicing evidence of their own (typology
category 7). This was flagged as an open question during planning and
resolved by *not* special-casing it -- the simpler, uniform rule was applied
consistently rather than adding resolution-specific logic on top of the
`splice_measure` fix.

Across all seven outputs, Data 5, the reclassification export, Data 6,
Figures 5 & 6, and the Extended Data figures now always agree row-for-row --
and now agree with Figure 4 too, for every category in the typology below.

## 8. A typology of splice-related variants, by output -- preprint vs. final

Seven categories, each anchored to a real row in the dataset, run through
all seven outputs. &#10003; = kept, &#10007; = excluded, &#9679; = always
present (the checkpoint is pre-filter by construction, in both preprint and
final versions). A single symbol means preprint and final agree; where they
differ, the cell reads **preprint &rarr; final**.

| # | Category | Real example | `splice_variant` | Checkpoint | Data 5 | Reclass. | Data 6 | Fig 4 | Fig 5&6 | Ext. Data |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | Canonical splice site, splice-aware assay | `BAP1_Waters_2024` c.37+1G>A, `splice_donor_variant`, DS_DG=0.99, nt | Yes (score + consequence) | &#9679; | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | **&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** |
| 2 | Canonical splice site, non-splice-aware assay | `ASPA_Grønbæk-Thygesen_2024` c.235_237delinsTAA, `splice_acceptor/donor_variant`, aa | Yes (consequence only) | &#9679; | &#10007; | &#10007; | &#10007; | &#10007; | &#10007; | &#10007; |
| 3 | Cryptic score hit (missense), splice-aware, aa-res. | `PALB2_Boonen_2026_SGE` c.2998G>A, p.Gly1000Cys, `missense_variant`, DS_AG=0.50, aa | Yes (score only) | &#9679; | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | **&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** |
| 4 | Cryptic score hit (missense), non-splice-aware, nt-res. | `BRCA2_Hu_2024` c.7712A>T, p.Glu2571Val, `missense_variant`, DS_AG=0.95, nt | Yes (score only) | &#9679; | &#10007; | &#10007; | &#10007; | &#10007; | &#10007; | &#10007; |
| 5 | Low-impact `splicing_variant` label, NaN score | `DDX3X_Radford_2023` c.679+3_679+4delinsTT, `splice_donor_region_variant`, DS=NaN, nt | No (neither leg fires) | &#9679; | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** |
| 6 | Low-impact `splicing_variant` label, scored but < 0.2 | `BAP1_Waters_2024` c.37+3G>A, `splice_donor_region_variant`, DS_DG=0.04, nt | No (neither leg fires) | &#9679; | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** | **&#10003;** |
| 7 | AA "guilt by association" -- sibling realization is Yes, this one isn't | `PALB2_Boonen_2026_SGE` c.2998_2999delinsTC, p.Gly1000Ser, `missense_variant`, DS=NaN (all four), aa | No (own evidence: none) | &#9679; | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** | **&#10003;** | &#10007;&rarr;**&#10003;** | &#10007;&rarr;**&#10003;** |

Only three of the seven categories actually change. Rows 1, 3, and 7 flip
from excluded-everywhere-but-Figure-4 (preprint) to kept everywhere (final),
now that every consumer checks `splice_measure`. Row 7 in particular
carries zero SpliceAI or consequence evidence of its own in either
version -- it's kept, in the final version, purely because `splice_var_amino`
broadcasts `Yes` across every nucleotide realization of `p.Gly1000Ser` once
one of them (`c.2998G>A`, row 3) scores &ge; 0.2, and
`PALB2_Boonen_2026_SGE`'s `splice_measure == 'Yes'` now spares that whole
broadcast group -- unchanged from the preprint in Figure 4, which never
depended on this broadcast to begin with (Section 6).

Rows 2 and 4 (non-splice-aware assays) and rows 5 and 6 (the
`"splicing_variant"`/NaN-score scope decision from Section 1) are identical
in preprint and final, as intended -- Figure 4 already agreed with the rest
of the pipeline on these four categories even before the fix, and still
does.

## Related docs

- [`docs/splice_variant_assay_consistency.md`](splice_variant_assay_consistency.md) -- the motivation for the `splice_measure` exception and its observed impact on each output.
- [`docs/annotate_simplified_consequence.md`](annotate_simplified_consequence.md) -- Step 17's consequence mapping (doesn't cover the splice threshold logic; see Section 2 above instead).
- [`docs/variant_annotation_pipeline.md`](variant_annotation_pipeline.md) -- full annotation pipeline architecture (same gap).
- [`docs/variant_classification.md`](variant_classification.md) -- the classification notebook's `VariantNotes` vocabulary and dedup-strategy parameters in full.

# Assay priority questions

## 1. VUS reclassification dedup doesn't use assay priority

### Summary

`src/lib/assay_priority.py`'s `ASSAY_PRIORITY_LIST` exists to resolve variants
scored by more than one assay for the same gene. It's only wired into the
`controls_aa`/`clingen_aa` calibration-diagnostic tables in
`notebooks/analysis/Variant_Classification_analysis.ipynb` and
`OddsPath_classifications.ipynb`. The actual VUS reclassification path — the
one that produces `Class_GeneSpecific_REVEL/AM/MP2`, i.e. the numbers that
matter — never consults the priority list, and its stand-in (a tag that
happens to correlate with which assay had the bigger functional-points
value) is fragile enough that it gets the wrong answer whenever the
colliding assays aren't the same type. See
[Effect on VUS classification](#effect-on-vus-classification) for exactly
when it does and doesn't track the strongest evidence.

### Where this lives

`Variant_Classification_analysis.ipynb`, cell 91:

```python
VUS_no_dup = (
    VUS
    .sort_values(by="VariantNotes", na_position="last")
    .drop_duplicates(subset=['Gene', 'hg38_start', 'ref_allele', 'alt_allele'], keep="first")
)
```

Compare to the `controls_aa`/`clingen_aa` path (cells 82–85, 107–108), which
does use the priority list:

```python
controls_aa["assay_priority"] = controls_aa["Dataset"].map(assay_priority_map).fillna(9999)
...
controls_aa_drop_REVEL_YP = (
    controls_aa[...]
    .sort_values("assay_priority")
    .drop_duplicates(subset=["Gene", "aa_pos", "aa_ref", "aa_alt", "Ref_seq_transcript_ID_stripped"], keep="first")
)
```

And to `catch_mis_2` (cell 88), which is the "keep whichever row actually has
evidence" version, correct in spirit but only ever applied to the
`controls`/`ClinGen` tables, never to `VUS`:

```python
def catch_mis_2(df, group_cols, points_col='Fxn_points'):
    """
    Handle duplicates by keeping the row with the highest functional points.
    """
    group_cols = ['Gene', 'Chrom', 'hg38_start', 'ref_allele', 'alt_allele']
    df_sorted = df.sort_values(by=points_col, ascending=False, na_position='last')
    cleaned = df_sorted.drop_duplicates(subset=group_cols, keep='first')
    return cleaned
```

`catch_mis_2` sorts on the signed `Fxn_points` value, not its absolute
value, even though its docstring and the pipeline README both describe it as
"highest [absolute] functional points." Concretely, `sort_values(ascending=False)`
on the signed value means: a positive `Fxn_points` always beats a negative
one no matter the magnitude (`+1` beats `-8`), and between two negatives it
keeps the one *closer to zero* — the weaker benign-direction evidence, not
the stronger one. Confirmed against real data: `BRCA1` chr17:43047691 T>C is
scored by `BRCA1_Findlay_2018` (`Fxn_points` −5) and
`BRCA1_Adamovich_2022_Cisplatin_Resistance` (`Fxn_points` −2) — `catch_mis_2`
keeps the Adamovich row, discarding the stronger (−5) benign signal for the
weaker (−2) one. See [category problems](#known-problems-by-category) below.

So there are three different dedup strategies for "same variant, multiple
assays" in this notebook, and the one applied to `VUS` is the one that
doesn't care about points or priority at all — it just keeps whichever row
sorts first alphabetically by `VariantNotes` (with nulls pushed last).

### `controls_aa`/`clingen_aa`: what "aa" means here, and who else uses it

"aa" is the amino-acid/protein-level half of `nucleotide_or_aa` — assays
that report a variant's effect via a protein-coding change (`aa_pos`,
`aa_ref`, `aa_alt`) rather than genomic coordinates directly. `controls_aa`
and `clingen_aa` are the aa-assay subset of the `controls` (cell 76: ClinVar
P/LP and B/LB) and `ClinGen_Repo` (cell 104) frames, split off from their
nt counterparts at cells 79 and 105.

`controls` and `ClinGen_Repo` are the only two categories that get this
split-by-assay-type treatment, and each goes through three dedup passes, not
one — nt and aa are deduped separately, then the two survivors are merged
and deduped again:

| Stage | Cells (`controls` / `ClinGen_Repo`) | Input | Sort key | Uses `assay_priority`? |
|---|---|---|---|---|
| 1. nt subset dedup | 80 / 106 | `controls_nuc` / `clingen_nuc` | `VariantNotes` (same weak sort as `VUS`) | no |
| 2. aa subset dedup | 85 / 108 | `controls_aa` / `clingen_aa` | `assay_priority` | **yes** — the one place `ASSAY_PRIORITY_LIST` is used |
| 3. merge nt+aa survivors, dedup again | 89 / 111 | `pd.concat([stage-1 output, stage-2 output])` | `catch_mis_2`, signed `Fxn_points` | no |

Stage 3 exists because stages 1 and 2 each only guarantee one row *per
assay type*. If the same physical variant was scored by both an nt assay and
an aa assay, both survive stages 1–2 and collide again on genomic
coordinates — so stage 3 dedups the combined nt+aa result a second time.
This is the same nt-vs-aa collision `VUS` has (see below); `controls`/
`ClinGen_Repo` at least make an explicit (if inconsistent — signed, not
`assay_priority` or absolute value) choice about how to resolve it. `VUS`
does not.

Stage 3's output feeds directly into `dfs` (cell 112) as
`controls_REVEL/AM/MP2_GeneSpecific` and `ClinGen_Repo_REVEL/AM/MP2_GeneSpecific`
— those Supplementary Data 5 tabs.

`VUS`, `gnomAD`, and `Unobserved` (`Unseen` in the notebook) skip stages 1–3
entirely — nt and aa rows are deduped together, in one pass, on genomic
coordinates only, via a single `VariantNotes` sort (cells 91, 102, 96). None
of the three ever look at `assay_priority`. Their outputs map straight to
the `VUS_*`, `gnomAD_*`, and `Unobserved_*` tabs.

| Category | Split by assay type? | Dedup stages | Uses `assay_priority`? | Supp Data 5 tabs |
|---|---|---|---|---|
| `controls` | yes | 3 (nt / aa / merge, table above) | yes (stage 2 only) | `controls_REVEL/AM/MP2_GeneSpecific` |
| `ClinGen_Repo` | yes | 3 (nt / aa / merge, table above) | yes (stage 2 only) | `ClinGen_Repo_REVEL/AM/MP2_GeneSpecific` |
| `VUS` | no | 1, `VariantNotes` over all rows | no | `VUS_REVEL/AM/MP2` |
| `gnomAD` | no | 1, `VariantNotes` over all rows | no | `gnomAD_REVEL/AM/MP2` |
| `Unobserved` | no | 1, `VariantNotes` over all rows | no | `Unobserved_REVEL/AM/MP2` |

None of this three-vs-one-stage difference is explained anywhere. The
pipeline README's "Deduplication Logic" section describes one
conceptual strategy ("select variants with the highest absolute functional
points") that matches neither `controls`/`ClinGen_Repo`'s actual three
sort keys (only one of which even looks at points, and that one uses the
signed value) nor `VUS`/`gnomAD`/`Unobserved`'s single alphabetical sort.

#### Known problems by category

Each category's dedup has a different, specific way of not doing "keep the
strongest evidence." None of them is a uniform coin-flip — each one is
wrong in one particular, describable direction.

| Category | Problem | Confirmed impact |
|---|---|---|
| `controls` / `ClinGen_Repo` | (1) `assay_priority` only covers the aa half | 170/1,147 nt-type multi-assay `controls` groups pick a different `Dataset` than `assay_priority` would; 0 REVEL/MP2, 1 AM classification flip |
| `controls` / `ClinGen_Repo` | (2) stage 3's sign bias | discards stronger evidence whenever the two survivors disagree in sign, or agree on a negative sign |
| `VUS` | systematic nt-over-aa bias | 921 nt/aa collisions, 46/84/68 REVEL/AM/MP2 flips |
| `gnomAD`, `Unobserved` | same nt-over-aa bias (same code as `VUS`) | not quantified |
| all five | string-sort fragility (tag rename/addition/tie) | none currently — see caveat below |

**(1) `assay_priority` only covers the aa half.** `ASSAY_PRIORITY_LIST` isn't
an aa-only list — it ranks nt assays too, and its own docstring flags
"BRCA2_Hu_2024 was and remains prioritized over BRCA2_IGVF and the Sahu
assays" as an open question. But stage 1 (nt, cells 80/106) never consults
it — only stage 2 (aa) does. Confirmed live: 981 `controls`-category BRCA2
variants are scored by more than one nt-level assay (`BRCA2_Hu_2024`, the
`Sahu_2023_exon13_*` family, `BRCA2_IGVF`, `BRCA2_Huang_2025_SGE`), and on
170 of the 1,147 nt-type multi-assay `controls` groups checked, the
`VariantNotes`-tag pick disagrees with what `assay_priority` would have
picked — e.g. `BRCA2_Sahu_2023_exon13_Cisplatin_Resistance` (priority 50)
beats `BRCA2_Sahu_2023_exon13_SGE` (priority 48) purely on `Fxn_points`
magnitude. This rarely changes classification today (0
`Class_GeneSpecific_REVEL`/MP2 flips, 1 AM flip, in this check), because the
magnitude-based pick and the curated ranking mostly agree by coincidence —
but the curated ranking isn't actually driving this decision, and BRCA2's
nt-assay count keeps growing. (See [section 2](#2-open-questions-on-assay-priority-order-for-oddspath-calibration)
below for the open question on whether `BRCA2_Hu_2024` should even still be
ranked #1.)

**(2) Stage 3's sign bias.** `catch_mis_2` sorts signed `Fxn_points`
descending, so a positive value always beats a negative one regardless of
magnitude (`+1` beats `-8`), and between two negatives it keeps the one
closer to zero. Confirmed with real data: `BRCA1` chr17:43047691 T>C keeps
`Adamovich_2022_Cisplatin_Resistance` (−2) over `Findlay_2018` (−5),
discarding the stronger benign-direction evidence. Unlike (1), this isn't
an nt-vs-aa bias — it can go either direction depending on which survivor
happens to be less negative.

**`VUS`: systematic nt-over-aa bias.** Single `VariantNotes` sort across
nt+aa rows together: `First_max_fxn_pts` (nt) always sorts ahead of
`max_fxn_pts` (aa), so whenever a variant is scored by both an nt and an aa
assay, the nt row wins regardless of which one has the bigger `Fxn_points`.
Confirmed: 921 VUS have this exact collision, and 46/84/68 get a different
REVEL/AM/MP2 classification than they would if the stronger row won. Within
a single assay type this doesn't happen — see
[Effect on VUS classification](#effect-on-vus-classification).

**`gnomAD`, `Unobserved`: same bias, unquantified.** Same code pattern as
`VUS` (cells 102, 96: single `VariantNotes` sort over nt+aa rows together),
so the same nt-over-aa bias applies whenever one of their variants is
scored by both an nt and an aa assay. Not separately checked here — nobody's
counted how many `gnomAD`/`Unobserved` variants actually hit that collision.

**All five: string-sort fragility.** See
[Alphabetical, not evidence-based](#alphabetical-not-evidence-based):
renaming a tag, adding a new one, or hitting an exact tie changes dedup
behavior with nothing to catch it. Doesn't currently misclassify anything
on its own — a real tie in `Fxn_points` means identical evidence either
way — but nothing protects against it going forward.

### `VariantNotes` values

Not an input column — nothing in `data/input/` or `Supplementary_Data_3/4.xlsx`
has it. Built up inside `Variant_Classification_analysis.ipynb`, on the full
per-assay-row dataset, before it's split into `controls`/`ClinGen_Repo`/
`VUS`/`gnomAD`/`Unobserved`:

| Cell(s) | Value | Meaning |
|---|---|---|
| 43 | `splice_variant_not_measured` | Splice variant with no functional measurement |
| 44 | `start_lost_variant_not_measured` | Start-lost aa variant with no functional measurement |
| 46 / 51 | `conflicting_fxn_data` | This coordinate group has rows with opposite-sign `Fxn_points` across assays |
| 47 | `First_max_fxn_pts` | Row has the max absolute `Fxn_points` in its nt-level coordinate group |
| 52 | `max_fxn_pts` | Same, aa-level coordinate group |
| — | `''` / `NaN` | None of the above — not flagged, not its group's max |
| 46/51 | e.g. `splice_variant_not_measured;conflicting_fxn_data` | Multiple conditions, joined with `;` |

These are written to
`data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`
at cell 67 and reloaded at cell 68. Cell 74 then drops every row tagged
`conflicting_fxn_data`, `splice_variant_not_measured`, or
`start_lost_variant_not_measured` from the shared `sankey_f` — before any
category is split off, so `controls`, `ClinGen_Repo`, `VUS`, `gnomAD`, and
`Unobserved` all inherit the same cleanup. By the time `VUS` is filtered out
at cell 90, the only `VariantNotes` values left are `First_max_fxn_pts`,
`max_fxn_pts`, or blank/`NaN`.

### Effect on VUS classification

Whichever row cell 91 keeps determines which assay's `Fxn_points` and
calibration feed `Total_Points_GeneSpecific_*`, and from there
`Class_GeneSpecific_REVEL/AM/MP2`. So does it pick the row with the
strongest evidence? Mostly, yes — and the failure mode is specific, not
random.

Checked against `integrated_variant_effect_dataset_analysis.csv.gz`,
reproducing cells 69/74/75 (drop SFPQ, drop conflicting/splice/start-lost,
drop `Flag='*'`) and grouping exactly as cell 91 does (`Gene`, `hg38_start`,
`ref_allele`, `alt_allele`, nt and aa rows together, ties broken by original
row order): 4,573 VUS coordinate-groups are scored by more than one assay —
3,652 by assays of the same type (all-nt or all-aa), 921 by a mix of both.

Comparing cell 91's pick to the row with the highest absolute `Fxn_points`
(what the `max_fxn_pts`/`First_max_fxn_pts` tags are themselves supposed to
mark): they disagree on 149 groups, and every one of the 149 is a
mixed-type group — zero disagreements among the 3,652 same-type groups.
Within a single assay type, the alphabetical sort works, because exactly
one row per group is tagged and that tag always sorts ahead of blank. It's
the mixed-type case that breaks it, and it breaks it the same way every
time: `First_max_fxn_pts` (nt) starts with a capital `F`, `max_fxn_pts` (aa)
starts with a lowercase `m`, and `F` (70) sorts before `m` (109) — so
whenever an nt assay and an aa assay both cover the same variant, **the nt
row wins regardless of which one actually has the larger `Fxn_points`.**
Checked directly: all 46 REVEL-classification flips below have `nucleotide_or_aa`
mismatched between the current pick and the abs-max pick, and the current
pick is the nt row in all 46.

Example: `BRCA1`, chr17:43057129 A>G is scored by `BRCA1_Findlay_2018` (nt,
`Fxn_points` 5) and `BRCA1_Adamovich_2022_HDR` (aa, `Fxn_points` 8, the
larger value). Cell 91 keeps `Findlay_2018` — `Likely Pathogenic` — purely
because `First_max_fxn_pts` outsorts `max_fxn_pts`; `Adamovich_2022_HDR`'s
`Fxn_points` of 8 would classify as `Pathogenic`.

Comparing cell 91's pick to the other two candidate sort keys:

| Alternative | Picks a different row | `Class_GeneSpecific_REVEL` flips | AM flips | MP2 flips |
|---|---|---|---|---|
| abs-max `Fxn_points` | 149 / 4,573 (all mixed-type) | 46 | 84 | 68 |
| `assay_priority` sort | 1,179 / 4,573 | 308 | 336 | 292 |
| signed max `Fxn_points` (`catch_mis_2` as coded) | 2,264 / 4,573 | 1,204 | 738 | 889 |

These two disagree far more than abs-max does, but that's a weaker claim
than it sounds: `assay_priority` and signed-`Fxn_points` are answering a
different question (which assay should we trust / which direction is the
evidence) than "which row has the biggest effect," so disagreement with
either one isn't itself evidence of a bug — the nt-beats-aa pattern above
is.

### Alphabetical, not evidence-based

`sort_values(by="VariantNotes", na_position="last")` sorts the values above
as strings. Today that order happens to be:

```
First_max_fxn_pts   (F = 70)
max_fxn_pts          (m = 109)
NaN                  (pushed last by na_position='last')
```

(`conflicting_fxn_data` and the splice/start-lost tags no longer appear by
the time `VUS` is filtered — see above — but the same sort is reused
verbatim for `gnomAD`/`Unobserved`/`ClinGen_Repo`'s nt subset, cells
96/102/106, where the same reasoning applies.)

This order has nothing to do with evidence strength; it's alphabetical.
It only avoids obviously wrong picks today because the two tags that matter
(`First_max_fxn_pts`, `max_fxn_pts`) happen to sort ahead of blank. Any of
the following would change dedup behavior with no warning and no test to
catch it:

- Renaming a tag (`max_fxn_pts` → `max_fxnpts` moves it past
  `splice_variant_not_measured` alphabetically, if that tag were ever
  reintroduced downstream of the cell-74 filter).
- Adding a new tag that happens to sort before `First_max_fxn_pts`/
  `max_fxn_pts`.
- Two rows in the same group sharing the same tag (both blank, or both
  `max_fxn_pts` on a numeric tie): `sort_values` defaults to `kind="quicksort"`,
  which isn't stable, so which row `drop_duplicates(keep="first")` keeps
  isn't even reliably "whichever came first upstream" — it's whatever
  numpy's quicksort does with equal keys. Doesn't currently change any
  classification (a real tie in `Fxn_points` means identical evidence, so
  either pick lands in the same points bucket), but it means the choice
  between tied rows isn't reproducible reasoning, just an implementation
  detail.

If a tag-based sort is kept at all, replace the tag strings with an explicit
numeric rank (tag → priority code), the same way `assay_priority` is a
number rather than a sort over `Dataset` name strings. That makes the
intended order something you read off a mapping, not something that falls
out of spelling and capitalization.

### Implications

`VUS_no_dup`'s sort isn't blind to evidence in general — within a single
assay type it agrees with "keep the biggest `Fxn_points`" every time, in
this dataset. Where it breaks is specifically the case `ASSAY_PRIORITY_LIST`
exists to handle: a variant scored by both an nt assay and an aa assay for
the same gene. There, `VUS_no_dup` doesn't pick the stronger of the two —
it picks the nt one, always, because `First_max_fxn_pts` sorts ahead of
`max_fxn_pts`. If the nt assay happens to be uncalibratable or weak and the
aa assay has the real signal, the aa evidence is discarded regardless.

### Current impact

An earlier version of this doc checked only SCN5A: none of its 7 ClinVar
VUS are scored by both `SCN5A_Glazer_2020` and `SCN5A_Ma_2024` (6 fall in
Glazer, 1 in Ma), so the gap wasn't suppressing evidence for that gene pair
— though 4 non-VUS SCN5A variants confirmed the cross-assay collision case
does happen for that gene in general.

The dataset-wide check in [Effect on VUS classification](#effect-on-vus-classification)
replaces that single-gene spot check: 921 VUS are scored by both an nt and
an aa assay, cell 91 keeps the nt row in every one of those groups by
construction, and 46/84/68 of them get a different REVEL/AM/MP2
classification than they would if the stronger-evidence row were kept
instead. That's happening today, not just in some future collision case —
it's just concentrated in nt/aa collisions rather than spread evenly across
all multi-assay VUS.

### Suggested fix

Replace the `VariantNotes`-sort dedup in cell 91 (and the identical pattern
at cells 96 and 102) with either:

- `catch_mis_2`, fixed to sort on `Fxn_points.abs()` rather than the signed
  value, applied the same way `controls`/`ClinGen_Repo` apply it at their
  stage-3 merge — this directly fixes the nt-beats-aa problem, since it
  compares actual point values instead of tag spelling, or
- an explicit `assay_priority` sort, consistent with `controls_aa`/`clingen_aa`'s
  stage 2 — but note from the table above that `assay_priority` disagrees
  with abs-max `Fxn_points` more often than the current sort does, so this
  changes more than just the nt/aa cases; it's a different policy choice
  (trust the assay ranking over trusting whichever assay happened to score
  higher), not strictly a superset fix.

Either way, stop sorting on `VariantNotes` strings for this purpose, and
consider giving `controls`/`ClinGen_Repo`'s nt stage (cells 80/106) and
their stage-3 merge (cells 89/111) the same fix — they inherit the identical
tag-sort problem for nt-subset ties, and stage 3's signed-value sort has the
same "signed, not absolute" inconsistency as `catch_mis_2` itself.

## 2. Open questions on assay priority order for OddsPath calibration

### What this order is for

`ASSAY_PRIORITY_LIST` is consulted in exactly one place in the pipeline
today: the amino-acid-level dedup for `controls` and `ClinGen_Repo` (stage 2
in [section 1](#controls_aaclingen_aa-what-aa-means-here-and-who-else-uses-it)).
When an aa-type variant there has functional scores from more than one MAVE
assay for the same gene, this list is the tie-breaker — the score from
whichever assay ranks highest is kept, the rest are dropped.

It is *not* consulted for `VUS`, `gnomAD`, or `Unobserved` at all, and not
even for `controls`/`ClinGen_Repo`'s own nucleotide-level dedup — despite
ranking nt-type assays too (see
[Known problems by category](#known-problems-by-category)). So whether a
given gene's ordering below is actually deciding anything in production
right now depends on whether its overlapping assays are aa-type or nt-type:
`LDLR`, `HMBS`, and `PAX6` are all aa-type, so their ordering is live.
`BRCA2` and `CBS`'s overlapping assays are nt-type, so their ordering
currently has no effect at all until the nt-side gap above is fixed — for
those two, question 1 and part of question 3 below are "which order should
this become," not "which order is currently deciding anything."

The current order reflects what was used in the preprint submission
(section 1 of the list below). It isn't specified for all genes, and new
datasets have been added since submission.

### Open questions

1. **BRCA2** (nt-type — not currently applied; see above). `BRCA2_Hu_2024`
   was ranked #1 for BRCA2 in the original submission and remains #1 in the
   current list, ahead of `BRCA2_IGVF` and all the Sahu assays
   (`BRCA2_Sahu_2025_SGE`, the four `BRCA2_Sahu_2023_exon13_*` sets, and
   `BRCA2_Huang_2025_SGE`). That ordering was kept by default since it was
   never revisited. Should `BRCA2_Hu_2024` still win, or should `BRCA2_IGVF`
   take priority? `BRCA2_IGVF` has ExCALIBR evidence but no functional
   classes/OddsPath, so the current order may still be correct.

2. **LDLR** (aa-type — live today). Three new assays
   (`LDLR_Tabet_2025_uptake`, `LDLR_Tabet_2025_abundance`,
   `LDLR_Tabet_2025_presence_VLDL`) have no established priority order and
   are currently listed in that order as a placeholder. Which one should be
   preferred?

3. **CBS** (nt-type — not currently applied), **HMBS, PAX6** (aa-type —
   live today). None of these three genes' assays were ever prioritized,
   and they're currently listed in an arbitrary order just to keep the
   pipeline's behavior well-defined. Better orderings for any of them are
   welcome — for HMBS and PAX6 that ordering is already deciding output
   today; for CBS it isn't yet, for the same reason as BRCA2.

Feedback on any other gene's ordering is also welcome.

### Current order (91 assays)

1 = highest priority.

**Section 1 — the original preprint submission's priority order:**

1. `BRCA1_Findlay_2018`
2. `BRCA2_Hu_2024`
3. `VHL_Buckley_2024`
4. `JAG1_Gilbert_2024`
5. `BARD1_IGVF`
6. `PALB2_IGVF`
7. `SFPQ_IGVF`
8. `RAD51D_IGVF`
9. `CTCF_IGVF`
10. `BAP1_Waters_2024`
11. `DDX3X_Radford_2023`
12. `RHO_Wan_2019`
13. `RAD51C_Olvera-León_2024`
14. `FKRP_Ma_2024`
15. `LARGE1_Ma_2024`
16. `CARD11_Meitlis_2020_SGE_Ibrutinib_GoF`
17. `CARD11_Meitlis_2020_SGE_LoF`
18. `ASPA_Grønbæk-Thygesen_2024_abundance`
19. `ASPA_Grønbæk-Thygesen_2024_toxicity`
20. `BRCA1_Adamovich_2022_Cisplatin_Resistance`
21. `BRCA1_Adamovich_2022_HDR`
22. `CHEK2_Gebbia_2024`
23. `CRX_Shepherdson_2024`
24. `F9_Popp_2025_model`
25. `G6PD_IGVF`
26. `GCK_Gersing_2023_complementation`
27. `GCK_Gersing_2024_abundance`
28. `KCNE1_Muhammad_2024_trafficking`
29. `KCNE1_Muhammad_2024_potassium_flux`
30. `KCNE1_Muhammad_2024_trafficking_WT_background_DN`
31. `KCNH2_Jiang_2022`
32. `KCNH2_Kozek_Glazer_2020`
33. `KCNH2_O_Neill_2024_surface_expression`
34. `MSH2_Jia_2021`
35. `NDUFAF6_Sung_2024`
36. `OTC_Lo_2023`
37. `PTEN_Matreyek_2018`
38. `PTEN_Mighell_2018`
39. `SCN5A_Glazer_2020`
40. `SCN5A_Ma_2024`
41. `SGCB_Li_2023`
42. `TP53_Fayer_2021_meta`
43. `TP53_Fortuno_2021`
44. `TSC2_IGVF`
45. `KCNQ4_Zheng_2022_current_homozygous`
46. `KCNQ4_Zheng_2022_v12_homozygous`

**Section 2 — assays that existed at the time of submission but weren't
explicitly ranked (best-guess order, unreviewed; added to make the output
well-defined):**

47. `BRCA2_IGVF`
48. `BRCA2_Sahu_2025_SGE`
49. `BRCA2_Sahu_2023_exon13_SGE`
50. `BRCA2_Sahu_2023_exon13_global_score`
51. `BRCA2_Sahu_2023_exon13_Cisplatin_Resistance`
52. `BRCA2_Sahu_2023_exon13_Olaparib_Resistance`
53. `CALM1_CALM2_CALM3_Weile_2017`
54. `CBS_Sun_2020_high_B6`
55. `CBS_Sun_2020_low_B6`
56. `F9_Popp_2025_heavy_chain`
57. `F9_Popp_2025_light_chain`
58. `F9_Popp_2025_carboxy_F9_specific`
59. `F9_Popp_2025_carboxy_gla_motif`
60. `F9_Popp_2025_strep_2`
61. `HMBS_van_Loggerenberg_2023_combined`
62. `HMBS_van_Loggerenberg_2023_erythroid`
63. `HMBS_van_Loggerenberg_2023_ubquitous`
64. `PAX6_McDonnell_2024_BLX_geneticin`
65. `PAX6_McDonnell_2024_BLX_no_geneticin`
66. `PAX6_McDonnell_2024_LE9_geneticin`
67. `PAX6_McDonnell_2024_LE9_no_geneticin`
68. `TARDBP_Bolognesi_Faure_2019`
69. `TP53_Kato_2003_AIP1nWT`
70. `TP53_Kato_2003_BAXnWT`
71. `TP53_Kato_2003_GADD45nWT`
72. `TP53_Kato_2003_h1433snWT`
73. `TP53_Kato_2003_MDM2nWT`
74. `TP53_Kato_2003_NOXAnWT`
75. `TP53_Kato_2003_P53R2nWT`
76. `TP53_Kato_2003_WAF1nWT`
77. `TP53_Giacomelli_2018_combined_score`
78. `TP53_Giacomelli_2018_p53WT_Nutlin3`
79. `TP53_Giacomelli_2018_p53null_Nutlin3`
80. `TP53_Giacomelli_2018_p53null_etoposide`
81. `TP53_Boettcher_2019`
82. `TPK1_Weile_2017`
83. `XRCC2_IGVF`

**Section 3 — assays added after submission (also a best guess,
unreviewed):**

84. `BRCA2_Huang_2025_SGE`
85. `CHEK2_McCarthy-Leo_2024`
86. `LDLR_Tabet_2025_uptake`
87. `LDLR_Tabet_2025_abundance`
88. `LDLR_Tabet_2025_presence_VLDL`
89. `PALB2_Boonen_2026`
90. `PALB2_Boonen_2026_SGE`
91. `TP53_Funk_2025`

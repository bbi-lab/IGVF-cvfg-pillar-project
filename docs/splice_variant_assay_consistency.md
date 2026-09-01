# Splice variant inclusion for splice-aware assays

Supplementary Data 5 and 6, the reclassification export, and the OddsPath
calibration now include a splice-affecting variant whenever its source assay
is curated as able to detect splicing effects, instead of excluding every
splice-affecting variant unconditionally.

## Motivation

The pipeline already includes any variant a functional assay actually
measures a consequence for -- synonymous, nonsense, frameshift, and every
other non-missense category pass through into classification and
calibration as long as the originating assay produced a functional score
for them. Splice-affecting variants were the one exception: they were
excluded from every downstream output regardless of whether the assay could
detect splicing disruption at all.

For SGE-based assays in particular, this was inconsistent. These assays
edit genomic DNA across intron/exon boundaries and read out a direct
functional consequence for a splice-disrupting variant the same way they do
for a missense or nonsense variant elsewhere in the gene. `Supplementary_
Data_3.xlsx` already records, per dataset, whether the assay is designed to
detect splicing effects ("Detects Splicing Variants?"). There is no
principled reason to treat splice-affecting variants differently from any
other consequence category once that's known -- so a variant should be
excluded for lacking assay coverage, not merely for being splice-related.

This also brings the rest of the pipeline in line with Figure 4's data
preparation (`figure_4/data_utils.py: Scoreset.splicing_filter`), which
already implemented exactly this exception; Supplementary Data 5/6, the
reclassification export, and the OddsPath notebooks simply hadn't been
updated to match.

## What changed

A variant flagged splice-affecting (`splice_var_amino == 'Yes'`) is now
excluded only when its dataset is *not* curated as splice-aware
(`splice_measure != 'Yes'`) -- mirroring how every other variant category is
already handled. This exception was applied everywhere a splice-affecting
variant was previously dropped unconditionally:

- `Variant_Classification_analysis.ipynb` -- the `VariantNotes` tagging
  step and the two filters feeding Supplementary Data 5.
- `OddsPath_calculations.ipynb` and `OddsPath_classifications.ipynb` -- the
  corresponding tagging step and filters feeding Supplementary Data 6.
- `src/build_variant_reclassification_dataset.py` -- the same filter,
  independently applied to the reclassification export.

See [`docs/splice_variant_filtering_pipeline.md`](splice_variant_filtering_pipeline.md) for the full,
cell-by-cell trace of how splice variants are annotated and filtered
end to end.

## Observed impact

Compared against a pipeline run from before this change (`data/output.
original-splice-variant-filter/`), across the outputs this change touches:

| Output | Distinct variants added | Distinct variants removed |
|---|---|---|
| Supplementary Data 5 | +3,196 | 0 |
| Supplementary Data 6 | +3,321 | 71 (see note below) |
| `integrated_variant_effect_reclassification.tsv.gz` | +4,551 | 0 |
| `integrated_variant_effect_dataset_analysis.csv.gz` (checkpoint) | 2,669 rows gained a corrected splice-related annotation tag; row count unchanged (975,049 both) -- this file is written before any exclusion filter runs, so it was never missing rows in the first place | -- |

Every added variant comes from a dataset curated `splice_measure == 'Yes'`
in `Supplementary_Data_3.xlsx` -- `BAP1_Waters_2024`, `DDX3X_Radford_2023`,
`BARD1_IGVF`, `BRCA2_Huang_2025_SGE`, `RAD51C_Olvera-León_2024`,
`BRCA1_Findlay_2018`, `TP53_Funk_2025`, `VHL_Buckley_2024`, and several
smaller SGE datasets.

**Supplementary Data 6's 71 "removed" rows are a reassignment, not a
coverage loss.** BRCA2 has five overlapping SGE studies covering the same
exon 13 region; where more than one independently measured the same
nucleotide variant, the pipeline's deduplication picks one representative
dataset per position. All 71 of these rows are cases where a different one
of those datasets is now the representative for the identical genomic
variant -- every one of the 71 positions is still fully present in
Supplementary Data 6, just attributed to a different qualifying dataset.
Checked directly against the checkpoint (which only reflects this change,
not that deduplication step): the annotation and priority-tag values at
these positions are identical before and after, confirming this
reassignment isn't caused by this change.

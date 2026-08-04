# MAVE Dataset Stats

`src/mave_dataset_stats.py` reports summary statistics for the integrated MAVE
variant effect dataset (`data/mave_data/integrated_variant_effect_dataset*.tsv.gz`),
split into three groupings of datasets:

- **Community (IGVF)** -- datasets Supplementary Data 3 marks `IGVF Produced? = Yes`
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

## Inputs

- **Condensed variant effect dataset**
  (default `data/mave_data/integrated_variant_effect_dataset.condensed.tsv.gz`):
  one row per variant effect measurement/score. For protein-resolution
  datasets, one row can correspond to several DNA-level changes (pipe-delimited
  within a cell) -- the *expanded* sibling file
  (`integrated_variant_effect_dataset.tsv.gz`) explodes those into one row per
  DNA variant, but this script intentionally works from the condensed file
  since the requested counts are measurement-level, not DNA-variant-level.
- **Dataset metadata** (default `data/mave_data/Supplementary_Data_3.xlsx`,
  `Curation` sheet): one row per dataset, keyed by `Dataset Name` (joined
  against the condensed file's `Dataset` column). Used for the
  `IGVF Produced?` and `Primary Score Set or Meta-analysis?` classifications.

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
  data/mave_data/integrated_variant_effect_dataset.condensed.tsv.gz \
  data/mave_data/Supplementary_Data_3.xlsx \
  [--output stats.csv]
```

Both positional arguments default to the paths above, so a bare invocation
works from the repo root.

Via Docker (same image as `flag_variants`, see `compose.yaml`):

```bash
src/scripts/run_mave_dataset_stats.sh [condensed-file] [metadata-file] [--output path]
```

Unlike `run_flag_variants.sh`, this wrapper doesn't map paths against a
`/work` staging mount -- the script only reads the repo's checked-in
`data/mave_data/` files, which are already available at `/usr/src/app` via the
whole-repo bind mount.

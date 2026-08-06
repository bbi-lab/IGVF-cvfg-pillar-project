# Annotate Simplified Consequence

`src/annotate_simplified_consequence.py` adds a `simplified_consequence`
column to a CVFG variants TSV, collapsing each DNA candidate's VEP
consequence term down to a single, coarser "SO summary term" (e.g.
`missense_variant` rather than `missense_variant^splice_region_variant`).

## Source logic

Ported from the `get_simplified_consequence` cell in
`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`, which
built a `VEP output term -> SO summary term` mapping from a curated
Ensembl-consequence table (`extended_ensembl_consequence.csv.gz`, extracted
from Ensembl's [predicted data
page](https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html)
on 2025-03-18) -- picking, for any VEP term listed more than once, the row
with the highest `Importance` -- and applied it to a single flat
`consequence` column (comma-separated VEP terms, one row per variant),
taking the first comma-separated term with a mapping-table entry.

The variant-annotation pipeline's `annotate_vep.py` (step 10) instead
produces a **pipe-delimited** `vep.most_severe_mutational_consequence`
column, aligned one segment per DNA candidate (the same convention
`flag_variants.py`'s `Flag` column and `recalculate_clingen_classification.py`
use), with each segment already reduced to VEP's own single most-severe term
for that candidate. This script maps each segment independently through the
same mapping table and re-joins the results with `|`, rather than picking
among multiple comma-separated terms itself the way the notebook did.

A term with no mapping-table entry (e.g. `no_change`, a
`variant-annotation`-specific value not in Ensembl's own vocabulary) maps to
an empty string, same as the notebook's `NaN` case for an unmapped term.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.annotate_simplified_consequence \
  data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
```

Via Docker (same image as `flag_variants`/`recalculate_clingen_classification`,
see `compose.yaml`):

```bash
src/scripts/run_annotate_simplified_consequence.sh \
  data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
```

Like `run_flag_variants.sh`, this wrapper maps its input/output paths against
the `/work` staging mount (`${VARIANT_DATA_DIR:-./data}`), since it reads
from the same staged pipeline data those files live in rather than a
committed repo file. `--consequence-map-file` is remapped the same way
`run_flag_variants.sh` remaps `--filtering-dir` -- its default resolves
against this project's own bind mount (`/usr/src/app`) instead, since
`extended_ensembl_consequence.csv.gz` is a committed repo file, not staged
pipeline data.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--consequence-col` | `vep.most_severe_mutational_consequence` | Column with pipe-delimited VEP consequence term(s) per DNA candidate |
| `--consequence-map-file` | `data/input/consequence/extended_ensembl_consequence.csv.gz` | CSV mapping `VEP output term` to `SO summary term`, with an `Importance` column for resolving duplicates |

Raises `click.ClickException` if `--consequence-col` isn't present in the
input file.

## Stdout report

After simplifying, the script prints per-candidate `simplified_consequence`
value counts (a candidate with an empty segment -- no VEP consequence, or an
unmapped term -- is reported as `(no consequence)`).

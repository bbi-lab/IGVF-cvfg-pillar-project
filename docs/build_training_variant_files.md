# Build Training Variant Files

`src/build_training_variant_files.py` regenerates
`revel_training_variants.tsv` and `mutpred2_training_variants.tsv` from
richer upstream training-variant sources, for use with
`annotate_predictors.py`'s `--revel-training-file` and
`--mutpred2-training-file` options (variant-annotation pipeline Step 12 --
REVEL and AlphaMissense; see `docs/variant_annotation_pipeline.md`). It is a
**preparatory step that runs before Step 12** -- `step_12` in
`scripts/variant_annotation_pipeline.sh` calls it automatically before
`run_annotate_predictors.sh`.

## Inputs and outputs

| | Path | Contents |
|---|---|---|
| Input | `data/input/predictors/clustering_variants_revel_training_overlap.csv.gz` | Genome-wide REVEL-training-overlap list with genomic coordinates and a 1-letter `protein_variant` column (e.g. `Q416K`). All rows are SNVs. |
| Input | `data/input/predictors/MutPred2_scores_1.csv.gz`, `data/input/predictors/MutPred2_scores_2.csv.gz` | Two independently-generated MutPred2-training exports (`Gene`/`AA`/`MP2_train` and `Gene`/`variant`/`is_mp2_train` respectively) that fully agree on every `(gene, variant)` pair they share; unioned for maximum coverage. |
| Output | `data/intermediate/variant_annotation/data/revel_training_variants.tsv` | REVEL training-variant list, keyed on hg38 genomic coordinates. |
| Output | `data/intermediate/variant_annotation/data/mutpred2_training_variants.tsv` | MutPred2 training-variant list, keyed on `(gene_symbol, unqualified_hgvs_p)`. |

The output directory is the gitignored staging directory
`scripts/run_variant_annotation_pipeline.sh` mounts as the
variant-annotation pipeline's `VARIANT_DATA_DIR` -- the script creates it if
missing, and its output lands in the same place `annotate_predictors.py`
already reads `data_frame_missense_variants_MP2_properties.csv.gz` from.

Both `.csv.gz` sources are read transparently (no need to decompress first).

## CALM1/CALM2/CALM3

CALM1, CALM2, and CALM3 encode an identical protein, so a MAVE assay of
"calmodulin" can't attribute a variant to one paralog over another. This
pipeline's `gene_symbol` column always records these variants as `CALM1`
(never a joint label), so MutPred2 training entries for CALM2/CALM3 are
aliased onto CALM1 -- otherwise they would never match the input's
`gene_symbol` and the training-set overlap would silently look empty for
every CALM1 variant. REVEL matching is genomic-coordinate-based (no gene
symbol involved), so no aliasing is needed there.

## Usage

Locally (with the Poetry environment):

```bash
poetry run python -m src.build_training_variant_files
```

Both directories default to the paths above, so a bare invocation works from
the repo root.

Via Docker:

```bash
src/scripts/run_build_training_variant_files.sh
```

Unlike `run_flag_variants.sh`/`run_recalculate_clingen_classification.sh`,
this wrapper doesn't map paths against a `/work` staging mount -- both the
input and output directories are fixed under this project's own tree, so
they're already reachable via the whole-repo bind mount (`/usr/src/app`),
the same as `run_mave_dataset_stats.sh`/`run_load_excalibr_calibrations.sh`.

## CLI options

| Option | Default | Description |
|---|---|---|
| `--input-dir` | `data/input/predictors` | Directory containing the upstream source files |
| `--output-dir` | `data/intermediate/variant_annotation/data` | Directory to write the training-variant TSVs into (created if missing) |
| `--revel-source` | `<input-dir>/clustering_variants_revel_training_overlap.csv.gz` | REVEL training-overlap source file |
| `--mutpred2-source-1` | `<input-dir>/MutPred2_scores_1.csv.gz` | First MutPred2-training export |
| `--mutpred2-source-2` | `<input-dir>/MutPred2_scores_2.csv.gz` | Second MutPred2-training export |
| `--revel-dest` | `<output-dir>/revel_training_variants.tsv` | REVEL training-variant output path |
| `--mutpred2-dest` | `<output-dir>/mutpred2_training_variants.tsv` | MutPred2 training-variant output path |

## Stdout report

Prints the number of rows written to each output file.

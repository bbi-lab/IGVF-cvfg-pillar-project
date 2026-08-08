# A scalable approach to resolving variants of uncertain significance

This repository contains the data-processing pipeline and figure-generation
code behind the paper *"A scalable approach to resolving variants of
uncertain significance."* It harmonizes variant-effect measurements from ~83
published multiplexed assays of variant effect (MAVE) datasets, maps them to
genomic coordinates, annotates them with ClinVar/gnomAD/predictor evidence,
runs the OddsPath/ACMG-AMP classification analysis described in the paper,
and generates the manuscript's main and extended data figures.

> **Preprint:** *A scalable approach to resolving variants of uncertain
> significance.* bioRxiv (2026).
> https://www.biorxiv.org/content/10.64898/2026.02.14.705848v2

The pipeline runs in three stages, each covered in its own section below:

1. [Data preparation and variant annotation](#1-data-preparation-and-variant-annotation)
2. [Data analysis / variant classification and table preparation](#2-data-analysis--variant-classification-and-table-preparation)
3. [Figure preparation](#3-figure-preparation)

For agent/contributor-facing conventions (repo layout, coding style, the
notebook-to-script migration in progress, roadmap) see [AGENTS.md](AGENTS.md)
— this README is the user-facing overview of what the pipeline does and how
to run it.

---

## 1. Data preparation and variant annotation

The first stage turns curated MAVE inputs into the integrated variant-effect
dataset (Supplementary Data 1) that everything downstream reads from. It used
to be a single notebook (`notebooks/analysis/Integrated_variant_effect_dataset_pipeline.ipynb`,
now archived); that notebook has been replaced by a Dockerized pipeline built
on the sibling **[variant-annotation](https://github.com/bbi-lab/variant-annotation)**
project.

### Relation to the variant-annotation pipeline

`variant-annotation` is a general-purpose, Dockerized variant-mapping and
annotation pipeline (reverse translation, ClinVar/gnomAD/SpliceAI/ClinGen
evidence repository/VEP/MaveDB/predictor annotation) that is developed and
versioned independently of this project. This repo depends on it rather than
reimplementing it:

- `vendor/variant-annotation` is a git submodule pinned to a known-good
  commit of `variant-annotation`.
- `scripts/variant_annotation_pipeline.sh` drives that pipeline's steps
  against this project's MAVE data, interleaving a handful of this project's
  own steps (`postprocess_mavedb_functional_classifications`,
  `add_mavedb_active_calibration_columns`, `derive_score_set_urn`,
  `annotate_simplified_consequence`, `recalculate_clingen_classification`,
  `flag_variants`) between `variant-annotation`'s numbered steps.
- `scripts/run_variant_annotation_pipeline.sh` is the orchestrator you
  actually run: it stages this project's inputs, points the
  `variant-annotation` checkout at them via `VARIANT_DATA_DIR`, runs the
  pipeline, and collects the final gzipped output.

Full architecture, the submodule/override convention, and a non-obvious
`VARIANT_DATA_DIR` path-mapping subtlety are documented in
[`docs/variant_annotation_pipeline.md`](docs/variant_annotation_pipeline.md)
— read that before touching `scripts/variant_annotation_pipeline.sh` or its
wrapper scripts.

### Environmental prerequisites

- **Docker and Docker Compose** — both `variant-annotation`'s steps and this
  project's own annotation steps run as Docker containers, not directly on
  the host.
- **A `variant-annotation` checkout** — either the vendored submodule
  (`git submodule update --init vendor/variant-annotation`) or your own
  existing checkout that already has the large reference caches downloaded
  (see below), pointed to via `export VARIANT_ANNOTATION_DIR=/path/to/checkout`.
- **Poetry** (see [Software requirements](#software-requirements) below) —
  needed to build this project's own Docker image (`Dockerfile`/`compose.yaml`)
  for the steps it contributes to the pipeline.
- **A gnomAD Hail table cache**, built once per `variant-annotation`
  checkout via `scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache`
  (can take 6-7 hours) — not needed if `VARIANT_ANNOTATION_DIR` already
  points at a checkout where this has been done.

### Data-file prerequisites

Committed to this repo (no action needed):

- `data/input/maves/cvfg_variants.0.tsv`, `data/input/maves/score_sets.tsv`,
  `data/input/maves/Supp_Data1_column_description.xlsx` — this project's own
  variant list, dataset-name mapping, and column documentation.
- `data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz` — RefSeq/Ensembl
  transcript mapping.
- `data/input/predictors/AlphaMissense_hg38.tsv.gz` (+ `.tbi`),
  `data/input/predictors/revel_hg38.tsv.gz` (+ `.tbi`),
  `data/input/predictors/data_frame_missense_variants_MP2_properties.csv.gz`
  — pinned predictor score files specific to this pipeline run.
- `data/input/mave_calibration/excalibr/json/` — per-dataset exCALIBR JSON
  calibrations (see [Section 2](#2-data-analysis--variant-classification-and-table-preparation)).

Not committed — you must populate or point at these yourself:

- `data/raw_mave_data/` — hand-populated MAVE study inputs; see its own
  `README.md` for exactly what belongs there before running the pipeline.
- A `variant-annotation` checkout with the large, generic reference caches
  already downloaded: SpliceAI VCFs, dbNSFP, `clinvar_cache/`, and the
  gnomAD Hail table cache above. These are tens of GB and are deliberately
  not duplicated into this repo.
- `IGVFFI3804AVJR.csv.gz`, downloaded separately from
  https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ — only needed later,
  for a subset of figures (see [Section 3](#3-figure-preparation)), not for
  this stage.

### Running it

```bash
git submodule update --init vendor/variant-annotation   # if using the vendored submodule
poetry install --all-extras

# One-time per variant-annotation checkout (skip if already done):
scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache

# Full run:
scripts/run_variant_annotation_pipeline.sh
```

This produces `data/output/maves/integrated_variant_effect_dataset.tsv.gz`
(expanded, one row per DNA-level variant) and its `.condensed.tsv.gz`
sibling (one row per measurement). Pass `--step N` to (re)run a single
numbered step for debugging — see
[`docs/variant_annotation_pipeline.md`](docs/variant_annotation_pipeline.md)
for the full step list and every per-step doc it links to.

---

## 2. Data analysis / variant classification and table preparation

The second stage takes the integrated variant-effect dataset from Stage 1
and produces the manuscript's classification tables: OddsPath calibrations,
ACMG/AMP evidence-point assignment, and the final per-gene, per-predictor
classification files (Supplementary Data 5).

This stage still lives in Jupyter notebooks under `notebooks/analysis/`
rather than `src/` scripts — converting them is a deliberately deferred item
on the [roadmap in AGENTS.md](AGENTS.md#roadmap). Each notebook has a
companion `README_*.md` in the same directory documenting its inputs,
methods, and outputs in detail — read those before running or modifying one.

Run in this order:

1. **`OddsPath_calculations.ipynb`** — OddsPath likelihood-ratio
   calculations per dataset, following Brnich et al. (2019).
2. **`OddsPath_classifications.ipynb`** — turns those OddsPath values into
   ACMG/AMP evidence points.
3. **`Variant_Classification_analysis.ipynb`** — combines functional
   evidence (via exCALIBR calibrations or OddsPath) with computational
   predictor evidence (REVEL, AlphaMissense, MutPred2; genome-wide and
   gene-specific), assigns final classifications, deduplicates variants
   seen in multiple assays, and exports the per-category classification
   files that become Supplementary Data 5. See
   [`notebooks/analysis/README_Variant_Classification_analysis.md`](notebooks/analysis/README_Variant_Classification_analysis.md).

Two `src/` scripts support this stage and are already converted:

- **`src/load_excalibr_calibrations.py`** — loads exCALIBR JSON calibrations
  (`data/input/mave_calibration/excalibr/json/`) into Supplementary Data 4
  so the workbook stays in sync with the latest calibration run. See
  [`docs/load_excalibr_calibrations.md`](docs/load_excalibr_calibrations.md).
- **`src/mave_dataset_stats.py`** — reports summary statistics (dataset,
  measurement, and variant counts; predictor score coverage; clinical
  attribute breakdowns; ExCALIBR calibration coverage; reclassification
  agreement) over the integrated dataset. See
  [`docs/mave_dataset_stats.md`](docs/mave_dataset_stats.md).

### Running it

```bash
poetry install --all-extras
poetry run python -m ipykernel install --user --name igvf-cvfg-pillar-project \
  --display-name "IGVF CVFG Pillar Project (Poetry)"

poetry run python -m src.load_excalibr_calibrations   # refresh Supplementary Data 4
poetry run jupyter lab                                # then run the three notebooks above, in order
```

Outputs land under `data/output/supplementary_data/` (`Supplementary_Data_4.xlsx`,
`Supplementary_Data_5.xlsx`) and `data/output/predictor_calibration/` (the
per-gene control files also used by Extended Data Figure 5).

---

## 3. Figure preparation

The third stage generates the manuscript's main and extended data figures
from Stage 2's outputs. Each figure directory under `Main_Figures/` and
`Extended_Data_Figures/` contains its own plotting notebook(s)/script(s) and
any small supporting inputs.

Two toolchains are in play, depending on the figure:

- **Python notebooks** (Altair/matplotlib), run with the Poetry environment
  and `jupyter nbconvert --execute` — no Docker needed.
- **R scripts/`.Rmd` files** (tidyverse, ggplot2 extensions like `ggsankey`,
  `patchwork`, `ggh4x`), run via the `r-figures` Docker Compose service
  (`Dockerfile.r`) rather than a local R install, since they depend on
  system fonts (Arial) and `cairo_pdf` output that are easiest to reproduce
  in a container.

A few figures additionally require `IGVFFI3804AVJR.csv.gz`, downloaded
separately from https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ (nothing
in this repo produces it) and placed in the figure's own directory.

As with Stage 2, notebooks here are still notebooks, not yet converted to
`src/` scripts (see [AGENTS.md](AGENTS.md#roadmap)).

### Running it

```bash
# Python figure notebooks
poetry install --all-extras
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --output executed_<notebook>.ipynb \
  Main_Figures/<figure_dir>/<notebook>.ipynb

# R figure scripts (one-time image build, then per script)
docker compose build r-figures
docker compose run --rm -w /usr/src/app/<figure_dir> r-figures <script>.R
docker compose run --rm -w /usr/src/app/<figure_dir> \
  r-figures -e 'rmarkdown::render("<script>.Rmd")'
```

See [`docs/figures.md`](docs/figures.md) for the exact commands, inputs, and
output paths for every individual figure and panel — including several
figure-specific gotchas (missing output directories, date-stamped filenames,
figures that only reconstruct from a cached intermediate).

---

## Data Availability

Large supporting datasets are hosted externally on Zenodo to comply with
GitHub file size limits.

**Zenodo record:** https://zenodo.org/records/18637474

The Zenodo archive includes:

- Supplementary Data
- Large intermediate files used for creation of the integrated variant
  effect dataset

Files included in this GitHub repository are:

- Analysis scripts and notebooks
- Figure-generation scripts and notebooks
- Small supporting inputs

## Software Requirements

- **Python 3.12** — dependencies are Poetry-managed (`pyproject.toml`); see
  [AGENTS.md](AGENTS.md#environment-poetry--ruff) for setup
  (`poetry install --all-extras`).
- **Docker and Docker Compose** — required for the variant-annotation
  pipeline (Stage 1) and the R figure scripts (Stage 3, `r-figures` service,
  R 4.4.2 via `rocker/r-ver`). A local R install is not required.
- Additional per-script package dependencies are declared in
  `pyproject.toml` (Python) or the top of each `.R`/`.Rmd` file (R).

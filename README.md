# A scalable approach to resolving variants of uncertain significance

This repository contains the data-processing pipeline and figure-generation
code behind the paper *"A scalable approach to resolving variants of
uncertain significance."* It harmonizes variant-effect measurements from 91
multiplexed assays of variant effect (MAVE) datasets, maps them to genomic
coordinates, annotates them with ClinVar/gnomAD/predictor evidence, runs the
OddsPath/ACMG-AMP classification analysis described in the paper, and
generates the manuscript's main and extended data figures.

> **Preprint:** *A scalable approach to resolving variants of uncertain
> significance.* bioRxiv (2026).
> https://www.biorxiv.org/content/10.64898/2026.02.14.705848v2

The pipeline runs in three stages, each covered in its own section below:

1. [Data preparation and variant annotation](#1-data-preparation-and-variant-annotation)
2. [Data analysis / variant classification and table preparation](#2-data-analysis--variant-classification-and-table-preparation)
3. [Figure preparation](#3-figure-preparation)

See [Roadmap](#roadmap) below for the current state of Dockerizing the
remaining notebook-based stages.

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

- A `variant-annotation` checkout with the large, generic reference caches
  already downloaded: SpliceAI VCFs, dbNSFP, `clinvar_cache/`, and the
  gnomAD Hail table cache above. These are tens of GB and are deliberately
  not duplicated into this repo.
- `IGVFFI3804AVJR.csv.gz`, downloaded separately from
  https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ and placed at
  `data/input/biobank/IGVFFI3804AVJR.csv.gz` — only needed later, for a
  subset of figures (see [Section 3](#3-figure-preparation)), not for this
  stage.

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

This stage still runs as Jupyter notebooks under `notebooks/analysis/`,
executed locally via Poetry rather than in Docker — Dockerizing them (as
notebooks, not converted to `src/` scripts) is the remaining item on the
[Roadmap](#roadmap) below. Each notebook has a companion `README_*.md` in
the same directory documenting its inputs, methods, and outputs in detail —
read those before running or modifying one.

Run in this order:

1. **`OddsPath_calculations.ipynb`** — OddsPath likelihood-ratio
   calculations per dataset, following Brnich et al. (2019).
2. **`src/load_oddspath_calibrations.py`** — loads that notebook's output
   (`data/output/mave_calibration/OddsPath_calibrations.csv.gz`) into the
   `OddsPath_calibrations` sheet of Supplementary Data 4, analogous to
   `src/load_excalibr_calibrations.py` below. See
   [`docs/load_oddspath_calibrations.md`](docs/load_oddspath_calibrations.md).
3. **`OddsPath_classifications.ipynb`** — turns those OddsPath values into
   ACMG/AMP evidence points.
4. **`Variant_Classification_analysis.ipynb`** — combines functional
   evidence (via exCALIBR calibrations or OddsPath) with computational
   predictor evidence (REVEL, AlphaMissense, MutPred2; genome-wide and
   gene-specific), assigns final classifications, deduplicates variants
   seen in multiple assays, and exports the per-category classification
   files that become Supplementary Data 5. See
   [`notebooks/analysis/README_Variant_Classification_analysis.md`](notebooks/analysis/README_Variant_Classification_analysis.md).
5. **`src/build_variant_reclassification_dataset.py`** — reads that
   notebook's own cell-69 checkpoint, re-applies its downstream exclusions
   (`SFPQ`, the CHEK2 QC flag, conflicting/unmeasured-splice
   `VariantNotes` tags), adds four ACMG/AMP points columns, and
   deduplicates to one row per DNA variant. Produces
   `data/output/reclassification/integrated_variant_effect_reclassification.tsv.gz`,
   which collaborators use for downstream biobank analysis.

Two more `src/` scripts support this stage (beyond `load_oddspath_calibrations.py`
above) and are already converted:

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

# Refresh Supplementary Data 4 from the latest exCALIBR calibrations
poetry run python -m src.load_excalibr_calibrations

# 1. OddsPath likelihood-ratio calculations
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_OddsPath_calculations.ipynb \
  notebooks/analysis/OddsPath_calculations.ipynb

# 2. Refresh Supplementary Data 4's OddsPath_calibrations sheet from that run
poetry run python -m src.load_oddspath_calibrations

# 3-4. Remaining notebooks above, in order
for nb in OddsPath_classifications Variant_Classification_analysis; do
  poetry run jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    notebooks/analysis/${nb}.ipynb
done

# 5. Build the biobank-analysis reclassification export from that notebook's checkpoint
poetry run python -m src.build_variant_reclassification_dataset
```

Each `nbconvert --execute` leaves a side-effect `executed_<name>.ipynb` next
to the original (nbconvert's copy with outputs attached) — gitignored, so no
need to delete them.

Outputs land under `data/output/supplementary_data/` (`Supplementary_Data_4.xlsx`,
`Supplementary_Data_5.xlsx`), `data/output/predictor_calibration/` (the
per-gene control files also used by Extended Data Figure 5), and
`data/output/reclassification/integrated_variant_effect_reclassification.tsv.gz`
(the biobank-analysis export from step 5 above).

---

## 3. Figure preparation

The third stage generates the manuscript's main and extended data figures
from Stage 2's outputs. Each figure directory under `notebooks/figures/`
contains its own plotting notebook(s)/script(s) and any small supporting
inputs.

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
in this repo produces it) and placed at
`data/input/biobank/IGVFFI3804AVJR.csv.gz`.

As with Stage 2, the Python notebooks here still run locally via Poetry
rather than in Docker (see [Roadmap](#roadmap) below) — the R
scripts/`.Rmd` files already do, via the `r-figures` service.

### Running it

One-time setup:

```bash
poetry install --all-extras
poetry run python -m ipykernel install --user --name igvf-cvfg-pillar-project \
  --display-name "IGVF CVFG Pillar Project (Poetry)"
docker compose build r-figures
```

#### Figure 2

```bash
mkdir -p data/output/figures/figure_2/Histogram_wStripplot

# 1. Prep notebook (run first; everything else in this figure depends on it)
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_PP_ProcessBigDataFrame.ipynb \
  notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb

# 2. Panel notebooks (independent of each other; all read step 1's output)
for nb in PP_ClinVarPrecisionRecall PP_Fig2_Heatmaps PP_ResolutionOverview PP_SeqFunctionMap PP_StackedHistograms; do
  poetry run jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    notebooks/figures/figure_2/${nb}.ipynb
done

# 3. Figure_2i.R -- reads data/input/biobank/IGVFFI3804AVJR.csv.gz (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_2 \
  r-figures Figure_2i.R
```

#### Figure 3

Not yet documented — see `notebooks/figures/figure_3/curation_summary_figure3.Rmd`
and its `Figure3a/c/d.csv.gz` inputs.

#### Figure 4

```bash
# 1. Rebuild figure4_data.json.gz, carrying forward fields that can't be
#    regenerated (see docs/build_figure4_data.md)
poetry run python -m src.build_figure4_data \
  --cached-json notebooks/figures/figure_4/old_figure4_data.json.gz

# 2. Execute the notebook to produce data/output/figures/figure_4.png
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_figure4.ipynb \
  notebooks/figures/figure_4/figure4.ipynb
```

#### Figure 5/6

```bash
# Figure_6b.R -- reads data/input/biobank/IGVFFI3804AVJR.csv.gz (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures Figure_6b.R

# Figure5_6.Rmd: writes executed_Figure5_6.html next to the .Rmd (gitignored)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures -e 'rmarkdown::render("Figure5_6.Rmd", output_file = "executed_Figure5_6.html")'
```

#### Extended Data Figure 2

```bash
# Extended_Data_Figure_2.R -- reads data/input/biobank/IGVFFI3804AVJR.csv.gz (see above)
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_2 \
  r-figures Extended_Data_Figure_2.R
```

#### Extended Data Figure 5

```bash
mkdir -p data/output/figures/extended_data_figure_5

poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_extended_data_figure_5.ipynb \
  notebooks/figures/extended_data_figure_5/Extended_Data_Figure_5.ipynb
```

#### Extended Data Figures 4-9 (`Extended_data_figures.Rmd`)

```bash
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_4_6_7_8_9 \
  r-figures -e 'rmarkdown::render("Extended_data_figures.Rmd")'
```

Every `nbconvert --execute` call above leaves a side-effect `executed_<name>.ipynb`
copy next to the source, and `Figure5_6.Rmd`'s `rmarkdown::render()` call
does the same with `executed_Figure5_6.html` — both gitignored, so no need
to delete them. `Extended_data_figures.Rmd` above is the one exception: it
still renders the plain, non-gitignored `Extended_data_figures.html` next to
itself — delete that manually if you don't want it in the working tree. See
[`docs/figures.md`](docs/figures.md) for output paths for every
individual panel and several figure-specific gotchas (date-stamped
filenames, figures that only reconstruct from a cached intermediate, output
directories that aren't created automatically).

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

- **Python 3.12** — dependencies are Poetry-managed (`pyproject.toml`):
  ```bash
  poetry env use python3.12   # one-time, only if `poetry env use` picks the wrong interpreter
  poetry install --all-extras
  ```
  This creates an in-project `.venv/` with the pipeline dependencies plus
  the `dev` (Ruff, nbstripout, pre-commit), `tests` (pytest), and
  `notebooks` (JupyterLab/ipykernel) extras.
- **Docker and Docker Compose** — required for the variant-annotation
  pipeline (Stage 1) and the R figure scripts (Stage 3, `r-figures` service,
  R 4.4.2 via `rocker/r-ver`). A local R install is not required.
- Additional per-script package dependencies are declared in
  `pyproject.toml` (Python) or the top of each `.R`/`.Rmd` file (R).

## Roadmap

Stage 1 (the `variant-annotation` pipeline) and Stage 3's R figures
(`r-figures` service) already run in Docker. The remaining piece is
Dockerizing Stage 2's three analysis notebooks
(`OddsPath_calculations.ipynb`, `OddsPath_classifications.ipynb`,
`Variant_Classification_analysis.ipynb`) and Stage 3's Python figure
notebooks, which still run locally via `poetry run jupyter nbconvert` — as
notebooks, not converted to `src/` scripts; the plan is to run them as-is
inside a container, the same way Stage 1 and the R figures already do.

# Figures

How to (re)generate each manuscript figure from `notebooks/figures/`.

## Figure 2

Directory: `notebooks/figures/figure_2/`. Six scripts: one prep notebook everything
else depends on, five independent Python panel notebooks, and one standalone
R script. The Python notebooks only need the Poetry environment (`pandas`,
`altair`, `matplotlib`, `vl-convert-python`) -- no Docker/R involved until
`Figure_2i.R`.

### 1. `PP_ProcessBigDataFrame.ipynb` (run this first)

Reads `data/output/maves/integrated_variant_effect_dataset.tsv.gz` (the
integrated MAVE dataset produced upstream of this figure) and splits out the
SGE genes (`BARD1, PALB2, BRCA2, RAD51D, XRCC2, CTCF, SFPQ`) and VAMP-seq
genes (`G6PD, TSC2, F9`; `F9` keeps only its heavy-chain-antibody dataset)
into three xlsx files under `data/intermediate/figures/figure_2/`:
`SGEsubset.xlsx`, `VAMPseqsubset_wDups.xlsx`, and `CAVAseqsubset.xlsx` (the
concatenation of the first two -- not currently read by anything else in
this directory). Every panel notebook below reads the first two of these.

```bash
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_PP_ProcessBigDataFrame.ipynb \
  notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb
```

### 2. Panel notebooks (independent of each other; all read step 1's output)

- `PP_ClinVarPrecisionRecall.ipynb` -- precision/recall of SGE and VAMP-seq
  against ClinVar, with Wilson-interval CIs. Saves
  `data/output/figures/figure_2/PillarProject_PRvsClinVar_wErrorBar_grey.svg`.
- `PP_Fig2_Heatmaps.ipynb` -- three charts in the figure's center column,
  despite its own docstring claiming "all" heatmaps: a RAD51D SGE
  amino-acid-position heatmap, a RAD51D SGE genomic-position map (a second,
  differently-rendered take on the same region `PP_SeqFunctionMap.ipynb`
  covers), and a G6PD VAMP-seq amino-acid-position heatmap. Each was only
  ever `.display()`ed inline -- no `.save()` call existed for any of them,
  so running the notebook wrote nothing to disk (hence its ~40 MB file size,
  all from embedded cell output). `.save()` calls have been added after each
  `.display()`, writing `RAD51D_sge_aa_heatmap.svg`,
  `RAD51D_sge_genomic_map.svg`, and `G6PD_vampseq_aa_heatmap.svg` to
  `data/output/figures/figure_2/`.
- `PP_ResolutionOverview.ipynb` -- the VAMP-seq vs. SGE genomic-position and
  amino-acid-change coverage bar chart. Saves
  `data/output/figures/figure_2/vampseq_sge_bars.svg`.
- `PP_SeqFunctionMap.ipynb` -- the RAD51D sequence-function map. Saves
  `data/output/figures/figure_2/RAD51D_X9_draft_SeqFunc_map_extended.svg`.
- `PP_StackedHistograms.ipynb` -- stacked score histograms (with a ClinVar
  density overlay) for both assays, plus per-gene SGE insets. Saves two SVGs
  under `data/output/figures/figure_2/Histogram_wStripplot/` and one
  `data/output/figures/figure_2/sge_histogram_inset_<gene>.svg` per SGE gene.

```bash
for nb in PP_ClinVarPrecisionRecall PP_Fig2_Heatmaps PP_ResolutionOverview PP_SeqFunctionMap PP_StackedHistograms; do
  poetry run jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    notebooks/figures/figure_2/${nb}.ipynb
done
```

`data/output/figures/figure_2/` (and its `Histogram_wStripplot/`
subdirectory) needs to already exist -- unlike the R scripts'
`save_my_plot()` helper, neither Altair's `Chart.save()` nor matplotlib's
`Figure.savefig()` create missing parent directories. On a fresh checkout,
run this first:

```bash
mkdir -p data/output/figures/figure_2/Histogram_wStripplot
```

Each `nbconvert --execute` above also leaves a side-effect
`executed_<name>.ipynb` next to the original (nbconvert's copy with outputs
attached, same as the Figure 4 step) -- gitignored, so no need to delete
them.

### 3. `Figure_2i.R` (separate from the notebook pipeline above)

Same pattern as `Figure_6b.R` / `Extended_Data_Figure_2.R`: a plain Rscript
run via the `r-figures` Docker service, reading
`data/input/biobank/IGVFFI3804AVJR.csv.gz` (downloaded separately from
https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ -- nothing in this repo
produces it). Builds panel i's odds-ratio plot for the IGVF functional
assays, excluding `TSC2_IGVF` since its RapGAP dataset is broken out
separately. Saves `data/output/figures/figure_2/figure_2i.pdf`.

```bash
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_2 \
  r-figures Figure_2i.R
```

## Figure 3

Directory: `notebooks/figures/figure_3/`. Not yet documented --
`curation_summary_figure3.Rmd` plus its `Figure3a/c/d.csv.gz` inputs.

## Figure 4

Directory: `notebooks/figures/figure_4/`. `figure4.ipynb` only *loads* a cached
`figure4_data.json.gz` and calls `plot_figure4()` (in `plot_utils.py`); it
does not compute anything itself. That cache used to be committed directly,
pre-built, with no in-repo script that produced it. It isn't tracked in git
at all right now -- it's moving to an intermediate data directory shortly --
so regenerate it locally with step 1 below before running the notebook.

`src/build_figure4_data.py` now rebuilds most of that cache from current
pipeline outputs. Three sub-panels (the panel a density-band fit, panel c's
out-of-bag confusion matrices, and panel d's illustrative cartoon) have no
source anywhere in this repo and can only be carried forward from a prior
cache -- see `docs/build_figure4_data.md` for exactly which fields and why.

### Steps

```bash
# One-time environment setup (skip if already done)
poetry install --all-extras
poetry run python -m ipykernel install --user --name igvf-cvfg-pillar-project \
  --display-name "IGVF CVFG Pillar Project (Poetry)"

# 1. Rebuild figure4_data.json.gz, carrying forward the fields that can't be
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

Step 1 writes to `notebooks/figures/figure_4/figure4_data.json.gz` by default
(override with `--output`); step 2 reads that same path. Step 2 writes
`data/output/figures/figure_4.png` and a side-effect `executed_figure4.ipynb`
(nbconvert's copy of the notebook with outputs attached, left in
`notebooks/figures/figure_4/`) -- gitignored, so no need to delete it.

`--cached-json` must point at a `figure4_data.json.gz`-shaped file;
`old_figure4_data.json.gz` (a byte-identical backup of the originally
committed cache) works. Without `--cached-json`, step 1 raises immediately
and names the fields it can't source rather than guessing.

See `docs/build_figure4_data.md` for what's actually reconstructed vs.
carried forward, and the known discrepancies (ClinVar/gnomAD drift since the
original figure, and `finalout_4f`'s row-count mismatch in panel f).

## Figure 5/6

Directory: `notebooks/figures/figure_5_6/`. Two independent scripts, run separately:

- `Figure5_6.Rmd` -- Fig 5a/b (controls and ClinGen-repo sankeys, confusion
  matrices, and REVEL/AlphaMissense/MutPred2 metric bar charts) plus Fig
  6a/c/d/e (VUS reclassification sankey + confusion matrix, three-ring donut
  plots by points bin, gnomAD sankey + confusion matrix, and
  unobserved-variant sankey + confusion matrix). Reads
  `Supplementary_Data_5.xlsx` under `data/output/supplementary_data/` -- the
  same workbook the Extended Data Figures use; see `docs/mave_dataset_stats.md`
  for how it's built.
- `Figure_6b.R` -- Fig 6b's per-gene odds-ratio forest plot, faceted by points
  bin and disease group. Reads `data/input/biobank/IGVFFI3804AVJR.csv.gz`,
  downloaded separately from
  https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ (nothing in this repo
  produces it).

Like Extended Data Figures, both need the `r-figures` Docker service
(`ggsankey`, `ggforce`, `patchwork`, `ggh4x`, `extrafont`/Arial, `cairo_pdf`)
rather than a local R install -- see the one-time `docker compose build
r-figures` in that section above.

### Steps

```bash
# Figure_6b.R: reads data/input/biobank/IGVFFI3804AVJR.csv.gz (see link
# above). Writes figure_6b.pdf/figure_6b.svg to data/output/figures/figure_6/.
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures Figure_6b.R

# Figure5_6.Rmd: writes executed_Figure5_6.html next to the .Rmd (gitignored)
docker compose run --rm -w /usr/src/app/notebooks/figures/figure_5_6 \
  r-figures -e 'rmarkdown::render("Figure5_6.Rmd", output_file = "executed_Figure5_6.html")'
```

`Figure5_6.Rmd`'s first chunk sets `DATA_DIR`/`OUT_DIR` to
`/usr/src/app/data/output[/figures]` -- the container's mount point, not a
path relative to the `.Rmd` like `Extended_data_figures.Rmd`'s
`DATA_DIR = "../../../data/output"`. Keep that in mind if you ever run this `.Rmd`
outside the `r-figures` container (e.g. RStudio on macOS): those two lines
will need to point at wherever `data/output` actually lives instead.

`rmarkdown::render()`'s `output_file` arg above renames the knitted HTML
byproduct from the default `Figure5_6.html` to `executed_Figure5_6.html`,
the same `executed_`-prefix convention the notebook steps elsewhere in this
doc use for their own nbconvert-with-outputs-attached byproducts. It still
lands next to the `.Rmd` in `notebooks/figures/figure_5_6/`, and like those
notebook byproducts is gitignored (`executed_*.html`/`executed_*.ipynb` in
`.gitignore`), so no need to delete it by hand.

Fig 5's plots save through `save_my_plot()` to `data/output/figures/figure_5/`;
Fig 6a/6c/6d/6e's plots (VUS, three-ring donut, gnomAD, and unobserved
sankeys/confusion matrices) save via direct `ggsave()` calls to
`data/output/figures/figure_6/`, the same directory `Figure_6b.R` (above)
writes `figure_6b.pdf`/`figure_6b.svg` to.

Every panel `Figure5_6.Rmd` itself builds (i.e. everything except
`Figure_6b.R`'s own output) can be restricted to missense-only variants via
the `consequence_filter` render param -- see [`README.md`](../README.md#figure-56)
for the command. This writes to `figure_5_missense/`/`figure_6_missense/`
instead of `figure_5/`/`figure_6/`, so it doesn't overwrite the default,
all-consequences run; `figure_6_missense/` will not include
`figure_6b.pdf`/`figure_6b.svg`, since those come from the separate,
non-parameterized `Figure_6b.R` script above.

The published Fig 5a/d and Fig 6a/c/d/e panels have counts hand-overlaid on
top of these raw plots in the assembled `PP_Final_Figure3_5_6.ai`. The
sankeys' node counts match the totals already printed on their companion
confusion-matrix-style tiles (`cm_*.png`), so no extra output is needed for
those. The `make_sankey_controls()`/`make_sankey_clingen()` functions (Fig
5a/d) and the three-ring donut functions (Fig 6c) now each save two variants
per plot: a `_lab` version with node/ring count labels baked in, and a
`_nolab`/plain version without (the Fig 6a/d/e VUS/gnomAD/unobserved sankeys
already had a labeled variant; only the Fig 5 controls/ClinGen sankeys and
Fig 6c donuts needed a labeled or unlabeled counterpart added). Pick whichever
matches how a given panel was assembled by hand.

## Extended Data Figures

Directories: `notebooks/figures/extended_data_figure_2/`,
`notebooks/figures/extended_data_figure_5/`, and
`notebooks/figures/extended_data_figure_4_6_7_8_9/`. Ext. Data Figs 2 and 5 are
standalone scripts, documented in their own subsections below;
`Extended_data_figures.Rmd` covers Figs 4-9 and is documented further down.
Ext. Data Fig 7alt and Ext. Data Fig 10 are standalone Python scripts under
`src/` instead (`src/make_extended_data_figure_7alt.py`,
`src/ablation_variant_reclassification.py`), also documented in their own
subsections below.

### Extended Data Figure 2

`Extended_Data_Figure_2.R` -- odds-ratio forest plots for gene-specific vs.
genome-wide predictor calibration (top) and per-assay ExCALIBR classification
(bottom), faceted by gene/disease group. Plain `Rscript`, not an `.Rmd` --
same `tidyverse`/`patchwork`/`ggh4x`/`extrafont` dependencies and `r-figures`
Docker service as `Figure_6b.R` (see Figure 5/6 above). Reads
`data/input/biobank/IGVFFI3804AVJR.csv.gz`, downloaded separately from
https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ (nothing in this repo
produces it).

```bash
# Reads data/input/biobank/IGVFFI3804AVJR.csv.gz (see link above).
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_2 \
  r-figures Extended_Data_Figure_2.R
```

Saves `data/output/figures/extended_data_figure_2.pdf`.

### Extended Data Figure 5

`Extended_Data_Figure_5.ipynb` -- a Python notebook (unlike everything else
in this directory), not R. Flags genes with excess ClinVar pathogenic/benign
discordance against REVEL gene-specific calls: computes each gene's
leave-one-out background discordance rate, tests observed vs. expected with a
binomial test + BH FDR correction, and scatter-plots discordance rate vs.
variant count, colored by significance.

Reads `data/output/predictor_calibration/gene_specific/controls_REVEL_GeneSpecific.csv`,
written by `notebooks/analysis/Variant_Classification_analysis.ipynb` -- see
`notebooks/analysis/README_Variant_Classification_analysis.md` for how to
produce it. No R/Docker dependency here; it only needs the Poetry environment
(`pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`).

```bash
# One-time environment setup (skip if already done -- see Figure 4 above)
poetry install --all-extras
poetry run python -m ipykernel install --user --name igvf-cvfg-pillar-project \
  --display-name "IGVF CVFG Pillar Project (Poetry)"

# The savefig() target directory isn't created for you -- make it first if
# it doesn't already exist.
mkdir -p data/output/figures/extended_data_figure_5

poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_extended_data_figure_5.ipynb \
  notebooks/figures/extended_data_figure_5/Extended_Data_Figure_5.ipynb
```

The notebook resolves its own `PROJECT_ROOT` as `../../..` relative to the
kernel's working directory; `nbconvert` sets that to the notebook's own
directory (`notebooks/figures/extended_data_figure_5/`) automatically, so
`PROJECT_ROOT` lands on the repo root without needing to set the env var it
also supports. Writes
`clinvar_discordance_per_gene.png` to
`data/output/figures/extended_data_figure_5/`, plus a side-effect
`executed_extended_data_figure_5.ipynb` (nbconvert's copy of the notebook
with outputs attached) in `notebooks/figures/extended_data_figure_5/` --
gitignored, so no need to delete it.

### Extended Data Figures 4-9 (`Extended_data_figures.Rmd`)

`Extended_data_figures.Rmd` builds Ext. Data Figs 4-9 (Sankey diagrams,
confusion matrices, VUS reclassification plots, and multi-ring donut plots)
from `Supplementary_Data_5.xlsx` / `Supplementary_Data_6.xlsx` under
`data/output/supplementary_data/` (its `DATA_DIR = "../../../data/output"`
constant, resolved relative to the `.Rmd`'s own directory) -- see
`docs/mave_dataset_stats.md` for how those two workbooks get built.

Unlike the Python `src/` steps, this script needs R packages
(`ggsankey` from GitHub, `ggforce`, Arial via `extrafont`, `cairo_pdf`
output) that aren't part of the Poetry environment. Run it with the
`r-figures` Docker service (`Dockerfile.r`) rather than a local R install.

### Steps

```bash
# One-time image build (skip if the image already exists)
docker compose build r-figures

# Render the whole .Rmd -- knits every chunk in order and writes each
# figure's PNGs via ggsave()/save_my_plot()
docker compose run --rm -w /usr/src/app/notebooks/figures/extended_data_figure_4_6_7_8_9 \
  r-figures -e 'rmarkdown::render("Extended_data_figures.Rmd")'
```

`-w` overrides the service's default working directory (see the `r-figures`
comment in `compose.yaml`) because the `.Rmd` reads/writes paths like
`../../../data/output/...` assuming its own directory is cwd. The repo is
bind-mounted into the container, so the rendered `Extended_data_figures.html`
and every chunk's saved PNGs land directly in your working tree, not just
inside the container.

Every chunk writes under `data/output/figures/`, but not to the same
subfolder -- the fig4/6 controls and ClinGen chunks save through
`save_my_plot()` to `OUT_DIR = "../../../data/output/figures/extended_data_figures"`
(i.e. `data/output/figures/extended_data_figures/Ext_Figure3_5/`), while the
later VUS/gnomAD/donut chunks call `ggsave()` directly with hardcoded
`../../../data/output/figures/extended_data_figure_6/` and
`../../../data/output/figures/extended_data_figures/` paths. Check both
`extended_data_figures/` and `extended_data_figure_6/` under
`data/output/figures/` if a figure you expect isn't where you thought.

### Extended Data Figure 4 missense

A missense-only companion to Extended Data Figure 4 (AlphaMissense/MutPred2 x
ClinVar controls/ClinGen Evidence Repository, via ExCALIBR/GeneSpecific
calibration), built the same way as Extended Data Figure 6alt below: a new,
isolated chunk right after the existing Ext. Data Fig 4/6 block in
`Extended_data_figures.Rmd` -- doesn't touch or reuse any of that block's own
objects. Adds a missense-only variant of every panel (`simplified_consequence
== "missense_variant"`, same column/idea as `Figure5_6.Rmd`'s own
`consequence_filter` param and Ext. Data Figure 6alt), doubling the current 8
sankeys+confusion-matrices to 16 -- meant to sit as an All-variants block and
a Missense-only block side by side. Every chart is shrunk 50% linearly from
Ext. Data Fig 4's own calibrated core dimensions (sankey 37x64mm ->
18.5x32mm, confusion matrix 31x27mm -> 15.5x13.5mm) with unchanged font
sizes, same reasoning as Ext. Data Figure 6alt. Sankey node labels are
abbreviated to P/LP/VUS/LB/B (`label_overrides`), including ClinGen's own
"No Classification" bucket, mapped to "VUS" like "Uncertain".

Writes 16 PDFs to `data/output/figures/extended_data_figure_4_missense/`,
named like `sankey_clinvar_ExOP_AM_GeneSpecific_missense_calibrated.pdf` /
`cm_clingen_ExOP_MP2_GeneSpecific_all_calibrated.pdf`.

### Extended Data Figure 6alt

An alternate version of Extended Data Figure 6 (REVEL/AlphaMissense/MutPred2
x ClinVar controls/ClinGen Evidence Repository, all via OddsPath/"Universal"
calibration), built as a new, isolated chunk right after the existing Ext.
Data Fig 6 block in `Extended_data_figures.Rmd` -- doesn't touch or reuse any
of that block's own objects. Adds a missense-only variant of every panel
(`simplified_consequence == "missense_variant"`, same column/idea as
`Figure5_6.Rmd`'s own `consequence_filter` param), doubling the current 12
sankeys+confusion-matrices to 24 -- meant to sit as an All-variants block and
a Missense-only block side by side. Every chart is shrunk 50% linearly from
Ext. Data Fig 6's own calibrated core dimensions (sankey 37x64mm ->
18.5x32mm, confusion matrix 26x23mm -> 13x11.5mm) with unchanged font sizes,
so each chart's own *returned* width/height (label/title overflow included,
same as every other calibrated chart here) ends up noticeably larger than
that core size -- the same-size text no longer shrinks with the box. Sankey
node labels are abbreviated to P/LP/VUS/LB/B (`label_overrides`), including
ClinGen's own "No Classification" bucket, mapped to "VUS" like "Uncertain".

Writes 24 PDFs to `data/output/figures/extended_data_figure_6alt/`, named
like `sankey_clinvar_OP_REVEL_Universal_missense_calibrated.pdf` /
`cm_clingen_OP_MP2_Universal_all_calibrated.pdf`. This also fixed a real bug
in `confusion_matrix_calibrated.R`: its canvas width previously assumed the
title/x-axis-label text would always be narrower than the matrix box itself
(true at the original 26x23/31x27mm sizes) -- at 13x11.5mm with the same
unchanged font size that stopped holding, clipping the title against the PDF
page edge. Fixed by widening `x_range` to also cover the title's/x-label's
own text width when it exceeds the matrix's; a no-op for every existing
(larger) calibrated confusion matrix in this project.

### Extended Data Figure 7alt (`src/make_extended_data_figure_7alt.py`)

A candidate replacement for Extended Data Figure 7: a single calibrated
(exact-print-size) heatmap superseding the classification info spread across
Figure 6 (VUS/gnomAD/unobserved x REVEL only) and Extended Data Figure 7 (VUS
x REVEL/AM/MP2 only). Nine rows -- REVEL/AM/MP2 for each of the VUS, gnomAD,
and unobserved variant sets, with a gap between the three groups -- x five
classification columns (P/LP/VUS/LB/B, reading each sheet's own
`Class_{REVEL,AM,MP2}` column), plus a Total column. Each cell shows its count
and row-wise percentage, colored on a single monochromatic 0-100% scale
(with a legend) shared across every cell.

Unlike the rest of Extended Data Figs 4-9, this one is a standalone Python
script (matplotlib) rather than an `Extended_data_figures.Rmd` chunk, and
reads `Supplementary_Data_5.xlsx` directly -- no R/Docker required:

```bash
poetry run python -m src.make_extended_data_figure_7alt
```

Writes `data/output/figures/extended_data_figure_7alt/new_classification_heatmap.pdf`
by default; see `--input`/`--output` (`--help`) to override either path.

`--consequence-filter missense` restricts every sheet to `simplified_consequence
== "missense_variant"` rows (the same column/idea as `Figure5_6.Rmd`'s own
`consequence_filter` param) and writes to
`extended_data_figure_7alt_missense/new_classification_heatmap_missense.pdf`
instead, so it doesn't overwrite the all-consequences run.

### Extended Data Figure 10 (`src/ablation_variant_reclassification.py`)

The functional-vs-predictor evidence ablation analysis (how much of the
reclassification pipeline's effect comes from functional evidence alone,
predictor evidence alone, or genuine added value from combining the two --
see [`docs/ablation_variant_reclassification.md`](ablation_variant_reclassification.md)
for the full method). Reads the same `Variant_Classification_analysis.ipynb`
checkpoint as `src/build_variant_reclassification_dataset.py`
(`data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz`
by default) directly -- no R/Docker required, and no dependency on
Supplementary Data 5/6.

Unlike every other Extended Data Figure here, this one is generated as many
small single-panel PDFs (one per chart type x variant-category scope,
following the same "individual panels assembled by hand" pattern as
Extended Data Figure 6alt's 24 PDFs above) rather than one PNG/PDF or one
`.Rmd`-rendered document -- `--document-split-dir` renders exactly the
per-section chart data `--document`'s combined multi-section image would,
but writes each `(chart type, scope)` block to its own file instead, each
already carrying its own legend (no shared or document-level legend to
reconstruct). See `docs/ablation_variant_reclassification.md`'s "Splitting
the document into individual chart files" section for the full mechanics.

All three predictors (REVEL, AlphaMissense, MutPred2 -- the default when
`--predictor` is omitted) and all four chart types (ablation, comparison,
concordance, gain -- the default when `--document-chart` is omitted):

```bash
poetry run python -m src.ablation_variant_reclassification \
  --document-split-dir data/output/figures/extended_data_figure_10
```

Writes 25 PDFs (`--document-split-dir` creates
`data/output/figures/extended_data_figure_10/` if it doesn't already exist):
one `ablation_<scope>.pdf`, `comparison_<scope>.pdf`, and `gain_<scope>.pdf`
per `--scope` category (`vus`, `gnomad`, `unobserved`, `clinvar_control`,
`clinvar_control_missense_only`, `clingen_control`,
`clingen_control_missense_only` -- all seven, regardless of `--scope`, same
as `--document`), plus one `concordance_<scope>.pdf` for just the four
control scopes (`clinvar_control`/`clinvar_control_missense_only`/
`clingen_control`/`clingen_control_missense_only` -- the other three scopes
carry no known truth to check concordance against, so that file is skipped
there rather than written empty). 7 + 7 + 4 + 7 = 25 PDFs total. The
`_missense_only`-suffixed scopes are strictly `condensed_consequence ==
"missense_variant"` (excluding start-loss, see `is_missense_only`).

Takes several minutes end to end on the real integrated dataset (dominated
by rendering all 25 panels' worth of charts, not the dedup pass) -- expected,
not a hang; the "Wrote N chart files" line at the end confirms completion.
Pass `--predictor`/`--document-chart` to narrow it down if only a subset of
panels is needed. `--consequence missense_only` (see
`docs/ablation_variant_reclassification.md`) restricts the whole run to
missense variants if needed for a supplementary check -- not part of this
figure's regular generation.

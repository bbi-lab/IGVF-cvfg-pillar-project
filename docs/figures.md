# Figures

How to (re)generate each manuscript figure from `Main_Figures/`.

## Figure 2

Directory: `Main_Figures/Figure_2/`. Six scripts: one prep notebook everything
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
`20260101_SGEsubset.xlsx`, `20260101_VAMPseqsubset_wDups.xlsx`, and
`20260101_CAVAseqsubset.xlsx` (the concatenation of the first two --
not currently read by anything else in this directory). Every panel
notebook below reads the first two of these.

```bash
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_PP_ProcessBigDataFrame.ipynb \
  Main_Figures/Figure_2/PP_ProcessBigDataFrame.ipynb
```

### 2. Panel notebooks (independent of each other; all read step 1's output)

- `PP_ClinVarPrecisionRecall.ipynb` -- precision/recall of SGE and VAMP-seq
  against ClinVar, with Wilson-interval CIs. Saves
  `data/output/figures/figure_2/<YYYYMMDD>_PillarProject_PRvsClinVar_wErrorBar_grey.svg`
  (filename is date-stamped by the run date, not fixed).
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
  `data/output/figures/figure_2/20260120_vampseq_sge_bars.svg`.
- `PP_SeqFunctionMap.ipynb` -- the RAD51D sequence-function map. Saves
  `data/output/figures/figure_2/20251204_RAD51D_X9_draft_SeqFunc_map_extended.svg`.
- `PP_StackedHistograms.ipynb` -- stacked score histograms (with a ClinVar
  density overlay) for both assays, plus per-gene SGE insets. Saves two SVGs
  under `data/output/figures/figure_2/Histogram_wStripplot/` and one
  `data/output/figures/figure_2/20260120_sge_histogram_inset_<gene>.svg` per
  SGE gene.

```bash
for nb in PP_ClinVarPrecisionRecall PP_Fig2_Heatmaps PP_ResolutionOverview PP_SeqFunctionMap PP_StackedHistograms; do
  poetry run jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
    --ExecutePreprocessor.timeout=600 \
    --output executed_${nb}.ipynb \
    Main_Figures/Figure_2/${nb}.ipynb
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
attached, same as the Figure 4 step) -- delete these if you don't want them
in the working tree.

### 3. `Figure_2i.R` (separate from the notebook pipeline above)

Same pattern as `Figure_6b.R` / `Extended_Data_Figure_2.R`: a plain Rscript
run via the `r-figures` Docker service, reading `IGVFFI3804AVJR.csv.gz` from
its own directory (downloaded separately from
https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ -- nothing in this repo
produces it). Builds panel i's odds-ratio plot for the IGVF functional
assays, excluding `TSC2_IGVF` since its RapGAP dataset is broken out
separately. Saves `fig2i.pdf` directly into `Main_Figures/Figure_2/`, not
under `data/output/figures/`.

```bash
docker compose run --rm -w /usr/src/app/Main_Figures/Figure_2 \
  r-figures Figure_2i.R
```

## Figure 3

Directory: `Main_Figures/Figure_3/`. Not yet documented --
`curation_summary_figure3.Rmd` plus its `Figure3a/c/d.csv.gz` inputs.

## Figure 4

Directory: `Main_Figures/Figure_4/`. `figure4.ipynb` only *loads* a cached
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
  --cached-json Main_Figures/Figure_4/old_figure4_data.json.gz

# 2. Execute the notebook to produce fig4.png
poetry run jupyter nbconvert --to notebook --execute \
  --ExecutePreprocessor.kernel_name=igvf-cvfg-pillar-project \
  --ExecutePreprocessor.timeout=600 \
  --output executed_figure4.ipynb \
  Main_Figures/Figure_4/figure4.ipynb
```

Step 1 writes to `Main_Figures/Figure_4/figure4_data.json.gz` by default
(override with `--output`); step 2 reads that same path. Step 2 writes
`Main_Figures/Figure_4/fig4.png` and a side-effect `executed_figure4.ipynb`
(nbconvert's copy of the notebook with outputs attached) -- delete the latter
afterward if you don't want it in the working tree.

`--cached-json` must point at a `figure4_data.json.gz`-shaped file;
`old_figure4_data.json.gz` (a byte-identical backup of the originally
committed cache) works. Without `--cached-json`, step 1 raises immediately
and names the fields it can't source rather than guessing.

See `docs/build_figure4_data.md` for what's actually reconstructed vs.
carried forward, and the known discrepancies (ClinVar/gnomAD drift since the
original figure, and `finalout_4f`'s row-count mismatch in panel f).

## Figure 5/6

Directory: `Main_Figures/Figure5_6/`. Two independent scripts, run separately:

- `Figure5_6.Rmd` -- Fig 5a/b (controls and ClinGen-repo sankeys, confusion
  matrices, and REVEL/AlphaMissense/MutPred2 metric bar charts) plus Fig
  6a/c/d/e (VUS reclassification sankey + confusion matrix, three-ring donut
  plots by points bin, gnomAD sankey + confusion matrix, and
  unobserved-variant sankey + confusion matrix). Reads
  `Supplementary_Data_5.xlsx` under `data/output/supplementary_data/` -- the
  same workbook the Extended Data Figures use; see `docs/mave_dataset_stats.md`
  for how it's built.
- `Figure_6b.R` -- Fig 6b's per-gene odds-ratio forest plot, faceted by points
  bin and disease group. Reads `IGVFFI3804AVJR.csv.gz`, downloaded separately
  from https://data.igvf.org/tabular-files/IGVFFI3804AVJR/ (nothing in this
  repo produces it) and expected in the same directory as the script.

Like Extended Data Figures, both need the `r-figures` Docker service
(`ggsankey`, `ggforce`, `patchwork`, `ggh4x`, `extrafont`/Arial, `cairo_pdf`)
rather than a local R install -- see the one-time `docker compose build
r-figures` in that section above.

### Steps

```bash
# Figure_6b.R: place IGVFFI3804AVJR.csv.gz in this directory first (see link
# above). Writes fig6b.pdf/fig6b.svg into the same directory.
docker compose run --rm -w /usr/src/app/Main_Figures/Figure5_6 \
  r-figures Figure_6b.R

# Figure5_6.Rmd
docker compose run --rm -w /usr/src/app/Main_Figures/Figure5_6 \
  r-figures -e 'rmarkdown::render("Figure5_6.Rmd")'
```

`Figure5_6.Rmd`'s first chunk sets `DATA_DIR`/`OUT_DIR` to
`/usr/src/app/data/output[/figures]` -- the container's mount point, not a
path relative to the `.Rmd` like `Extended_data_figures.Rmd`'s
`DATA_DIR = "../data/output"`. Keep that in mind if you ever run this `.Rmd`
outside the `r-figures` container (e.g. RStudio on macOS): those two lines
will need to point at wherever `data/output` actually lives instead.

Fig 5's plots save through `save_my_plot()` to `data/output/figures/figure_5/`;
Fig 6a/6c/6d/6e's plots (VUS, three-ring donut, gnomAD, and unobserved
sankeys/confusion matrices) save via direct `ggsave()` calls to
`data/output/figures/figure_6/`.


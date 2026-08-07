# Figures

How to (re)generate each manuscript figure from `Main_Figures/`. Per
AGENTS.md's roadmap, these are still notebooks/scripts rather than converted
`src/` pipeline steps, so "how do I run this" isn't always obvious from the
directory alone -- this doc is where that gets written down, one figure at a
time, starting with Figure 4.

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


# Data

Overview of every file and directory under `data/`: what's bundled in this
repo, what you need to add or download yourself, the structure of the
gitignored intermediate/staging area, and what each pipeline stage writes out.
See [`docs/variant_annotation_pipeline.md`](variant_annotation_pipeline.md),
[`docs/build_training_variant_files.md`](build_training_variant_files.md),
[`docs/fetch_mavedb_scores.md`](fetch_mavedb_scores.md),
[`docs/load_excalibr_calibrations.md`](load_excalibr_calibrations.md), and
[`docs/mave_dataset_stats.md`](mave_dataset_stats.md) for the operational
detail behind each script mentioned here; this page is the map, not the
manual.

## Input files

### Bundled in this repo (no action needed)

| Path | Contents |
|---|---|
| `data/input/maves/score_sets.tsv` | Dataset-name / MaveDB score-set URN mapping. Read by variant-annotation pipeline Steps 11 and 15. |
| `data/input/maves/Supplementary_Data_3.xlsx` | Assay/dataset curation metadata (the paper's Supplementary Data 3): one row per dataset, `IGVF Produced?`, `Primary Score Set or Meta-analysis?`, gene, HGNC/Entrez IDs, etc. Read by variant-annotation pipeline Step 16, `build_training_variant_files`, and `mave_dataset_stats`. |
| `data/input/maves/Supp_Data1_column_description.xlsx` | Column documentation for Supplementary Data 1 (the integrated variant-effect dataset). |
| `data/input/maves/CHEK2_Gebbia_2024.xlsx` | CHEK2 variant filtering annotations, referenced by the OddsPath notebooks. |
| `data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz` | RefSeq/Ensembl transcript accession mapping. Small enough to commit despite being generic reference data; read by variant-annotation pipeline Step 2. |
| `data/input/reference/extended_ensembl_consequence.csv.gz` | CVFG-specific VEP-term-to-SO-summary-term mapping, read by `annotate_simplified_consequence` (Step 17). |
| `data/input/predictors/clustering_variants_revel_training_overlap.csv.gz`, `data/input/predictors/mp2_actual_training_data.txt.gz`, `data/input/predictors/supplementary_data_3_gene_entrez_ids.tsv` | Upstream training-variant sources consumed by `build_training_variant_files` (a preparatory step Step 13 runs automatically) — see [`docs/build_training_variant_files.md`](build_training_variant_files.md). |
| `data/input/predictors/mp2_annotations.csv.gz`, `data/input/predictors/mp2_gene_symbol_map.tsv` | MutPred2 scores (`gene_symbol`/`AA`/`MutPred2 score`) and the `CALM1, CALM2, CALM3` → `CALM1` gene-symbol remap, read directly by Step 16 (`--mutpred2-properties-file`/`--mutpred2-gene-symbol-map-file`). Small enough to commit, unlike the AlphaMissense/REVEL files below. See [Where `mp2_annotations.csv.gz` comes from](#where-mp2_annotationscsvgz-comes-from) below. |
| `data/filtering/*.tsv` | Per-dataset variant-filtering annotations (`BAP1_Waters_2024.tsv`, `CHEK2_Gebbia_2024.tsv`, `DDX3X_Radford_2023.tsv`, `OTC_Lo_2023.tsv`, `RAD51C_Olvera-León_2024.tsv`), read by `flag_variants` (`--filtering-dir`, default `data/filtering`). |
| `data/input/genes/gencc-submissions.csv.gz` | GenCC gene-disease validity submissions, downloaded 2026-09-10 from `https://thegencc.org/download/action/submissions-export-csv` (legacy UUID-based export; CC0 1.0 Public Domain, attribution requested). ⚠️ GenCC is discontinuing this legacy format on 2026-09-30 in favor of an SGC-ID-based export at `?format=new` — re-fetching after that date needs `build_figure3_data.py`'s GenCC-column handling updated to match. Read by `build_figure3_data` (Figure 3a). |
| `data/input/genes/uniprotkb_9606_reviewed.tsv.gz` | Reviewed (Swiss-Prot) human proteome, downloaded 2026-09-10 from UniProt's REST API (`rest.uniprot.org/uniprotkb/stream`, `reviewed:true AND organism_id:9606`, CC BY 4.0 — attribution required). Read by `build_figure3_data` (Figure 3a, protein length → possible SNVs). |
| `data/input/genes/test_condition_gene.txt.gz` | NCBI Genetic Testing Registry's public bulk export, downloaded 2026-09-10 from `https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt` (US government work, public domain, attribution requested). Read by `build_figure3_data` (Figure 3a, clinical test counts per gene). |

The three `data/input/genes/` files above are periodic snapshots, not
one-time downloads — see
[`docs/build_figure3_data.md`](build_figure3_data.md#refreshing-these-files)
for the exact commands to refresh them, and
[`README.md`](../README.md#third-party-data--licenses) for the attribution
their licenses require.

`data/input/predictors/` also holds `data_frame_missense_variants_MP2_properties.csv.gz`,
`mp2_annotated.tsv`, and `ddx3x_mutpred2_annotated.tsv` — the three raw
sources `mp2_annotations.csv.gz` is assembled from (see [Where
`mp2_annotations.csv.gz` comes from](#where-mp2_annotationscsvgz-comes-from)
below) — plus `ddx3x_mutpred2_scores_hgvsp.tsv`, which isn't read by
anything in the current pipeline or by `combine_mutpred2_annotations.py`
(it looks like an upstream, pre-gene/AA-annotation predecessor of
`ddx3x_mutpred2_annotated.tsv`) — see [Roadmap](#roadmap) below.

### Not bundled — you must add or download or provide these yourself

Not committed to git, either because the file is too large or because it's
fetched/regenerated per checkout:

| What | Where it goes | How to get it |
|---|---|---|
| `cvfg_variants.0.tsv` | `data/input/maves/` | Produced via `src/mavedb_scores.sql` (DB access) or `src/fetch_mavedb_scores.py` (public API) — see [Where `cvfg_variants.0.tsv` comes from](#where-cvfg_variants0tsv-comes-from) below. |
| `AlphaMissense_hg38.tsv.gz` (+ `.tbi`), `revel_hg38.tsv.gz` (+ `.tbi`) | `data/input/predictors/` | Pinned predictor-score extracts specific to this pipeline run, read by Step 16. Hundreds of MB each — provide your own copies rather than committing them. |
| `data/input/mave_calibration/excalibr/json/`, `data/input/mave_calibration/excalibr/png/` | `data/input/mave_calibration/excalibr/` | Per-dataset exCALIBR point-value calibrations (JSON, read by `load_excalibr_calibrations`) and calibration plots (PNG, visual reference only). Gitignored (`data/input/mave_calibration/excalibr/*` in `.gitignore`) — large and regenerated/hand-supplied per checkout; only `.gitkeep` is committed. |
| Large, generic reference caches: SpliceAI VCFs, dbNSFP, `clinvar_cache/`, and the gnomAD Hail table cache | Inside your `variant-annotation` checkout (not this repo) | Follow `variant-annotation`'s own setup; point `VARIANT_ANNOTATION_DIR` at a checkout that already has them. Tens of GB — deliberately not duplicated here. Build/refresh the gnomAD cache once per checkout with `scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache` (~6-7 hours). |
| `IGVFFI3804AVJR.csv.gz` | `data/input/biobank/` | Download from https://data.igvf.org/tabular-files/IGVFFI3804AVJR/. Only needed for a subset of figures (Section 3 of the [top-level README](../README.md)), not for pipeline Stages 1-2. Gitignored (`data/input/biobank/*` in `.gitignore`) -- only `.gitkeep` is committed. |
| A `variant-annotation` checkout itself | Wherever `VARIANT_ANNOTATION_DIR` points, or `vendor/variant-annotation/` (vendored submodule) | `git submodule update --init vendor/variant-annotation`, or point `VARIANT_ANNOTATION_DIR` at your own existing checkout. |

### Where `cvfg_variants.0.tsv` comes from

`data/input/maves/cvfg_variants.0.tsv` is *not* committed to the repo — it's
regenerated locally from one of two different data sources:

- **`src/mavedb_scores.sql`** — the source of truth. Run directly against a
  local Postgres mirror of the MaveDB database. Requires DB access most
  contributors won't have. Also holds the score-set URN list and five
  dataset-specific manual overrides (CHEK2 preferred transcript, BRCA1_Findlay
  MANE Select version lift, NDUFAF6/PTEN `raw_hgvs_nt` blanking, and
  JAG1/TARDBP/SCN5A/KCNH2 target-sequence corrections).
- **`src/fetch_mavedb_scores.py`** — an API-based equivalent, for contributors
  without DB access. Hits MaveDB's public API (`https://api.mavedb.org`)
  instead of a DB mirror, parsing the same score-set list and manual overrides
  directly out of `src/mavedb_scores.sql` at run time so the two can't drift
  apart. **Not yet validated against the SQL version's output** — this
  equivalence is in progress and will be tested soon. Known differences
  (no internal `variant_id`, at most one of `hgvs_g`/`_c`/`_p` populated per
  row) are documented in [`docs/fetch_mavedb_scores.md`](fetch_mavedb_scores.md).

  ```bash
  poetry run python -m src.fetch_mavedb_scores \
    --output data/input/maves/cvfg_variants.0.api.tsv
  ```

### Where `mp2_annotations.csv.gz` comes from

Assembled from three sources, each reduced to `gene_symbol`/`AA`/`MutPred2
score` and concatenated row-wise (not joined): MaveDB's own per-variant
MP2-properties export (`data_frame_missense_variants_MP2_properties.csv.gz`,
1,120,798 rows), a small curated gene/amino-acid-change-keyed supplement
(`mp2_annotated.tsv`, 985 rows), and a DDX3X-specific supplement
(`ddx3x_mutpred2_annotated.tsv`, 3,867 rows -- DDX3X appears in neither of
the other two sources, confirmed directly by grepping both).

`src/combine_mutpred2_annotations.py` only accepts two positional arguments
(`properties_file`, `annotated_file`), so it can reproduce the first two
sources but not the DDX3X supplement baked into the committed
`mp2_annotations.csv.gz` -- confirmed by row-count arithmetic (1,120,798 +
985 + 3,867 matches the committed file's row count) and by the DDX3X gene
check above. Regenerating `mp2_annotations.csv.gz` from scratch today needs
an extra manual step to fold `ddx3x_mutpred2_annotated.tsv` in -- e.g. run
the command below, then separately concatenate `ddx3x_mutpred2_annotated.tsv`
(reduced to the same three columns, with its `MutPred2_score` column renamed
to `MutPred2 score`) onto the result -- or extend the script to accept a
third file. `mp2_annotations.csv.gz` and `mp2_gene_symbol_map.tsv` are
committed directly (see [Bundled](#bundled-in-this-repo-no-action-needed)
above) rather than regenerated per run, so this only matters when refreshing
them after a predictor-data update:

```bash
poetry run python -m src.combine_mutpred2_annotations \
  data/input/predictors/data_frame_missense_variants_MP2_properties.csv.gz \
  data/input/predictors/mp2_annotated.tsv \
  data/input/predictors/mp2_annotations.csv.gz
```

### Variant-annotation pipeline: data dependencies / prerequisites

Everything the Dockerized variant-annotation pipeline
(`scripts/run_variant_annotation_pipeline.sh`) needs before a run, beyond
what's listed above:

| Prerequisite | Where to get it | Where it goes |
|---|---|---|
| Docker and Docker Compose | Your package manager / docker.com | N/A |
| A `variant-annotation` checkout with large reference caches already downloaded | See [Not bundled](#not-bundled--you-must-add-or-download-or-provide-these-yourself) above | `export VARIANT_ANNOTATION_DIR=/path/to/checkout`, or the vendored submodule at `vendor/variant-annotation` |
| gnomAD Hail table cache (Step 7 prerequisite) | Built locally, once per `variant-annotation` checkout | `scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache`; cache lands in the `variant-annotation-gnomad-cache` Docker volume scoped to that checkout |
| Poetry environment | This repo | `poetry install --all-extras` — needed to build this project's own Docker image (`Dockerfile`/`compose.yaml`) for the steps it contributes to the pipeline |
| All of the [Bundled](#bundled-in-this-repo-no-action-needed) input files, plus `cvfg_variants.0.tsv` and `AlphaMissense_hg38.tsv.gz`/`revel_hg38.tsv.gz` from [Not bundled](#not-bundled--you-must-add-or-download-or-provide-these-yourself) above | This repo (bundled ones) / provided by you (the rest) | `scripts/run_variant_annotation_pipeline.sh` stages the relevant ones into `data/intermediate/variant_annotation/data/` automatically per run |

See [`docs/variant_annotation_pipeline.md`](variant_annotation_pipeline.md)
for the full step-by-step data flow and the `VARIANT_DATA_DIR` path-mapping
subtlety.

## `data/intermediate/`: structure

Gitignored entirely (`data/intermediate/` in `.gitignore`). Regenerated by
`scripts/run_variant_annotation_pipeline.sh` on every run — nothing here
should be hand-edited or relied on across runs.

```
data/intermediate/
├── variant_annotation/
│   └── data/                              # Mounted as VARIANT_DATA_DIR for the
│       │                                   # variant-annotation pipeline run.
│       ├── cvfg_variants.0.tsv             # Staged from data/input/maves/ (Step 1 only)
│       ├── cvfg_variants.1.tsv ... .19.tsv # One per numbered pipeline step
│       ├── cvfg_variants.20*.tsv           # Step 20's flattening intermediates
│       ├── integrated_variant_effect_dataset*.tsv  # Step 21's final outputs, pre-gzip
│       ├── score_sets.tsv, Supplementary_Data_3.xlsx,
│       │   MANE.GRCh38.v1.5.summary.txt.gz # Staged from data/input/ every run
│       ├── AlphaMissense_hg38.tsv.gz, revel_hg38.tsv.gz (+ .tbi),
│       │   mp2_annotations.csv.gz, mp2_gene_symbol_map.tsv
│       │                                   # Staged from data/input/predictors/ for Step 16
│       └── revel_training_variants.tsv, mutpred2_training_variants.tsv
│                                           # Written by build_training_variant_files
│                                           # ahead of Step 13
└── figures/
    └── figure_2/                          # PP_ProcessBigDataFrame.ipynb's split-out
        ├── SGEsubset.xlsx                 # subsets for Figure 2's panel notebooks
        ├── VAMPseqsubset_wDups.xlsx
        └── CAVAseqsubset.xlsx
```

## Output files

Gitignored except per-directory `.gitkeep` placeholders (see `.gitignore`).

| Path | Contents | Produced by |
|---|---|---|
| `data/output/maves/integrated_variant_effect_dataset.tsv.gz` | Expanded integrated variant-effect dataset — one row per DNA-level variant (Supplementary Data 1) | `scripts/run_variant_annotation_pipeline.sh` (Stage 1) |
| `data/output/maves/integrated_variant_effect_dataset.condensed.tsv.gz` | Condensed sibling — one row per variant-effect measurement/score | same |
| `data/output/mave_calibration/OddsPath_calibrations.csv.gz` | Per-dataset OddsPath likelihood ratios | `OddsPath_calculations.ipynb` |
| `data/output/mave_calibration/oddspath/*.csv` | Genome-wide OddsPath-based classification files, one row per variant (representative-only), one per category (`controls`, `VUS`, `ClinGen_Repo`, `gnomAD`, `Unobserved`) × predictor (REVEL/AM/MP2) | `OddsPath_classifications.ipynb` |
| `data/output/mave_calibration/oddspath_with_secondary_variants/*.csv` | Same, but every candidate NT variant kept (`Variant_Role` column) instead of just the aa-stage tie-break winner | `OddsPath_classifications.ipynb` |
| `data/output/predictor_calibration/gene_specific/*.csv` | Same category × predictor breakdown as above, but using gene-specific rather than genome-wide predictor thresholds; one row per variant (representative-only) | `Variant_Classification_analysis.ipynb` |
| `data/output/predictor_calibration/gene_specific_with_secondary_variants/*.csv` | Same, but every candidate NT variant kept (`Variant_Role` column) | `Variant_Classification_analysis.ipynb` |
| `data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz` | Full per-variant reclassification analysis dataset | `Variant_Classification_analysis.ipynb` |
| `data/output/supplementary_data/Supplementary_Data_4.xlsx` | exCALIBR calibrations workbook | `load_excalibr_calibrations` (`ExCALIBR_calibrations` sheet), consumed by `OddsPath_classifications.ipynb` |
| `data/output/supplementary_data/Supplementary_Data_5.xlsx` | Gene-specific calibration/reclassification results, combined into Excel -- one row per variant | `Variant_Classification_analysis.ipynb` |
| `data/output/supplementary_data/Supplementary_Data_5.with_secondary_variants.xlsx(.gz)` | Same, but every candidate NT variant kept (`Variant_Role` column) instead of just the aa-stage tie-break winner -- see [`docs/variant_classification.md`](variant_classification.md#aa-stage-tie-break-clinvarclingen-evidence-quality) | `Variant_Classification_analysis.ipynb` |
| `data/output/supplementary_data/Supplementary_Data_6.xlsx(.gz)` | Genome-wide OddsPath-based classification results, combined into Excel -- one row per variant | `OddsPath_classifications.ipynb` |
| `data/output/supplementary_data/Supplementary_Data_6.with_secondary_variants.xlsx(.gz)` | Same, but every candidate NT variant kept (`Variant_Role` column) | `OddsPath_classifications.ipynb` |
| `data/output/figures/` | Manuscript figure outputs (SVG/PNG), one subdirectory per figure | `notebooks/figures/*` notebooks/scripts |
| `stats.txt` (or wherever `--output` points) | Dataset/measurement/variant counts, predictor score coverage, clinical-attribute breakdowns, ExCALIBR calibration coverage, reclassification agreement | `src/mave_dataset_stats.py` |

## Roadmap

Suggestions for keeping the data layout internally consistent, in rough
priority order:

1. **Add explicit `.gitignore` entries for the "never commit" predictor/MAVE
   inputs.** `cvfg_variants.0.tsv`, `AlphaMissense_hg38.tsv.gz`/`.tbi`, and
   `revel_hg38.tsv.gz`/`.tbi` are now a deliberate policy of "provide
   locally, don't commit" (this page's [Not bundled](#not-bundled--you-must-add-or-download-or-provide-these-yourself)
   section), the same way `data/input/mave_calibration/excalibr/*` already
   is. Without matching `.gitignore` patterns, a future `git add -A` could
   commit them by accident.
2. **Extend `combine_mutpred2_annotations.py` to also fold in
   `ddx3x_mutpred2_annotated.tsv`** (see [Where `mp2_annotations.csv.gz`
   comes from](#where-mp2_annotationscsvgz-comes-from) above), so
   regenerating `mp2_annotations.csv.gz` after a predictor-data refresh is a
   single command again instead of needing a manual extra concatenation
   step. Also decide whether `ddx3x_mutpred2_scores_hgvsp.tsv` (not read by
   anything) should be removed.
3. **Update stale `data/raw_mave_data/` references now that it's gone.**
   `scripts/run_variant_annotation_pipeline.sh`'s header comment and data-flow
   diagram, `docs/variant_annotation_pipeline.md`'s "Data flow" section and
   several per-file exception notes, `README.md`'s Stage 1 description, and
   `AGENTS.md`'s repository map and "See also" section all still describe or
   link to `data/raw_mave_data/README.md`. The directory's actual staging
   step (`rsync -a "$CVFG_PROJECT_DIR/data/raw_mave_data/" ...`) was already
   commented out in `scripts/run_variant_annotation_pipeline.sh` before this
   removal, so the directory was dead weight in practice — but the prose
   describing it is now actively wrong and should be updated to point at
   `data/input/maves/`, `data/input/predictors/`, and `data/input/reference/`
   instead. The `.gitignore` comment above `data/intermediate/` ("Staging
   copy of `data/raw_mave_data/`...") should be corrected too.
4. **Remove `data/figures/`.** It's an untracked, byte-for-byte duplicate of
   `data/output/figures/extended_data_figures/Ext_Figure3_5/` and isn't
   referenced by any script, notebook, or doc — looks like a stray copy left
   over from a manual run.
5. **Clean up `data/to organize/`.** Contains a stale `stats.txt` (referencing
   an old `data/mave_data/` output path that predates the `data/output/maves/`
   rename) and a duplicate `raw_mave_data/README.md`. Worth triaging and
   deleting now that the real `README.md` is gone.
6. **Relocate the stray `data/input/mave_calibration/DDX3X_Radford_2023.json`.**
   It sits directly under `data/input/mave_calibration/` rather than
   `data/input/mave_calibration/excalibr/json/` alongside every other
   calibration file, and its content is a slightly older revision of
   `data/input/mave_calibration/excalibr/json/DDX3X_Radford_2023.json`
   (fewer `point_ranges` entries) rather than an identical copy — worth
   confirming which is current before deleting the other.
7. **Consolidate the root-level `Supplementary_Data_*.v1.xlsx` files** (and
   the `data/output/supplementary_data/work/` scratch copies) into
   `data/output/supplementary_data/` or a dedicated `data/output/submission/`
   directory, rather than sitting at the repo root — keeps generated
   manuscript-submission artifacts out of the top-level listing alongside
   source files.
8. **Refresh the per-notebook `README_*.md` files** under
   `notebooks/analysis/` (`README_OddsPath_calculations.md`,
   `README_OddsPath_classifications.md`, `README_Variant_Classification_analysis.md`)
   — they still describe outputs landing in `outputs/`, `outputs/GeneSpecific_calibrations/`,
   and inputs at `data/Supplemental_Data/` and `data/raw_mave_data/`, all of
   which predate the current `data/input/`/`data/output/` layout this page
   documents. `AGENTS.md`'s notebook-to-script conversion conventions already
   call for folding these into `docs/` — doing so would resolve the drift.
9. **Fix the stale `exCALIBR/` entry in `AGENTS.md`'s repository map.** It
   describes a top-level `exCALIBR/` directory that no longer exists;
   calibration data now lives entirely under
   `data/input/mave_calibration/excalibr/`.

# build_figure3_data

Rebuilds `data/intermediate/figures/figure_3/Figure3a.csv.gz`,
`Figure3c.csv.gz`, `Figure3c_measurements.csv.gz`, and `Figure3d.csv.gz` —
all four data files `notebooks/figures/figure_3/curation_summary_figure3.Rmd`
reads. All four were originally produced by an ad hoc, uncommitted notebook
(`data/output/Curation_summary_V5_cleaned.ipynb`) pointed at a stale personal
directory layout; this script reproduces their logic against this repo's
actual current outputs, writing to the gitignored `data/intermediate/`
staging area (the same convention `notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb`
uses for its `data/intermediate/figures/figure_2/` outputs) so they're
regenerated on demand rather than committed as stale snapshots.

## What gets rebuilt

| Output | Contents | Source |
|---|---|---|
| `Figure3a.csv.gz` | per-gene ClinVar classification counts, GenCC gene-disease validity, UniProt protein length → possible SNVs, clinical test counts, IGVF flag | `data/output/maves/integrated_variant_effect_dataset.tsv.gz` joined against the three files under `data/input/genes/` (see [`docs/data.md`](data.md)) and `Supplementary_Data_3.xlsx`'s `Curation` sheet (`IGVF Produced?`) |
| `Figure3c.csv.gz` | unique-variant counts (a variant tested by several sibling datasets for the same gene counts once, globally), tagged SGE/Vamp-seq/IGVF | `data/output/maves/integrated_variant_effect_dataset.tsv.gz` joined against `data/input/maves/Supplementary_Data_3.xlsx`'s `Curation` sheet (`Assay Name`, `IGVF Produced?`, `Primary Score Set or Meta-analysis?`) |
| `Figure3c_measurements.csv.gz` | the same, but per-dataset: a variant tested by N sibling datasets counts N times | same as `Figure3c.csv.gz` |
| `Figure3d.csv.gz` | ClinVar control (Benign/Likely benign/Pathogenic/Likely pathogenic and combined labels), gnomAD, unreported-SNV, and VUS counts | `data/output/maves/integrated_variant_effect_dataset.tsv.gz` alone |

### `Figure3c` vs `Figure3c_measurements`

Several genes have multiple non-meta-analysis datasets that assay largely
the same variant library under different conditions or antibody tags -- e.g.
F9's 5 `Popp_2025` datasets (one per epitope tag) all assay ~9,700 of the
same variants; CBS's two `Sun_2020` selection conditions and CARD11's two
`Meitlis_2020` conditions are similar. `Figure3c.csv.gz` counts each variant
once across all of a gene's sibling datasets (a true count of distinct
variants -- the default/current methodology); `Figure3c_measurements.csv.gz`
counts it once *per dataset* that tested it, so the same variant tested by 5
sibling datasets contributes 5x to the gene's total (a count of
measurements, not distinct variants -- this matches the methodology behind
the committed manuscript snapshot, `git show
d4d7770:Main_Figures/Figure_3/Figure3c.csv.gz`, which predates this script).
Both files record which scope produced them in a `variant_dedup_scope`
column ("global"/"per-dataset"), which `curation_summary_figure3.Rmd` uses
to pick each rendered bar chart's y-axis title ("Total unique variants" vs
"Total Variant effect measurements") automatically.

`Figure3a` and `Figure3c`'s IGVF flags, and `Figure3c`'s SGE/Vamp-seq/
meta-analysis flags, are all joined from `Supplementary_Data_3.xlsx`'s
`Curation` sheet rather than ported as hardcoded per-dataset/per-gene sets
(the original notebook's approach) — confirmed the curation sheet's
`Assay Name`/`IGVF Produced?`/`Primary Score Set or Meta-analysis?` columns
reproduce the old hardcoded sets exactly for every dataset/gene that was
already curated when those sets were written, so nothing changes for old
data; new datasets/genes are now classified automatically instead of
silently defaulting to "Other"/"No".

`Figure3d`'s priority genes (`BRCA1`, `PTEN`, `MSH2`, `TP53`) use each
variant's 2018 ClinVar snapshot (`clinvar_sig_2018`) rather than the current
one; every other gene uses `clinvar_sig_2025`. This matches the original
notebook and isn't something this script tries to change.

## The three `data/input/genes/` references `Figure3a` needs

- **`gencc-submissions.csv.gz`** — GenCC gene-disease validity submissions
  (legacy UUID-based export). Filtered to Definitive/Strong/Moderate
  classifications, deduplicated by (gene, disease). ⚠️ This legacy CSV format
  is scheduled for removal by GenCC on 2026-09-30; after that, re-fetching
  needs the newer SGC-ID-based export and a matching update to this script's
  GenCC column handling (`gene_symbol`, `classification_title`,
  `disease_curie` are the only columns actually read).
- **`uniprotkb_9606_reviewed.tsv.gz`** — reviewed (Swiss-Prot) human
  proteome. `Length * 9` approximates the possible single-nucleotide
  missense/nonsense changes per gene. Deduplicated on `Gene Names (primary)`
  before joining, since a small number of gene symbols (not any gene
  currently in this pipeline, but a real possibility for a future one) have
  multiple reviewed UniProt entries (isoforms/paralogs sharing a primary gene
  name).
- **`test_condition_gene.txt.gz`** — NCBI Genetic Testing Registry's public
  bulk export. Filtered to `test_type == "Clinical"` and `object == "gene"`
  rows, counted per gene for `gene_test_count`.

See [`docs/data.md`](data.md) for exact download dates and licenses, and
[`README.md`](../README.md#third-party-data--licenses) for the attribution
these licenses require.

### Refreshing these files

All three are committed snapshots (unlike `Figure3a/c/d.csv.gz` themselves,
which are gitignored and always rebuilt) since they change slowly and the
pipeline needs *a* consistent snapshot of each to build from. To pull a fresh
copy of each:

```bash
# GenCC (legacy UUID-based export -- see the removal warning above)
curl -sL "https://thegencc.org/download/action/submissions-export-csv" \
  | gzip > data/input/genes/gencc-submissions.csv.gz

# UniProt reviewed human proteome
curl -sL "https://rest.uniprot.org/uniprotkb/stream?query=%28reviewed%3Atrue%29+AND+%28organism_id%3A9606%29&format=tsv&fields=accession,reviewed,id,protein_name,gene_names,gene_primary,organism_name,length" \
  | gzip > data/input/genes/uniprotkb_9606_reviewed.tsv.gz

# NCBI Genetic Testing Registry bulk export
curl -sL "https://ftp.ncbi.nlm.nih.gov/pub/GTR/data/test_condition_gene.txt" \
  | gzip > data/input/genes/test_condition_gene.txt.gz
```

Update the download date in `docs/data.md` after refreshing. The `gene_primary`
field in the UniProt query is required — it's what this script joins on
(`Gene Names (primary)` in the downloaded TSV); the plain `gene_names` field
alone is a space-separated list of all names (primary + synonyms) and isn't a
usable join key.

## Known discrepancies vs. the original figure

Rebuilding from the current pipeline output (~6.5 months of pipeline changes
since the committed files were generated, plus a live GenCC/UniProt/GTR
snapshot rather than whatever was pulled originally) changes contents, not
just freshness:

- `Figure3a` gains a 41st gene, `LDLR` (absent from the previously-committed
  snapshot).
- `Figure3c` gains 15 datasets absent from the committed snapshot, including
  `LDLR`'s 3 datasets (`LDLR_Tabet_2025_abundance`,
  `LDLR_Tabet_2025_presence_VLDL`, `LDLR_Tabet_2025_uptake`) plus
  `BRCA2_Huang_2025_SGE`, `CHEK2_McCarthy-Leo_2024`, `PALB2_Boonen_2026`
  (and `_SGE`), `TP53_Funk_2025`, and 7 more.
- `Figure3d`'s control/gnomAD/VUS counts are all higher than the committed
  snapshot (e.g. `Pathogenic` 2,318 → 2,828), consistent with more variants
  having flowed through the pipeline since February.

## Usage

```bash
poetry run python -m src.build_figure3_data
```

Optional flags: `--integrated-dataset`, `--curation-sheet`, `--gencc`,
`--uniprot`, `--testing-registry`, `--figure3a-output`, `--figure3c-output`,
`--figure3c-measurements-output`, `--figure3d-output` (the last four default
under `data/intermediate/figures/figure_3/`).

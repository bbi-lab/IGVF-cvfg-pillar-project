# Smoke Tests

Two scripts let you verify your environment setup against tiny, real-data
fixtures instead of the full pipeline's multi-hour/multi-GB run. Neither
touches your real `data/intermediate/`, `data/output/`, or
`data/output/supplementary_data/` -- each runs against an isolated scratch
directory that's removed on success (left in place on failure, for
debugging).

Fixtures for both live under `tests/fixtures/smoke_test/` -- see "Fixture
provenance" below for how they were built and why they're shaped the way
they are.

## `scripts/smoke_test_variant_annotation.sh` -- Stage 1

```bash
scripts/smoke_test_variant_annotation.sh [--keep]
```

Runs the real, Dockerized `scripts/run_variant_annotation_pipeline.sh` (all
21 steps) against `tests/fixtures/smoke_test/cvfg_variants.0.smoke.tsv` (15
real `BARD1` rows) instead of the full ~550k-row
`data/input/maves/cvfg_variants.0.tsv`.

This is the **full pipeline**, just with a tiny input -- it still requires
every prerequisite in
[`docs/variant_annotation_pipeline.md`](variant_annotation_pipeline.md):
Docker/Compose, a `variant-annotation` checkout with its large reference
caches already downloaded (SpliceAI, dbNSFP, `clinvar_cache/`, etc.), and the
gnomAD Hail table cache prepared once via `--prepare-gnomad-cache`. Set those
up first; this script only checks that the pipeline wiring itself (staging,
`VARIANT_DATA_DIR` mounts, this project's own Dockerized steps) works, and
does so much faster than a real run since there's almost nothing to process.

It runs against scratch staging/output directories rather than
`data/intermediate/variant_annotation/` and `data/output/maves/`, via three
env var overrides added to `run_variant_annotation_pipeline.sh` for this
purpose (`CVFG_STAGE_DIR`, `CVFG_OUTPUT_DIR`, `CVFG_VARIANTS_0_FILE` -- see
that script's own header comment). Pass `--keep` to keep the scratch
directory after a successful run too (a failing run always leaves it in
place, with its path printed, for debugging).

### Known risk: this needs live network access to `ftp.ncbi.nlm.nih.gov`

Every fixture row is a "Case 1" variant in `map_variants.py` (a raw
nucleotide HGVS string given directly, e.g. `NM_000465.4:c.*43T>G`), which
triggers a ClinGen-normalization attempt via `dcd_mapping`/`cool_seq_tool`.
Importing `cool_seq_tool.paths` does something surprising: unless a
`MANE_SUMMARY_PATH` env var is already set, it unconditionally opens a live
FTP connection to `ftp.ncbi.nlm.nih.gov` to check/fetch the MANE summary --
every single invocation, with no offline fallback. `map_variants.py`'s
`_try_import_dcd_mapping` only catches `ImportError`, so any other failure
here (e.g. `socket.gaierror` from a transient DNS/network hiccup) propagates
uncaught and crashes the whole Step 1 container.

This isn't specific to the smoke test's fixture -- any real MAVE dataset
with at least one directly-provided nucleotide HGVS variant hits the same
code path in a real pipeline run. It's just guaranteed to be *exercised* by
this fixture (all 15 rows are Case 1), so a transient network issue at the
moment Step 1 runs will surface here reliably. `scripts/variant_annotation_pipeline.sh`
now runs under `set -euo pipefail`, so a crash here aborts the whole run
loudly instead of silently continuing into later steps against a
truncated/empty `cvfg_variants.1.tsv` (the previous, harder-to-diagnose
behavior). As of this writing there's no
env-var pre-seeding in place to avoid the live FTP dependency -- if Step 1
fails with a `socket.gaierror`/`ftplib` traceback, first check that the
Docker container running `map-variants` actually has outbound network
access (e.g. `docker compose --profile tools run --rm --entrypoint sh
map-variants -c "python3 -c \"import socket;
print(socket.gethostbyname('ftp.ncbi.nlm.nih.gov'))\""` from the
`variant-annotation` checkout) and retry; it's very likely transient.

## `scripts/smoke_test_analysis.sh` -- Stage 2

```bash
scripts/smoke_test_analysis.sh [--keep]
```

Runs the full Stage 2 sequence documented in the repo README's "2. Data
analysis / variant classification and table preparation" section --
`src.load_excalibr_calibrations`, `OddsPath_calculations.ipynb`,
`src.load_oddspath_calibrations`, `Variant_Classification_analysis.ipynb`,
`OddsPath_classifications.ipynb`, `src.build_variant_reclassification_dataset`,
in that order -- against
`tests/fixtures/smoke_test/integrated_variant_effect_dataset.smoke.tsv.gz`
(14 real rows from two real datasets, `BARD1_IGVF` and `G6PD_IGVF`) instead
of the real `data/output/maves/integrated_variant_effect_dataset.tsv.gz`.

Unlike Stage 1, this needs no Docker -- just the Poetry environment and the
`igvf-cvfg-pillar-project` Jupyter kernel (see the repo README's
"Environment: Poetry + Ruff" section for how to register it; the script
checks for it up front and prints the registration command if it's
missing). It validates that notebook execution itself works end to end:
dependencies, the kernel, and the whole calibration -> classification ->
reclassification chain.

All three notebooks resolve their own data directory from a `PROJECT_ROOT`
env var (default `../..`, i.e. two directories up from
`notebooks/analysis/`; see the first cell of each notebook). This script
points `PROJECT_ROOT` at a scratch directory containing a symlink to this
project's own `src/` (so `from src.lib... import ...` still resolves) plus
the fixture files, so the real notebooks run unmodified against isolated
data. `data/input/maves/Supplementary_Data_3.xlsx` is committed to the repo
and used as-is (real dataset metadata) rather than faked, since the
fixture's `Dataset` values are real entries in it.

## Fixture provenance

`tests/fixtures/smoke_test/` holds four files, all derived from real
project data rather than hand-written, so they exercise the pipeline's
actual column schema and value conventions rather than a synthetic
approximation of it:

- **`cvfg_variants.0.smoke.tsv`** -- the first 15 `BARD1` rows (by row
  order) sliced verbatim from `data/input/maves/cvfg_variants.0.tsv`, same
  19 columns.
- **`integrated_variant_effect_dataset.smoke.tsv.gz`** -- 14 rows sliced
  from `data/output/maves/integrated_variant_effect_dataset.tsv.gz`, same
  95 columns. 12 rows are from `BARD1_IGVF` (nucleotide-resolution),
  chosen to span several ClinVar significances (Pathogenic, Likely
  pathogenic, Pathogenic/Likely pathogenic, Uncertain significance,
  Conflicting classifications of pathogenicity, Likely benign, Benign, and
  no ClinVar record) and consequence types (missense, synonymous,
  stop-gained, splice-region, intron, 3' UTR, in-frame deletion); 2 rows
  are from `G6PD_IGVF` (protein-resolution, `nucleotide_or_aa == "aa"`) --
  `BARD1_IGVF` alone is entirely nucleotide-resolution, and
  `Variant_Classification_analysis.ipynb`'s protein-resolution branch
  (`sankey_aa`) raises `ValueError: cannot set a frame with no defined
  index and a scalar` if that subset is completely empty. Both datasets are
  real entries in `data/input/maves/Supplementary_Data_3.xlsx`'s `Curation`
  sheet (`BARD1_IGVF` has `Score Intervals Reported? = Reported` and
  `Functional Classification Provided? = Yes`, which
  `OddsPath_calculations.ipynb` requires of at least one dataset in its
  input to avoid an empty-groupby `KeyError`).
- **`CHEK2_Gebbia_2024.smoke.xlsx`** -- header-only stub (0 rows, same 9
  columns as the real `data/input/maves/CHEK2_Gebbia_2024.xlsx`).
  `OddsPath_calculations.ipynb` and `Variant_Classification_analysis.ipynb`
  both read this file unconditionally; since the fixture has no CHEK2 rows,
  its CHEK2-specific merge/flagging logic simply matches nothing.
- **`Supplementary_Data_4.smoke.xlsx`** -- header-only stub (0 rows in each
  sheet) of the 5 data sheets `Variant_Classification_analysis.ipynb` reads
  (`ExCALIBR_calibrations`, `REVEL_gene_specific_calibration`,
  `MP2_gene_specific_calibrations`, `AM_gene_specific_calibrations`,
  `OddsPath_calibrations`) and `src.load_excalibr_calibrations` /
  `src.load_oddspath_calibrations` overwrite in place.
  `data/output/supplementary_data/` is gitignored, so nothing else builds
  this workbook from a fresh checkout -- gene-specific calibration lookups
  simply find nothing for `BARD1`/`G6PD`, the same as any other gene without
  a gene-specific calibration in the real data.

If you need to regenerate these fixtures (e.g. after a schema change to the
integrated dataset), rebuild them from the real files the same way -- slice
real rows rather than hand-authoring values, since several of the crash
risks above only surface when a value falls outside a hardcoded set of
recognized real-data strings (e.g. `summarize_clnstar`'s ClinVar
review-status check in `notebooks/analysis/Variant_Classification_analysis.ipynb`
and `notebooks/analysis/OddsPath_classifications.ipynb`).

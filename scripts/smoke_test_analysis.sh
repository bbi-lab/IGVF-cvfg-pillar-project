#!/usr/bin/env bash
set -euo pipefail

########################################################################################################################
# Smoke test for Stage 2 (data analysis / variant classification) -- see
# README.md's "2. Data analysis / variant classification and table
# preparation" section, whose exact script/notebook order this mirrors.
#
# Runs the full Stage 2 sequence -- load-excalibr-calibrations,
# OddsPath_calculations.ipynb, load-oddspath-calibrations,
# Variant_Classification_analysis.ipynb, OddsPath_classifications.ipynb,
# build-variant-reclassification-dataset -- against a tiny fixture
# (tests/fixtures/smoke_test/integrated_variant_effect_dataset.smoke.tsv.gz,
# 14 real rows from two real datasets: BARD1_IGVF and G6PD_IGVF) instead of
# the real data/output/maves/integrated_variant_effect_dataset.tsv.gz. This
# validates the Docker environment setup (image build, the "notebooks"
# Poetry extra, notebook execution) end to end -- Stage 2 now runs the same
# way Stage 1 and the R figures already do, no local Poetry/Jupyter needed.
#
# All three notebooks resolve their own data directory from a PROJECT_ROOT
# env var (default "../.."; see the first cell of each notebook). This
# script points PROJECT_ROOT at an isolated scratch directory (via
# src/scripts/run_notebook.sh's --env passthrough) containing a symlink to
# this project's own src/ (so `from src.lib... import ...` still resolves)
# plus the fixture files below, so the smoke test never reads or writes
# your real data/output/ or data/input/maves/CHEK2_Gebbia_2024.xlsx -- safe
# to run alongside real Stage 2 output you already have.
#
# The scratch directory lives under data/intermediate/ (gitignored) rather
# than /tmp: none of the Stage 2 Docker services (analysis-notebooks,
# load-excalibr-calibrations, load-oddspath-calibrations,
# build-variant-reclassification-dataset) mount anything but this repo's own
# tree, so scratch data has to live inside the repo for them to see it --
# same reasoning as scripts/smoke_test_variant_annotation.sh's scratch
# directory and Step 16's --output-dir fix in
# scripts/variant_annotation_pipeline.sh.
#
# Fixtures used (see tests/fixtures/smoke_test/):
#   - integrated_variant_effect_dataset.smoke.tsv.gz: the Stage 2 input.
#   - CHEK2_Gebbia_2024.smoke.xlsx: header-only stub -- both notebooks that
#     touch CHEK2 read this file unconditionally, but the fixture doesn't
#     exercise CHEK2-specific logic (no CHEK2 rows in the fixture).
#   - Supplementary_Data_4.smoke.xlsx: header-only stub of the workbook
#     Variant_Classification_analysis.ipynb reads and
#     load-oddspath-calibrations/load-excalibr-calibrations update in
#     place -- data/output/supplementary_data/ is gitignored, so nothing
#     builds this workbook from scratch otherwise.
# data/input/maves/Supplementary_Data_3.xlsx is committed to the repo and
# used as-is (real dataset metadata -- the fixture's Dataset values,
# BARD1_IGVF and G6PD_IGVF, are real entries in it).
#
# Usage:
#   scripts/smoke_test_analysis.sh [--keep]
#
# --keep leaves the scratch directory in place afterwards for inspection
# (normally removed on success, left in place on failure regardless of this
# flag).
########################################################################################################################

keep=0
if [[ "${1:-}" == "--keep" ]]; then
  keep=1
elif [[ $# -gt 0 ]]; then
  echo "error: unrecognized argument '$1'" >&2
  echo "usage: $0 [--keep]" >&2
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd "$script_dir/.." && pwd)"
fixture_dir="$project_dir/tests/fixtures/smoke_test"

for f in integrated_variant_effect_dataset.smoke.tsv.gz CHEK2_Gebbia_2024.smoke.xlsx Supplementary_Data_4.smoke.xlsx; do
  if [[ ! -f "$fixture_dir/$f" ]]; then
    echo "error: smoke-test fixture not found at $fixture_dir/$f" >&2
    exit 1
  fi
done

scratch_root="$project_dir/data/intermediate"
mkdir -p "$scratch_root"
scratch_dir="$(mktemp -d "$scratch_root/cvfg_smoke_analysis.XXXXXX")"
rel="${scratch_dir#"$project_dir"/}"
container_project_root="/usr/src/app/$rel"

cleanup() {
  local exit_code=$?
  if [[ "$keep" -eq 1 || $exit_code -ne 0 ]]; then
    echo "Scratch directory left at: $scratch_dir" >&2
  else
    rm -rf "$scratch_dir"
  fi
}
trap cleanup EXIT

echo "Staging fixtures into scratch PROJECT_ROOT: $scratch_dir ..."
# Target is the *container* path (/usr/src/app/src), not a host path: this
# symlink is only ever resolved from inside the analysis-notebooks
# container, where /usr/src/app is the repo root bind mount.
ln -s /usr/src/app/src "$scratch_dir/src"
mkdir -p \
  "$scratch_dir/data/input/maves" \
  "$scratch_dir/data/output/maves" \
  "$scratch_dir/data/output/supplementary_data" \
  "$scratch_dir/data/input/mave_calibration/excalibr/json" \
  "$scratch_dir/executed"
cp "$project_dir/data/input/maves/Supplementary_Data_3.xlsx" "$scratch_dir/data/input/maves/Supplementary_Data_3.xlsx"
cp "$fixture_dir/CHEK2_Gebbia_2024.smoke.xlsx" "$scratch_dir/data/input/maves/CHEK2_Gebbia_2024.xlsx"
cp "$fixture_dir/integrated_variant_effect_dataset.smoke.tsv.gz" "$scratch_dir/data/output/maves/integrated_variant_effect_dataset.tsv.gz"
cp "$fixture_dir/Supplementary_Data_4.smoke.xlsx" "$scratch_dir/data/output/supplementary_data/Supplementary_Data_4.xlsx"

run_notebook() {
  local name="$1"
  echo
  echo "Running notebooks/analysis/${name}.ipynb ..."
  "$project_dir/src/scripts/run_notebook.sh" \
    --env "PROJECT_ROOT=$container_project_root" \
    --to notebook --execute \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=300 \
    --output-dir "$container_project_root/executed" \
    --output "executed_${name}.ipynb" \
    "notebooks/analysis/${name}.ipynb"
}

cd "$project_dir"

echo "Refreshing ExCALIBR_calibrations sheet (empty fixture JSON dir -- expected 0 rows) ..."
src/scripts/run_load_excalibr_calibrations.sh \
  "$rel/data/input/mave_calibration/excalibr/json" \
  "$rel/data/output/supplementary_data/Supplementary_Data_4.xlsx"

run_notebook "OddsPath_calculations"

odds_path_csv="$scratch_dir/data/output/mave_calibration/OddsPath_calibrations.csv.gz"
if [[ ! -s "$odds_path_csv" ]]; then
  echo "FAIL: expected output missing or empty: $odds_path_csv" >&2
  exit 1
fi

echo
echo "Refreshing OddsPath_calibrations sheet ..."
src/scripts/run_load_oddspath_calibrations.sh \
  "$rel/data/output/mave_calibration/OddsPath_calibrations.csv.gz" \
  "$rel/data/output/supplementary_data/Supplementary_Data_4.xlsx"

run_notebook "Variant_Classification_analysis"

checkpoint="$scratch_dir/data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz"
sd5="$scratch_dir/data/output/supplementary_data/Supplementary_Data_5.xlsx"
for f in "$checkpoint" "$sd5"; do
  if [[ ! -s "$f" ]]; then
    echo "FAIL: expected output missing or empty: $f" >&2
    exit 1
  fi
done

run_notebook "OddsPath_classifications"

sd6="$scratch_dir/data/output/supplementary_data/Supplementary_Data_6.xlsx"
if [[ ! -s "$sd6" ]]; then
  echo "FAIL: expected output missing or empty: $sd6" >&2
  exit 1
fi

echo
echo "Building the biobank-analysis reclassification export ..."
reclassification_output_rel="$rel/data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz"
src/scripts/run_build_variant_reclassification_dataset.sh \
  "$rel/data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz" \
  --chek2-file "$rel/data/input/maves/CHEK2_Gebbia_2024.xlsx" \
  --output "$reclassification_output_rel"

reclassification_output="$scratch_dir/data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz"
if [[ ! -s "$reclassification_output" ]]; then
  echo "FAIL: expected output missing or empty: $reclassification_output" >&2
  exit 1
fi

rows=$(($(gzip -dc "$reclassification_output" | wc -l) - 1))
if [[ "$rows" -lt 1 ]]; then
  echo "FAIL: expected at least 1 row in $reclassification_output, got $rows" >&2
  exit 1
fi

echo
echo "PASS: Stage 2 analysis smoke test succeeded ($rows row(s) in the final reclassification export)."

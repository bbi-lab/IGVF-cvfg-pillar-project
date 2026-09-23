#!/usr/bin/env bash
set -euo pipefail

########################################################################################################################
# Full end-to-end smoke test for scripts/run_variant_annotation_pipeline.sh.
#
# Runs the real Dockerized variant-annotation pipeline (all 21 steps) against
# a tiny fixture (tests/fixtures/smoke_test/cvfg_variants.0.smoke.tsv -- 15
# real BARD1 rows sliced from data/input/maves/cvfg_variants.0.tsv) instead
# of the full ~550k-row cvfg_variants.0.tsv. This exercises the whole
# environment setup end-to-end -- Docker, the variant-annotation submodule/
# VARIANT_ANNOTATION_DIR checkout, VEP/ClinVar/gnomAD/SpliceAI/predictor
# annotation, this project's own Dockerized steps -- fast, since the input is
# tiny. It does NOT shortcut any of the heavy prerequisites documented in
# docs/variant_annotation_pipeline.md: Docker/Compose, a variant-annotation
# checkout with its large reference caches already downloaded, and the
# gnomAD Hail table cache prepared via --prepare-gnomad-cache. Set those up
# first (see docs/variant_annotation_pipeline.md and the repo README); this
# script only checks that the pipeline wiring itself works end to end.
#
# Runs against scratch staging/output directories (via the CVFG_STAGE_DIR/
# CVFG_OUTPUT_DIR/CVFG_VARIANTS_0_FILE overrides added to
# run_variant_annotation_pipeline.sh for this purpose), so it never touches
# your real data/intermediate/variant_annotation/ or data/output/maves/ --
# safe to run alongside (or in between) real pipeline runs. The scratch
# directory lives under data/intermediate/ (gitignored) rather than /tmp:
# some Dockerized steps (e.g. build-training-variant-files, invoked from
# Step 16) only bind-mount this repo's own tree, not an arbitrary
# VARIANT_DATA_DIR path, so scratch data has to live inside the repo for
# those steps to see it -- see Step 16's --output-dir handling in
# scripts/variant_annotation_pipeline.sh.
#
# Usage:
#   scripts/smoke_test_variant_annotation.sh [--keep]
#
# --keep leaves the scratch staging/output directories in place afterwards
# for inspection (normally removed on success, left in place on failure
# regardless of this flag).
#
# Env vars: VARIANT_ANNOTATION_DIR, same meaning as for
# run_variant_annotation_pipeline.sh -- see docs/variant_annotation_pipeline.md.
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
fixture="$project_dir/tests/fixtures/smoke_test/cvfg_variants.0.smoke.tsv"

if [[ ! -f "$fixture" ]]; then
  echo "error: smoke-test fixture not found at $fixture" >&2
  exit 1
fi

scratch_root="$project_dir/data/intermediate"
mkdir -p "$scratch_root"
scratch_dir="$(mktemp -d "$scratch_root/cvfg_smoke_variant_annotation.XXXXXX")"
stage_dir="$scratch_dir/intermediate"
output_dir="$scratch_dir/output"
mkdir -p "$stage_dir" "$output_dir"

cleanup() {
  local exit_code=$?
  if [[ "$keep" -eq 1 || $exit_code -ne 0 ]]; then
    echo "Scratch staging/output directories left at: $scratch_dir" >&2
  else
    rm -rf "$scratch_dir"
  fi
}
trap cleanup EXIT

echo "Running scripts/run_variant_annotation_pipeline.sh against $fixture"
echo "  (staging: $stage_dir, output: $output_dir) ..."
echo

CVFG_STAGE_DIR="$stage_dir" \
CVFG_OUTPUT_DIR="$output_dir" \
CVFG_VARIANTS_0_FILE="$fixture" \
  "$script_dir/run_variant_annotation_pipeline.sh"

echo
echo "Checking final outputs ..."
expanded="$output_dir/integrated_variant_effect_dataset.tsv.gz"
condensed="$output_dir/integrated_variant_effect_dataset.condensed.tsv.gz"

for f in "$expanded" "$condensed"; do
  if [[ ! -s "$f" ]]; then
    echo "FAIL: expected output missing or empty: $f" >&2
    exit 1
  fi
done

expanded_rows=$(($(gzip -dc "$expanded" | wc -l) - 1))
condensed_rows=$(($(gzip -dc "$condensed" | wc -l) - 1))
if [[ "$expanded_rows" -lt 1 || "$condensed_rows" -lt 1 ]]; then
  echo "FAIL: expected at least 1 output row, got $expanded_rows expanded / $condensed_rows condensed" >&2
  exit 1
fi

echo
echo "PASS: variant-annotation pipeline smoke test succeeded."
echo "  $expanded_rows expanded row(s), $condensed_rows condensed row(s)."

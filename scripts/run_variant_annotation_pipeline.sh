#!/usr/bin/env bash
set -euo pipefail

########################################################################################################################
# Runs scripts/variant_annotation_pipeline.sh (the Dockerized variant-annotation
# pipeline steps) against this project's own data, without writing into the
# variant-annotation checkout's own working tree.
#
# Data flow (see data/raw_mave_data/README.md and
# docs/variant_annotation_pipeline.md for the why):
#   data/raw_mave_data/  --(staged, non-destructively)-->
#   data/intermediate/variant_annotation/data/  --(mounted as VARIANT_DATA_DIR)-->
#   variant_annotation_pipeline.sh Steps 1-21  -->
#   data/mave_data/*.tsv.gz
#
# Usage:
#   scripts/run_variant_annotation_pipeline.sh                     # full pipeline
#   scripts/run_variant_annotation_pipeline.sh --step N            # just step N
#   scripts/run_variant_annotation_pipeline.sh --prepare-gnomad-cache
#
# --step N still stages data/raw_mave_data/ and sets up VARIANT_DATA_DIR/
# CVFG_PROJECT_DIR exactly as a full run does, but runs only step N of
# scripts/variant_annotation_pipeline.sh (N reads
# data/cvfg_variants.<N-1>.tsv and writes data/cvfg_variants.<N>.tsv, except
# step 21, which writes the final integrated MAVE dataset files directly; see
# that script's header for the full step list) and skips the final gzip into
# data/mave_data/ either way, since even --step 21 only writes its outputs to
# the staged intermediate directory.
#
# --prepare-gnomad-cache builds/refreshes Step 7's local gnomAD Hail table
# cache and exits -- it skips all of the CVFG-specific data staging below,
# since that data is irrelevant to a cache build. This is a prerequisite for
# Step 7, not something a normal run needs to do -- see the "gnomAD Hail
# table cache" section of docs/variant_annotation_pipeline.md for when this
# is (and isn't) necessary.
#
# Env vars:
#   VARIANT_ANNOTATION_DIR  Path to a variant-annotation checkout. Defaults to
#                           the vendored submodule at vendor/variant-annotation.
#                           Set this to an existing checkout that already has
#                           the large reference files / caches downloaded
#                           (AlphaMissense, REVEL, SpliceAI, clinvar_cache,
#                           etc.) to avoid re-fetching multi-GB data.
########################################################################################################################

step=""
prepare_gnomad_cache=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --step)
      step="${2:?--step requires a step number}"
      shift 2
      ;;
    --prepare-gnomad-cache)
      prepare_gnomad_cache=1
      shift
      ;;
    *)
      echo "error: unrecognized argument '$1'" >&2
      echo "usage: $0 [--step N | --prepare-gnomad-cache]" >&2
      exit 1
      ;;
  esac
done

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export CVFG_PROJECT_DIR
CVFG_PROJECT_DIR="$(cd "$script_dir/.." && pwd)"

va_dir="${VARIANT_ANNOTATION_DIR:-$CVFG_PROJECT_DIR/vendor/variant-annotation}"
if [[ ! -d "$va_dir" || ! -f "$va_dir/compose.yaml" ]]; then
  cat >&2 <<EOF
error: no variant-annotation checkout found at $va_dir

Set VARIANT_ANNOTATION_DIR to an existing checkout (recommended if you
already have the large reference files downloaded there), or fetch the
vendored submodule with:
  git submodule update --init vendor/variant-annotation
EOF
  exit 1
fi

stage_dir="$CVFG_PROJECT_DIR/data/intermediate/variant_annotation"
mkdir -p "$stage_dir/data"

if [[ "$prepare_gnomad_cache" -eq 1 ]]; then
  export VARIANT_DATA_DIR="$stage_dir"
  echo "Running prepare-gnomad-cache of scripts/variant_annotation_pipeline.sh in $va_dir (VARIANT_DATA_DIR=$VARIANT_DATA_DIR) ..."
  ( cd "$va_dir" && bash "$CVFG_PROJECT_DIR/scripts/variant_annotation_pipeline.sh" prepare-gnomad-cache )
  echo "Done: gnomAD Hail table cache built/refreshed in the variant-annotation-gnomad-cache Docker volume."
  exit 0
fi

echo "Staging data/raw_mave_data/ -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/ ..."
rsync -a "$CVFG_PROJECT_DIR/data/raw_mave_data/" "$stage_dir/data/"

# score_sets.tsv and Supplementary_Data_3.xlsx live in data/input/maves/
# (not data/raw_mave_data/) so they're staged here explicitly rather than
# picked up by the rsync above. Steps 11/15 and Step 16 of
# scripts/variant_annotation_pipeline.sh read them as
# /work/data/score_sets.tsv and /work/data/Supplementary_Data_3.xlsx -- see
# docs/variant_annotation_pipeline.md.
echo "Staging data/input/maves/ -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/ ..."
cp "$CVFG_PROJECT_DIR/data/input/maves/score_sets.tsv" "$stage_dir/data/score_sets.tsv"
cp "$CVFG_PROJECT_DIR/data/input/maves/Supplementary_Data_3.xlsx" "$stage_dir/data/Supplementary_Data_3.xlsx"

# MANE.GRCh38.v1.5.summary.txt.gz lives in data/input/reference/ (not
# data/raw_mave_data/) so it's staged here explicitly too. Step 2 of
# scripts/variant_annotation_pipeline.sh reads it as
# /work/data/MANE.GRCh38.v1.5.summary.txt.gz -- see
# docs/variant_annotation_pipeline.md.
echo "Staging data/input/reference/ -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/ ..."
cp "$CVFG_PROJECT_DIR/data/input/reference/MANE.GRCh38.v1.5.summary.txt.gz" "$stage_dir/data/MANE.GRCh38.v1.5.summary.txt.gz"

# cvfg_variants.0.tsv (Step 1's input) only needs staging when Step 1 is
# actually about to run -- a full run, or --step 1. Later --step N runs read
# an already-staged cvfg_variants.<N-1>.tsv and shouldn't require this file.
if [[ -z "$step" || "$step" == "1" ]]; then
  echo "Staging data/input/maves/cvfg_variants.0.tsv -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/cvfg_variants.0.tsv ..."
  cp "$CVFG_PROJECT_DIR/data/input/maves/cvfg_variants.0.tsv" "$stage_dir/data/cvfg_variants.0.tsv"
fi

# data_frame_missense_variants_MP2_properties.csv.gz (Step 13's MutPred2
# scores input) is ~350MB, so it's only staged when Step 13 is actually about
# to run -- a full run, or --step 13 -- same gating as cvfg_variants.0.tsv
# above. Step 13 reads it as /work/data/data_frame_missense_variants_MP2_properties.csv.gz;
# see docs/variant_annotation_pipeline.md.
if [[ -z "$step" || "$step" == "13" ]]; then
  echo "Staging data/input/predictors/data_frame_missense_variants_MP2_properties.csv.gz -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/ ..."
  cp "$CVFG_PROJECT_DIR/data/input/predictors/data_frame_missense_variants_MP2_properties.csv.gz" "$stage_dir/data/data_frame_missense_variants_MP2_properties.csv.gz"
fi

export VARIANT_DATA_DIR="$stage_dir"

if [[ -n "$step" ]]; then
  echo "Running step $step of scripts/variant_annotation_pipeline.sh in $va_dir (VARIANT_DATA_DIR=$VARIANT_DATA_DIR) ..."
  ( cd "$va_dir" && bash "$CVFG_PROJECT_DIR/scripts/variant_annotation_pipeline.sh" "$step" )
  if [[ "$step" == "21" ]]; then
    echo "Done:"
    echo "  $stage_dir/data/integrated_variant_effect_dataset.tsv"
    echo "  $stage_dir/data/integrated_variant_effect_dataset.condensed.tsv"
  else
    echo "Done: $stage_dir/data/cvfg_variants.$step.tsv"
  fi
  exit 0
fi

echo "Running scripts/variant_annotation_pipeline.sh in $va_dir (VARIANT_DATA_DIR=$VARIANT_DATA_DIR) ..."
( cd "$va_dir" && bash "$CVFG_PROJECT_DIR/scripts/variant_annotation_pipeline.sh" )

expanded="$stage_dir/data/integrated_variant_effect_dataset.tsv"
condensed="$stage_dir/data/integrated_variant_effect_dataset.condensed.tsv"
out_dir="$CVFG_PROJECT_DIR/data/mave_data"
mkdir -p "$out_dir"

echo "Gzipping final outputs into data/mave_data/ ..."
gzip -kf -c "$expanded" > "$out_dir/integrated_variant_effect_dataset.tsv.gz"
gzip -kf -c "$condensed" > "$out_dir/integrated_variant_effect_dataset.condensed.tsv.gz"

echo "Done:"
echo "  $out_dir/integrated_variant_effect_dataset.tsv.gz"
echo "  $out_dir/integrated_variant_effect_dataset.condensed.tsv.gz"

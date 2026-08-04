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
#   variant_annotation_pipeline.sh Steps 1-17 + condensed/expanded frame assembly  -->
#   data/mave_data/*.tsv.gz
#
# Env vars:
#   VARIANT_ANNOTATION_DIR  Path to a variant-annotation checkout. Defaults to
#                           the vendored submodule at vendor/variant-annotation.
#                           Set this to an existing checkout that already has
#                           the large reference files / caches downloaded
#                           (AlphaMissense, REVEL, SpliceAI, MANE, clinvar_cache,
#                           etc.) to avoid re-fetching multi-GB data.
########################################################################################################################

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

echo "Staging data/raw_mave_data/ -> ${stage_dir#"$CVFG_PROJECT_DIR"/}/data/ ..."
rsync -a "$CVFG_PROJECT_DIR/data/raw_mave_data/" "$stage_dir/data/"

export VARIANT_DATA_DIR="$stage_dir"

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

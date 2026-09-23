#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_build_variant_reclassification_dataset.sh [checkpoint-file] [--chek2-file path] [--output path] [--dedup] [flags]

Examples:
  src/scripts/run_build_variant_reclassification_dataset.sh
  src/scripts/run_build_variant_reclassification_dataset.sh --dedup

Notes:
  - Defaults to
    data/output/reclassification/integrated_variant_effect_dataset_analysis.csv.gz
    (Variant_Classification_analysis.ipynb's own checkpoint) and
    data/input/maves/CHEK2_Gebbia_2024.xlsx (both bind-mounted at
    /usr/src/app since this service only reads/writes this repo's own
    tree, unlike flag-variants).
  - Paths are interpreted relative to /usr/src/app in the container.
  - Writes data/output/reclassification/integrated_variant_effect_biobank_input_data.tsv.gz
    by default.
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF

compose_build_flag=""
compose_no_cache_flag=""
args=()

for arg in "$@"; do
  case "$arg" in
    --rebuild-image)
      compose_build_flag="--build"
      ;;
    --no-build-cache)
      compose_no_cache_flag="--no-cache"
      ;;
    *)
      args+=("$arg")
      ;;
  esac
done

cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
cmd+=(--rm build-variant-reclassification-dataset)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

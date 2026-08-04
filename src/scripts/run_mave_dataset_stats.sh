#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_mave_dataset_stats.sh [condensed-file] [metadata-file] [--output path] [flags]

Examples:
  src/scripts/run_mave_dataset_stats.sh
  src/scripts/run_mave_dataset_stats.sh --output data/mave_data/mave_dataset_stats.csv

Notes:
  - Defaults to data/mave_data/integrated_variant_effect_dataset.condensed.tsv.gz
    and data/mave_data/Supplementary_Data_3.xlsx (both bind-mounted at
    /usr/src/app since this service only reads the repo's checked-in
    data/mave_data/ files, unlike flag-variants).
  - Paths are interpreted relative to /usr/src/app in the container.
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
cmd+=(--rm mave-dataset-stats)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

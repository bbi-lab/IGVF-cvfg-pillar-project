#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_build_training_variant_files.sh [--input-dir dir] [--output-dir dir] [flags]

Examples:
  src/scripts/run_build_training_variant_files.sh
  src/scripts/run_build_training_variant_files.sh \
    --input-dir data/input/predictors \
    --output-dir data/intermediate/variant_annotation/data

Notes:
  - Defaults to data/input/predictors and
    data/intermediate/variant_annotation/data, both bind-mounted at
    /usr/src/app since this service only reads/writes this project's own
    tree, unlike flag-variants/recalculate-clingen-classification.
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
cmd+=(--rm build-training-variant-files)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi

# cd to the project root (where compose.yaml lives) so this still works when
# called with a different cwd, e.g. from scripts/variant_annotation_pipeline.sh's
# preparatory step before step_13, which runs from a variant-annotation checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir/../.."
exec "${cmd[@]}"

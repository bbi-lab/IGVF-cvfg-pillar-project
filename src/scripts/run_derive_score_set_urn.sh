#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_derive_score_set_urn.sh <input-file> <output-file> [derive_score_set_urn options...]

Examples:
  src/scripts/run_derive_score_set_urn.sh data/cvfg/v13/cvfg_variants.13.tsv data/cvfg/v13/cvfg_variants.14.temp.tsv
  src/scripts/run_derive_score_set_urn.sh input.tsv output.tsv --variant-urn-column variant_urn

Notes:
  - Paths are interpreted relative to /work in the container, unless they
    already exist in this project's working tree, in which case they're
    interpreted relative to /usr/src/app (the whole-repo bind mount).
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path -- set this
    to the same value used for the variant-annotation pipeline steps so both
    containers see the same staged data.
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF
  exit 1
fi

input_path="$1"
output_path="$2"
shift 2

compose_build_flag=""
compose_no_cache_flag=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rebuild-image)
      compose_build_flag="--build"
      shift
      ;;
    --no-build-cache)
      compose_no_cache_flag="--no-cache"
      shift
      ;;
    --)
      shift
      break
      ;;
    *)
      break
      ;;
  esac
done

map_to_container_path() {
  local path="$1"

  if [[ "$path" == /work/* || "$path" == /usr/src/app/* ]]; then
    printf '%s\n' "$path"
    return
  fi

  if [[ -f "$path" ]]; then
    printf '/usr/src/app/%s\n' "$path"
  else
    printf '/work/%s\n' "$path"
  fi
}

input_in_container="$(map_to_container_path "$input_path")"

if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
  output_in_container="$output_path"
elif [[ "$input_in_container" == /usr/src/app/* ]]; then
  output_in_container="/usr/src/app/$output_path"
else
  output_in_container="/work/$output_path"
fi

cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
cmd+=(--rm derive-score-set-urn "$input_in_container" "$output_in_container")
if [[ $# -gt 0 ]]; then
  cmd+=("$@")
fi

# cd to the project root (where compose.yaml lives) so this still works when
# called with a different cwd, e.g. from scripts/variant_annotation_pipeline.sh's
# step_14, which runs from a variant-annotation checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir/../.."
exec "${cmd[@]}"

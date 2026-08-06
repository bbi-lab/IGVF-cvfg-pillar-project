#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_annotate_simplified_consequence.sh <input-file> <output-file> [annotate_simplified_consequence options...]

Examples:
  src/scripts/run_annotate_simplified_consequence.sh data/cvfg_variants.16.tsv data/cvfg_variants.17.tsv
  src/scripts/run_annotate_simplified_consequence.sh input.tsv output.tsv --consequence-col vep.most_severe_mutational_consequence

Notes:
  - Paths are interpreted relative to /work in the container, unless they
    already exist in this project's working tree, in which case they're
    interpreted relative to /usr/src/app (the whole-repo bind mount).
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path -- set this
    to the same value used for the variant-annotation pipeline steps so both
    containers see the same staged data.
  - --consequence-map-file defaults to
    data/input/consequence/extended_ensembl_consequence.csv.gz, which
    resolves against this project's own bind mount (/usr/src/app), not /work.
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

# Remap a --consequence-map-file path argument to a container path, same as
# the positional input/output paths above.
mapped_args=()
skip_next=0
for arg in "$@"; do
  if [[ "$skip_next" -eq 1 ]]; then
    mapped_args+=("$(map_to_container_path "$arg")")
    skip_next=0
  elif [[ "$arg" == --consequence-map-file ]]; then
    mapped_args+=("$arg")
    skip_next=1
  else
    mapped_args+=("$arg")
  fi
done

cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
cmd+=(--rm annotate-simplified-consequence "$input_in_container" "$output_in_container")
if [[ ${#mapped_args[@]} -gt 0 ]]; then
  cmd+=("${mapped_args[@]}")
fi

# cd to the project root (where compose.yaml lives) so this still works when
# called with a different cwd, e.g. from scripts/variant_annotation_pipeline.sh's
# step_17, which runs from a variant-annotation checkout.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir/../.."
exec "${cmd[@]}"

#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_make_extended_data_figure_7alt_op.sh [flags]

OddsPath-calibrated variant of run_make_extended_data_figure_7alt.sh -- see
docs/figures.md.

Examples:
  src/scripts/run_make_extended_data_figure_7alt_op.sh
  src/scripts/run_make_extended_data_figure_7alt_op.sh --consequence-filter missense

Notes:
  - All CLI args default to committed/generated repo files (bind-mounted at
    /usr/src/app since this service only reads/writes this repo's own
    tree, unlike flag-variants).
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
cmd+=(--rm make-extended-data-figure-7alt-op)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

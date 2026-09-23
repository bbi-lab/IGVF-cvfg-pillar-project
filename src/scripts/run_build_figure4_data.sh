#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_build_figure4_data.sh --cached-json path [flags]

Examples:
  src/scripts/run_build_figure4_data.sh \
    --cached-json notebooks/figures/figure_4/old_figure4_data.json.gz

Notes:
  - --cached-json has no usable default: three sub-panels can't be
    regenerated from scratch, only carried forward from a prior run (see
    docs/build_figure4_data.md), so it must be supplied explicitly.
  - Paths are interpreted relative to /usr/src/app in the container
    (bind-mounted there since this service only reads/writes this repo's
    own tree, unlike flag-variants).
  - Writes notebooks/figures/figure_4/figure4_data.json.gz by default.
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
cmd+=(--rm build-figure4-data)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

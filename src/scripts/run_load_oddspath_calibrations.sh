#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_load_oddspath_calibrations.sh [csv] [workbook] [flags]

Examples:
  src/scripts/run_load_oddspath_calibrations.sh
  src/scripts/run_load_oddspath_calibrations.sh \
    data/output/mave_calibration/OddsPath_calibrations.csv.gz \
    data/output/supplementary_data/Supplementary_Data_4.xlsx

Notes:
  - Defaults to
    data/output/mave_calibration/OddsPath_calibrations.csv.gz
    and data/output/supplementary_data/Supplementary_Data_4.xlsx (both
    bind-mounted at /usr/src/app since this service only reads/writes
    data/output/ files, unlike flag-variants).
  - Paths are interpreted relative to /usr/src/app in the container.
  - Overwrites the OddsPath_calibrations sheet of the workbook in place; all
    other sheets are left untouched.
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
cmd+=(--rm load-oddspath-calibrations)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

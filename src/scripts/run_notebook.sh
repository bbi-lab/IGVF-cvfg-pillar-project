#!/usr/bin/env bash
set -euo pipefail

cat <<'EOF' >&2
Usage: src/scripts/run_notebook.sh [nbconvert args...] <notebook-path>

Generic runner for the Stage 2 analysis notebooks and Stage 3 Python figure
notebooks, via the analysis-notebooks service (entrypoint: jupyter nbconvert).

Examples:
  src/scripts/run_notebook.sh --to notebook --execute \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=600 \
    --output executed_OddsPath_calculations.ipynb \
    notebooks/analysis/OddsPath_calculations.ipynb

  src/scripts/run_notebook.sh --to notebook --execute \
    --ExecutePreprocessor.kernel_name=python3 \
    --ExecutePreprocessor.timeout=600 \
    --output executed_PP_ProcessBigDataFrame.ipynb \
    notebooks/figures/figure_2/PP_ProcessBigDataFrame.ipynb

Notes:
  - The notebook path and all nbconvert args are passed through verbatim --
    this wrapper doesn't default or remap anything.
  - Paths are interpreted relative to /usr/src/app in the container, i.e.
    repo-relative, same as running `poetry run jupyter nbconvert ...` from
    the repo root.
  - No --workdir override is needed (unlike r-figures): nbconvert's
    ExecutePreprocessor already defaults the kernel's cwd to the notebook's
    own directory, regardless of the container's own cwd.
  - kernel_name=python3 is ipykernel's default kernel, not a registered
    project-specific one -- nothing in the repo depends on a particular
    kernel name.
  - Pass --env KEY=VALUE (repeatable, wrapper-only, must come before any
    nbconvert args) to set a container env var, e.g. --env
    PROJECT_ROOT=/usr/src/app/data/intermediate/some/scratch/dir to point a
    notebook's PROJECT_ROOT override at an isolated directory instead of
    the real repo tree (see scripts/smoke_test_analysis.sh).
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF

compose_build_flag=""
compose_no_cache_flag=""
env_flags=()
args=()

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
    --env)
      env_flags+=(-e "${2:?--env requires a KEY=VALUE argument}")
      shift 2
      ;;
    --env=*)
      env_flags+=(-e "${1#--env=}")
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
if [[ ${#env_flags[@]} -gt 0 ]]; then
  cmd+=("${env_flags[@]}")
fi
cmd+=(--rm analysis-notebooks)
if [[ ${#args[@]} -gt 0 ]]; then
  cmd+=("${args[@]}")
fi
exec "${cmd[@]}"

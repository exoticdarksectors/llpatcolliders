#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

PYTHON_BIN="${HNL_PYTHON:-python}"
PRODUCTION_WORKERS="${HNL_PRODUCTION_WORKERS:-6}"
PROMPT_TAU_CORES="${HNL_PROMPT_TAU_CORES:-12}"
WZ_JOBS="${HNL_WZ_JOBS:-8}"
WZ_CORES_PER_JOB="${HNL_WZ_CORES_PER_JOB:-1}"
ANALYSIS_WORKERS="${HNL_ANALYSIS_WORKERS:-3}"
TEMPLATE_PYTHON="${HNL_TEMPLATE_PYTHON:-}"

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "ERROR: Python executable not found: $PYTHON_BIN" >&2
  echo "Activate the hnl environment or set HNL_PYTHON." >&2
  exit 1
fi

export HNL_RUN_TAG="${HNL_RUN_TAG:-full_$(date +%Y%m%d_%H%M%S)}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-$PWD/tmp/.mplconfig}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-$PWD/tmp/.cache}"
export PYTHONUNBUFFERED=1

mkdir -p "$MPLCONFIGDIR" "$XDG_CACHE_HOME"

"$PYTHON_BIN" -c \
  "import matplotlib, mpmath, numba, numpy, pandas, particle, rtree, scipy, sympy, tqdm, trimesh"

printf 'HNL full run\n'
printf '  Python: %s\n' "$("$PYTHON_BIN" -c 'import sys; print(sys.executable)')"
printf '  Run tag: %s\n' "$HNL_RUN_TAG"
printf '  Production workers: %s\n' "$PRODUCTION_WORKERS"
printf '  Prompt-tau MG5 cores: %s\n' "$PROMPT_TAU_CORES"
printf '  W/Z shards: %s x %s core(s)\n' "$WZ_JOBS" "$WZ_CORES_PER_JOB"
printf '  Analysis workers: %s\n\n' "$ANALYSIS_WORKERS"

if [[ -z "$TEMPLATE_PYTHON" ]] && "$PYTHON_BIN" -c \
  'import ROOT; assert hasattr(ROOT, "TPythia8")' >/dev/null 2>&1; then
  TEMPLATE_PYTHON="$PYTHON_BIN"
fi

if [[ -n "$TEMPLATE_PYTHON" ]]; then
  if ! command -v "$TEMPLATE_PYTHON" >/dev/null 2>&1; then
    echo "ERROR: template Python executable not found: $TEMPLATE_PYTHON" >&2
    exit 1
  fi
  "$TEMPLATE_PYTHON" -c \
    'import numpy, scipy, ROOT; assert hasattr(ROOT, "TPythia8")'
  printf '  Template Python: %s\n\n' \
    "$("$TEMPLATE_PYTHON" -c 'import sys; print(sys.executable)')"
else
  printf '%s\n\n' \
    '  Template generation: skipped (pre-generate templates or set HNL_TEMPLATE_PYTHON)'
fi

"$PYTHON_BIN" -u run_all.py \
  --flavor Ue Umu Utau \
  --no-wz \
  --skip-combine \
  --workers "$PRODUCTION_WORKERS" \
  --prompt-tau-nb-core "$PROMPT_TAU_CORES" \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_prod_nonwz.log"

for flavor in Ue Umu Utau; do
  "$PYTHON_BIN" -u production/madgraph/run_wz_sharded.py \
    --flavor "$flavor" \
    --jobs "$WZ_JOBS" \
    --nb-core "$WZ_CORES_PER_JOB" \
    --skip-existing \
    2>&1 | tee "tmp/${HNL_RUN_TAG}_prod_wz_${flavor}.log"
done

"$PYTHON_BIN" -m production.combine_channels --flavor Ue Umu Utau \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_combine.log"

# Rest-frame templates require a ROOT/Pythia-enabled interpreter, which is
# normally distinct from the production/analysis conda environment.
if [[ -n "$TEMPLATE_PYTHON" ]]; then
  "$TEMPLATE_PYTHON" -u -m analysis.generate_decay_templates \
    --flavor Ue Umu Utau \
    --skip-existing \
    2>&1 | tee "tmp/${HNL_RUN_TAG}_templates.log"
fi

"$PYTHON_BIN" -u run_analysis.py \
  --flavor Ue Umu Utau \
  --workers "$ANALYSIS_WORKERS" \
  2>&1 | tee "tmp/${HNL_RUN_TAG}_analysis.log"

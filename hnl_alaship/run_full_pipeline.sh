#!/usr/bin/env bash
set -euo pipefail

# ──────────────────────────────────────────────────────────────────────
# GARGOYLE HNL full pipeline — optimized parallel execution
#
# Machine budget: 12 CPUs, 18 GB RAM, 2 GB swap
#
# Tier 1 (concurrent, ~2-4h bottleneck = WZ):
#   - generate_meson_csvs  (1 CPU, ~300 MB)
#   - run_tau_production    (1-2 CPU via Docker, ~500 MB)
#   - run_wz_production ×3  (3 CPU each via Docker --nb-core 3, ~800 MB each)
#   - generate_ctau_tables  (1 CPU, ~200 MB)
#   Total: 12 CPU, ~3.5 GB peak
#
# Tier 2 (after Tier 1): combine_channels (~30s)
# Tier 3 (after Tier 2): run_sensitivity --workers 3 (~45-90 min)
# ──────────────────────────────────────────────────────────────────────

cd "$(dirname "$0")"
PROJECT_DIR="$(pwd)"
CONDA_ENV="${CONDA_ENV:-llpatcolliders}"
LOG_DIR="$PROJECT_DIR/output/logs"
STATUS_FILE="$PROJECT_DIR/output/pipeline_status.json"

mkdir -p "$LOG_DIR"

# ── helper: update status JSON ───────────────────────────────────────
update_status() {
    local step="$1" status="$2"
    local ts
    ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    echo "[$ts] $step: $status" >> "$LOG_DIR/pipeline.log"
    python3 -c "
import json
path = '$STATUS_FILE'
try:
    data = json.load(open(path))
except Exception:
    data = {'steps': {}, 'start_ts': '$ts'}
data['steps']['$step'] = {'status': '$status', 'ts': '$ts'}
data['last_update'] = '$ts'
json.dump(data, open(path, 'w'), indent=2)
"
}

# ── helper: launch a step in the background, store PID ──────────────
declare -A PIDS
launch_step() {
    local label="$1" logfile="$2"; shift 2
    update_status "$label" "running"
    (
        conda run -n "$CONDA_ENV" python "$@" > "$LOG_DIR/$logfile" 2>&1
        echo "EXIT=$?" >> "$LOG_DIR/$logfile"
    ) &
    PIDS[$label]=$!
    printf "  %-15s PID=%5d  -> %s\n" "$label" "${PIDS[$label]}" "$LOG_DIR/$logfile"
}

# ── helper: run a step serially, exit on failure ────────────────────
run_step() {
    local label="$1" logfile="$2"; shift 2
    update_status "$label" "running"
    if conda run -n "$CONDA_ENV" python "$@" > "$LOG_DIR/$logfile" 2>&1; then
        update_status "$label" "done"
        echo "  ✓ $label completed"
    else
        update_status "$label" "FAILED"
        echo "  ✗ $label FAILED — see $LOG_DIR/$logfile"
        return 1
    fi
}

echo "================================================================"
echo "  GARGOYLE HNL Full Pipeline — Parallel Execution"
echo "  $(date)"
echo "  CPUs: $(sysctl -n hw.ncpu)  RAM: $(( $(sysctl -n hw.memsize) / 1073741824 )) GB"
echo "================================================================"
echo ""

PIPELINE_START=$(date +%s)

# ══════════════════════════════════════════════════════════════════════
# TIER 1: Launch all production steps in parallel
# ══════════════════════════════════════════════════════════════════════

echo "[TIER 1] Launching 6 parallel production processes..."
echo ""

launch_step meson_csvs     step1_meson_csvs.log      production/decay_engine/generate_meson_csvs.py --flavor Ue Umu Utau --channel all
launch_step tau_production step2_tau_production.log  production/madgraph/run_tau_production.py      --flavor Ue Umu Utau
launch_step wz_Ue          step3_wz_Ue.log           production/madgraph/run_wz_production.py       --flavor Ue   --nb-core 3
launch_step wz_Umu         step3_wz_Umu.log          production/madgraph/run_wz_production.py       --flavor Umu  --nb-core 3
launch_step wz_Utau        step3_wz_Utau.log         production/madgraph/run_wz_production.py       --flavor Utau --nb-core 3
launch_step ctau_tables    step5_ctau.log            production/generate_ctau_tables.py             --flavor Ue Umu Utau

echo ""
echo "  All ${#PIDS[@]} processes launched. Waiting for completion..."
echo "  Monitor: tail -f $LOG_DIR/pipeline.log"
echo ""

# ── Wait for all Tier 1 processes ────────────────────────────────────
TIER1_OK=true
for label in "${!PIDS[@]}"; do
    pid="${PIDS[$label]}"
    if wait "$pid"; then
        update_status "$label" "done"
        echo "  ✓ $label (PID $pid) completed successfully"
    else
        update_status "$label" "FAILED"
        echo "  ✗ $label (PID $pid) FAILED — check $LOG_DIR/"
        TIER1_OK=false
    fi
done

echo ""
echo "  Tier 1 elapsed: $(( $(date +%s) - PIPELINE_START ))s"

if [ "$TIER1_OK" = false ]; then
    echo ""
    echo "WARNING: Some Tier 1 steps failed. Proceeding with available data..."
fi

# ══════════════════════════════════════════════════════════════════════
# TIER 2: Combine channels
# ══════════════════════════════════════════════════════════════════════

echo ""
echo "[TIER 2] Combining channels..."
run_step combine_channels step4_combine.log production/combine_channels.py --flavor Ue Umu Utau || exit 1

# ══════════════════════════════════════════════════════════════════════
# TIER 3: Sensitivity scan (3 workers, ~9 CPUs, ~4.5 GB)
# ══════════════════════════════════════════════════════════════════════

echo ""
echo "[TIER 3] Running sensitivity scan (3 workers, all 3 flavors)..."
update_status "sensitivity" "running"

conda run -n "$CONDA_ENV" python analysis/run_sensitivity.py \
    --flavor Ue Umu Utau --workers 3 \
    2>&1 | tee "$LOG_DIR/step6_sensitivity.log"

if [ ${PIPESTATUS[0]} -eq 0 ]; then
    update_status "sensitivity" "done"
    echo "  ✓ sensitivity scan completed"
else
    update_status "sensitivity" "FAILED"
    echo "  ✗ sensitivity scan FAILED"
fi

# ══════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════

TOTAL_S=$(( $(date +%s) - PIPELINE_START ))
echo ""
echo "================================================================"
echo "  Pipeline complete in ${TOTAL_S}s ($(( TOTAL_S / 60 ))m $(( TOTAL_S % 60 ))s)"
echo "  Results: $PROJECT_DIR/output/analysis/"
echo "  Status:  $STATUS_FILE"
echo "================================================================"
update_status "pipeline" "complete:${TOTAL_S}s"

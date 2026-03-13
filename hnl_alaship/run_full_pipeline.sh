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
CONDA_ENV="llpatcolliders"
LOG_DIR="$PROJECT_DIR/output/logs"
STATUS_FILE="$PROJECT_DIR/output/pipeline_status.json"

mkdir -p "$LOG_DIR"

# ── helper: update status JSON ───────────────────────────────────────
update_status() {
    local step="$1" status="$2"
    local ts
    ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    # Append to a simple log too
    echo "[$ts] $step: $status" >> "$LOG_DIR/pipeline.log"
    # Write machine-readable status
    python3 -c "
import json, os, time
path = '$STATUS_FILE'
try:
    data = json.load(open(path))
except:
    data = {'steps': {}, 'start_ts': '$ts'}
data['steps']['$step'] = {'status': '$status', 'ts': '$ts'}
data['last_update'] = '$ts'
json.dump(data, open(path, 'w'), indent=2)
"
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

# --- Step 1: Generate meson CSVs (FONLL) ---
update_status "meson_csvs" "running"
(
    conda run -n "$CONDA_ENV" python production/decay_engine/generate_meson_csvs.py \
        --flavor Ue Umu Utau --channel all \
        > "$LOG_DIR/step1_meson_csvs.log" 2>&1
    echo "STEP1_EXIT=$?" >> "$LOG_DIR/step1_meson_csvs.log"
) &
PID_MESON=$!
echo "  Step 1 [meson CSVs]:     PID=$PID_MESON  → $LOG_DIR/step1_meson_csvs.log"

# --- Step 2: Tau production (Docker + decay) ---
update_status "tau_production" "running"
(
    conda run -n "$CONDA_ENV" python production/madgraph/run_tau_production.py \
        --flavor Ue Umu Utau \
        > "$LOG_DIR/step2_tau_production.log" 2>&1
    echo "STEP2_EXIT=$?" >> "$LOG_DIR/step2_tau_production.log"
) &
PID_TAU=$!
echo "  Step 2 [tau production]: PID=$PID_TAU  → $LOG_DIR/step2_tau_production.log"

# --- Step 3a/b/c: WZ production — 3 flavors in parallel, 3 cores each ---
update_status "wz_Ue" "running"
(
    conda run -n "$CONDA_ENV" python production/madgraph/run_wz_production.py \
        --flavor Ue --nb-core 3 \
        > "$LOG_DIR/step3_wz_Ue.log" 2>&1
    echo "STEP3A_EXIT=$?" >> "$LOG_DIR/step3_wz_Ue.log"
) &
PID_WZ_UE=$!
echo "  Step 3a [WZ Ue]:         PID=$PID_WZ_UE  → $LOG_DIR/step3_wz_Ue.log"

update_status "wz_Umu" "running"
(
    conda run -n "$CONDA_ENV" python production/madgraph/run_wz_production.py \
        --flavor Umu --nb-core 3 \
        > "$LOG_DIR/step3_wz_Umu.log" 2>&1
    echo "STEP3B_EXIT=$?" >> "$LOG_DIR/step3_wz_Umu.log"
) &
PID_WZ_UMU=$!
echo "  Step 3b [WZ Umu]:        PID=$PID_WZ_UMU → $LOG_DIR/step3_wz_Umu.log"

update_status "wz_Utau" "running"
(
    conda run -n "$CONDA_ENV" python production/madgraph/run_wz_production.py \
        --flavor Utau --nb-core 3 \
        > "$LOG_DIR/step3_wz_Utau.log" 2>&1
    echo "STEP3C_EXIT=$?" >> "$LOG_DIR/step3_wz_Utau.log"
) &
PID_WZ_UTAU=$!
echo "  Step 3c [WZ Utau]:       PID=$PID_WZ_UTAU → $LOG_DIR/step3_wz_Utau.log"

# --- Step 5: cτ tables (independent of production) ---
update_status "ctau_tables" "running"
(
    conda run -n "$CONDA_ENV" python production/generate_ctau_tables.py \
        --flavor Ue Umu Utau \
        > "$LOG_DIR/step5_ctau.log" 2>&1
    echo "STEP5_EXIT=$?" >> "$LOG_DIR/step5_ctau.log"
) &
PID_CTAU=$!
echo "  Step 5 [cτ tables]:      PID=$PID_CTAU  → $LOG_DIR/step5_ctau.log"

echo ""
echo "  All 6 processes launched. Waiting for completion..."
echo "  Monitor: tail -f $LOG_DIR/pipeline.log"
echo "  Dashboard: watch -n5 'cat $STATUS_FILE 2>/dev/null; echo; for f in $LOG_DIR/step*.log; do echo \"=== \$(basename \$f) ===\"; tail -1 \$f; done'"
echo ""

# ── Wait for all Tier 1 processes ────────────────────────────────────
TIER1_OK=true

for label_pid in "meson_csvs:$PID_MESON" "tau_production:$PID_TAU" \
                 "wz_Ue:$PID_WZ_UE" "wz_Umu:$PID_WZ_UMU" "wz_Utau:$PID_WZ_UTAU" \
                 "ctau_tables:$PID_CTAU"; do
    label="${label_pid%%:*}"
    pid="${label_pid##*:}"
    if wait "$pid"; then
        update_status "$label" "done"
        echo "  ✓ $label (PID $pid) completed successfully"
    else
        update_status "$label" "FAILED"
        echo "  ✗ $label (PID $pid) FAILED — check $LOG_DIR/"
        TIER1_OK=false
    fi
done

TIER1_END=$(date +%s)
echo ""
echo "  Tier 1 elapsed: $(( TIER1_END - PIPELINE_START ))s"

if [ "$TIER1_OK" = false ]; then
    echo ""
    echo "WARNING: Some Tier 1 steps failed. Proceeding with available data..."
fi

# ══════════════════════════════════════════════════════════════════════
# TIER 2: Combine channels
# ══════════════════════════════════════════════════════════════════════

echo ""
echo "[TIER 2] Combining channels..."
update_status "combine_channels" "running"

conda run -n "$CONDA_ENV" python production/combine_channels.py \
    --flavor Ue Umu Utau \
    > "$LOG_DIR/step4_combine.log" 2>&1

if [ $? -eq 0 ]; then
    update_status "combine_channels" "done"
    echo "  ✓ combine_channels completed"
else
    update_status "combine_channels" "FAILED"
    echo "  ✗ combine_channels FAILED"
    exit 1
fi

# ══════════════════════════════════════════════════════════════════════
# TIER 3: Sensitivity scan (3 workers, ~9 CPUs, ~4.5 GB)
# ══════════════════════════════════════════════════════════════════════

echo ""
echo "[TIER 3] Running sensitivity scan (3 workers, all 3 flavors)..."
update_status "sensitivity" "running"

conda run -n "$CONDA_ENV" python analysis/run_sensitivity.py \
    --flavor Ue Umu Utau \
    --workers 3 \
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

PIPELINE_END=$(date +%s)
TOTAL_S=$(( PIPELINE_END - PIPELINE_START ))
echo ""
echo "================================================================"
echo "  Pipeline complete in ${TOTAL_S}s ($(( TOTAL_S / 60 ))m $(( TOTAL_S % 60 ))s)"
echo "  Results: $PROJECT_DIR/output/analysis/"
echo "  Status:  $STATUS_FILE"
echo "================================================================"
update_status "pipeline" "complete:${TOTAL_S}s"

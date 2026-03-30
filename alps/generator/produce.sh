#!/usr/bin/env bash
set -euo pipefail

# Usage: bash produce.sh <cmnd_file> <target_events> <output_name> [n_jobs] [events_per_job]
#
# Event-targeted production driver for ALP samples.
#
# - Stops when accumulated generated events (from batch _meta.json) reaches target_events.
# - Keeps parallel batching for throughput.
# - Merges LLP CSVs with event-id renumbering.
#
# Notes:
# - For h -> aa runs with 25:onIfAny = 6000113, each generated event should
#   carry two LLPs (up to generator failures).
# - BR normalization is done in analysis; generation stays at BR=1.
#
# Examples:
#   bash produce.sh heavy_alp.cmnd 10000 alp_heavy_m15
#   bash produce.sh heavy_alp.cmnd 10000 alp_heavy_m15 4
#   bash produce.sh light_alp.cmnd 10000 alp_light_m1 4 50000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CMND="${1:?Usage: $0 <cmnd_file> <target_events> <output_name> [n_jobs] [events_per_job]}"
TARGET_EVENTS="${2:?}"
NAME="$(basename "${3:?}")"

# Host CPU info (best-effort; useful to interpret n_jobs setting).
SYSCTL_BIN="$(command -v sysctl || true)"
if [[ -x /usr/sbin/sysctl ]]; then
  SYSCTL_BIN="/usr/sbin/sysctl"
fi
if [[ -n "$SYSCTL_BIN" ]]; then
  PHYS_CPU="$("$SYSCTL_BIN" -n hw.physicalcpu 2>/dev/null || echo unknown)"
  LOGICAL_CPU="$("$SYSCTL_BIN" -n hw.logicalcpu 2>/dev/null || echo unknown)"
else
  PHYS_CPU="unknown"
  LOGICAL_CPU="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo unknown)"
fi

if [[ ! "$LOGICAL_CPU" =~ ^[0-9]+$ ]]; then
  LOGICAL_CPU="$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo unknown)"
fi

DEFAULT_JOBS=1
if [[ "$PHYS_CPU" =~ ^[0-9]+$ ]]; then
  DEFAULT_JOBS=$((PHYS_CPU - 1))
elif [[ "$LOGICAL_CPU" =~ ^[0-9]+$ ]]; then
  DEFAULT_JOBS=$((LOGICAL_CPU - 1))
fi
if (( DEFAULT_JOBS < 1 )); then
  DEFAULT_JOBS=1
fi

NJOBS="${4:-$DEFAULT_JOBS}"
EVENTS_PER_JOB="${5:-0}"

# Validate integer-like inputs.
for val_name in TARGET_EVENTS NJOBS EVENTS_PER_JOB; do
  val="${!val_name}"
  if [[ ! "$val" =~ ^[0-9]+$ ]]; then
    echo "ERROR: $val_name must be a non-negative integer (got '$val')." >&2
    exit 1
  fi
done
if (( TARGET_EVENTS <= 0 )); then
  echo "ERROR: target_events must be > 0." >&2
  exit 1
fi
if (( NJOBS <= 0 )); then
  echo "ERROR: n_jobs must be > 0." >&2
  exit 1
fi

# Resolve command file:
# - absolute path: use as-is
# - relative path: first try current working directory, then script directory
if [[ "$CMND" != /* ]]; then
  if [[ -f "$CMND" ]]; then
    CMND="$(cd "$(dirname "$CMND")" && pwd)/$(basename "$CMND")"
  else
    CMND="${SCRIPT_DIR}/${CMND}"
  fi
fi
if [[ ! -f "$CMND" ]]; then
  echo "ERROR: command file not found: $CMND" >&2
  exit 1
fi

OUTDIR="${SCRIPT_DIR}/../output"
DATADIR="${OUTDIR}/data"
IMAGEDIR="${OUTDIR}/images"
mkdir -p "$DATADIR" "$IMAGEDIR"
CSV="${DATADIR}/${NAME}.csv"
TMP="${DATADIR}/.${NAME}_part"

# Initialize output CSV with header
printf 'event,id,pt,eta,phi,momentum,mass\n' > "$CSV"

seed=0
offset=0
round=0
MAX_ROUNDS=500

total_gen=0
total_with_llp=0
total_llp=0
llp_pdg_id=""

if (( EVENTS_PER_JOB > 0 )); then
  SPLIT_DESC="$EVENTS_PER_JOB events/job"
else
  SPLIT_DESC="auto split events/job"
fi

echo "Target: $TARGET_EVENTS generated events | up to $NJOBS parallel jobs | $SPLIT_DESC"
echo "CPU context: physical=$PHYS_CPU logical=$LOGICAL_CPU | generator uses process-level parallelism (1 process/job)"
echo "Output root: $OUTDIR"
echo "Output LLP CSV: $CSV"

while (( total_gen < TARGET_EVENTS )); do
  round=$((round + 1))
  if (( round > MAX_ROUNDS )); then
    echo "ERROR: hit $MAX_ROUNDS rounds without reaching target_events. Aborting." >&2
    exit 1
  fi

  remaining=$((TARGET_EVENTS - total_gen))
  if (( EVENTS_PER_JOB > 0 )); then
    per_job=$EVENTS_PER_JOB
  else
    # Auto split: use up to NJOBS jobs this round.
    per_job=$(( (remaining + NJOBS - 1) / NJOBS ))
  fi
  (( per_job <= 0 )) && per_job=1

  n_launch=$(( (remaining + per_job - 1) / per_job ))
  (( n_launch > NJOBS )) && n_launch=$NJOBS

  echo "--- Round $round plan: remaining=$remaining, launch=$n_launch job(s), baseline events/job=$per_job ---"

  pids=()
  job_seeds=()
  req_list=()
  req_total=0

  for j in $(seq 1 "$n_launch"); do
    req=$per_job
    left=$((remaining - per_job * (j - 1)))
    (( left < req )) && req=$left
    (( req <= 0 )) && continue

    seed=$((seed + 1))
    job_seeds+=("$seed")
    req_list+=("$req")
    req_total=$((req_total + req))

    "${SCRIPT_DIR}/generator" -c "$CMND" -n "$req" -s "$seed" -o "${TMP}_${seed}" &
    pids+=($!)
  done

  echo "  Requested events per job this round: ${req_list[*]} (sum=$req_total)"

  failed=0
  for pid in "${pids[@]}"; do
    wait "$pid" || ((failed++)) || true
  done
  if (( failed > 0 )); then
    echo "  WARNING: $failed job(s) failed in this round."
  fi

  gen_before=$total_gen
  round_gen=0
  round_with_llp=0
  round_llp=0

  for idx in "${!job_seeds[@]}"; do
    s="${job_seeds[$idx]}"
    f="${TMP}_${s}.csv"

    if [[ -f "$f" ]]; then
      awk -F',[[:space:]]*' -v OFS=',' -v off="$offset" \
        'NR>1{$1=$1+off; print $1,$2,$3,$4,$5,$6,$7}' "$f" >> "$CSV"

      last_evt=$(awk -F',[[:space:]]*' 'NR>1{last=$1} END{print (NR>1 ? last : -1)+0}' "$f")
      offset=$((offset + last_evt + 1))

      # Only count generated events when the LLP CSV was produced successfully.
      mf="${TMP}_${s}_meta.json"
      if [[ -f "$mf" ]]; then
        ng=$(awk -F': *' '/"n_generated"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")
        nw=$(awk -F': *' '/"n_with_llp"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")
        nl=$(awk -F': *' '/"n_total_llp"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")

        total_gen=$((total_gen + ${ng:-0}))
        total_with_llp=$((total_with_llp + ${nw:-0}))
        total_llp=$((total_llp + ${nl:-0}))
        round_gen=$((round_gen + ${ng:-0}))
        round_with_llp=$((round_with_llp + ${nw:-0}))
        round_llp=$((round_llp + ${nl:-0}))

        if [[ -z "$llp_pdg_id" ]]; then
          llp_pdg_id=$(awk -F': *' '/"llp_pdg_id"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")
        fi
        rm -f "$mf"
      fi
    fi

    # Clean up temp files unconditionally (handles failed batches).
    rm -f "$f" "${TMP}_${s}.log" "${TMP}_${s}_meta.json"
  done

  if (( total_gen == gen_before )); then
    echo "ERROR: no new generated events collected in round $round. Check generator/cmd setup." >&2
    exit 1
  fi

  done_pct=$(awk -v a="$total_gen" -v b="$TARGET_EVENTS" 'BEGIN{printf "%.1f", 100.0*a/b}')
  echo "  Round $round result: +$round_gen generated, +$round_with_llp events with LLP, +$round_llp LLP rows"
  echo "  Progress: $total_gen / $TARGET_EVENTS generated events ($done_pct%)"
done

# Final counts from merged CSV
count_llp=$(tail -n +2 "$CSV" | wc -l | tr -d ' ')

if (( total_llp > 0 && count_llp != total_llp )); then
  echo "WARNING: merged LLP rows ($count_llp) != aggregated n_total_llp ($total_llp)."
fi

# Write aggregated metadata sidecar
META="${DATADIR}/${NAME}_meta.json"
cat > "$META" <<EOF
{
  "n_generated": $total_gen,
  "n_with_llp": $total_with_llp,
  "n_total_llp": $count_llp,
  "llp_pdg_id": ${llp_pdg_id:-null}
}
EOF

if (( total_gen > 0 )); then
  llp_per_event=$(awk -v a="$count_llp" -v b="$total_gen" 'BEGIN{printf "%.3f", a/b}')
else
  llp_per_event="nan"
fi

echo "Done: generated_events=$total_gen (target=$TARGET_EVENTS)"
echo "  LLP rows: $count_llp | events with LLP: $total_with_llp"
echo "  Effective LLP/event: $llp_per_event"
echo "  Metadata: $META"

#!/usr/bin/env bash
set -euo pipefail

# Usage: bash produce.sh <cmnd_file> <target_llps> <output_name> [n_jobs] [batch_size]
#
# Runs generator in parallel batches until the CSV has >= target_llps rows.
# Output CSV goes to ../output/<output_name>.csv
#
# Examples:
#   bash produce.sh heavy_alp.cmnd  10000 alp_heavy_m15
#   bash produce.sh light_alp.cmnd 10000 alp_light_m1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CMND="${1:?Usage: $0 <cmnd_file> <target_llps> <output_name> [n_jobs] [batch_size]}"
TARGET="${2:?}"
NAME="$(basename "${3:?}")"
NJOBS="${4:-$(( $(sysctl -n hw.physicalcpu) - 1 ))}"
BATCH="${5:-10000}"         # Pythia events per job per round

# Resolve cmnd file relative to script dir if not an absolute path
[[ "$CMND" != /* ]] && CMND="${SCRIPT_DIR}/${CMND}"

OUTDIR="${SCRIPT_DIR}/../output"
mkdir -p "$OUTDIR"
CSV="${OUTDIR}/${NAME}.csv"
TMP="${OUTDIR}/.${NAME}_part"

# Initialize output CSV with header
printf 'event,id,pt,eta,phi,momentum,mass\n' > "$CSV"

count=0
round=0
seed=0
offset=0
MAX_ROUNDS=200
consec_empty=0
total_gen=0
total_with_llp=0

echo "Target: $TARGET LLP rows | $NJOBS parallel jobs | $BATCH events/job/round"
echo "Output: $CSV"

while [ "$count" -lt "$TARGET" ]; do
    round=$((round + 1))
    if [ "$round" -gt "$MAX_ROUNDS" ]; then
        echo "ERROR: hit $MAX_ROUNDS rounds without reaching target. Aborting." >&2
        exit 1
    fi
    echo "--- Round $round (have $count / $TARGET LLPs) ---"

    # Launch parallel batch
    pids=()
    for i in $(seq 1 "$NJOBS"); do
        seed=$((seed + 1))
        "${SCRIPT_DIR}/generator" -c "$CMND" -n "$BATCH" -s "$seed" -o "${TMP}_${seed}" &
        pids+=($!)
    done

    # Wait
    failed=0
    for pid in "${pids[@]}"; do wait "$pid" || ((failed++)) || true; done
    if [ "$failed" -gt 0 ]; then echo "  WARNING: $failed job(s) failed"; fi

    # Append results with renumbered event IDs; aggregate metadata
    for i in $(seq $((seed - NJOBS + 1)) "$seed"); do
        f="${TMP}_${i}.csv"
        [ -f "$f" ] || continue
        awk -F',[[:space:]]*' -v OFS=',' -v off="$offset" \
            'NR>1{$1=$1+off; print $1,$2,$3,$4,$5,$6,$7}' "$f" >> "$CSV"
        last_evt=$(awk -F',[[:space:]]*' 'NR>1{last=$1} END{print (NR>1 ? last : -1)+0}' "$f")
        offset=$((offset + last_evt + 1))
        # Accumulate metadata from batch sidecar
        mf="${TMP}_${i}_meta.json"
        if [ -f "$mf" ]; then
            ng=$(awk -F': *' '/"n_generated"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")
            nw=$(awk -F': *' '/"n_with_llp"/{gsub(/[^0-9]/,"",$2); print $2}' "$mf")
            total_gen=$((total_gen + ${ng:-0}))
            total_with_llp=$((total_with_llp + ${nw:-0}))
            rm -f "$mf"
        fi
        rm -f "$f" "${TMP}_${i}.log"
    done

    prev_count=$count
    count=$(tail -n +2 "$CSV" | wc -l | tr -d ' ')
    if [ "$count" -eq "$prev_count" ]; then
        consec_empty=$((consec_empty + 1))
        if [ "$consec_empty" -ge 3 ]; then
            echo "ERROR: 3 consecutive rounds with zero new LLPs. Check cmnd file or binary." >&2
            exit 1
        fi
    else
        consec_empty=0
    fi
done

# Write aggregated metadata sidecar
META="${OUTDIR}/${NAME}_meta.json"
cat > "$META" <<EOF
{
  "n_generated": $total_gen,
  "n_with_llp": $total_with_llp,
  "n_total_llp": $count
}
EOF

echo "Done: $count LLP rows in $CSV (after $round round(s), $seed total Pythia jobs)"
echo "  Generated events: $total_gen | Events with LLP: $total_with_llp"
echo "  Metadata: $META"

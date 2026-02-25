#!/usr/bin/env bash
set -euo pipefail

# Usage: bash parallel_produce.sh <cmnd_file> <target_llps> <output_name> [n_jobs] [batch_size]
#
# Runs main144 in parallel batches until the CSV has >= target_llps rows.
# Output CSV goes to ../output/<output_name>.csv
#
# Examples:
#   bash parallel_produce.sh higgsLL.cmnd  10000 alp_heavy_m15
#   bash parallel_produce.sh alp_meson.cmnd 10000 alp_light_m1

CMND="${1:?Usage: $0 <cmnd_file> <target_llps> <output_name> [n_jobs] [batch_size]}"
TARGET="${2:?}"
NAME="$(basename "${3:?}")"
NJOBS="${4:-$(( $(sysctl -n hw.physicalcpu) - 1 ))}"
BATCH="${5:-10000}"         # Pythia events per job per round

OUTDIR="../output"
mkdir -p "$OUTDIR"
CSV="${OUTDIR}/${NAME}.csv"
TMP="${OUTDIR}/.${NAME}_part"

# Initialize output CSV with header
printf 'event,\tid,\tpt,\teta,\tphi,\tmomentum,\tmass\n' > "$CSV"

count=0
round=0
seed=0
offset=0

echo "Target: $TARGET LLP rows | $NJOBS parallel jobs | $BATCH events/job/round"
echo "Output: $CSV"

while [ "$count" -lt "$TARGET" ]; do
    round=$((round + 1))
    echo "--- Round $round (have $count / $TARGET LLPs) ---"

    # Launch parallel batch
    pids=()
    for i in $(seq 1 "$NJOBS"); do
        seed=$((seed + 1))
        ./main144 -c "$CMND" -n "$BATCH" -s "$seed" -o "${TMP}_${seed}" &
        pids+=($!)
    done

    # Wait
    failed=0
    for pid in "${pids[@]}"; do wait "$pid" || ((failed++)); done
    if [ "$failed" -gt 0 ]; then echo "  WARNING: $failed job(s) failed"; fi

    # Append results with renumbered event IDs
    for i in $(seq $((seed - NJOBS + 1)) "$seed"); do
        f="${TMP}_${i}.csv"
        [ -f "$f" ] || continue
        awk -F',\t' -v OFS=',\t' -v off="$offset" 'NR>1{$1=$1+off; print}' "$f" >> "$CSV"
        last_evt=$(awk -F',\t' 'END{print $1+0}' "$f")
        offset=$((offset + last_evt + 1))
        rm -f "$f" "${TMP}_${i}.log"
    done

    count=$(tail -n +2 "$CSV" | wc -l | tr -d ' ')
done

echo "Done: $count LLP rows in $CSV (after $round round(s), $seed total Pythia jobs)"

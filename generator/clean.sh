#!/usr/bin/env bash
set -euo pipefail

# Cleanup helper for produce.sh outputs.
#
# Usage:
#   bash clean.sh <output_name> [--outdir DIR]
#   bash clean.sh --all [--outdir DIR]
#
# Examples:
#   bash clean.sh alp_heavy_m15 --outdir output/alps
#   bash clean.sh --all --outdir output/dark_photon

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse --outdir
CUSTOM_OUTDIR=""
ARGS=()
for arg in "$@"; do
  if [[ "$arg" == "--outdir" ]]; then
    shift_next=1; continue
  fi
  if [[ "${shift_next:-}" == "1" ]]; then
    CUSTOM_OUTDIR="$arg"; shift_next=0; continue
  fi
  ARGS+=("$arg")
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

OUTDIR="${CUSTOM_OUTDIR:-${SCRIPT_DIR}/../output}"
DATADIR="${OUTDIR}/data"

usage() {
  cat <<'EOF'
Usage:
  bash clean.sh <output_name>
  bash clean.sh --all

Removes files produced by produce.sh:
  - output/data/<output_name>.csv
  - output/data/<output_name>_daughters.csv
  - output/data/<output_name>_meta.json
  - output/data/.<output_name>_part_*.csv
  - output/data/.<output_name>_part_*.log
  - output/data/.<output_name>_part_*.root

With --all, removes all matching production outputs in output/.
EOF
}

if [[ ! -d "$OUTDIR" ]]; then
  echo "Output directory not found: $OUTDIR"
  exit 1
fi
mkdir -p "$DATADIR"

if [[ $# -ne 1 ]]; then
  usage
  exit 1
fi

target="$1"
removed=0

remove_if_exists() {
  local f
  for f in "$@"; do
    if [[ -e "$f" ]]; then
      rm -f -- "$f"
      removed=$((removed + 1))
      echo "removed: $f"
    fi
  done
}

if [[ "$target" == "--all" ]]; then
  shopt -s nullglob
  all_files=(
    "$DATADIR"/*.csv
    "$DATADIR"/*_meta.json
    "$DATADIR"/.*_part_*.csv
    "$DATADIR"/.*_part_*.log
    "$DATADIR"/.*_part_*.root
    "$OUTDIR"/*.csv
    "$OUTDIR"/*_meta.json
    "$OUTDIR"/.*_part_*.csv
    "$OUTDIR"/.*_part_*.log
    "$OUTDIR"/.*_part_*.root
  )
  shopt -u nullglob

  if [[ ${#all_files[@]} -eq 0 ]]; then
    echo "No production files found in $OUTDIR"
    exit 0
  fi

  remove_if_exists "${all_files[@]}"
  echo "Done. Removed $removed file(s)."
  exit 0
fi

name="$(basename "$target")"
files=()
if [[ -e "$DATADIR/$name.csv" ]]; then
  files+=("$DATADIR/$name.csv")
fi
if [[ -e "$DATADIR/${name}_daughters.csv" ]]; then
  files+=("$DATADIR/${name}_daughters.csv")
fi
if [[ -e "$DATADIR/${name}_meta.json" ]]; then
  files+=("$DATADIR/${name}_meta.json")
fi

# Backward compatibility with legacy flat output layout.
if [[ -e "$OUTDIR/$name.csv" ]]; then
  files+=("$OUTDIR/$name.csv")
fi
if [[ -e "$OUTDIR/${name}_meta.json" ]]; then
  files+=("$OUTDIR/${name}_meta.json")
fi

shopt -s nullglob
files+=(
  "$DATADIR/.${name}_part_"*.csv
  "$DATADIR/.${name}_part_"*.log
  "$DATADIR/.${name}_part_"*.root
  "$OUTDIR/.${name}_part_"*.csv
  "$OUTDIR/.${name}_part_"*.log
  "$OUTDIR/.${name}_part_"*.root
)
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
  echo "No files found for output_name='$name' in $OUTDIR"
  exit 0
fi

remove_if_exists "${files[@]}"
echo "Done. Removed $removed file(s) for '$name'."

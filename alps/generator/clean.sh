#!/usr/bin/env bash
set -euo pipefail

# Cleanup helper for alps/generator/produce.sh outputs.
#
# Usage:
#   bash clean.sh <output_name>
#   bash clean.sh --all
#
# Examples:
#   bash clean.sh alp_light_m1
#   bash clean.sh --all

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTDIR="${SCRIPT_DIR}/../output"

usage() {
  cat <<'EOF'
Usage:
  bash clean.sh <output_name>
  bash clean.sh --all

Removes files produced by produce.sh:
  - output/<output_name>.csv
  - output/.<output_name>_part_*.csv
  - output/.<output_name>_part_*.log
  - output/.<output_name>_part_*.root

With --all, removes all matching production outputs in output/.
EOF
}

if [[ ! -d "$OUTDIR" ]]; then
  echo "Output directory not found: $OUTDIR"
  exit 1
fi

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
if [[ -e "$OUTDIR/$name.csv" ]]; then
  files+=("$OUTDIR/$name.csv")
fi
if [[ -e "$OUTDIR/${name}_meta.json" ]]; then
  files+=("$OUTDIR/${name}_meta.json")
fi

shopt -s nullglob
files+=(
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

#!/usr/bin/env bash
# Compile generator against Pythia8 only (no ROOT needed — output is LLP.csv).
# ROOT output is ifdef-guarded in generator.cc; omitting -DPY8ROOT disables it cleanly.
#
# Pythia8 path resolution: $PYTHIA8_DIR > <repo-root>/../pythia8315 (sibling of repo) > /usr/local
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PYTHIA="${PYTHIA8_DIR:-$(cd "${SCRIPT_DIR}/../../../pythia8315" 2>/dev/null && pwd || echo /usr/local)}"

g++ "${SCRIPT_DIR}/generator.cc" -o "${SCRIPT_DIR}/generator" -w \
  -I${PYTHIA}/include \
  -O2 -std=c++17 -fPIC -pthread \
  -L${PYTHIA}/lib -Wl,-rpath,${PYTHIA}/lib \
  -lpythia8 -ldl

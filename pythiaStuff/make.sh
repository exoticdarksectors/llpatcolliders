#!/usr/bin/env bash
# Compile main144 against Pythia8 only (no ROOT needed — output is LLP.csv).
# ROOT output is ifdef-guarded in main144.cc; omitting -DPY8ROOT disables it cleanly.
#
# Pythia8 path resolution: $PYTHIA8_DIR > ../pythia8315 (sibling dir) > /usr/local
PYTHIA="${PYTHIA8_DIR:-$(cd "$(dirname "$0")/../../pythia8315" 2>/dev/null && pwd || echo /usr/local)}"

g++ main144.cc -o main144 -w \
  -I${PYTHIA}/include \
  -O2 -std=c++17 -fPIC -pthread \
  -L${PYTHIA}/lib -Wl,-rpath,${PYTHIA}/lib \
  -lpythia8 -ldl

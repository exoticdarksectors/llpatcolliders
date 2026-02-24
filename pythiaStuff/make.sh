#!/usr/bin/env bash
# Compile main144 against Pythia8 only (no ROOT needed — output is LLP.csv).
# ROOT output is ifdef-guarded in main144.cc; omitting -DPY8ROOT disables it cleanly.

PYTHIA=/Users/fredi/sandbox-offline/pythia8315

g++ main144.cc -o main144 -w \
  -I${PYTHIA}/include \
  -O2 -std=c++17 -fPIC -pthread \
  -L${PYTHIA}/lib -Wl,-rpath,${PYTHIA}/lib \
  -lpythia8 -ldl

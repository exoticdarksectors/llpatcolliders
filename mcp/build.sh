#!/bin/bash
set -e
P8=/Users/mcitron/pythia8315
g++ mcp_dy.cc -o mcp_dy -O2 -std=c++11 -pthread \
    -I$P8/include -L$P8/lib -Wl,-rpath,$P8/lib -lpythia8 -ldl

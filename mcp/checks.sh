#!/bin/bash
# Phase-space-cut and LO-scale checks at the low-mass end (m = 10 GeV), where
# the photon propagator is steepest.  Reproduces the numbers in README.md.
set -e
printf 'PhaseSpace:mHatMin = 0.\n'   > c_mh0.cmnd
printf 'PhaseSpace:mHatMin = 20.\n'  > c_mh20.cmnd
printf 'PhaseSpace:mHatMin = 24.\n'  > c_mh24.cmnd
printf 'PhaseSpace:pTHatMinDiverge = 0.5\nPhaseSpace:pTHatMin = 0.\n' > c_pt.cmnd
printf 'SigmaProcess:renormMultFac = 0.25\nSigmaProcess:factorMultFac = 0.25\n' > c_lo.cmnd
printf 'SigmaProcess:renormMultFac = 4.0\nSigmaProcess:factorMultFac = 4.0\n'   > c_hi.cmnd
for t in mh0 mh20 mh24 pt; do
  ./mcp_dy 10 200000 1010 14000 out/chk_${t}.txt c_${t}.cmnd > out/chk_${t}.log 2>&1 &
done
for m in 10 100; do for t in lo hi; do
  ./mcp_dy $m 200000 $((1010+m)) 14000 out/scl_${m}_${t}.txt c_${t}.cmnd \
    > out/scl_${m}_${t}.log 2>&1 &
done; done
wait

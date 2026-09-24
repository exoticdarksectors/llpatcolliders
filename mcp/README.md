# DY production of a heavy fractionally charged particle, folded through A(eta)

Answers the handoff spec: generate `q qbar -> gamma*/Z -> chi chibar` at 14 TeV for
m in {10, 20, 30, 40, 50, 100, 200, 300, 500, 700, 1000} GeV, histogram single-particle eta in the
bins of `A_eta.npy`, and return `R(m) = int A f deta / 0.02384`, plus beta*gamma.

## Generator

MadGraph is not installed here, and Pythia's built-in `ffbar2ffbar(s:gmZ)` cannot be
used for a heavy final state (it treats the outgoing fermions as massless in the
phase space even though its matrix element carries the mass terms).  So the process
is implemented as a Pythia **semi-internal process** in `mcp_dy.cc`:

    dsigma_hat/dt_hat = (pi alpha^2 / s_hat^2) * e_q^2 * Q_chi^2
                        * [ (1 + cos^2 th) + (1 - beta^2)(1 - cos^2 th) ] / 3

photon-mediated, and with `MCP:includeZ = on` the Z is added as well:

    coef = e_q^2 gamProp Q_chi^2  +  e_q v_q intProp (Q_chi v_chi)
                                  +  (v_q^2 + a_q^2) resProp v_chi^2
    dsigma_hat/dt_hat = coef * [ (1 + cos^2 th) + (1 - beta^2)(1 - cos^2 th) ] / 3

**The Z must be included for these limits.**  Kinetic mixing has to be written with
hypercharge, `X_{mu nu} B^{mu nu}` - a mixing with the photon field strength is not
SU(2)xU(1) gauge invariant above the EW scale.  Shifting `X -> X - eps B` to remove
the mixing leaves the Higgs mass term untouched (it involves only `B`), so `X` stays
massless and unmixed, and the chi is left coupling to `B` with strength `-eps g_X`.
Writing `B = cos(thetaW) A - sin(thetaW) Z`, the chi couples to BOTH:

    photon:  eps g_X cos(thetaW)  =  Q_chi e
    Z:       eps g_X sin(thetaW)  =  Q_chi e tan(thetaW)

i.e. the chi behaves exactly as a fermion with `T3 = 0` and electric charge
`Q_chi`, giving `v_chi = -4 sin^2(thetaW) Q_chi` and `a_chi = 0` in Pythia's
coupling convention.  The Z coupling does not vanish for a massless dark photon -
anything carrying effective hypercharge couples to the Z.  (The same conclusion
follows in the other basis, where the shift instead induces `Z0`-`X` mass mixing and
the massive eigenstate picks up an `eps sin(thetaW) X` component.)

Because `a_chi = 0` the chi vertex is pure vector for both mediators, so the
final-state Lorentz structure is identical to the photon-only case and there is no
forward-backward asymmetry.  The Z enters *only* by reweighting the `sHat`
spectrum - which is precisely why it matters so much below the pole.

`chi` is `id 6000015`, spin 1/2, colour singlet, stable.  Its matrix-element charge
is `Q_chi = 1`; the cross section for charge `eps` is `eps^2 x sigma` and the eta
spectrum is independent of `eps`.  Its particle-data charge is set to 0 so Pythia
does not shower QED off it (that radiation is `O(eps^2)` and negligible).

Settings: 14 TeV pp, Pythia 8.315 default PDF (NNPDF2.3 LO), `mu_R = mu_F = mHat`,
ISR and FSR on, MPI and hadronisation off (they do not touch the chi kinematics),
400k events per mass.  LO only - no K-factor applied.

### Validation

* **Normalisation.** Setting `m_chi = 0.5 GeV` and cutting `mHat > 200 (1000) GeV`
  reproduces Pythia's own pure-gamma* Drell-Yan to muons: 1.4661e-9 vs 1.466e-9 mb
  (and 3.951e-12 vs 3.934e-12 mb), i.e. agreement at the sub-percent MC level.
  This checks the prefactor, the colour factor, the parton flux and the `1+cos^2`
  term.
* **Mass and angular terms.** `COSNEAR` / `COSFAST` in each output file are the
  `cos(theta*)` distributions for `beta < 0.5` and `beta > 0.9`.  The first is flat
  (isotropic at threshold, as `1+c^2+(1-beta^2)(1-c^2) -> 2`), the second is
  U-shaped with edge/centre = 2 (the `1+cos^2` limit).

* **s-channel sampling.** `resonanceA()` returns 23 when the Z is on, so that
  `PhaseSpace` adds a Breit-Wigner term at `sHat = mZ^2` to the tau sampling.
  Without it the sampling is only `1/tau + 1/tau^2` and the pole is badly
  under-sampled whenever it sits well above the `2 m_chi` threshold: the total
  cross section stays right (the integral is unbiased) but the accepted event
  sample under-represents the pole, which is the more central population, so
  acceptances come out too low.  At m = 10 GeV this biased `P(>=1 chi at
  |eta|<1)` from 0.336 down to 0.267, and made the 13.6 and 14 TeV results
  disagree by 12% when they must agree to under 1%.  That energy-consistency
  check is the cheap way to catch it.  Pythia drops the term by itself once
  `mHatMin` rises above the pole, so masses above `m_Z/2` are unaffected.

* **Z normalisation.** In a narrow window `86 < mHat < 96 GeV` the m = 10 GeV
  `gamma*/Z` cross section is 1.145e-6 mb, against 1.3315e-6 mb for Pythia's own
  Z-only Drell-Yan to muons in the same window.  Subtracting the photon continuum
  measured in that window (5.51e-9 mb) gives a ratio of 0.856, against the 0.8497
  predicted from `(v_chi^2 + a_chi^2)/(v_mu^2 + a_mu^2)` times the threshold factor
  `beta(3-beta^2)/2`.  The residual sub-percent difference is the gamma-Z
  interference, which is not separated out.  `MHATHIST` puts the pole at the right
  place with the right width.

### Phase-space cuts and LO scale dependence

The m = 10 GeV point sits at 808 pb for Q = e, so it is worth confirming that no
generator cut is shaping the low-mass points differently from the rest.  It is not:
the threshold is entirely physical.  `mHatMin` is set by Pythia as
`max(2*m_chi, PhaseSpace:mHatMin)`, and the 4 GeV default is far below `2*m_chi = 20
GeV` for even the lightest point, so it never binds.  `pTHatMinDiverge` (default 1
GeV) applies only when a produced particle has mass below 1 GeV, so it never
applies either.  The `MHATHIST` record in each output file starts in exactly the
bin containing `2*m_chi` for every mass.  Explicitly (m = 10 GeV, 200k events):

| variation | sigma [pb] | R |
|---|---|---|
| default (`mHatMin = 4`) | 808.4 | 0.1527 +- 0.0011 |
| `mHatMin = 0` | 808.4 | 0.1527 +- 0.0011 |
| `mHatMin = 20` (= 2m) | 808.4 | 0.1527 +- 0.0011 |
| `pTHatMinDiverge = 0.5`, `pTHatMin = 0` | 808.4 | 0.1527 +- 0.0011 |
| `mHatMin = 24` (= 1.2 x 2m) | 607.6 | 0.1711 +- 0.0012 |

The first four are identical; the last is a deliberate cut, included to show the
test has teeth - a stray cut at 1.2x threshold would move sigma by 25% and R by 12%.

There is also no propagator divergence to regulate: the `1/sHat^2` photon
propagator is cut off physically at `sHat = 4 m_chi^2`.  The low-mass points are in
fact *less* threshold-dominated than the high-mass ones - the fraction of the cross
section below `1.2 x 2m` rises monotonically from 0.257 at m = 10 GeV to 0.520 at
m = 1000 GeV.

The LO scale dependence is a real caveat, but for sigma rather than for R.  Varying
`mu_R = mu_F` by a factor 2 either way:

| m [GeV] | mu | sigma [pb] | ratio | R |
|---|---|---|---|---|
| 10 | mHat/2 | 597.8 | 0.740 | 0.1500 +- 0.0011 |
| 10 | mHat | 807.7 | 1.000 | 0.1527 +- 0.0008 |
| 10 | 2 mHat | 1039 | 1.286 | 0.1568 +- 0.0011 |
| 100 | mHat/2 | 0.999 | 0.945 | 0.2265 +- 0.0014 |
| 100 | mHat | 1.057 | 1.000 | 0.2282 +- 0.0010 |
| 100 | 2 mHat | 1.106 | 1.047 | 0.2275 +- 0.0014 |

So the m = 10 GeV cross section carries a -26%/+29% LO scale uncertainty (and the
low-mass DY K-factor is not small), while R moves by only +-2%.  R is an
acceptance ratio and the scale largely cancels in it.  Run `./checks.sh` to
reproduce.

## The isotropic normalisation

The 0.02384 that arrived with the original handoff is exactly
`int A(eta) * (1/2)sech^2(eta) deta` for the hand-made A(eta) of the time.  The
reference f is therefore an **isotropic single particle (flat in cos theta),
normalised over all eta**, and `R(m)` uses f normalised over all eta too, so it
is the per-particle acceptance of the real DY spectrum relative to that
isotropic assumption.

`analyze.py` now derives that closure from whichever A(eta) is loaded rather
than hard-coding it (0.02525 for the current full-fiducial array), and records
it as `a_iso` in the results.  `R` is therefore a pure shape ratio, unchanged
to a few parts in a thousand if the acceptance normalisation changes, and
comparable across geometry definitions.  (For contrast, if f is instead renormalised inside |eta| < 1.2 the
ratio is a nearly mass-independent 0.87-0.94 - column `<A>|window` in
`analyze.py`'s output.)

## Results

| m [GeV] | sigma gmZ [fb] | **R gmZ** | sigma gam [fb] | R gam | sigma ratio | R ratio |
|---|---|---|---|---|---|---|
| 10 | 2.227e+06 | 0.2347 +- 0.0010 | 8.077e+05 | 0.1535 +- 0.0008 | 2.757 | 1.529 |
| 20 | 1.475e+06 | 0.2415 +- 0.0010 | 1.307e+05 | 0.1672 +- 0.0008 | 11.287 | 1.444 |
| 30 | 1.256e+06 | 0.2042 +- 0.0009 | 4.19e+04 | 0.1743 +- 0.0008 | 29.977 | 1.171 |
| 35 | 1.103e+06 | 0.1758 +- 0.0008 | 2.685e+04 | 0.1823 +- 0.0008 | 41.067 | 0.965 |
| 36 | 1.06e+06 | 0.1694 +- 0.0008 | 2.473e+04 | 0.1833 +- 0.0008 | 42.877 | 0.924 |
| 37 | 1.018e+06 | 0.1604 +- 0.0008 | 2.282e+04 | 0.1820 +- 0.0008 | 44.583 | 0.881 |
| 38 | 9.7e+05 | 0.1552 +- 0.0008 | 2.109e+04 | 0.1842 +- 0.0008 | 45.987 | 0.843 |
| 39 | 9.137e+05 | 0.1472 +- 0.0008 | 1.953e+04 | 0.1854 +- 0.0008 | 46.781 | 0.794 |
| 40 | 8.506e+05 | 0.1410 +- 0.0007 | 1.814e+04 | 0.1853 +- 0.0008 | 46.889 | 0.761 |
| 41 | 7.769e+05 | 0.1328 +- 0.0007 | 1.685e+04 | 0.1849 +- 0.0008 | 46.112 | 0.718 |
| 42 | 6.927e+05 | 0.1227 +- 0.0007 | 1.568e+04 | 0.1876 +- 0.0009 | 44.170 | 0.654 |
| 43 | 5.883e+05 | 0.1145 +- 0.0007 | 1.463e+04 | 0.1866 +- 0.0009 | 40.200 | 0.614 |
| 44 | 4.602e+05 | 0.1053 +- 0.0006 | 1.367e+04 | 0.1866 +- 0.0008 | 33.663 | 0.564 |
| 45 | 2.86e+05 | 0.0994 +- 0.0006 | 1.276e+04 | 0.1916 +- 0.0009 | 22.420 | 0.519 |
| 46 | 1.113e+05 | 0.1120 +- 0.0007 | 1.195e+04 | 0.1908 +- 0.0009 | 9.314 | 0.587 |
| 47 | 5.561e+04 | 0.1298 +- 0.0007 | 1.119e+04 | 0.1918 +- 0.0009 | 4.968 | 0.677 |
| 48 | 3.672e+04 | 0.1416 +- 0.0007 | 1.053e+04 | 0.1897 +- 0.0009 | 3.487 | 0.746 |
| 49 | 2.754e+04 | 0.1515 +- 0.0008 | 9906 | 0.1929 +- 0.0009 | 2.780 | 0.785 |
| 50 | 2.201e+04 | 0.1585 +- 0.0008 | 9301 | 0.1936 +- 0.0009 | 2.366 | 0.819 |
| 100 | 1039 | 0.2276 +- 0.0009 | 1057 | 0.2277 +- 0.0009 | 0.983 | 0.999 |
| 200 | 92.21 | 0.2831 +- 0.0010 | 98.83 | 0.2839 +- 0.0011 | 0.933 | 0.997 |
| 300 | 20.15 | 0.3276 +- 0.0011 | 21.75 | 0.3254 +- 0.0011 | 0.927 | 1.007 |
| 500 | 2.433 | 0.3953 +- 0.0012 | 2.633 | 0.3967 +- 0.0012 | 0.924 | 0.997 |
| 700 | 0.5057 | 0.4553 +- 0.0013 | 0.5497 | 0.4518 +- 0.0013 | 0.920 | 1.008 |
| 1000 | 0.07685 | 0.5271 +- 0.0014 | 0.0835 | 0.5262 +- 0.0014 | 0.920 | 1.002 |

`gmZ` is the physical case (hypercharge kinetic mixing); `gam` is photon-only, kept
for comparison.  Uncertainties are MC statistical only.  Above the pole the Z costs
7-8% in sigma and leaves R alone; below it the Z dominates completely - sigma grows
by up to a factor 47 and R moves by -24% to +45%, non-monotonically, because the
sHat spectrum is pinned at m_Z rather than falling away from `2 m_chi`.

The 1 GeV scan over 35-50 GeV resolves that structure, which is the chi chi
threshold walking up through the Z pole.  R falls monotonically from 0.1747 at
35 GeV to a minimum of 0.0993 at 45 GeV, then turns over and climbs back:
0.1116, 0.1276, 0.1426, 0.1510, 0.1583 at 46-50 GeV.  The turning point is
`m_chi = m_Z/2 = 45.59 GeV`, where `2 m_chi = m_Z` exactly and the chi are produced
at rest on the pole - median `beta*gamma` bottoms out at 0.44 there.  Above it the
pole is closed and only the Breit-Wigner tail survives, so sigma collapses by an
order of magnitude over 45-50 GeV (2.86e5 -> 2.20e4 fb) while R recovers.  The
photon-only R is flat at ~0.186 across the whole window, so every bit of this
structure is the Z.  `plot_zoom.py` draws the region.

## Per-particle vs per-event

(gamma*/Z case.)  `R(m)` above is **per particle**: `int A f deta` with `f` the single-particle eta
pdf, matching the single-particle isotropic reference behind 0.02384.

`per_event.py` gives the per-event numbers, treating `A(eta)` as a per-particle
probability and using the real `eta1`-`eta2` correlation from the pair boost (which
is why `P(both)` is not `R_pp^2`).  Reference row = isotropic and uncorrelated.

| m [GeV] | R per particle | <N_acc>/event | P(>=1 chi) | P(both chi) | P(both)/R_pp^2 |
|---|---|---|---|---|---|
| 10 | 0.2486 | 0.01185 | 0.01178 | 7.63e-05 | 2.17 |
| 20 | 0.2557 | 0.01219 | 0.01211 | 7.94e-05 | 2.14 |
| 30 | 0.2162 | 0.01031 | 0.01023 | 7.57e-05 | 2.85 |
| 35 | 0.1862 | 0.00888 | 0.00881 | 7.07e-05 | 3.58 |
| 36 | 0.1794 | 0.00855 | 0.00849 | 6.57e-05 | 3.59 |
| 37 | 0.1699 | 0.00810 | 0.00804 | 6.41e-05 | 3.91 |
| 38 | 0.1644 | 0.00784 | 0.00778 | 6.25e-05 | 4.07 |
| 39 | 0.1559 | 0.00743 | 0.00737 | 6.1e-05 | 4.42 |
| 40 | 0.1494 | 0.00712 | 0.00706 | 5.94e-05 | 4.68 |
| 41 | 0.1407 | 0.00671 | 0.00665 | 5.73e-05 | 5.09 |
| 42 | 0.1299 | 0.00620 | 0.00614 | 5.38e-05 | 5.61 |
| 43 | 0.1213 | 0.00578 | 0.00573 | 5.41e-05 | 6.47 |
| 44 | 0.1115 | 0.00532 | 0.00526 | 5.3e-05 | 7.51 |
| 45 | 0.1053 | 0.00502 | 0.00497 | 5.48e-05 | 8.70 |
| 46 | 0.1186 | 0.00565 | 0.00560 | 5.81e-05 | 7.26 |
| 47 | 0.1374 | 0.00655 | 0.00649 | 6.02e-05 | 5.61 |
| 48 | 0.1500 | 0.00715 | 0.00709 | 6.1e-05 | 4.77 |
| 49 | 0.1605 | 0.00765 | 0.00759 | 6.3e-05 | 4.30 |
| 50 | 0.1679 | 0.00800 | 0.00794 | 6.53e-05 | 4.08 |
| 100 | 0.2410 | 0.01149 | 0.01141 | 8.22e-05 | 2.49 |
| 200 | 0.2999 | 0.01430 | 0.01419 | 0.000104 | 2.03 |
| 300 | 0.3469 | 0.01654 | 0.01642 | 0.000123 | 1.79 |
| 500 | 0.4187 | 0.01996 | 0.01981 | 0.000154 | 1.54 |
| 700 | 0.4823 | 0.02299 | 0.02281 | 0.000189 | 1.43 |
| 1000 | 0.5583 | 0.02662 | 0.02638 | 0.000236 | 1.34 |
| *isotropic* | 1.0000 | 0.04768 | 0.04711 | 0.000568 | 1.00 |

The last column is the pair-correlation enhancement of the both-hit probability over
the naive square of the per-particle acceptance: at low mass the pair is boosted
together and lands on the same side, so requiring both is ~3x less costly than
independence would suggest.

## Velocity / time of flight

gamma*/Z case.  Acceptance-weighted (weights `A(eta)`, i.e. the chi that actually reach the
detector).  `dt` is the delay relative to a massless particle over 10 m of path.

| m [GeV] | bg 2% | bg 16% | bg med | bg 84% | beta med | f(bg<1) | f(bg<2) | dt/10m med | dt/10m 98% |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.46 | 1.30 | 3.33 | 4.81 | 0.958 | 0.11 | 0.27 | 1.5 ns | 47 ns |
| 20 | 0.49 | 1.11 | 1.82 | 2.39 | 0.876 | 0.12 | 0.63 | 4.7 ns | 42 ns |
| 30 | 0.35 | 0.72 | 1.08 | 1.46 | 0.735 | 0.40 | 0.95 | 12.0 ns | 67 ns |
| 35 | 0.27 | 0.56 | 0.82 | 1.16 | 0.636 | 0.74 | 0.97 | 19.1 ns | 93 ns |
| 36 | 0.26 | 0.53 | 0.78 | 1.12 | 0.614 | 0.77 | 0.97 | 20.9 ns | 98 ns |
| 37 | 0.25 | 0.49 | 0.73 | 1.07 | 0.589 | 0.81 | 0.97 | 23.3 ns | 106 ns |
| 38 | 0.23 | 0.46 | 0.68 | 1.03 | 0.562 | 0.83 | 0.98 | 26.0 ns | 114 ns |
| 39 | 0.21 | 0.43 | 0.64 | 0.98 | 0.536 | 0.85 | 0.98 | 28.9 ns | 131 ns |
| 40 | 0.20 | 0.39 | 0.59 | 0.95 | 0.508 | 0.86 | 0.98 | 32.4 ns | 140 ns |
| 41 | 0.18 | 0.36 | 0.55 | 0.92 | 0.479 | 0.87 | 0.98 | 36.3 ns | 156 ns |
| 42 | 0.16 | 0.32 | 0.50 | 0.89 | 0.450 | 0.88 | 0.98 | 40.8 ns | 175 ns |
| 43 | 0.14 | 0.29 | 0.47 | 0.90 | 0.425 | 0.87 | 0.98 | 45.2 ns | 205 ns |
| 44 | 0.12 | 0.25 | 0.44 | 0.90 | 0.405 | 0.87 | 0.98 | 49.1 ns | 243 ns |
| 45 | 0.10 | 0.22 | 0.45 | 0.95 | 0.408 | 0.85 | 0.97 | 48.4 ns | 296 ns |
| 46 | 0.11 | 0.26 | 0.53 | 1.12 | 0.471 | 0.80 | 0.96 | 37.5 ns | 271 ns |
| 47 | 0.13 | 0.30 | 0.63 | 1.31 | 0.531 | 0.74 | 0.94 | 29.4 ns | 228 ns |
| 48 | 0.15 | 0.34 | 0.70 | 1.44 | 0.573 | 0.69 | 0.92 | 24.9 ns | 195 ns |
| 49 | 0.16 | 0.37 | 0.75 | 1.52 | 0.599 | 0.66 | 0.91 | 22.3 ns | 176 ns |
| 50 | 0.17 | 0.39 | 0.79 | 1.58 | 0.620 | 0.64 | 0.91 | 20.5 ns | 166 ns |
| 100 | 0.24 | 0.53 | 1.00 | 1.86 | 0.706 | 0.50 | 0.86 | 13.9 ns | 109 ns |
| 200 | 0.23 | 0.51 | 0.94 | 1.69 | 0.684 | 0.54 | 0.90 | 15.4 ns | 117 ns |
| 300 | 0.22 | 0.49 | 0.89 | 1.57 | 0.665 | 0.58 | 0.92 | 16.8 ns | 121 ns |
| 500 | 0.21 | 0.45 | 0.82 | 1.39 | 0.635 | 0.64 | 0.96 | 19.2 ns | 132 ns |
| 700 | 0.19 | 0.42 | 0.76 | 1.26 | 0.606 | 0.70 | 0.98 | 21.7 ns | 143 ns |
| 1000 | 0.18 | 0.39 | 0.69 | 1.12 | 0.570 | 0.77 | 0.99 | 25.1 ns | 155 ns |

## Comparison with an external reference

`~/Downloads/dy_xsecs_eta1_RUN3.txt` (Run 3, 13.6 TeV) has three columns: mass,
inclusive cross section in pb for `Q_chi = e`, and an acceptance fraction.
`plot_ref.py` compares it with this generation.  It has the `m_Z/2` cusp, so it
includes the Z.

* **Column 3 is `P(at least one chi within |eta| < 1)`, per event.**  It matches
  the quantity computed here to 0.6-3% over 10-70 GeV: 0.3377/0.3356, 0.3373/0.3454,
  0.2859/0.2896, 0.1935/0.1994, 0.2265/0.2258, 0.2777/0.2752, 0.2887/0.2942 at
  m = 10, 20, 30, 40, 50, 60, 70.  It is not a per-particle fraction (that is 0.21
  at m = 10, not 0.34) and not the both-chi probability (0.089).
* **Column 2 is inclusive** and sits 1.20-1.26 above the LO cross section here at
  matched beam energy.  Roughly half of that is PDF order - switching to
  NNPDF2.3 NLO (`PDF:pSet = 15`) moves the ratio to 1.10-1.11 at m = 20, 40 and 50
  alike - and the remaining ~11% is consistent with an NLO matrix element against
  this LO one.  It is not a per-particle/per-event factor, which would be exactly 2.
* The acceptance ratio dips to 0.90 at m = 45-46.  That is interpolation of the
  reference's own coarse mass grid across the cusp, not a disagreement.

## Sensitivity: the (mass, charge) reach plot

`make_mcp_plot.py` turns the numbers above into the GRENDEL MCP reach at the
HL-LHC.  It reads `out/results_gmZ.json` and `out/per_event_gmZ.json` directly,
so the plot always reflects the current generation - there is no second copy of
the cross sections to keep in sync.

    python3 make_mcp_plot.py           # closed contours, with the muon upper cut
    python3 make_mcp_plot.py --open    # one-sided Qmin curves, no upper cut

Writes `out/grendel_mcp_{closed,open}.png` and one CSV per configuration.
Signal model:

    N(m, Q) = L * sigma(m) * Q^2 * pAny(m) * eps(Q)

| term | source |
|---|---|
| `sigma(m)` | `sigma_fb` - pure `gamma*/Z` pair production at `Q = e`, 14 TeV, from `mcp_dy.cc`.  Scales as `Q^2`: the chi is an SU(2) singlet, so both its photon and its Z coupling are proportional to Q (see the Generator section). |
| `pAny(m)` | `pAny` - per-EVENT probability that at least one of the two chi crosses the tracker. |
| `eps(Q)` | `[ P(NMIN <= N <= NMAX | lambda = NPE * Q^2) ]^4` - Poisson photon count in each of the four scintillator planes. |

**Use the per-event number, not the per-particle one.**  `pAny = Nacc - pBoth`,
and `Nacc = 2 * perPart`, so the per-particle acceptance (`perPart`, identical
to `I_all` in the results file) is roughly half of `pAny`.  Feeding it in would
undercount the rate by about a factor 2.  The script asserts that the two JSON
files describe the same sample, and `check_a_iso()` re-derives the isotropic
closure from `external/A_eta.npy` and checks it against the hard-coded
`A_ISO = 0.02384`, so the acceptance normalisation cannot drift.

### The upper edge

A muon is a `Q = e` particle, so its mean photon count per plane is `NPE`.
Requiring `[P(N <= NMAX | NPE)]^4 <= 1 / 2e9` sets `NMAX`:

| NPE at Q = e | NMAX | per-plane survival | NMAX / NPE |
|---|---|---|---|
| 100 | 74 | 4.0e-3 | 0.74 (2.60 sigma) |
| 200 | 163 | 4.0e-3 | 0.815 (2.62 sigma) |

The same cut bounds the signal from above, which is what closes the contour: a
`Q -> e` particle is indistinguishable from a muon by ionisation alone.  This
treats dE/dx as the *only* muon discriminant and so understates the real
rejection - the floor and right-wall veto panels, timing and pointing all
contribute.  If those carry most of the rejection, `NMAX` relaxes and the upper
edge moves toward `Q = e`.

### Steepness, and what softens it

Near the lower edge `lambda = NPE * Q^2` is a few, deep in the Poisson tail:

    d ln(eps) / d ln(Q) = 2 * NLAYER * lambda * p(NMIN-1; lambda) / P(N >= NMIN)

At `NMIN = 10`, 100 PE this is ~45, so `N ~ Q^47` including the `Q^2` from
production.  That is why two orders of magnitude in cross section move `Qmin` by
only ~10%, and why the curves turn up at high mass (the index falls to 2 once
`eps -> 1`).

That exponent assumes a hard photon-count threshold with no smearing.  SiPM gain
spread, path-length variation through the bars and the beta-dependent dE/dx
enhancement would all soften it - and the samples here are slow: median beta is
0.57-0.96 and 22-86% of the accepted chi have `beta*gamma < 1` (see the velocity
table), so that last term is not a small correction.

Interpolation in mass is log-log linear.  The mass grid is dense enough
everywhere except across the `m_Z/2` cusp, where sigma falls by 13x between the
45 and 46 GeV samples; a straight line in log-log across that gap is the weakest
part of the interpolation.

### External curves (`external/`)

| file | contents |
|---|---|
| `Collider.csv` | LEP bound.  Excluded region is ABOVE the curve, for `m < m_Z/2 = 45.6 GeV`, where it terminates. |
| `milliQanRun3Fix.csv` | milliQan Run 3 LOWER boundary only; mass-dependent, 0.1-22.7 GeV.  Its upper boundary is the LEP line (the two meet at 0.235 vs 0.237 at 22.7 GeV). |
| `CMS.csv` | CMS FCP search, closed contour, m = 50-641 GeV, Q = 0.33-0.91.  Derived at `Q = e/3`; the sigma limit cannot be rescaled to other charges because the efficiency is dE/dx dependent. |
| `FORMOSA_nominal2ab.csv` | FORMOSA projection at 2 ab^-1 - a projection, not an exclusion, so it is drawn as an open curve.  Covers 0.01-97.5 GeV, so only its high-mass tail appears on a 10-1000 GeV axis. |

### Geometry provenance

`external/A_eta.npy` is the azimuthal fraction, at each pseudorapidity, of
directions from the CMS IP that cross the GRENDEL fiducial volume.  It is
**built from `../higgs/grendel_geometry.py`** by `build_aeta.py` - the geometry
is not taken on trust and not frozen:

    python3 build_aeta.py -n 40000      # what produced the current file

48 bins of width 0.05, whole fiducial volume (`--surface full`), 24 cm inset,
40k rays per bin.  Isotropic closure 0.02525, i.e. **0.3173 sr** - matching the
0.3161 sr quoted for the 673.4 m^3 fiducial volume.

The superseded hand-made array is kept as
`external/A_eta_2026-09-02_preflip.npy` (and `.csv`).  Two things were wrong
with it, both now moot:

* **It was mirrored in eta.**  It is dated 2026-09-02; commit `31f8512` ("flip
  beam Z to CMS +z/-z convention") landed 2026-09-19, so its eta axis ran the
  other way - peak at eta = -0.43 where today's geometry puts +0.43.  Testing it
  against the ray-cast `frac(eta)` gives chi2 = 13396 over 24 bins; against
  `frac(-eta)`, chi2 = 288.  It never mattered, because the DY eta spectrum is
  symmetric, but a rebuild comes out in today's convention automatically.
* **Its surface definition matched neither live option.**  Its closure sat at
  0.95 of the full fiducial solid angle, where `points_on_tracker()`
  (arch/ceiling + left wall) keeps 0.77.  Whatever "tracker" meant then covered
  nearly the whole acceptance.

Switching to the live full-volume array moved `pAny` up ~6% (0.01109 -> 0.01178
at m = 10), left `R` alone (0.2340 -> 0.2347), and strengthened `Qmin` on the
reach plot by 0.1-3.3%.

### Choosing a different surface

`build_aeta.py` can also produce the tracker-only acceptance, or drop the
fiducial inset:

    python3 build_aeta.py --surface tracker   # arch/ceiling + left wall only
    python3 build_aeta.py --inset 0           # the un-inset tunnel wall

Every script reads `external/A_eta.npy` unless `MCP_AETA` names another file,
so an alternative can be tried without disturbing the default:

    MCP_AETA=external/A_eta_live_tracker.npy python3 analyze.py gmZ
    MCP_AETA=external/A_eta_live_tracker.npy python3 per_event.py gmZ

Isotropic closures:

| A(eta) | closure | solid angle | vs default |
|---|---|---|---|
| `A_eta.npy` - full fiducial (default) | 0.02525 | 0.3173 sr | 1.000 |
| `A_eta_live_tracker.npy` | 0.01924 | 0.2418 sr | 0.762 |
| `A_eta_2026-09-02_preflip.npy` | 0.02384 | 0.2996 sr | 0.944 |

The tracker-only definition costs 24% of the rate and weakens `Qmin` by 0.5% at
low mass, up to 12% at high mass - small at low mass because of the steep
Poisson exponent above, growing once `eps -> 1` and the index falls to 2.  `R`
is insensitive to the choice, since it divides by the closure of the same array.

`analyze.py` records `a_iso` and `a_eta` in the results, and
`make_mcp_plot.py` re-derives the closure from that file and refuses to run if
it has changed - so the acceptances and the reach plot cannot be built from
different geometries.


## Files

    mcp_dy.cc     generator (semi-internal Pythia process)
    builtin_dy.cc reference pure-gamma* DY to muons, used for the xsec validation
    build.sh      compile against /Users/mcitron/pythia8315
    analyze.py    folds the eta spectra through A(eta); prints the tables.
                  Variant argument: `gam` (photon only) or `gmZ`
    per_event.py  per-event acceptance from the PAIR records; same argument
    gmZ.cmnd      one line, `MCP:includeZ = on`
    plot.py       out/mcp_dy_summary.png
    plot_zoom.py  out/mcp_dy_zoom.png, the 35-50 GeV region
    plot_ref.py   out/mcp_dy_ref_compare.png, against the external reference
    zmumu.cc      independent Z -> mu mu eta cross-check
    make_mcp_plot.py  the (mass, charge) reach plot; reads out/*.json
    build_aeta.py     builds external/A_eta.npy by ray-casting
                      ../higgs/grendel_geometry.py from the CMS IP
    external/     A_eta.npy (built by build_aeta.py), the superseded
                  A_eta_2026-09-02_preflip.{npy,csv}, the digitised curves from
                  LEP / milliQan / CMS / FORMOSA, and the reference Run-3 cross
                  sections used by plot_ref.py
    make_readme_tables.py  regenerates the tables in this file
    checks.sh     cut and scale variations behind the checks section above
    out/mcp_m*.txt  per-mass output: header, full-eta histogram, cos/mHat
                    validation histograms, PAIR records (both chi eta per event,
                    for events with >=1 chi at |eta|<1.2), then
                    (eta, betagamma, pt) for |eta|<1.5
    out/results_{gam,gmZ}.{csv,json}, out/per_event_{gam,gmZ}.json

Reproduce:

    ./build.sh
    python3 build_aeta.py -n 40000        # A(eta) from ../higgs/, if not present
    M="10 20 30 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 100 200 300 500 700 1000"
    for m in $M; do
      ./mcp_dy $m 400000 $((1000+m)) 14000 out/mcp_m${m}.txt > out/log_m${m}.txt
      ./mcp_dy $m 400000 $((1000+m)) 14000 out/z_m${m}.txt gmZ.cmnd > out/zlog_m${m}.txt
    done
    for t in gam gmZ; do python3 analyze.py $t; python3 per_event.py $t; done
    python3 plot.py && python3 plot_zoom.py && python3 plot_ref.py
    python3 make_mcp_plot.py && python3 make_mcp_plot.py --open
    python3 make_readme_tables.py

Masses are discovered from `out/`, so adding a mass only means running
`mcp_dy` for it and rerunning the analysis; nothing hard-codes the grid.

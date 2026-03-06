# HNL Portal: Heavy Neutral Leptons

Sensitivity projections for GARGOYLE to Heavy Neutral Leptons (HNLs) at HL-LHC
(14 TeV, 3 ab^-1), using the FONLL NLO+NLL approach.

## Physics

HNLs mix with active neutrinos via |U_alpha|^2 (alpha = e, mu, tau).
Production channels:

| Channel | Mass range | Method |
|---------|-----------|--------|
| B meson (B -> D(*) l N) | 0.5-5.3 GeV | FONLL + decay kinematics |
| D meson (D -> K(*) l N) | 0.2-1.8 GeV | FONLL + decay kinematics |
| Bc (Bc -> X l N) | 0.2-6.3 GeV | FONLL + decay kinematics |
| tau (W -> tau nu -> N X) | 0.2-1.78 GeV | MadGraph + tau decay |
| W/Z (W/Z -> l N) | 5-80 GeV | MadGraph (SM_HeavyN model) |

PBC benchmarks: BC6 (|Ue|^2), BC7 (|Umu|^2), BC8 (|Utau|^2).

## Quick Start

```bash
# 1. Generate ctau tables (requires HNLCalc)
python hnl/production/generate_ctau_tables.py --br-vis

# 2. Generate meson -> HNL 4-vectors (FONLL channels)
python hnl/production/decay_engine/generate_meson_csvs.py --flavor Ue --masses 1.0

# 3. Generate tau -> HNL 4-vectors (requires MadGraph)
python hnl/production/madgraph/run_tau_production.py --flavor Ue

# 4. Generate W/Z -> l N (requires MadGraph + SM_HeavyN model)
python hnl/production/madgraph/run_wz_production.py --flavor Ue

# 5. Combine all channels
python hnl/production/combine_channels.py --flavor Ue

# 6. Run sensitivity analysis
python hnl/analysis/run_sensitivity.py --flavor Ue --workers 6

# 7. Run tests
python -m pytest hnl/tests/ -v
```

## Directory Structure

```
hnl/
  config_mass_grid.py         Mass grid (0.2-10 GeV, 110 points)
  production/
    constants.py              FONLL cross-sections, fragmentation fractions
    combine_channels.py       Merge all channel CSVs per (flavor, mass)
    generate_ctau_tables.py   HNL lifetime tables via HNLCalc
    fonll/                    FONLL grid parser + inverse-CDF meson sampler
    decay_engine/             Vectorized 2-body + 3-body decay kinematics
    madgraph/                 W/Z and tau production via MG5_aMC@NLO
  analysis/
    constants.py              Luminosity, cuts, U^2 scan range
    sensitivity.py            Acceptance-weighted decay probability
    exclusion.py              U^2 exclusion band extraction
    run_sensitivity.py        CLI driver (parallel mass scan)
    plot_exclusion.py         1x3 panel exclusion plot
    format_bridge.py          CSV -> geometry-ready arrays
    reference_curves.py       Load comparison experiment contours
    plot_diagnostics.py       Diagnostic plots (BR_vis, ctau, weight sums)
  tests/                      Unit tests (format_bridge, kinematics, sensitivity)
  vendored/
    HNLCalc/                  HNL branching ratio / lifetime calculator
    fonll_grids/              FONLL meson-level d sigma/dpT/dy tables (14 TeV)
    reference_curves/         Processed comparison curves (mass u2_min u2_max)
  external/                   Raw digitized curves from other experiments
```

## Dependencies

- numpy, scipy, pandas, matplotlib, particle
- trimesh (for GARGOYLE geometry ray-casting)
- MadGraph5_aMC@NLO (for W/Z and tau channels; set MG5_DIR env var)

## References

- FONLL: Cacciari, Greco, Nason (NLO+NLL heavy-quark production)
- HNL BRs: Bondarenko et al., arXiv:1805.08567
- PBC benchmarks: arXiv:1901.09966 (BC6/BC7/BC8)
- MATHUSLA: arXiv:1901.04040
- CODEX-b: arXiv:1911.00481

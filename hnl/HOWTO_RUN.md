# HNL Analysis: Full Run Guide

Copy-paste commands to run the complete GARGOYLE HNL sensitivity analysis.
All commands run from the repo root (`llpatcolliders/`).

## Prerequisites

```bash
# 1. Activate conda environment
conda activate llpatcolliders
# or prefix every command with: conda run -n llpatcolliders

# 2. MadGraph (needed for tau and W/Z channels only)
export MG5_DIR=/path/to/MG5_aMC_v3_6_6
# The SM_HeavyN_CKM_AllMasses_LO model must be installed in MG5:
#   python $MG5_DIR/bin/mg5_aMC
#   MG5> install model SM_HeavyN_CKM_AllMasses_LO

# 3. Shared geometry must be in place
ls geometry/gargoyle_geometry.py   # must exist
```

## Quick single-mass test (5 min)

Verify the pipeline end-to-end with one mass point before running the full scan.

```bash
# ctau table (single mass)
python hnl/production/generate_ctau_tables.py --flavor Umu --masses 1.0

# FONLL meson channels
python hnl/production/decay_engine/generate_meson_csvs.py \
    --flavor Umu --channel all --masses 1.0 --n-pool 10000

# Combine
python hnl/production/combine_channels.py --flavor Umu --masses 1.0

# Sensitivity (single worker)
python hnl/analysis/run_sensitivity.py --flavor Umu --mass 1.0 --workers 1

# Check output
cat hnl/output/analysis/gargoyle_hnl_sensitivity.csv
```

## Full Analysis Cycle

### Step 1: Generate ctau and BR_vis tables

```bash
python hnl/production/generate_ctau_tables.py
```

This computes HNL proper lifetime (ctau at U^2=1) and visible branching ratio
for all flavors (Ue, Umu, Utau) over the full mass grid (0.2-10 GeV).

Output:
- `hnl/output/ctau/ctau_{Ue,Umu,Utau}.dat`
- `hnl/output/ctau/br_vis_{Ue,Umu,Utau}.dat`

To skip BR_vis (faster): add `--no-br-vis`

### Step 2: Generate FONLL meson -> HNL 4-vectors

```bash
# All channels (B, D, Bc) and all flavors — ~30-60 min
python hnl/production/decay_engine/generate_meson_csvs.py

# Or run per-channel for more control:
python hnl/production/decay_engine/generate_meson_csvs.py --channel Bmeson --flavor Ue Umu
python hnl/production/decay_engine/generate_meson_csvs.py --channel Dmeson --flavor Ue Umu
python hnl/production/decay_engine/generate_meson_csvs.py --channel Bc --flavor Ue Umu
```

Options:
- `--n-pool 100000` (default) — number of parent meson 4-vectors per pool
- `--seed 42` (default) — RNG seed for reproducibility
- `--masses 0.5 1.0 2.0` — run specific masses only

Output: `hnl/output/llp_4vectors/{flavor}/{Bmeson,Dmeson,Bc}/mN_{mass}.csv`

### Step 3: Generate tau -> HNL 4-vectors (requires MadGraph)

```bash
# Full run: MG5 tau pool + decay for all flavors
python hnl/production/madgraph/run_tau_production.py

# Or skip MG5 if tau pool already cached:
python hnl/production/madgraph/run_tau_production.py --skip-mg5

# Quick test (1k events, Utau only):
python hnl/production/madgraph/run_tau_production.py --test
```

Output: `hnl/output/llp_4vectors/{flavor}/tau/mN_{mass}.csv`

### Step 4: Generate W/Z -> l N (requires MadGraph + SM_HeavyN model)

This channel is relevant for m_N > 5 GeV. Each (flavor, mass) point runs a
separate MG5 job, so this is the slowest step.

```bash
# Full run: all flavors, all masses
python hnl/production/madgraph/run_wz_production.py

# Run specific masses (e.g., above FONLL range):
python hnl/production/madgraph/run_wz_production.py \
    --flavor Ue Umu --masses 5.0 6.0 7.0 8.0 9.0 10.0

# Quick test (1k events, single point):
python hnl/production/madgraph/run_wz_production.py --test
```

Options:
- `--nb-core 4` — CPU cores per MG5 job
- `--min-mass 5.0` — skip masses below this

Output: `hnl/output/llp_4vectors/{flavor}/WZ/mN_{mass}.csv`

### Step 5: Combine all channels

```bash
python hnl/production/combine_channels.py

# Or specific flavors:
python hnl/production/combine_channels.py --flavor Ue Umu
```

Merges Bmeson + Dmeson + Bc + tau + WZ into a single CSV per (flavor, mass).
Missing channels are silently skipped.

Output: `hnl/output/llp_4vectors/{flavor}/combined/mN_{mass}.csv`

### Step 6: Run sensitivity analysis

```bash
# Default: Ue + Umu, 6 workers, full FONLL mass grid (m_N <= 5 GeV)
python hnl/analysis/run_sensitivity.py

# All options:
python hnl/analysis/run_sensitivity.py \
    --flavor Ue Umu \
    --workers 6 \
    --diagnostics           # save N_signal vs U^2 plots per mass point

# Single mass for debugging:
python hnl/analysis/run_sensitivity.py --flavor Umu --mass 1.0 --workers 1

# Re-plot from existing results (no recomputation):
python hnl/analysis/run_sensitivity.py --plot-only
```

Output:
- `hnl/output/analysis/gargoyle_hnl_sensitivity.csv` — results table
- `hnl/output/analysis/gargoyle_hnl_exclusion.pdf` — 1x3 panel exclusion plot
- `hnl/output/analysis/gargoyle_hnl_exclusion.png`
- `hnl/output/analysis/run_metadata.json`
- `hnl/output/analysis/scan_status.json` — live progress (see Monitoring below)

### Step 7: Diagnostic plots

```bash
python hnl/analysis/plot_diagnostics.py

# Single flavor:
python hnl/analysis/plot_diagnostics.py --flavor Umu
```

Produces: BR_vis, ctau, event counts, weight sums, sensitivity details.

Output: `hnl/output/analysis/diagnostics/diag_*.png`

### Step 8: Run tests

```bash
python -m pytest hnl/tests/ -v
```

## Monitoring a running scan

The sensitivity scan writes a JSON status file that updates as mass points complete.

```bash
# Live dashboard (updates every 5s)
watch -n 5 cat hnl/output/analysis/scan_status.json

# One-shot check
cat hnl/output/analysis/scan_status.json

# Follow latest log output
# (run_sensitivity.py prints progress to stdout)
```

Status file fields: `ts`, `done`, `total`, `n_sensitive`, `elapsed_s`, `eta_s`.

## Directory layout after a full run

```
hnl/output/
  ctau/
    ctau_Ue.dat, ctau_Umu.dat, ctau_Utau.dat
    br_vis_Ue.dat, br_vis_Umu.dat, br_vis_Utau.dat
  llp_4vectors/
    Ue/  (and Umu/, Utau/)
      Bmeson/mN_*.csv
      Dmeson/mN_*.csv
      Bc/mN_*.csv
      tau/mN_*.csv
      WZ/mN_*.csv
      combined/mN_*.csv
  analysis/
    gargoyle_hnl_sensitivity.csv
    gargoyle_hnl_exclusion.pdf
    gargoyle_hnl_exclusion.png
    run_metadata.json
    scan_status.json
    geometry_cache/{flavor}/geom_*.npz
    diagnostics/diag_*.png
```

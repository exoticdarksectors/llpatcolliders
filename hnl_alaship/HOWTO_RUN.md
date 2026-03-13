# HOWTO: Run the FairShip-based HNL Pipeline

## Prerequisites

```bash
conda activate llpatcolliders
```

ROOT with Pythia8 must be available (`python -c "import ROOT; print(ROOT.gROOT.GetVersion())"`).

All commands below assume you are inside `hnl_alaship/`.
All paths are `__file__`-based, so running from the repo root works too — just
prefix each script path with `hnl_alaship/`.

## Step 1: Generate meson -> HNL 4-vectors (FONLL channels)

```bash
conda run -n llpatcolliders python production/decay_engine/generate_meson_csvs.py \
    --flavor Ue Umu Utau --channel all
```

## Step 2: Generate tau -> HNL (requires Docker + MG5)

```bash
conda run -n llpatcolliders python production/madgraph/run_tau_production.py \
    --flavor Ue Umu Utau
```

## Step 3: Generate W/Z -> HNL (requires Docker + MG5)

```bash
conda run -n llpatcolliders python production/madgraph/run_wz_production.py \
    --flavor Ue Umu Utau
```

## Step 4: Combine channels

```bash
conda run -n llpatcolliders python production/combine_channels.py \
    --flavor Ue Umu Utau
```

## Step 5: Generate ctau tables (FairShip backend)

```bash
conda run -n llpatcolliders python production/generate_ctau_tables.py \
    --flavor Ue Umu Utau
```

## Step 6: Run sensitivity analysis

```bash
# Single flavor (fast)
conda run -n llpatcolliders python analysis/run_sensitivity.py --flavor Umu --workers 1

# All flavors
conda run -n llpatcolliders python analysis/run_sensitivity.py --flavor Ue Umu Utau --workers 1

# Reduced sampling (faster, good for validation)
conda run -n llpatcolliders python analysis/run_sensitivity.py --flavor Ue Umu Utau --workers 1 \
    --mother-samples 1000 --position-bins 6 --decays-per-bin 3

# Re-plot only (from existing CSV)
conda run -n llpatcolliders python analysis/run_sensitivity.py --plot-only
```

**Performance notes:**
- `--workers 1` is recommended on laptops. Multi-worker mode spawns
  subprocesses that each load ROOT + trimesh, causing high memory usage.
- First run builds geometry and decay-acceptance caches (~5-10 min per flavor
  with reduced sampling). Subsequent runs with the same data hit cache and
  complete in ~1 min for all 3 flavors.
- Reduced sampling (`--mother-samples 1000 --position-bins 6 --decays-per-bin 3`)
  is ~20x faster and sufficient for validation. Use default settings for
  publication-quality results.

## Step 7: Run tests

```bash
conda run -n llpatcolliders python -m unittest discover -s tests -v
```

## Monitoring

Live status during sensitivity scan:
```bash
watch -n2 'cat output/analysis/scan_status.json 2>/dev/null || echo "not started"'
```

## Output

- `output/llp_4vectors/{flavor}/combined/` — combined 4-vector CSVs
- `output/ctau/` — ctau tables
- `output/analysis/gargoyle_hnl_sensitivity.csv` — sensitivity results
- `output/analysis/gargoyle_hnl_exclusion.{pdf,png}` — exclusion plots

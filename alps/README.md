# ALP Portal

Two production channels:

- **Heavy ALP** (PDG 6000113): h -> aa, mass ~ 1-60 GeV, decay a -> mu+ mu-
- **Light ALP** (PDG 9000001): B -> K(*) a via FCNC, mass ~ 0.2-4 GeV

Both cmnd files use BR=1 (efficiency maps); rescale by true BR via `--xsec`.

## Cross-sections

- Heavy ALP: `--xsec 54700` (sigma(gg->h) = 54.7 pb, N3LO, 14 TeV; LHCXSWG)
- Light ALP: `--xsec 373000000` (sigma(pp->bb) ~ 373 ub at LO, Pythia8 HardQCD)

**Note on bb̄ cross-section:** The 373 μb value is the Pythia8 leading-order
`HardQCD:gg2bbbar + qqbar2bbbar` cross-section at 14 TeV. The FONLL NLO+NLL
value is ~495 μb (used by the HNL portal in `hnl_alaship/`). The ~33% difference
is the expected NLO/LO K-factor. Since the light ALP events are generated with
Pythia8 LO kinematics, using the Pythia8 LO cross-section is self-consistent.
For publication, this LO systematic should be noted; applying a K-factor of
~1.33 to the light ALP exclusion would strengthen the limit proportionally.

## Production

```bash
bash generator/produce.sh alps/cmnd/heavy_alp.cmnd 50000 alp_heavy_m15 4 --outdir output/alps
bash generator/produce.sh alps/cmnd/light_alp.cmnd 50000 alp_light_m1 4 50000 --outdir output/alps
```

## Analysis

```bash
python analysis/decayProbPerEvent_2body.py output/alps/data/alp_heavy_m15.csv \
    --xsec 54700 --outdir output/alps --external-dir higgs/external

python analysis/decayProbPerEvent_2body.py output/alps/data/alp_light_m1.csv \
    --xsec 373000000 --outdir output/alps

python analysis/signal_surface_hitmap.py output/alps/data/alp_heavy_m15.csv \
    --outdir output/alps
```

**Note:** `alps/external/` is intentionally empty. The heavy ALP (h->aa) uses the
same PBC benchmark curves (h->SS) as the Higgs portal, located in `higgs/external/`.

## Note on acceptance model

The current analysis uses the analytic 2-body acceptance formula (opening angle
from Lorentz kinematics, separation = theta × d_remaining). This is a good
approximation for 2-body decays at high gamma but does not ray-cast individual
daughters to the detector mesh. For the HNL portal (multi-body decays), a full
3D boost + daughter ray-cast was found to give 15-30% tighter exclusion above
m > 1.5 GeV (see `hnl_alaship/`). For 2-body ALP decays the effect is expected
to be small, but a cross-check with the 3D method at low gamma (high mass) is
recommended as a future validation.

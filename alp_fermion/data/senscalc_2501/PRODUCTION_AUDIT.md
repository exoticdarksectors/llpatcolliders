# BC10 production audit against arXiv:2501.04525

_Audit date: 2026-07-13; status reconciled 2026-07-17._ This is an implementation audit for the GRENDEL
geometry, not a claim that every production mode discussed in the paper must
be included in every experiment. GRENDEL is a transverse detector covering
approximately `|eta| < 0.5`; both the total production probability and the
accepted angular/energy distribution matter.

## Pinned sources

- SensCalc `v.1.3.3`, commit
  `0bca050633aae16e148d47f21840fa07ff4b8724`.
- Official arXiv:2501.04525 companion repository
  <https://github.com/maksymovchynnikov/ALPs-phenomenology>, commit
  `e1443ac1a07ed8a702173e7ead97aa231eedb563` (2026-02-03).
- GKOZ arXiv:2310.03524 for the LHC channel hierarchy and the exclusive
  `B -> K^(i) a` calculation already used by GRENDEL.

The relevant SensCalc source hashes are:

| File | SHA-256 |
|---|---|
| `3. ALP-fermion sensitivity.nb` | `48883fb7ecce956bd73965709ecf44d0cae286e637067c672c1566d736072ef3` |
| `codes/generic.nb` | `95cc528268e41b273182377275bf9eeb00d8cc1cb92d149ddeaa9e72967f0d4a` |
| `codes/LLP distribution/prod-pheno-ALP-fermion.nb` | `e1a7b91a4e82421081d165561a1b7f23e90aefe451efca32a5a3459b1ac0d0f7` |
| `codes/Acceptances/ALP-fermion.nb` | `8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4` |
| `codes/experiments.nb` | `3e60809105f49d8c5274de9a352cd0063bac9cb8feba31188e5e49ca6e468bdf` |
| `Coefficients-ALP-fermion-Lambda=1.-TeV.mx` | `e6885569ee8d33fb5310723ba2f020bf8b4e3506c2b5ad258b5230f9aa981f43` |
| `ProductionProbability-Fragmentation-ALP-fermion.mx` | `4804205cfb14c0d40322fcabb7ec65843023ee44e91f6918603bca0e8835e9eb` |
| `DoubleDistr_ALP-fermion_Fragmentation_LHC.m` | `43116651d526358c446af54a16d43ca9b2233fb42ff7ab4bae6946fe21477a79` |
| `BrRatios-Msquared-LightMesonDecays-ALP-fermion-Lambda=1.-TeV.mx` | `c7468dcf8863a8ac182c2720d34d12b2dd12693700aec032460cd5d678985f71` |
| `ProdProb_ALP-fermion_Bremsstrahlung-AP_LHC.m` | `435022c25a2bb908ffd5c4b12a5645a0d5aa7f6914c49aab5e33fe0b8b60e3d2` |
| `sigmaDrellYan_ALP-fermion_LHC.txt` | `12e013f6df9a65c07514e75d8f1c47d53b53ebc72b21624ab13eb70afaed185d` |

## Channel decisions

| Production mode | GRENDEL status | Decision and evidence |
|---|---|---|
| Exclusive `B` decays | Included | The current code samples `B+` and `B0` from the 14 TeV FONLL spectrum and sums the nine-state kaon tower. This is the GKOZ channel set retained by the 2501 implementation and is more differential than SensCalc's merged `Bcharged -> PiCharged + a` event topology. Keep it. |
| Drell-Yan / gluon fusion | Not included | GKOZ's LHC production plot places this channel roughly five to six orders below the `B` tower over the current GRENDEL island. The tabulation starts at 1.5 GeV, while the canonical island closes at 3.293 GeV and `B` production remains open to about 4.8 GeV. It cannot affect the published island at current precision; document the omission. |
| Old flux-times-mixing | Deliberately excluded | ArXiv:2501.04525 identifies this approximation as chiral-rotation dependent and kinematically ambiguous. SensCalc still exposes `Old-Mixing-Pi0/Eta/EtaPr`; they must not be enabled as a shortcut. |
| Proton bremsstrahlung | Not included | The 2501 calculation is a forward quasi-real approximation with a large theory uncertainty. For a transverse detector it requires a dedicated angular acceptance calculation. It is not justified to add a total-rate reweighting to the FONLL `B` sample. |
| Quark fragmentation | Omitted off pole; pole windows excluded | The decoded generalized-mixing probability and LHC four-vector grid give a central cross section of `7.36e6 pb` at 0.96 GeV, about 52% of the current B-tower value. This point lies inside the eta-prime window where both arXiv:2501.04525 and SensCalc exclude the mixing description. Outside the excluded eta/eta-prime windows, the largest tested central contribution is 0.15% of the B tower (1.0 GeV); omit it at current precision. |
| Light-meson decays | Omitted | Exact decoded BNT branching coefficients are many orders below the rate needed to compete with the B tower. At 0.22 GeV, even the sum of eta-prime three-body modes is below `2e-9` at the reference coupling; multiplying by the entire parent yield before any transverse or decay-kinematic loss gives less than `70 pb`, versus about `1.5e7 pb` accepted from B decays. Pole-window values are not used. |
| `B_s -> phi a` | Not added | This is absent from GKOZ, ALPINIST, SensCalc's channel set, and the cited Boiarska kaon-tower form factors. It is an optional new calculation, not a missing implementation of arXiv:2501.04525. |

## Transverse parent-spectrum check

The five readable SensCalc LHC parent spectra are shape tables normalized by
the sampler; their physical multiplicities per collision are separate entries in
`codes/experiments.nb`. The source-pinned sampler in
`alp_fermion/production_spectra.py` gives the following parent-level fractions
inside `|eta| < 0.5`:

| Parent | Yield / collision | Parent fraction in `|eta| < 0.5` |
|---|---:|---:|
| eta | 3.64 | 0.0596 |
| eta-prime | 0.46 | 0.0653 |
| omega | 4.42 | 0.0555 |
| charged rho | 8.68181 | 0.0546 |
| K-short | 3.1 | 0.0606 |

These are not ALP acceptances: the daughter decay kinematics and the
mass-dependent branching ratios still have to be applied. The source tables
also use a 6.5 TeV beam energy (13 TeV collisions), whereas the GRENDEL FONLL
baseline is 14 TeV. They are suitable for deciding whether a dedicated 14 TeV
calculation is necessary, but not for silently mixing final event weights.

### Two-body daughter check

The revised SensCalc process list contains three two-body light-parent modes:
`Omega-to-ALP-gamma`, `RhoCh-to-ALP-PiCh`, and `KS-to-ALP-Pi0`. A 600,000
event sample per open channel and mass was generated from the pinned parent
spectra with the existing `production.decay_engine.kinematics.decay_2body`
Lorentz sampler. No branching coefficient was applied.

| Process | `m_a` [GeV] | ALP fraction in `|eta| < 0.5` | Central ALPs / collision / unit BR |
|---|---:|---:|---:|
| `omega -> a gamma` | 0.22 | 0.0731 | 0.3232 |
|  | 0.40 | 0.0627 | 0.2773 |
|  | 0.60 | 0.0568 | 0.2511 |
| `rho_charged -> a pi_charged` | 0.22 | 0.0703 | 0.6104 |
|  | 0.40 | 0.0603 | 0.5239 |
|  | 0.60 | 0.0546 | 0.4743 |
| `K_S -> a pi0` | 0.22 | 0.0660 | 0.2045 |
|  | 0.30 | 0.0622 | 0.1928 |

The binomial component of the MC uncertainty on each listed fraction is at
most `0.00034`; interpolation and source-model systematics are separate.
The daughter fractions are comparable to, and at low mass slightly larger
than, the parent fractions. Geometry alone therefore does not reject these
modes. The decoded 2501 branching coefficients do: at `m_a = 0.22 GeV`, the
largest coefficients are the eta-prime three-body modes, whose physical
branching ratios at `1/f = 1e-3 GeV^-1` are `7.10e-10` and `1.29e-9`.
Even assigning them 100% transverse acceptance gives less than `70 pb` after
multiplying by the full eta-prime yield and `72 mb` inelastic cross section,
over five orders below the accepted B-tower cross section. The remaining
light-parent modes are smaller. Their three-body matrix elements therefore do
not need to enter the GRENDEL event generator at present precision.

For an indicative rate threshold, the earlier 120,000-pool audit sample has
a central (`|eta| < 0.5`) reference cross section of `15.0--15.8 microbarn`
over `m_a = 0.22--0.60 GeV` at `1/f_BNT = 1e-3 GeV^-1`. Dividing by the
SensCalc LHC inelastic cross section (`72 mb`) and the per-unit-BR central
yields above, a single light-parent channel would need approximately:

- `BR(omega -> a gamma) = 6.8e-4--8.3e-4`;
- `BR(rho_charged -> a pi_charged) = 3.6e-4--4.4e-4`;
- `BR(K_S -> a pi0) = 1.1e-3` at `m_a = 0.22 GeV`;

to equal the accepted B-tower yield. Both the B branching ratios and revised
light-meson branching ratios scale as `(1/f_BNT)^2`, so this comparison is
coupling independent. It is only an orientation threshold because it compares
the 14 TeV FONLL B baseline with SensCalc's 13 TeV light-parent shapes. If a
decoded coefficient approaches the threshold, a common-energy production
calculation is required.

Reproduce the source check and parent-level summary with:

```bash
python alp_fermion/tools/export_senscalc_2501_production.py \
  --check-only /path/to/SensCalc
python -m alp_fermion.production_spectra /path/to/SensCalc \
  --two-body-events 600000
```

After Wolfram Engine is activated, omit `--check-only` to decode the binary
probability, matrix-element, bremsstrahlung, and fragmentation-grid assets into
the gitignored `alp_fermion/tmp/senscalc_2501_production/` directory.

SensCalc applies the raw meson-decay, fragmentation, and effective-production
coefficients as `(Fpi*g_Y/(2*v_h))^2 * coefficient_raw`, with
`Fpi = 0.093 GeV`. On the GRENDEL BNT axis, `1/f_BNT = g_Y/v_h`; consequently
the exported probability or branching coefficient is
`coefficient_bnt = Fpi^2 * coefficient_raw / 4`, and the physical probability
is `(1/f_BNT)^2 * coefficient_bnt`. The exporter writes only explicitly named
BNT coefficient tables; raw Drell-Yan cross-section coefficients remain marked
as raw and carry their separate normalization formula in metadata.

The pinned binary sources record `Windows-x86-64` in their MX headers. Wolfram
MX is system-dependent, so the source system ID is reported during preflight
and the decoder `$SystemID` is retained in the export manifest. A successful
cross-platform import must be demonstrated before any generated table is
adopted.

### Fragmentation acceptance and light-meson poles

The decoded fragmentation distribution is a normalized 232 by 183
`(theta, energy)` grid at each of 34 mass nodes. Integrating the same
log-linear interpolant used by SensCalc gives:

| `m_a` [GeV] | Central fragmentation [pb] | Central B tower [pb] | Ratio |
|---:|---:|---:|---:|
| 0.50 | `1.12e3` | `1.52e7` | `7.4e-5` |
| 0.60 | `6.70e2` | `1.50e7` | `4.5e-5` |
| 0.90 | `6.01e3` | `1.43e7` | `4.2e-4` |
| 0.96 | `7.36e6` | `1.41e7` | `0.52` |
| 1.00 | `2.10e4` | `1.44e7` | `1.5e-3` |
| 1.50 | `4.98e2` | `1.29e7` | `3.9e-5` |
| 2.00 | `1.89e2` | `1.16e7` | `1.6e-5` |

The apparent 0.96 GeV exception is not a physical point to retain. The paper
states that its ALP-meson diagonalization breaks down near the light-meson
poles, and the pinned SensCalc analysis implements
`0.538 < m_a < 0.555 GeV` and `0.94 < m_a < 0.974 GeV` as excluded windows.
GRENDEL now uses the same windows in `model.LIGHT_MESON_RESONANCE_WINDOWS`;
production and template generation skip them, sensitivity writes explicit
non-sensitive marker rows, and the plot is split rather than interpolated
through them. The lower pion window is below the BC10 scan threshold.

Reproduce the fragmentation table with:

```bash
python -m alp_fermion.production_spectra /path/to/SensCalc \
  --production-export-dir alp_fermion/tmp/senscalc_2501_production \
  --fragmentation-mass 0.5 --fragmentation-mass 0.6 \
  --fragmentation-mass 0.9 --fragmentation-mass 0.96 \
  --fragmentation-mass 1.0 --fragmentation-mass 1.5 \
  --fragmentation-mass 2.0
```

## SensCalc integration boundary

ArXiv:2305.13383 is the SensCalc methodology paper. SensCalc is beneficial
here as a pinned source of decay tables, branching ratios, matrix elements,
and production cross-checks. It should not replace the GRENDEL detector and
reconstruction simulation.

The `v.1.3.3` ALP-fermion sensitivity notebook exposes meson decays,
bremsstrahlung, Drell-Yan, and fragmentation when the revised `2501.04525`
production description is selected. Legacy neutral-meson mixing is retained
only as the alternate old description and must not be combined with the
revised channels. GRENDEL imports verified physics inputs and applies its own
geometry, decay sampling, stable-particle handling, and four-hit track
selection rather than substituting the SensCalc detector acceptance.

## Final-campaign gate and outcome

The production/template/sensitivity campaign was gated on:

1. the 2501 decay widths, branching ratios, and squared matrix elements have
   been exported and validated together, and the full branching mixture is
   sampled in stable-particle templates;
2. the excluded light-meson pole windows have explicit marker rows and are not
   bridged by the plotting/publication path;
3. any retained production mode has a real four-vector sample rather than a
   total-rate-only correction.

All three gates were satisfied before the canonical campaign. The published
result uses a 1.2-million-event high-pT importance production pool, 20,000
Pythia templates at each supported mass, a 145-point scan, and the exact
exported three-body matrix-element reweighting. Flat primary phase space is no
longer the central three-body model. The remaining two-gluon light-quark
surrogate is explicit and propagated through named template variations; it
does not require another central production campaign.

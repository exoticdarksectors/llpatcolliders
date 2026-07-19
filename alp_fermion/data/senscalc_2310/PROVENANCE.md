# SensCalc 2310 structural decay model

This directory is the exact fermionic-ALP decay model shipped by SensCalc
v1.3.3 for arXiv:2310.03524. It is retained as a **one-sided structural model
comparison** to the arXiv:2501.04525 central model, not as a calibrated
`+/-1 sigma` uncertainty and not as part of the pointwise FONLL/decay-acceptance/
`C_bs` halo.

## Pinned upstream inputs

- Repository: `https://github.com/maksymovchynnikov/SensCalc`
- Tag: `v.1.3.3`
- Commit: `0bca050633aae16e148d47f21840fa07ff4b8724`
- Directory: `phenomenology/ALP-fermion/decay widths/`

| Input | SHA-256 |
|---|---|
| `Widths-model-ALP-fermion-scale-1000.-GeV-2310.03524.m` | `40b51e8d35a8edbb3027cfcc3a2bb369faff8fc009a958863d6264f227410e93` |
| `Br-ratios-SensCalc-model-ALP-fermion-scale-1000.-GeV-2310.03524.m` | `130aaa7e95c50a3fea8f43a537fcddfe2bf2d4f676a6110bd3f88097da2826cd` |
| `Matrix-elements-squared-model-ALP-fermion-scale-1000.-GeV-2310.03524.m` | `89ef54601fc995418e9967d8269862e27cb86852a847bff1669e8948aa12dfca` |
| `codes/Acceptances/ALP-fermion.nb` | `8d205b65da456fb1a8fac8d1803eaa73ae839128a60d0624eeea17cacee5d4b4` |

The MX files were decoded by Wolfram Engine 15.0 on `MacOSX-ARM64` through the
same hash-verifying and schema-validating exporter used for the 2501 central
model:

```bash
python alp_fermion/tools/export_senscalc_2501.py /path/to/SensCalc \
  --decay-model 2310_structural
```

The raw upstream width convention is
`Gamma = (gY/(2 vH))^2 coefficient`. Every non-mass column in
`widths_bnt.csv` is exactly the raw coefficient divided by four, giving
`Gamma = (1/f_BNT)^2 coefficient` for `1/f_BNT = gY/vH`.

## Exact content and checks

The export contains all 32 exclusive branching channels and all 18 shipped
three-body squared matrix elements. Channel and matrix-element process order is
identical to the 2501 bundle, but the numerical widths, branching functions,
and matrix-element expressions are read from the 2310 MX files independently.
No matrix element or branching channel is approximated or borrowed from 2501.

Numerical anchors in the BNT coefficient table are:

| mass [GeV] | `Gamma_total/(1/f)^2` [GeV^3] | `BR(a -> mu mu)` |
|---:|---:|---:|
| 1.00 | `1.4807487381396991e-3` | `0.32305858705921875` |
| 1.30 | `1.2317291370015708e-2` | `0.05096158623008926` |
| 1.45 | `2.572263805341004e-2` | `0.027290142885435482` |
| 2.17 | `1.874645871303943e-2` | `0.05637073587242273` |

At 1.30 GeV the 2501 central total-width coefficient is 88.79 times the 2310
value. From the common perturbative switch at 2.17 GeV through 10 GeV, the
shipped 2310 and 2501 total widths and exclusive branching ratios agree
numerically. This makes the alternative specifically diagnostic of the
heavy-pseudoscalar structure that changes the lower-mass island topology; it
does not define a new high-mass perturbative-width error.

The exact exported 2310 `matrix_element_004` evaluates to
`2.774159983255917` at `mLLP=1.2`, `E1=0.4`, `E3=0.3` GeV in Wolfram Engine.
The Python evaluator is tested against this anchor.

As for the central templates, external Pythia cannot hadronize an isolated
color-singlet gluon pair. The `a -> gg` mode therefore retains the equal
`u/d/s` surrogate for this structural contour. The separate `gg_u`, `gg_d`,
and `gg_s` variations remain acceptance diagnostics in the ordinary halo and
are not folded into this one-sided model comparison.

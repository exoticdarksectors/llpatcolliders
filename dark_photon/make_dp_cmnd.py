#!/usr/bin/env python3
"""
make_dp_cmnd.py  --  generate Pythia8 cmnd files for dark photon A' production.

Three production modes (--production):
  heavy      (default) exotic Higgs decay  h → A' A'
  drell_yan  Drell-Yan  qq̄ → γ* → A'  (s-channel, kinetic mixing)
  meson      meson decay  η → A' γ  or  ω → A' π⁰

Usage (from dark_photon/generator/):

  # Heavy DP h→A'A':
  python make_dp_cmnd.py --mass 0.5  --outfile heavy_dp_m05.cmnd
  python make_dp_cmnd.py --mass 1.0  --outfile heavy_dp_m1.cmnd
  python make_dp_cmnd.py --mass 15.0 --outfile heavy_dp_m15.cmnd

  # Drell-Yan (meson channels closed at m=1 GeV):
  python make_dp_cmnd.py --production drell_yan --mass 1.0 --outfile dy_dp_m1.cmnd

  # Meson decay (sub-GeV; kinematic limits: eta < 548 MeV, omega < 648 MeV):
  python make_dp_cmnd.py --production meson --mass 0.5 --parent omega --outfile meson_dp_omega_m05.cmnd
  python make_dp_cmnd.py --production meson --mass 0.5 --parent eta   --outfile meson_dp_eta_m05.cmnd

BR source: br_tables/dp_brs_deliver.csv  (DeLiVeR, VMD + PDG R-ratio, <= 1.7 GeV)
           Above 1.7 GeV: perturbative QCD (R-ratio with alpha_s correction).
           Regenerate the table: python generate_br_table.py --deliver-path <clone>
"""

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
DEFAULT_BR_TABLE = HERE / "br_tables" / "dp_brs_deliver.csv"

# Particle thresholds (GeV): channel is kinematically open only above 2*m_daughter
M_PI  = 0.13957   # charged pion
M_PI0 = 0.13498   # neutral pion
M_MU  = 0.10566   # muon
M_TAU = 1.77686   # tau
M_K   = 0.49368   # charged kaon
M_K0  = 0.49765   # neutral kaon
M_ETA = 0.54785   # eta
M_P   = 0.93827   # proton

# Decay channel mapping: (csv_column, [Pythia8 PDG IDs], threshold_GeV, comment)
# Channels whose BR falls below MIN_BR are merged into pi+pi-.
CHANNELS = [
    ("ee",       [11, -11],                  2 * 0.000511, "e+ e-"),
    ("mumu",     [13, -13],                  2 * M_MU,     "mu+ mu-"),
    ("tau",      [15, -15],                  2 * M_TAU,    "tau+ tau-"),
    ("pipi",     [211, -211],                2 * M_PI,     "pi+ pi-"),
    ("pi3",      [111, 211, -211],           M_PI0 + 2 * M_PI,  "pi0 pi+ pi-"),
    ("pi4c",     [211, 211, -211, -211],     4 * M_PI,          "pi+ pi+ pi- pi-"),
    ("pi4n",     [111, 111, 211, -211],      2 * M_PI0 + 2 * M_PI, "pi0 pi0 pi+ pi-"),
    ("KKc",      [321, -321],               2 * M_K,       "K+ K-"),
    ("KKn",      [311, -311],               2 * M_K0,      "K0 K0bar"),
    ("PiGamma",  [111, 22],                  M_PI0,        "pi0 gamma"),
    ("EtaGamma", [221, 22],                  M_ETA,        "eta gamma"),
    ("ppbar",    [2212, -2212],             2 * M_P,        "p pbar"),
]

# Minor hadronic residual column in CSV always merged into pipi
CSV_MERGE_INTO_PIPI = ["had_other"]

# Channels with BR < this fraction are also merged into pipi
MIN_BR = 0.005


# --------------------------------------------------------------------------
# Drell-Yan coupling constants
# --------------------------------------------------------------------------
SIN2TW = 0.2312          # sin^2(theta_W) at M_Z (Pythia8 default)
COS2TW = 1.0 - SIN2TW
# Correction factor: sigma_real = sigma_Pythia8 * K_NORM
# Pythia8 NewGaugeBoson Z' uses internal coupling g_Z = e/(2*sin*cos).
# For a dark photon, the coupling at each vertex is eps*e*Q_f.
# K_NORM maps Pythia8's Z' normalization to the EM coupling convention.
# Since DY mode has zero PX56 sensitivity, this affects only the
# diagnostic cross-section printout, not the final exclusion.
K_NORM = (16.0 * SIN2TW * COS2TW) ** 2   # ~ 8.09

# EM charges (vector couplings for dark photon, axial = 0)
VU  =  2.0 / 3.0    # up-type quark
VD  = -1.0 / 3.0    # down-type quark
VE  = -1.0           # charged lepton
VNU =  0.0           # neutrino (no EM coupling)


# --------------------------------------------------------------------------
# Meson decay parent definitions
# --------------------------------------------------------------------------
PARENT_MESONS = {
    "eta": {
        "pdg": 221,
        "name": "η",
        "decay_desc": "η → A' γ",
        "daughters_pdg": "6000115 22",   # A' γ
        "daughters_comment": "A' gamma",
        "mass": 0.54785,
        "threshold_comment": "m_A' < m_η = 548 MeV",
    },
    "omega": {
        "pdg": 223,
        "name": "ω",
        "decay_desc": "ω → π⁰ A'",
        "daughters_pdg": "111 6000115",  # π⁰ A'
        "daughters_comment": "pi0 A'",
        "mass": 0.78266,
        "threshold_comment": "m_A' < m_ω - m_π⁰ = 648 MeV",
        "max_dp_mass": 0.78266 - 0.13498,  # m_ω - m_π⁰
    },
}


# --------------------------------------------------------------------------
# Perturbative QCD BRs for masses above the DeLiVeR table (> ~1.7 GeV)
# --------------------------------------------------------------------------
from dp_meson_brs import perturbative_brs as _perturbative_brs


# --------------------------------------------------------------------------
# CSV table loading + smoke checks
# --------------------------------------------------------------------------
def _load_table(path):
    cols = None
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if cols is None:
                cols = [c.strip() for c in line.split(",")]
                continue
            rows.append([float(x) for x in line.split(",")])
    data = np.array(rows)
    return cols, data


def _smoke_check_table(cols, data, tol=0.01):
    """
    Sanity checks on the loaded BR table.  Raises ValueError on failure.

    Checks:
      1. Required columns are present.
      2. Each row: BR(ee) + BR(mumu) + BR(tau) + BR(BRqcd) in [1-tol, 1+tol].
      3. Mass column is strictly increasing.
    """
    required = ["mass_GeV", "ee", "mumu", "tau", "BRqcd"]
    for col in required:
        if col not in cols:
            raise ValueError(
                f"BR table missing required column {col!r}. "
                f"Regenerate with generate_br_table.py."
            )
    ci = {c: i for i, c in enumerate(cols)}
    for row in data:
        m   = row[ci["mass_GeV"]]
        tot = row[ci["ee"]] + row[ci["mumu"]] + row[ci["tau"]] + row[ci["BRqcd"]]
        if not (1.0 - tol <= tot <= 1.0 + tol):
            raise ValueError(
                f"BR table row at m={m:.3f} GeV: "
                f"ee+mumu+tau+BRqcd = {tot:.5f}  (expected 1 ± {tol}). "
                f"Table may be corrupted; regenerate with generate_br_table.py."
            )
    masses = data[:, ci["mass_GeV"]]
    if not np.all(np.diff(masses) > 0):
        raise ValueError("BR table mass column is not strictly increasing.")
    print(f"  [check] BR table OK: {len(data)} rows, "
          f"m={masses[0]:.2f}–{masses[-1]:.2f} GeV, "
          f"BR sums 1.0 ± {tol} ✓")


def _interp(mass, cols, data):
    masses = data[:, 0]
    return {cols[j]: float(np.interp(mass, masses, data[:, j]))
            for j in range(len(cols))}


# --------------------------------------------------------------------------
# Channel assembly
# --------------------------------------------------------------------------
def _assemble_channels(mass, br_dict, use_perturbative=False):
    """
    Build list of (br, [pdg_ids], comment) from br_dict, applying thresholds
    and merging minor channels into pi+pi-.
    Returns normalised list.
    """
    if use_perturbative:
        return _assemble_perturbative(mass, br_dict)

    pipi_extra = 0.0
    out = []

    # Absorb explicit merge columns
    for col in CSV_MERGE_INTO_PIPI:
        pipi_extra += br_dict.get(col, 0.0)

    for col, pdgs, thresh, comment in CHANNELS:
        if mass < thresh:
            continue
        br = br_dict.get(col, 0.0)
        if br < MIN_BR:
            pipi_extra += br
        else:
            out.append([br, pdgs, comment])

    # Add pipi_extra to pi+pi- channel if present, else add a new one
    added = False
    for entry in out:
        if entry[1] == [211, -211]:
            entry[0] += pipi_extra
            added = True
            break
    if not added and pipi_extra > 0.0 and mass >= 2 * M_PI:
        out.append([pipi_extra, [211, -211], "pi+ pi- (minor had. residual)"])

    total = sum(e[0] for e in out)
    if total <= 0:
        raise ValueError(f"No open channels at m = {mass} GeV")

    return [(br / total, pdgs, comment) for br, pdgs, comment in out]


def _assemble_perturbative(mass, br_dict):
    """Build channel list from perturbative QCD BRs (quark-level)."""
    # br_dict is from _perturbative_brs()
    out = []

    lep_map = [
        ("ee",   [11, -11],    "e+ e-",   0.0),
        ("mumu", [13, -13],    "mu+ mu-", 2 * M_MU),
        ("tau",  [15, -15],    "tau+ tau-", 2 * M_TAU),
    ]
    for col, pdgs, comment, thresh in lep_map:
        br = br_dict.get(col, 0.0)
        if mass > thresh and br > 0:
            out.append([br, pdgs, comment, 0])   # meMode=0: phase-space

    quark_map = [
        ("uu_q", [2, -2],  "u ubar",  0.0),
        ("dd_q", [1, -1],  "d dbar",  0.0),
        ("ss_q", [3, -3],  "s sbar",  0.0),
        ("cc_q", [4, -4],  "c cbar",  2 * 1.27),
        ("bb_q", [5, -5],  "b bbar",  2 * 4.18),
    ]
    for col, pdgs, comment, thresh in quark_map:
        br = br_dict.get(col, 0.0)
        if mass > thresh and br > MIN_BR:
            out.append([br, pdgs, comment, 101])   # meMode=101: string-frag qq-bar

    total = sum(e[0] for e in out)
    if total <= 0:
        raise ValueError(f"No perturbative channels at m = {mass} GeV")

    return [(br / total, pdgs, comment, memode)
            for br, pdgs, comment, memode in out]


# --------------------------------------------------------------------------
# Shared BR channel loader (used by all production modes)
# --------------------------------------------------------------------------
def _load_br_channels(mass, br_table, perturb_above):
    """Load A' decay BRs and return (channels, br_source)."""
    if mass > perturb_above:
        br_dict = _perturbative_brs(mass)
        br_source = (f"perturbative QCD R-ratio with NLO alpha_s correction "
                     f"(m={mass} GeV >> resonance region)")
        channels = _assemble_perturbative(mass, br_dict)
    else:
        if not Path(br_table).exists():
            raise FileNotFoundError(
                f"BR table not found: {br_table}\n"
                f"Regenerate with: python generate_br_table.py "
                f"--deliver-path <DeLiVeR_clone> --outfile {br_table}"
            )
        cols, data = _load_table(br_table)
        _smoke_check_table(cols, data)
        table_mmax = float(data[-1, 0])
        if mass > table_mmax:
            raise ValueError(
                f"Mass {mass} GeV exceeds table range ({table_mmax} GeV). "
                f"Use --perturb-above {table_mmax} or regenerate the table with larger mmax."
            )
        br_dict = _interp(mass, cols, data)
        br_source = (f"DeLiVeR VMD + PDG R-ratio (arXiv:2201.01788), "
                     f"interpolated at m={mass} GeV")
        channels = _assemble_channels(mass, br_dict, use_perturbative=False)

    br_sum = sum(entry[0] for entry in channels)
    if abs(br_sum - 1.0) > 1e-6:
        raise RuntimeError(
            f"Assembled channel BRs sum to {br_sum:.8f}, expected 1.0. Bug in normalisation."
        )
    return channels, br_source


# --------------------------------------------------------------------------
# Heavy DP cmnd template and writer
# --------------------------------------------------------------------------
_CMND_TEMPLATE = """\
! Command file: heavy dark photon production via exotic Higgs decay  h -> A' A'
! A' = dark photon (massive vector boson, spin-1, neutral, uncoloured)
! Kinetic mixing portal: epsilon coupling to SM photon.
! Regime: few GeV to ~60 GeV  (Higgs portal coupling)
! Generated by make_dp_cmnd.py — do not edit manually.
! BR source: {br_source}

! ---- 1) Beam & collision -----------------------------------------------
Beams:idA = 2212                   ! proton
Beams:idB = 2212                   ! proton
Beams:eCM = {ecm:.0f}.                  ! CM energy in GeV ({ecm_label})

! ---- 2) SM Higgs production (all inclusive channels) -------------------
HiggsSM:all = on
25:onMode = off                    ! switch off all SM Higgs decays ...
25:m0 = 125                        ! ... and replace with h -> A' A' below
25:mMin = 125
25:mWidth = 0.004                  ! SM-like narrow width

! Add exotic decay  h -> A' A'  (BR set to 1; rescale yield in analysis)
25:addChannel = 1 1.0 100 {pdgid} {pdgid}
25:onIfAny    = {pdgid}            ! keep event only if an A' is present

! ---- 3) Dark photon particle definition --------------------------------
! Self-conjugate massive vector: "Aprime" for both particle and anti-particle.
! spinType=3 (spin-1 vector), chargeType=0 (neutral), colType=0 (no colour).
{pdgid}:new       = Aprime Aprime 3 0 0
{pdgid}:m0        = {mass}         ! A' mass in GeV  [scan: 0.2 -- 60 GeV]
{pdgid}:tau0      = {tau0}         ! proper decay length c*tau_0  [mm]
{pdgid}:isResonance = off           ! off = decay after event generation (LLP treatment)
{pdgid}:mayDecay  = on

! Override generator default LLP PDG ID (hardcoded as 6000113 in generator.cc)
LLP:pdgId = {pdgid}

! ---- A' decay channels (kinetic mixing, m={mass} GeV) ----------------------
! BRs from {br_source}
! Channels below {min_br_pct:.1f}% of total merged into pi+pi-; sum normalised to 1.
{channel_lines}
! ---- 4) Run settings ---------------------------------------------------
Init:showChangedSettings    = on
Init:showAllSettings        = off
Init:showChangedParticleData = on
Init:showAllParticleData    = off
Next:numberCount    = 1000
Next:numberShowLHA  = 1
Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent   = 1
Stat:showPartonLevel = on

! ---- 5) Output ---------------------------------------------------------
Main:numberOfEvents = 10000
Main:writeLog = on
"""


def write_cmnd(mass, tau0, pdgid, channels, outfile, br_source, ecm=13600.):
    lines = []
    for entry in channels:
        if len(entry) == 3:
            br, pdgs, comment = entry
            memode = 0
        else:
            br, pdgs, comment, memode = entry
        pdg_str = "  ".join(str(p) for p in pdgs)
        lines.append(
            f"{pdgid}:addChannel = 1 {br:.6f}  {memode}  {pdg_str}    ! {comment}"
        )

    ecm_label = "HL-LHC" if ecm >= 13900 else "LHC Run 3"
    content = _CMND_TEMPLATE.format(
        pdgid=pdgid,
        mass=f"{mass:.4g}",
        tau0=f"{tau0:.4g}",
        br_source=br_source,
        min_br_pct=MIN_BR * 100,
        channel_lines="\n".join(lines),
        ecm=ecm,
        ecm_label=ecm_label,
    )
    with open(outfile, "w") as fh:
        fh.write(content)

    print(f"Written {outfile}  (m={mass} GeV, {len(channels)} channels)")
    for entry in channels:
        br = entry[0]
        pdgs = entry[1]
        comment = entry[2]
        print(f"  {br:.4f}  {' '.join(str(p) for p in pdgs)}  ! {comment}")


# --------------------------------------------------------------------------
# Drell-Yan cmnd template and writer
# --------------------------------------------------------------------------
_CMND_TEMPLATE_DY = """\
! Command file: light dark photon production via Drell-Yan  qq-bar -> A'
! A' = dark photon (massive vector boson, spin-1, neutral, uncoloured)
! Kinetic mixing portal: epsilon coupling to SM photon.
! Production: pp -> qq-bar -> gamma* -> A' (Drell-Yan, s-channel)
! Cross-section at eps^2 = 1; multiply by eps^2 in analysis.
! Generated by make_dp_cmnd.py — do not edit manually.
! BR source (A' decay): {br_source}
!
! COUPLING NORMALISATION:
!   Pythia8 NewGaugeBoson uses internal g_Z = e/(sin*cos) normalisation.
!   sigma_Pythia = sigma_real / {knorm:.4f}.
!   Apply correction: --xsec = sigma_from_log * {knorm:.4f}

! ---- 1) Beam & collision -----------------------------------------------
Beams:idA = 2212                   ! proton
Beams:idB = 2212                   ! proton
Beams:eCM = {ecm:.0f}.                  ! CM energy in GeV ({ecm_label})

! ---- 2) Drell-Yan production via NewGaugeBoson framework -----------------
! Uses ffbar2gmZZprime (Sigma1, 2->1 process producing PDG 32 = Z').
! gmZmode = 3: pure Z' contribution only (no gamma*/Z interference).
NewGaugeBoson:ffbar2gmZZprime = on
Zprime:gmZmode = 3
! Phase-space cut: default mHatMin = 4 GeV would miss low-mass A'.
PhaseSpace:mHatMin = {mhatmin}

! ---- 3) Coupling setup: EM charges (dark photon) ------------------------
! v_f = Q_f (EM charge), a_f = 0.  Universality copies 1st-gen to all.
Zprime:universality = on
Zprime:vu   = {vu:.6f}               ! +2/3  (up-type EM charge)
Zprime:au   = 0.0
Zprime:vd   = {vd:.6f}               ! -1/3  (down-type EM charge)
Zprime:ad   = 0.0
Zprime:ve   = {ve:.6f}               ! -1    (charged lepton EM charge)
Zprime:ae   = 0.0
Zprime:vnue = {vnu:.6f}              ! 0     (neutrino: no EM coupling)
Zprime:anue = 0.0

! ---- 4) Dark photon (Z' = PDG 32) properties ----------------------------
32:m0         = {mass}             ! A' mass in GeV
32:mWidth     = {mwidth}           ! narrow (placeholder; isResonance=off)
32:mMin       = {mmin}             ! BW sampling lower bound
32:mMax       = {mmax}             ! BW sampling upper bound
32:isResonance = off                ! decay after event generation (LLP treatment)
32:mayDecay   = on
32:tau0       = {tau0}             ! proper decay length c*tau_0 [mm] (placeholder)

! Override generator default LLP PDG ID
LLP:pdgId = 32

! ---- 5) A' decay channels (kinetic mixing, m={mass} GeV) ----------------
! Turn off all default channels, then add VMD BRs.
32:onMode = off
! BRs from {br_source}
! Channels below {min_br_pct:.1f}% of total merged into pi+pi-; sum normalised to 1.
{channel_lines}

! ---- 6) Run settings ----------------------------------------------------
Init:showChangedSettings    = on
Init:showAllSettings        = off
Init:showChangedParticleData = on
Init:showAllParticleData    = off
Next:numberCount    = 1000
Next:numberShowLHA  = 1
Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent   = 1
Stat:showPartonLevel = on

! ---- 7) Output -----------------------------------------------------------
Main:numberOfEvents = 10000        ! overridden by produce.sh batch size
Main:writeLog = on
"""


def _write_cmnd_dy(args, channels, br_source):
    """Write Drell-Yan A' cmnd file (PDG 32 = Z' in NewGaugeBoson framework)."""
    # Build channel lines (Z' = PDG 32)
    lines = []
    for entry in channels:
        if len(entry) == 3:
            br, pdgs, comment = entry
            memode = 0
        else:
            br, pdgs, comment, memode = entry
        pdg_str = "  ".join(str(p) for p in pdgs)
        lines.append(
            f"32:addChannel = 1 {br:.6f}  {memode}  {pdg_str}    ! {comment}"
        )

    # Narrow A' width at eps=1 for BW sampling bounds
    sys.path.insert(0, str(HERE))
    from dp_meson_brs import dp_width_eps1
    gamma0 = dp_width_eps1(args.mass)
    mwidth_str = f"{gamma0:.4e}" if gamma0 > 0 else "0.001"
    mmin_str = f"{max(0.1, args.mass - 10 * gamma0):.4g}"
    mmax_str = f"{args.mass + 10 * gamma0:.4g}"
    mhatmin = max(0.1, args.mass - 50 * gamma0)

    ecm_label = "HL-LHC" if args.ecm >= 13900 else "LHC Run 3"
    content = _CMND_TEMPLATE_DY.format(
        mass=f"{args.mass:.4g}",
        tau0=f"{args.tau0:.4g}",
        mwidth=mwidth_str,
        mmin=mmin_str,
        mmax=mmax_str,
        mhatmin=f"{mhatmin:.4g}",
        br_source=br_source,
        min_br_pct=MIN_BR * 100,
        channel_lines="\n".join(lines),
        ecm=args.ecm,
        ecm_label=ecm_label,
        knorm=K_NORM,
        vu=VU,
        vd=VD,
        ve=VE,
        vnu=VNU,
    )
    with open(args.outfile, "w") as fh:
        fh.write(content)

    print(f"Written {args.outfile}")
    print(f"  Production: Drell-Yan qq-bar -> A' "
          f"(NewGaugeBoson:ffbar2gmZZprime, PDG 32)")
    print(f"  m_A' = {args.mass} GeV, tau0 = {args.tau0} mm")
    print(f"  eCM = {args.ecm:.0f} GeV ({ecm_label})")
    print(f"  Gamma_0(eps=1) = {gamma0:.4e} GeV")
    print(f"  Couplings: vu={VU:.4f}, vd={VD:.4f}, ve={VE:.4f}, vnu={VNU:.4f}")
    print(f"  Normalisation: sigma_real = sigma_Pythia * {K_NORM:.4f}")
    print(f"  A' decay channels ({len(channels)}):")
    for entry in channels:
        br = entry[0]
        comment = entry[2]
        print(f"    {br:.4f}  {comment}")


# --------------------------------------------------------------------------
# Meson decay cmnd template and writer
# --------------------------------------------------------------------------
_CMND_TEMPLATE_MESON = """\
! Command file: light dark photon production via meson decay
! A' = dark photon (massive vector boson, spin-1, neutral, uncoloured)
! Kinetic mixing portal: epsilon coupling to SM photon.
! Production: pp → X (SoftQCD) → {decay_desc}
! Generated by make_dp_cmnd.py — do not edit manually.
! BR source (A' decay): {br_source}

! ---- 1) Beam & collision -----------------------------------------------
Beams:idA = 2212                   ! proton
Beams:idB = 2212                   ! proton
Beams:eCM = {ecm:.0f}.                  ! CM energy in GeV ({ecm_label})

! ---- 2) Soft QCD (minimum-bias, inclusive hadron production) ------------
SoftQCD:nonDiffractive = on

! ---- 3) Dark photon particle definition --------------------------------
! Self-conjugate massive vector: spinType=3 (spin-1), chargeType=0, colType=0.
{pdgid}:new       = Aprime Aprime 3 0 0
{pdgid}:m0        = {mass}         ! A' mass in GeV
{pdgid}:tau0      = {tau0}         ! proper decay length c*tau_0 [mm] (placeholder)
{pdgid}:isResonance = off           ! decay after event generation (LLP treatment)
{pdgid}:mayDecay  = on

! Override generator default LLP PDG ID
LLP:pdgId = {pdgid}

! ---- 4) Force parent meson decay to A' (BR=1 efficiency map) ----------
! {threshold_comment}
! oneChannel wipes ALL SM decays of {meson_name}, forcing {decay_desc} with BR=1.
! The physical BR(ε², m_A') is applied in the analysis.
{meson_pdg}:oneChannel = 1 1.0 0 {decay_daughters}    ! {decay_comment}

! ---- A' decay channels (kinetic mixing, m={mass} GeV) ------------------
! BRs from {br_source}
! Channels below {min_br_pct:.1f}% of total merged into pi+pi-; sum normalised to 1.
{channel_lines}
! ---- 5) Run settings ---------------------------------------------------
Init:showChangedSettings    = on
Init:showAllSettings        = off
Init:showChangedParticleData = on
Init:showAllParticleData    = off
Next:numberCount    = 1000
Next:numberShowLHA  = 1
Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent   = 1
Stat:showPartonLevel = on

! ---- 6) Output ---------------------------------------------------------
Main:numberOfEvents = 100000        ! overridden by produce.sh batch size
Main:writeLog = on
"""


def _write_cmnd_meson(args, channels, br_source):
    """Write meson-decay A' cmnd file."""
    meson = PARENT_MESONS[args.parent]

    # Kinematic threshold check
    max_mass = meson.get("max_dp_mass", meson["mass"])
    if args.mass >= max_mass:
        print(f"ERROR: m_A' = {args.mass} GeV exceeds kinematic threshold "
              f"for {meson['decay_desc']} ({meson['threshold_comment']})")
        sys.exit(1)

    # Build channel lines
    lines = []
    for entry in channels:
        if len(entry) == 3:
            br, pdgs, comment = entry
            memode = 0
        else:
            br, pdgs, comment, memode = entry
        pdg_str = "  ".join(str(p) for p in pdgs)
        lines.append(
            f"{args.pdgid}:addChannel = 1 {br:.6f}  {memode}  {pdg_str}    ! {comment}"
        )

    ecm_label = "HL-LHC" if args.ecm >= 13900 else "LHC Run 3"
    content = _CMND_TEMPLATE_MESON.format(
        pdgid=args.pdgid,
        mass=f"{args.mass:.4g}",
        tau0=f"{args.tau0:.4g}",
        br_source=br_source,
        min_br_pct=MIN_BR * 100,
        channel_lines="\n".join(lines),
        ecm=args.ecm,
        ecm_label=ecm_label,
        meson_name=meson["name"],
        decay_desc=meson["decay_desc"],
        meson_pdg=meson["pdg"],
        decay_daughters=meson["daughters_pdg"],
        decay_comment=meson["daughters_comment"],
        threshold_comment=meson["threshold_comment"],
    )
    with open(args.outfile, "w") as fh:
        fh.write(content)

    print(f"Written {args.outfile}")
    print(f"  Parent: {meson['name']} ({meson['decay_desc']})")
    print(f"  m_A' = {args.mass} GeV, tau0 = {args.tau0} mm")
    print(f"  eCM = {args.ecm:.0f} GeV ({ecm_label})")
    print(f"  A' decay channels ({len(channels)}):")
    for entry in channels:
        br = entry[0]
        comment = entry[2]
        print(f"    {br:.4f}  {comment}")


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--production", choices=["heavy", "drell_yan", "meson"], default="heavy",
                    help="Production mode: heavy (h→A'A'), drell_yan (qq̄→A'), "
                         "meson (η/ω→A'X). Default: heavy (backward-compatible).")
    ap.add_argument("--parent",    choices=["eta", "omega"], default=None,
                    help="Parent meson for --production meson (required for that mode)")
    ap.add_argument("--mass",      type=float, required=True,
                    help="A' mass in GeV")
    ap.add_argument("--tau0",      type=float, default=1e4,
                    help="c*tau_0 in mm  (default 1e4)")
    ap.add_argument("--pdgid",     type=int,   default=6000115,
                    help="PDG ID for A'  (default 6000115; DY mode uses PDG 32 regardless)")
    ap.add_argument("--outfile",   required=True,
                    help="Output .cmnd filename")
    ap.add_argument("--br-table",  default=str(DEFAULT_BR_TABLE),
                    help="Path to BR table CSV  (default: br_tables/dp_brs_deliver.csv)")
    ap.add_argument("--perturb-above", type=float, default=1.7,
                    help="Use perturbative QCD above this mass in GeV  (default 1.7)")
    ap.add_argument("--ecm",          type=float, default=14000.,
                    help="Centre-of-mass energy in GeV (default 14000 = HL-LHC)")
    args = ap.parse_args()

    if args.production == "meson" and args.parent is None:
        ap.error("--parent {eta,omega} is required for --production meson")

    channels, br_source = _load_br_channels(args.mass, args.br_table, args.perturb_above)

    if args.production == "heavy":
        write_cmnd(args.mass, args.tau0, args.pdgid, channels, args.outfile, br_source,
                   ecm=args.ecm)
    elif args.production == "drell_yan":
        _write_cmnd_dy(args, channels, br_source)
    else:  # meson
        _write_cmnd_meson(args, channels, br_source)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
production/decay_engine/generate_meson_csvs.py

Driver for B/D/Bc → HNL CSV generation.

For each (flavor, mass_point):
  1. Load meson 4-vector pool (generated once per quark type)
  2. For each meson species, check kinematic threshold: m_N < m_meson - m_lepton
  3. Compute BR from HNLCalc (get_2body_br + integrate_3body_br) at U²=1
  4. Decay meson → HNL via 2-body or 3-body kinematics
  5. Assign weight: w_i = 2 × σ_FONLL × f_species × BR / N_sample

Factor of 2: FONLL gives (quark + antiquark)/2 → multiply by 2 for total.

Output: output/llp_4vectors/{Ue,Umu,Utau}/{Bmeson,Dmeson,Bc}/mN_{mass}.csv
Format: headerless, 5 columns: weight,E,px,py,pz
"""

import sys
import numpy as np
from pathlib import Path

# Project root
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.constants import (
    MESON_MASSES, LEPTON_MASSES, FLAVOR_TO_LEPTON_PDG,
    FRAG_B, FRAG_C, SIGMA_BC_PB,
)
from production.fonll.fonll_parser import get_sigma_total
from production.fonll.meson_sampler import sample_meson_4vectors
from production.decay_engine.kinematics import decay_2body, decay_3body_flat

# Output base directory
OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"

# Number of meson 4-vectors per pool
N_POOL = 100_000

# Channel labels for output directory
CHANNEL_LABELS = {
    "bottom": "Bmeson",
    "charm": "Dmeson",
    "bc": "Bc",
}

# HNLCalc 3-body channels for production (parent_pdg, daughter_pdg, type)
# type: "pseudo" for pseudoscalar, "vector" for vector meson daughters
THREEBODY_CHANNELS = {
    # D mesons (charm)
    421: [  # D0
        (-321, "pseudo"),   # D0 → K- ℓ+ N
        (-323, "vector"),   # D0 → K*- ℓ+ N
        (-211, "pseudo"),   # D0 → π- ℓ+ N
        (-213, "vector"),   # D0 → ρ- ℓ+ N
    ],
    411: [  # D+
        (-311, "pseudo"),   # D+ → K̄0 ℓ+ N
        (-313, "vector"),   # D+ → K̄*0 ℓ+ N
        (111, "pseudo"),    # D+ → π0 ℓ+ N
        (113, "vector"),    # D+ → ρ0 ℓ+ N
        (221, "pseudo"),    # D+ → η ℓ+ N
        (331, "pseudo"),    # D+ → η' ℓ+ N
    ],
    431: [  # Ds
        (221, "pseudo"),    # Ds → η ℓ+ N
        (331, "pseudo"),    # Ds → η' ℓ+ N
        (311, "pseudo"),    # Ds → K0 ℓ+ N
        (313, "vector"),    # Ds → K*0 ℓ+ N
        (333, "vector"),    # Ds → φ ℓ+ N
    ],
    # B mesons (bottom)
    521: [  # B+
        (-421, "pseudo"),   # B+ → D̄0 ℓ+ N
        (-423, "vector"),   # B+ → D̄*0 ℓ+ N
        (111, "pseudo"),    # B+ → π0 ℓ+ N
        (113, "vector"),    # B+ → ρ0 ℓ+ N
        (221, "pseudo"),    # B+ → η ℓ+ N
        (223, "vector"),    # B+ → ω ℓ+ N
        (331, "pseudo"),    # B+ → η' ℓ+ N
    ],
    511: [  # B0
        (-411, "pseudo"),   # B0 → D- ℓ+ N
        (-413, "vector"),   # B0 → D*- ℓ+ N
        (-211, "pseudo"),   # B0 → π- ℓ+ N
        (-213, "vector"),   # B0 → ρ- ℓ+ N
    ],
    531: [  # Bs
        (-431, "pseudo"),   # Bs → Ds- ℓ+ N
        (-433, "vector"),   # Bs → Ds*- ℓ+ N
        (-321, "pseudo"),   # Bs → K- ℓ+ N
        (-323, "vector"),   # Bs → K*- ℓ+ N
    ],
    541: [  # Bc
        (511, "pseudo"),    # Bc → B0 ℓ+ N
        (531, "pseudo"),    # Bc → Bs ℓ+ N
        (513, "vector"),    # Bc → B*0 ℓ+ N
        (533, "vector"),    # Bc → Bs* ℓ+ N
        (421, "pseudo"),    # Bc → D0 ℓ+ N
        (441, "pseudo"),    # Bc → ηc ℓ+ N
        (423, "vector"),    # Bc → D*0 ℓ+ N
        (443, "vector"),    # Bc → J/ψ ℓ+ N
    ],
}


def _init_hnlcalc(flavor):
    """Initialize HNLCalc with unit coupling for the given flavor."""
    from HNLCalc import HNLCalc
    if flavor == "Ue":
        return HNLCalc(ve=1, vmu=0, vtau=0)
    elif flavor == "Umu":
        return HNLCalc(ve=0, vmu=1, vtau=0)
    elif flavor == "Utau":
        return HNLCalc(ve=0, vmu=0, vtau=1)
    else:
        raise ValueError(f"Unknown flavor: {flavor}")


def _eval_2body_br(hnl, parent_pdg, lepton_pdg, m_N):
    """Evaluate 2-body BR: parent → ℓ N at given HNL mass, U²=1."""
    # Sign convention: parent+ → ℓ+ N → lepton is anti-lepton
    sign = "-" if parent_pdg > 0 else ""
    pid_lep = f"{sign}{abs(lepton_pdg)}"
    try:
        br_expr = hnl.get_2body_br(str(parent_pdg), pid_lep)
        # Restricted namespace covers all plausible variable names used by HNLCalc
        ns = {
            "mass": m_N, "mN": m_N, "mHNL": m_N, "m_N": m_N,
            "coupling": 1.0, "Usq": 1.0, "U2": 1.0,
            "np": np, "math": __import__("math"),
        }
        br_val = eval(br_expr, {"__builtins__": {}}, ns)  # noqa: S307
        if not np.isfinite(br_val) or br_val < 0:
            return 0.0
        return float(br_val)
    except Exception as e:
        import warnings
        warnings.warn(f"BR eval failed for PDG {parent_pdg}→{pid_lep}: {e}")
        return 0.0


def _eval_3body_br(hnl, parent_pdg, daughter_pdg, lepton_pdg, m_N, ch_type):
    """Evaluate 3-body BR by numerical integration at U²=1."""
    m_parent = MESON_MASSES.get(abs(parent_pdg), hnl.masses(parent_pdg))
    m_daughter = hnl.masses(daughter_pdg)
    m_lepton = hnl.masses(lepton_pdg)

    if m_N >= m_parent - m_daughter - m_lepton:
        return 0.0

    # Get differential BR expression
    sign_lep = "-" if parent_pdg > 0 else ""
    pid_lep_str = f"{sign_lep}{abs(lepton_pdg)}"

    try:
        if ch_type == "pseudo":
            dbr = hnl.get_3body_dbr_pseudoscalar(str(parent_pdg), str(daughter_pdg), pid_lep_str)
        elif ch_type == "vector":
            dbr = hnl.get_3body_dbr_vector(str(parent_pdg), str(daughter_pdg), pid_lep_str)
        else:
            return 0.0

        br_val = hnl.integrate_3body_br(
            dbr, m_N, m_parent, m_daughter, m_lepton,
            coupling=1.0, nsample=500
        )
        if br_val is None or np.isnan(br_val) or br_val < 0:
            return 0.0
        return float(br_val)
    except Exception:
        return 0.0


def compute_total_production_br(hnl, parent_pdg, lepton_pdg, m_N):
    """
    Compute total BR(parent → N + anything) summing 2-body and 3-body channels.

    Returns
    -------
    float
        Total branching ratio at U²=1.
    """
    _, _, br_total = compute_production_br_components(hnl, parent_pdg, lepton_pdg, m_N)
    return br_total


def compute_production_br_components(hnl, parent_pdg, lepton_pdg, m_N):
    """Return BR components: 2-body BR, 3-body channel list, and total BR."""
    m_parent = MESON_MASSES.get(abs(parent_pdg), hnl.masses(parent_pdg))
    m_lepton = hnl.masses(lepton_pdg)

    br_2body = 0.0
    if m_N < m_parent - m_lepton:
        br_2body = _eval_2body_br(hnl, parent_pdg, lepton_pdg, m_N)

    br_3body_channels = []
    abs_pdg = abs(parent_pdg)
    if abs_pdg in THREEBODY_CHANNELS:
        for daughter_pdg, ch_type in THREEBODY_CHANNELS[abs_pdg]:
            d_pdg = -daughter_pdg if parent_pdg < 0 else daughter_pdg
            br_3b = _eval_3body_br(hnl, parent_pdg, d_pdg, lepton_pdg, m_N, ch_type)
            if br_3b > 0:
                br_3body_channels.append((d_pdg, br_3b))

    br_3body_total = sum(br for _, br in br_3body_channels)
    br_total = br_2body + br_3body_total
    return br_2body, br_3body_channels, br_total


def _sample_hnl_from_mesons(parent_E, parent_px, parent_py, parent_pz, m_parent, m_lepton, m_N,
                            br_3body_channels, br_total, hnl, rng):
    """Sample HNL 4-vectors from 2-body/3-body decays with BR-weighted channel selection."""
    n_events = len(parent_E)
    hnl_4v = np.empty((n_events, 4))

    br_3body_total = sum(br for _, br in br_3body_channels)
    use_3body = np.zeros(n_events, dtype=bool)
    if br_3body_total > 0:
        p_3body = br_3body_total / br_total
        use_3body = rng.random(n_events) < p_3body

    use_2body = ~use_3body
    if use_2body.any():
        _, hnl_2b = decay_2body(
            parent_E[use_2body], parent_px[use_2body], parent_py[use_2body], parent_pz[use_2body],
            m_parent, m_lepton, m_N, rng=rng,
        )
        hnl_4v[use_2body] = hnl_2b

    if use_3body.any():
        br_arr = np.array([br for _, br in br_3body_channels], dtype=float)
        prob_arr = br_arr / br_arr.sum()
        channel_idx = rng.choice(len(br_3body_channels), size=use_3body.sum(), p=prob_arr)
        evt_idx = np.where(use_3body)[0]
        for idx_ch in np.unique(channel_idx):
            daughter_pdg, _ = br_3body_channels[int(idx_ch)]
            sel = evt_idx[channel_idx == idx_ch]
            m_daughter = hnl.masses(daughter_pdg)
            _, _, hnl_3b = decay_3body_flat(
                parent_E[sel], parent_px[sel], parent_py[sel], parent_pz[sel],
                m_parent, m_daughter, m_lepton, m_N, rng=rng,
            )
            hnl_4v[sel] = hnl_3b

    return hnl_4v


def generate_pool(quark, n_pool, rng):
    """Generate a reusable meson 4-vector pool for a given quark type."""
    print(f"  Generating {quark} meson pool ({n_pool} events)...")
    pool = sample_meson_4vectors(n_pool, quark, rng=rng)
    sigma = get_sigma_total(quark)
    print(f"  σ_FONLL({quark}) = {sigma:.3e} pb")
    return pool, sigma


def process_channel(flavor, quark, pool, sigma_fonll, masses, rng):
    """
    Process all mass points for a (flavor, quark) channel.

    For B/D mesons, uses fragmentation fractions from the pool species assignment.
    """
    hnl = _init_hnlcalc(flavor)
    lepton_pdg = FLAVOR_TO_LEPTON_PDG[flavor]
    m_lepton = LEPTON_MASSES[flavor]
    channel_label = CHANNEL_LABELS[quark]
    out_dir = OUTPUT_BASE / flavor / channel_label
    out_dir.mkdir(parents=True, exist_ok=True)

    if quark in ("bottom", "charm"):
        frag_map = FRAG_B if quark == "bottom" else FRAG_C
    else:
        frag_map = None  # Bc handled separately

    n_pool = len(pool['E'])

    for m_N in masses:
        mass_label = format_mass_for_filename(m_N)
        csv_path = out_dir / f"mN_{mass_label}.csv"

        all_weights = []
        all_E = []
        all_px = []
        all_py = []
        all_pz = []

        if quark == "bc":
            # Bc enriched channel: all events are Bc
            parent_pdg = 541
            m_parent = MESON_MASSES[541]
            if m_N >= m_parent - m_lepton:
                _write_empty_csv(csv_path)
                continue

            _, br_3body_channels, br = compute_production_br_components(
                hnl, parent_pdg, lepton_pdg, m_N
            )
            if br <= 0:
                _write_empty_csv(csv_path)
                continue

            # Weight: σ(pp → Bc + Bc⁻) × BR / N_pool  (no factor 2; SIGMA_BC_PB is total)
            w = SIGMA_BC_PB * br / n_pool

            # Decay each Bc → HNL with BR-weighted 2-body/3-body kinematics
            mask = np.ones(n_pool, dtype=bool)
            hnl_4v = _sample_hnl_from_mesons(
                pool['E'][mask], pool['px'][mask], pool['py'][mask], pool['pz'][mask],
                m_parent, m_lepton, m_N,
                br_3body_channels, br, hnl, rng,
            )
            all_weights.append(np.full(mask.sum(), w))
            all_E.append(hnl_4v[:, 0])
            all_px.append(hnl_4v[:, 1])
            all_py.append(hnl_4v[:, 2])
            all_pz.append(hnl_4v[:, 3])

        else:
            # B or D mesons: loop over species in the pool
            unique_species = np.unique(pool['species_pdg'])
            for species_pdg in unique_species:
                m_parent = MESON_MASSES[species_pdg]
                if m_N >= m_parent - m_lepton:
                    continue

                _, br_3body_channels, br = compute_production_br_components(
                    hnl, int(species_pdg), lepton_pdg, m_N
                )
                if br <= 0:
                    continue

                frag = frag_map.get(species_pdg, 0.0)
                if frag <= 0:
                    continue

                # Select events of this species
                mask = pool['species_pdg'] == species_pdg
                n_species = mask.sum()
                if n_species == 0:
                    continue

                # Weight: 2 × σ_FONLL × f_species × BR / N_species
                w = 2.0 * sigma_fonll * frag * br / n_species

                # Decay meson → HNL using BR-weighted 2-body/3-body kinematics
                hnl_4v = _sample_hnl_from_mesons(
                    pool['E'][mask], pool['px'][mask], pool['py'][mask], pool['pz'][mask],
                    m_parent, m_lepton, m_N,
                    br_3body_channels, br, hnl, rng,
                )

                all_weights.append(np.full(n_species, w))
                all_E.append(hnl_4v[:, 0])
                all_px.append(hnl_4v[:, 1])
                all_py.append(hnl_4v[:, 2])
                all_pz.append(hnl_4v[:, 3])

        if all_weights:
            weights = np.concatenate(all_weights)
            E = np.concatenate(all_E)
            px = np.concatenate(all_px)
            py = np.concatenate(all_py)
            pz = np.concatenate(all_pz)
            _write_csv(csv_path, weights, E, px, py, pz)
            print(f"    {csv_path.name}: {len(weights)} events, w_sum={weights.sum():.3e}")
        else:
            _write_empty_csv(csv_path)
            print(f"    {csv_path.name}: 0 events (below threshold)")


def _write_csv(path, weights, E, px, py, pz):
    """Write headerless CSV: weight,E,px,py,pz"""
    data = np.column_stack([weights, E, px, py, pz])
    np.savetxt(path, data, delimiter=",", fmt="%.8e")


def _write_empty_csv(path):
    """Write an empty file for mass points below threshold."""
    path.write_text("")


def generate_bc_pool(n_pool, rng):
    """
    Generate a Bc-enriched pool: FONLL bottom meson kinematics forced to Bc mass.

    We sample from the bottom FONLL grid but override the species to Bc
    and recompute 4-vectors with Bc mass.
    """
    from production.fonll.fonll_parser import parse_fonll_file, FONLL_FILES
    from production.fonll.meson_sampler import _build_cdf

    path = FONLL_FILES["bottom"]
    pt_arr, y_arr, dsigma_2d = parse_fonll_file(path)
    cdf, _, _, pt_widths, y_widths = _build_cdf(pt_arr, y_arr, dsigma_2d)

    n_pt = len(pt_arr)
    n_y = len(y_arr)

    u = rng.random(n_pool)
    flat_idx = np.searchsorted(cdf, u)
    flat_idx = np.clip(flat_idx, 0, n_pt * n_y - 1)

    i_pt = flat_idx // n_y
    i_y = flat_idx % n_y

    pt = pt_arr[i_pt] + rng.uniform(-0.5, 0.5, n_pool) * pt_widths[i_pt]
    y = y_arr[i_y] + rng.uniform(-0.5, 0.5, n_pool) * y_widths[i_y]
    pt = np.maximum(pt, 0.0)
    phi = rng.uniform(0, 2 * np.pi, n_pool)

    m = MESON_MASSES[541]  # Bc mass
    mt = np.sqrt(pt**2 + m**2)
    pz = mt * np.sinh(y)
    E = mt * np.cosh(y)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)

    return {
        'E': E, 'px': px, 'py': py, 'pz': pz,
        'species_pdg': np.full(n_pool, 541),
        'pt': pt, 'y': y, 'phi': phi,
    }


def main():
    """Main entry point: generate meson → HNL CSVs for all channels."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate meson → HNL 4-vector CSVs")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--channel", choices=["Bmeson", "Dmeson", "Bc", "all"],
                        default="all")
    parser.add_argument("--n-pool", type=int, default=N_POOL)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--masses", type=float, nargs="+", default=None,
                        help="Custom mass list (default: full grid)")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    masses = args.masses if args.masses else MASS_GRID

    # Determine which channels to run
    if args.channel == "all":
        channels = ["bottom", "charm", "bc"]
    else:
        label_to_quark = {"Bmeson": "bottom", "Dmeson": "charm", "Bc": "bc"}
        channels = [label_to_quark[args.channel]]

    # Generate pools (reused across flavors)
    pools = {}
    sigmas = {}
    for ch in channels:
        if ch == "bc":
            pools["bc"] = generate_bc_pool(args.n_pool, rng)
            sigmas["bc"] = 0.0  # Not used; Bc uses SIGMA_BC_PB directly
        else:
            pools[ch], sigmas[ch] = generate_pool(ch, args.n_pool, rng)

    # Process each (flavor, channel) combination
    for flavor in args.flavor:
        print(f"\n{'='*60}")
        print(f"Flavor: {flavor}")
        print(f"{'='*60}")
        for ch in channels:
            print(f"\n  Channel: {CHANNEL_LABELS[ch]}")
            process_channel(flavor, ch, pools[ch], sigmas.get(ch, 0.0), masses, rng)

    print("\nDone.")


if __name__ == "__main__":
    main()

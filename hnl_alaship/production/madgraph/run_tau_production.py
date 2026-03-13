#!/usr/bin/env python3
"""
production/madgraph/run_tau_production.py

Generate tau pool via MadGraph (pp → W → τν, pp → τ+τ-), then analytically
decay τ → N X using HNLCalc BRs for each (flavor, mass) point.

Two stages:
  Stage 1: MG5 → tau LHE → tau 4-vector pool (flavor-independent, run once)
  Stage 2: For each (flavor, m_N < m_τ):
           - Compute BR(τ → N X) from HNLCalc at U²=1
           - Decay τ → N via 2-body kinematics
           - Weight: w_i = σ_tau(MG5) × K_FACTOR_EW × BR(τ→NX) / N_sample

Output: output/llp_4vectors/{Ue,Umu,Utau}/tau/mN_{mass}.csv
Format: headerless, 5 columns: weight,E,px,py,pz

Usage:
    python run_tau_production.py                    # full scan
    python run_tau_production.py --skip-mg5         # use cached tau pool
    python run_tau_production.py --flavor Utau
    python run_tau_production.py --test
"""

import sys
import shutil
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "vendored" / "HNLCalc"))

from config_mass_grid import MASS_GRID, format_mass_for_filename
from production.constants import (
    K_FACTOR_EW, M_TAU,
)
from production.decay_engine.kinematics import decay_2body
from production.madgraph.docker_mg5 import (
    MG5_IMAGE,
    docker_repo_path,
    run_generate_events,
    run_mg5_command_file,
)

CARDS_DIR = Path(__file__).parent / "cards"
WORK_DIR = Path(__file__).parent / "work"
TAU_POOL_CSV = Path(__file__).parent / "work" / "tau_pool.csv"
OUTPUT_BASE = PROJECT_ROOT / "output" / "llp_4vectors"


def _write_tau_proc_card(cmd_file, proc_lines, output_dir):
    """Write the MG5 command file for tau production."""
    with open(cmd_file, 'w') as f:
        # Use SM model for tau production (no BSM)
        f.write("import model sm\n\n")
        f.write("set automatic_html_opening False\n")
        for line in proc_lines:
            stripped = line.strip()
            if ('generate' in stripped or 'add process' in stripped) and not stripped.startswith('#'):
                f.write(line)
        f.write(f"\noutput {output_dir} -nojpeg\n")
        f.write("quit\n")


def _disable_mg5_html_opening(cards_dir):
    """Disable MG5 browser auto-open for this process directory."""
    cfg = cards_dir / "me5_configuration.txt"
    if not cfg.exists():
        return
    text = cfg.read_text()
    new_text = text.replace("# automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# automatic_html_opening = False", "automatic_html_opening = False")
    new_text = new_text.replace("automatic_html_opening = True", "automatic_html_opening = False")
    new_text = new_text.replace("# web_browser = None", "web_browser = None")
    if "automatic_html_opening" not in new_text:
        new_text += "\nautomatic_html_opening = False\n"
    if new_text != text:
        cfg.write_text(new_text)


def generate_tau_pool(n_events=100_000):
    """Stage 1: Generate tau pool via MadGraph."""
    print("Stage 1: Generating tau pool via MadGraph...")
    print(f"  Using Docker image: {MG5_IMAGE}")

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    work_subdir = WORK_DIR / "tau_pool"

    # Generate process
    proc_card = CARDS_DIR / "proc_card_tau_production.dat"
    cmd_file = WORK_DIR / "mg5_gen_tau.txt"

    with open(proc_card) as f:
        proc_lines = f.readlines()

    log_file = WORK_DIR / "mg5_gen_tau.log"
    _write_tau_proc_card(cmd_file, proc_lines, docker_repo_path(PROJECT_ROOT, work_subdir))
    result = run_mg5_command_file(PROJECT_ROOT, cmd_file, log_file, timeout=300)

    if result.returncode != 0:
        print(f"  FAILED: process generation (see {log_file})")
        return None

    # Write run card
    cards_dir = work_subdir / "Cards"
    run_content = (CARDS_DIR / "run_card_template.dat").read_text()
    run_content = run_content.replace("N_EVENTS_PLACEHOLDER", str(n_events))
    # Add pT cut for leptons (10 GeV)
    run_content = run_content.replace("0.0  = ptl", "10.0  = ptl")
    (cards_dir / "run_card.dat").write_text(run_content)
    _disable_mg5_html_opening(cards_dir)

    # Run event generation
    gen_log = work_subdir / "generate_events.log"
    result = run_generate_events(PROJECT_ROOT, work_subdir, gen_log, timeout=3600)

    if result.returncode != 0:
        print(f"  FAILED: event generation (see {gen_log})")
        return None

    # Find LHE and convert to tau CSV
    events_dir = work_subdir / "Events"
    lhe_path = None
    if events_dir.exists():
        for run_dir in sorted(events_dir.glob("run_*")):
            for lhe in [run_dir / "unweighted_events.lhe.gz",
                        run_dir / "unweighted_events.lhe"]:
                if lhe.exists():
                    lhe_path = lhe
                    break
            if lhe_path:
                break

    if lhe_path is None:
        print("  FAILED: no LHE file found")
        return None

    # Extract tau 4-vectors
    from production.madgraph.lhe_to_csv import LHEParser
    parser = LHEParser(lhe_path)
    n_tau = parser.write_tau_csv(TAU_POOL_CSV)
    print(f"  Extracted {n_tau} tau 4-vectors → {TAU_POOL_CSV}")

    # Cleanup MG5 work directory
    shutil.rmtree(work_subdir, ignore_errors=True)
    cmd_file.unlink(missing_ok=True)

    return TAU_POOL_CSV


def _get_tau_production_components(hnl, m_N):
    """
    Compute total BR(τ → N + anything) and hadronic 2-body channel components.

    Uses HNLCalc's get_2body_br_tau for τ → π/K/ρ/K* + N channels.
    """
    br_total = 0.0
    hadronic_2body = []

    # 2-body leptonic: τ → N (via mixing) — this is the dominant channel
    # Actually τ → ℓ ν N is 3-body; the 2-body channels are τ → meson + N
    # HNLCalc handles this via get_2body_br_tau(tau_pdg, meson_pdg)

    # τ → π N
    for meson_pdg in [211, 321, 213, 323]:
        try:
            br_expr = hnl.get_2body_br_tau(15, meson_pdg)
            mass = m_N
            coupling = 1.0
            br_val = eval(br_expr)
            if not (np.isnan(br_val) or br_val < 0):
                br_total += br_val
                hadronic_2body.append((hnl.masses(meson_pdg), float(br_val)))
        except Exception:
            pass

    # 3-body leptonic: τ → ℓ ν N  (via get_3body_dbr_tau)
    # Channel 1: τ⁻ → ℓ⁻ ν_τ N  (neutrino is ν_τ, pid=16)
    for lep_pid in [11, 13]:
        try:
            dbr = hnl.get_3body_dbr_tau(15, -lep_pid, 16)
            m_lep = hnl.masses(lep_pid)
            br_val = hnl.integrate_3body_br(
                dbr, m_N, M_TAU, m_lep, 0.0,
                coupling=1.0, nsample=500, integration="dE",
            )
            if br_val is not None and not np.isnan(br_val) and br_val > 0:
                br_total += br_val
        except Exception:
            pass

    # Channel 2: τ⁻ → ℓ⁻ ν̄_ℓ N  (neutrino is ν̄_ℓ)
    for lep_pid in [11, 13]:
        nu_pid = lep_pid + 1  # 12 for e, 14 for mu
        try:
            dbr = hnl.get_3body_dbr_tau(15, -lep_pid, nu_pid)
            m_lep = hnl.masses(lep_pid)
            br_val = hnl.integrate_3body_br(
                dbr, m_N, M_TAU, m_lep, 0.0,
                coupling=1.0, nsample=500, integration="dE",
            )
            if br_val is not None and not np.isnan(br_val) and br_val > 0:
                br_total += br_val
        except Exception:
            pass

    return br_total, hadronic_2body


def _get_tau_production_br(hnl, m_N):
    """Backward-compatible wrapper returning only total BR."""
    br_total, _ = _get_tau_production_components(hnl, m_N)
    return br_total


def process_tau_decay(tau_pool_path, flavor, masses, rng):
    """
    Stage 2: Decay tau pool → HNL for each (flavor, mass).

    Weight: w_i × K_FACTOR_EW × BR(τ→NX) / 1.0
    (w_i from MG5 already encodes σ/N_events)
    """
    from HNLCalc import HNLCalc

    # Initialize HNLCalc
    if flavor == "Ue":
        hnl = HNLCalc(ve=1, vmu=0, vtau=0)
    elif flavor == "Umu":
        hnl = HNLCalc(ve=0, vmu=1, vtau=0)
    else:
        hnl = HNLCalc(ve=0, vmu=0, vtau=1)

    # Load tau pool
    try:
        data = np.loadtxt(tau_pool_path, delimiter=",")
    except ValueError:
        data = np.empty((0, 5))
    if data.size == 0:
        for m_N in masses:
            mass_label = format_mass_for_filename(m_N)
            csv_path = OUTPUT_BASE / flavor / "tau" / f"mN_{mass_label}.csv"
            csv_path.parent.mkdir(parents=True, exist_ok=True)
            csv_path.write_text("")
        print(f"    No tau events in pool for {flavor}; wrote empty tau CSVs.")
        return
    if data.ndim == 1:
        data = data.reshape(1, -1)
    tau_w = data[:, 0]
    tau_E = data[:, 1]
    tau_px = data[:, 2]
    tau_py = data[:, 3]
    tau_pz = data[:, 4]
    n_tau = len(tau_E)

    out_dir = OUTPUT_BASE / flavor / "tau"
    out_dir.mkdir(parents=True, exist_ok=True)

    for m_N in masses:
        mass_label = format_mass_for_filename(m_N)
        csv_path = out_dir / f"mN_{mass_label}.csv"

        if m_N >= M_TAU:
            csv_path.write_text("")
            print(f"    {csv_path.name}: 0 events (m_N >= m_tau)")
            continue

        # Compute BR(τ → N X)
        br, hadronic_2body = _get_tau_production_components(hnl, m_N)
        if br <= 0:
            csv_path.write_text("")
            print(f"    {csv_path.name}: 0 events (BR=0)")
            continue

        hnl_4v = np.empty((n_tau, 4))
        br_had_2body = sum(br_ch for _, br_ch in hadronic_2body)
        use_had_2body = np.zeros(n_tau, dtype=bool)
        if br_had_2body > 0:
            use_had_2body = rng.random(n_tau) < (br_had_2body / br)

        if use_had_2body.any():
            br_arr = np.array([br_ch for _, br_ch in hadronic_2body], dtype=float)
            prob_arr = br_arr / br_arr.sum()
            ch_idx = rng.choice(len(hadronic_2body), size=use_had_2body.sum(), p=prob_arr)
            evt_idx = np.where(use_had_2body)[0]
            for idx_ch in np.unique(ch_idx):
                m_daughter, _ = hadronic_2body[int(idx_ch)]
                sel = evt_idx[ch_idx == idx_ch]
                _, hnl_2b = decay_2body(
                    tau_E[sel], tau_px[sel], tau_py[sel], tau_pz[sel],
                    M_TAU, m_daughter, m_N, rng=rng,
                )
                hnl_4v[sel] = hnl_2b

        use_eff = ~use_had_2body
        if use_eff.any():
            _, hnl_eff = decay_2body(
                tau_E[use_eff], tau_px[use_eff], tau_py[use_eff], tau_pz[use_eff],
                M_TAU, 0.0, m_N, rng=rng,
            )
            hnl_4v[use_eff] = hnl_eff

        # Weight: MG5_weight × K_FACTOR_EW × BR
        weights = tau_w * K_FACTOR_EW * br

        out_data = np.column_stack([weights, hnl_4v[:, 0], hnl_4v[:, 1],
                                     hnl_4v[:, 2], hnl_4v[:, 3]])
        np.savetxt(csv_path, out_data, delimiter=",", fmt="%.8e")
        print(f"    {csv_path.name}: {n_tau} events, BR={br:.3e}, w_sum={weights.sum():.3e}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Tau → N production pipeline")
    parser.add_argument("--flavor", choices=["Ue", "Umu", "Utau"], nargs="+",
                        default=["Ue", "Umu", "Utau"])
    parser.add_argument("--masses", type=float, nargs="+", default=None)
    parser.add_argument("--nevents", type=int, default=100_000,
                        help="Events for tau pool generation")
    parser.add_argument("--skip-mg5", action="store_true",
                        help="Skip MG5 step; use cached tau pool")
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    if args.test:
        flavors = ["Utau"]
        masses = [0.5, 1.0]
        n_events_mg5 = 1000
    else:
        flavors = args.flavor
        masses = args.masses if args.masses else [m for m in MASS_GRID if m < M_TAU]
        n_events_mg5 = args.nevents

    # Stage 1: Generate tau pool (or use cached)
    if args.skip_mg5 and TAU_POOL_CSV.exists():
        print(f"Using cached tau pool: {TAU_POOL_CSV}")
    else:
        if shutil.which("docker") is None:
            print(f"ERROR: docker not found in PATH; tau production requires the {MG5_IMAGE} image.")
            return 1
        result = generate_tau_pool(n_events_mg5)
        if result is None:
            return 1

    # Stage 2: Decay taus → HNL
    for flavor in flavors:
        print(f"\n  Flavor: {flavor}")
        process_tau_decay(TAU_POOL_CSV, flavor, masses, rng)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

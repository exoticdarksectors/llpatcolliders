"""BC4 signal acceptance: S decays -> best-two-track -> shared GRENDEL reco.

Mirrors ``hnl/analysis/decay_reco_acceptance.py`` (the PR #15 template) but
replaces the FairShip rest-frame decay templates with an analytic scalar decay
engine.  Everything geometric is the *shared single source*: the GRENDEL mesh
(``higgs/grendel_geometry``) and the bounded 4-hit two-track reconstruction +
PR #13 selection (``higgs/reco_common``, re-exported by the hnl template) are
imported -- never copied.

Per S four-vector that hits the fiducial air volume:
  1. sample a decay vertex along the flight path (uniform; the decay-density
     weight is folded in by the lifetime reweighting in ``scan_u2``);
  2. sample a decay MODE by its branching ratio (mu mu / pi pi / K K / s s /
     c c / g g / tau tau / e e), generate the two charged daughters as a
     two-body S -> x xbar in the rest frame and boost to the lab;
  3. neutral sub-modes (pi0 pi0, K0 K0bar) yield no charged tracks and fail, so
     the visible (>= 2 charged track) branching fraction is folded into the
     reconstruction outcome exactly as the HNL folds it via its templates -- no
     separate BR_vis multiplier;
  4. take the two daughter directions, ray-cast to the wall, build the 4 hits,
     reconstruct and apply ``selection_mask``.

The decay geometry/outcome is independent of ``sin^2 theta`` (every partial
width scales as ``sin^2 theta`` so the mode composition is coupling-independent),
so the MC is built ONCE per mass and ``scan_u2`` reweights it to every coupling
by the production (linear) and lifetime (exponential) factors.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

_SCALAR_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _SCALAR_ROOT.parent
_HNL_ROOT = _REPO_ROOT / "hnl"
_HIGGS_DIR = _REPO_ROOT / "higgs"
for _p in (str(_REPO_ROOT), str(_HNL_ROOT), str(_HIGGS_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Shared reconstruction + selection + the coupling scan, imported from the hnl
# template (which itself imports the higgs/reco_common single source).
from analysis.decay_reco_acceptance import (          # noqa: E402
    reconstruct_decays, selection_mask, scan_u2,
    signal_contribution_diagnostics,
    HIT_RESOLUTION, P_CUT, CMS_ORIGIN)
from analysis._engine import (                          # noqa: E402
    compute_geometry, _eta_phi_to_directions_batch)
from analysis.exclusion import find_exclusion_band_refined  # noqa: E402
from analysis.format_bridge import load_combined_csv     # noqa: E402
from analysis.constants import L_INT_PB, N_THRESHOLD     # noqa: E402
from production.decay_engine.kinematics import _boost_to_lab  # noqa: E402
from reco_common import SIGMA_T_DEFAULT                  # noqa: E402

from scalar import model                                 # noqa: E402

# Decay mode -> (daughter mass for the 2-body proxy, charged-track fraction).
# The charged fraction encodes the neutral isospin sub-modes that give no
# tracks (pi0 pi0 = 1/3 of pi pi; K0 K0bar = 1/2 of K K).  The quark/gluon
# continuum (s s, c c, g g) and 4 pi are modelled by their two LEADING charged
# hadrons (pion-mass proxy) -- a documented approximation; these channels are
# sub-dominant in the < 2 GeV region that drives the reach.
_MODE_DAUGHTER = {
    "ee":     (model.M_ELECTRON, 1.0),
    "mumu":   (model.M_MUON,     1.0),
    "tautau": (model.M_TAU,      1.0),     # tau directions proxy the leading tracks
    "pipi":   (model.M_PIPLUS,   2.0 / 3.0),
    "KK":     (model.M_KPLUS,    0.5),
    "4pi":    (model.M_PIPLUS,   1.0),
    "ss":     (model.M_PIPLUS,   1.0),
    "cc":     (model.M_PIPLUS,   1.0),
    "gg":     (model.M_PIPLUS,   1.0),
}

# Bound the daughter-to-wall ray-cast allocation. Decays are independent; this
# changes only RNG ordering relative to the old monolithic call, not the MC
# estimator. It prevents a large production-pool x multi-sample campaign from
# exhausting memory during a single reconstruction call.
EVENT_RECO_CHUNK = 25_000


def _mode_table(m_S, width_scheme="winkler"):
    """Return (modes, probs, daughter_mass, charged_frac) arrays for ``m_S``."""
    br = model.branching_ratios(m_S, scheme=width_scheme)
    modes = [k for k in br if br[k] > 0.0]
    probs = np.array([br[k] for k in modes], float)
    probs = probs / probs.sum()
    m_d = np.array([_MODE_DAUGHTER[k][0] for k in modes], float)
    chf = np.array([_MODE_DAUGHTER[k][1] for k in modes], float)
    return modes, probs, m_d, chf


def build_event_mc(p4, direction, entry_d, exit_d, m_S, n_samples, rng,
                   sigma_hit=HIT_RESOLUTION, sigma_t=SIGMA_T_DEFAULT,
                   origin=CMS_ORIGIN, width_scheme="winkler",
                   reco_chunk=EVENT_RECO_CHUNK):
    """Sample decay vertices + scalar decays and reconstruct, for a batch of S
    four-vectors.  Returns ``(d, passed)`` each ``(n_events, n_samples)``; both
    are ``sin^2 theta``-independent (``scan_u2`` applies the coupling)."""
    origin = np.asarray(origin, float)
    p4 = np.asarray(p4, float)
    direction = np.asarray(direction, float)
    n_ev = len(entry_d)

    d = rng.uniform(np.asarray(entry_d)[:, None], np.asarray(exit_d)[:, None],
                    size=(n_ev, n_samples))
    M = n_ev * n_samples
    vtx = origin[None, :] + d.reshape(M, 1) * np.repeat(direction, n_samples, axis=0)
    parent = np.repeat(p4, n_samples, axis=0)            # (M, 4) S four-vectors

    # Sample a decay mode per decay, then its charged/neutral sub-mode.
    modes, probs, m_d_tab, chf_tab = _mode_table(
        m_S, width_scheme=width_scheme)
    mi = rng.choice(len(modes), size=M, p=probs)
    m_d = m_d_tab[mi]
    charged = rng.random(M) < chf_tab[mi]
    pstar = np.sqrt(np.clip(0.25 * m_S ** 2 - m_d ** 2, 0.0, None))
    valid = charged & (pstar > 0)

    passed = np.zeros(M, dtype=bool)
    if valid.any():
        idx = np.where(valid)[0]
        # Isotropic back-to-back daughters in the S rest frame.
        cos_t = rng.uniform(-1.0, 1.0, len(idx))
        sin_t = np.sqrt(np.clip(1.0 - cos_t ** 2, 0.0, None))
        phi = rng.uniform(0.0, 2 * np.pi, len(idx))
        nhat = np.column_stack([sin_t * np.cos(phi), sin_t * np.sin(phi), cos_t])
        ps = pstar[idx][:, None]
        E_d = np.full(len(idx), m_S / 2.0)
        d1_rest = np.column_stack([E_d, ps * nhat])
        d2_rest = np.column_stack([E_d, -ps * nhat])
        pe, ppx, ppy, ppz = (parent[idx, 0], parent[idx, 1],
                             parent[idx, 2], parent[idx, 3])
        d1 = _boost_to_lab(d1_rest, pe, ppx, ppy, ppz)
        d2 = _boost_to_lab(d2_rest, pe, ppx, ppy, ppz)
        p1 = np.linalg.norm(d1[:, 1:], axis=1)
        p2 = np.linalg.norm(d2[:, 1:], axis=1)
        dir1 = d1[:, 1:] / p1[:, None]
        dir2 = d2[:, 1:] / p2[:, None]
        p_soft = np.minimum(p1, p2)
        # Daughter speeds beta = |p|/E feed the timing model (soft pions/kaons
        # are genuinely slow; beta = 1 would overestimate the timing-cut
        # acceptance -- same fix as the HNL chain). Guard against E < |p|
        # roundoff after the boost.
        b1 = p1 / np.maximum(d1[:, 0], p1)
        b2 = p2 / np.maximum(d2[:, 0], p2)
        # Both daughters must clear the track momentum floor (mirrors best-two).
        ok = (p_soft > P_CUT)
        if ok.any():
            sub = idx[ok]
            ok_idx = np.where(ok)[0]
            for start in range(0, len(sub), reco_chunk):
                stop = min(start + reco_chunk, len(sub))
                take = ok_idx[start:stop]
                mc = reconstruct_decays(
                    vtx[sub[start:stop]], dir1[take], dir2[take], p_soft[take],
                    sigma_hit, sigma_t, rng, beta1=b1[take], beta2=b2[take])
                passed[sub[start:stop]] = selection_mask(mc)

    return d, passed.reshape(n_ev, n_samples)


def _geometry(csv_path, m_S, eta, phi, mesh):
    """Ray-cast with an on-disk cache per mass (mirrors the BC10/HNL pattern):
    the npz sits next to the four-vector CSVs and is invalidated whenever the
    CSV is newer (a regenerated pool must force a fresh ray-cast)."""
    cache_dir = csv_path.parent.parent / "geometry_cache"
    cache = cache_dir / f"geom_{csv_path.stem}.npz"
    if cache.exists() and cache.stat().st_mtime >= csv_path.stat().st_mtime:
        d = np.load(cache)
        return d["hits"].astype(bool), d["entry_d"], d["exit_d"]
    hits, entry_d, exit_d = compute_geometry(
        eta, phi, mesh, CMS_ORIGIN, batch_label=f"[mS={m_S:.3f}]")
    cache_dir.mkdir(parents=True, exist_ok=True)
    tmp = cache.with_suffix(cache.suffix + ".tmp")
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, hits=hits, entry_d=entry_d, exit_d=exit_d)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, cache)
    return hits, entry_d, exit_d


def process_mass_point(m_S, mesh, csv_path, sin2theta_grid, n_samples=100,
                       rng=None, width_scheme="winkler"):
    """Acceptance + coupling scan for one ``m_S``.  Returns its exclusion band
    dict (``find_exclusion_band``) augmented with mass / bookkeeping, or None if
    production is closed / no events hit."""
    csv_path = Path(csv_path)
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return None
    data = load_combined_csv(csv_path, m_S)
    if len(data["weight"]) == 0:
        return None
    if rng is None:
        rng = np.random.default_rng(12345)

    hits, entry_d, exit_d = _geometry(csv_path, m_S, data["eta"], data["phi"],
                                      mesh)
    idx = np.where(hits & np.isfinite(entry_d) & np.isfinite(exit_d))[0]
    n_hits = len(idx)
    base = {"mass_GeV": m_S, "n_events": int(len(data["weight"])), "n_hits": n_hits}
    if n_hits == 0:
        return {**base, "u2_min": np.nan, "u2_max": np.nan,
                "u2_min_open": False, "u2_max_open": False,
                "peak_N": 0.0, "peak_u2": np.nan, "has_sensitivity": False}

    direction = _eta_phi_to_directions_batch(data["eta"][idx], data["phi"][idx])
    p_mag = data["beta_gamma"][idx] * m_S
    energy = data["gamma"][idx] * m_S
    p4 = np.column_stack([energy, p_mag[:, None] * direction])

    d, passed = build_event_mc(p4, direction, entry_d[idx], exit_d[idx],
                               m_S, n_samples, rng,
                               width_scheme=width_scheme)
    ctau1 = model.ctau_sin2theta1(m_S, scheme=width_scheme)
    _, N_grid = scan_u2(
        d, passed, exit_d[idx] - entry_d[idx], data["weight"][idx],
        data["beta_gamma"][idx], ctau1, L_INT_PB, sin2theta_grid)

    def evaluate(sin2theta):
        _, signal = scan_u2(
            d, passed, exit_d[idx] - entry_d[idx], data["weight"][idx],
            data["beta_gamma"][idx], ctau1, L_INT_PB,
            np.asarray([sin2theta]),
        )
        return signal[0]

    result = find_exclusion_band_refined(
        sin2theta_grid, N_grid, evaluate, N_THRESHOLD,
    )
    result.update(base)
    result["ctau_sin2th1_m"] = ctau1
    for label in ("u2_min", "peak_u2", "u2_max"):
        coupling = result[label]
        if not np.isfinite(coupling):
            for field in ("sample_ess", "event_ess", "max_event_fraction"):
                result[f"{label}_{field}"] = np.nan
            continue
        diagnostics = signal_contribution_diagnostics(
            d,
            passed,
            exit_d[idx] - entry_d[idx],
            data["weight"][idx],
            data["beta_gamma"][idx],
            ctau1,
            coupling,
        )
        for field in ("sample_ess", "event_ess", "max_event_fraction"):
            result[f"{label}_{field}"] = diagnostics[field]
    return result, sin2theta_grid, N_grid

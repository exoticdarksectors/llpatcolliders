"""E14-16 -- the GRENDEL HNL money plot + A2 metadata + machine-readable bundle.

Built ON `analysis.plot_exclusion` (the proper renderer: excluded-region fill,
open-edge markers, ylim to the scan ceiling 1e-1, segment-aware drawing, and the
FONLL theory ribbon via `band_csv`). On top of that base it layers the additional
in-scope boundary bands:
  * lower edge: FONLL production band (orange, from `combine_band.py`),
                Bc normalization band (teal, `bc_nuisance.py`), and charged-kaon
                transport band (purple, `kaon_desc_band.csv`);
  * upper edge: decay-model width band (blue, `decay_model_band.py`).
Plus the F17 provenance line, the Umu low-mass `N->l pi` callout, the A2
hypothesis metadata (`run_metadata.json`), and a bundle dir with the tables.

Production + HNL->SM decay physics only; reconstruction / detector / background /
statistics are the partner layer, idealized (background-free, N>=3).

  python -m analysis.plot_money --run central_newgrids_20260623 --out-dir tmp/runs/v1_bundle
"""
from __future__ import annotations

import argparse
import json
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis.plot_exclusion import _plot_single_panel  # noqa: E402

HNL_ROOT = Path(__file__).resolve().parent.parent
RUNS = HNL_ROOT / "tmp" / "runs"
PUBLISHED = HNL_ROOT / "data" / "published"   # version-controlled, survives a clean clone
BUNDLE = PUBLISHED / "bundle"                  # tracked money-plot input CSVs (see its README)
LAB = {"Ue": r"$|U_e|^2$", "Umu": r"$|U_\mu|^2$", "Utau": r"$|U_\tau|^2$"}


def _provenance(have_fonll, have_kaon):
    fonll = ("FONLL production (orange)" if have_fonll
             else "FONLL production (variation PENDING -- not yet on this figure)")
    kaon = ("charged-kaon transport d_esc=[1,3] m (purple)" if have_kaon
            else "charged-kaon transport d_esc variation PENDING")
    return (
        f"Theory/model variations: {fonll}; {kaon}; Bc direct-norm ±40% "
        "(teal, LHCb arXiv:1910.13404; induced-tau Bc 0.06% negligible); "
        "HNL total-width/lifetime duality δ(m) (blue, cap 20% conservative vs ~10% post-QCD residual).  "
        "NOT banded (named limitations): production form factors, absolute visible-BR norm, neutral kaons, "
        "detailed material/magnetic transport beyond the d_esc proxy, FONLL αs.  "
        "Reconstruction / detector / background / statistics IDEALIZED at the partner handoff (background-free, N≥3).  "
        "Island closes on real refined grid points where peak_N crosses 3 (config_mass_grid 3.62-3.70): "
        "~3.62 GeV (Ue/Umu), ~3.69 GeV (Utau) -- no extrapolation/synthetic pinch.  "
        "Nuisance-induced topology changes are recorded in the bundle but not drawn as ordinary ribbons."
    )


def _dex_densify(bm, b_cen, b_lo, b_hi, cm, c_edge):
    """Densify one band edge onto the fine central grid.

    The ~18 band anchors are too coarse to plot as absolute edges -- straight
    log-segments between anchors detach from the wiggly dense central line. So we
    interpolate the band *width in dex* (relative to the band-run central) and
    re-apply it to the dense plotted central ``c_edge``, so the ribbon hugs every
    wiggle of the red curve and ends where the central does.
    """
    ok = (np.isfinite(b_cen) & (b_cen > 0) & np.isfinite(b_lo) & (b_lo > 0)
          & np.isfinite(b_hi) & (b_hi > 0))
    lo = np.full(len(cm), np.nan)
    hi = np.full(len(cm), np.nan)
    # Interpolate each contiguous finite anchor segment independently.  In
    # particular, do not extrapolate the last ordinary ribbon through a point
    # where a nuisance variation creates or destroys the exclusion island.
    valid = np.flatnonzero(ok)
    for run in np.split(valid, np.flatnonzero(np.diff(valid) > 1) + 1):
        if len(run) < 2:
            continue
        support = (cm >= bm[run[0]]) & (cm <= bm[run[-1]])
        dex_lo = np.log10(b_cen[run] / b_lo[run])
        dex_hi = np.log10(b_hi[run] / b_cen[run])
        lo[support] = c_edge[support] * 10.0 ** (
            -np.interp(cm[support], bm[run], dex_lo))
        hi[support] = c_edge[support] * 10.0 ** (
            +np.interp(cm[support], bm[run], dex_hi))
    return cm, lo, hi


def _ribbon(ax, m, lo, hi, color, label, zorder):
    good = np.isfinite(lo) & np.isfinite(hi) & (lo > 0) & (hi > 0)
    if good.any():
        ax.fill_between(m, np.where(good, lo, np.nan), np.where(good, hi, np.nan),
                        alpha=0.5, color=color, linewidth=0, zorder=zorder, label=label)




def _metadata(run, l_int_fb, p_cut_mev, have_fonll, have_kaon):
    in_band = ["HNL total-width/lifetime duality (decay-model band)",
               "direct Bc normalization (B4)"]
    if have_fonll:
        in_band.insert(0, "FONLL heavy-flavor production (scale/PDF/m_Q)")
    if have_kaon:
        in_band.append("charged-kaon transport escape length d_esc in [1,3] m")
    lim = ["production form factors (B3)", "absolute visible-BR normalization (B5)",
           "neutral-kaon contribution (B7)",
           "detailed material/magnetic transport beyond the d_esc proxy (B7)",
           "FONLL alpha_s (sub-dominant)"]
    if not have_fonll:
        lim.insert(0, "FONLL scale/PDF/m_Q band -- PENDING, not on this figure")
    if not have_kaon:
        lim.insert(0, "charged-kaon d_esc transport band -- PENDING, not on this figure")
    band_sources = {
        "fonll_lower_edge": ("hnl_band_fonll.csv (run_variation_band + combine_band)"
                             if have_fonll else "PENDING"),
        "decay_model_upper_edge": "decay_model_band.csv (width_band delta(m))",
        "bc_lower_edge": "bc_nuisance_band.csv (SIGMA_BC_REL_UNCERT=0.40)",
    }
    if have_kaon:
        band_sources["kaon_lower_edge"] = (
            "kaon_desc_band.csv (d_esc=1/1.5/3 m; nominal 1.5 m)")
    return {
        "result": "GRENDEL HNL theory/model variation diagnostic (single-flavor)",
        "hypothesis": {"mixing": "single-flavor", "flavors": ["Ue", "Umu", "Utau"],
                       "observable": "|U_alpha|^2", "nature": "Majorana",
                       "charge_conjugate_counting": "included (NDecayWidth x2/channel; production both charges)",
                       "production_mixing_equals_decay_mixing": True},
        "beam": {"collision": "pp", "sqrt_s_TeV": 14, "L_int_fb": l_int_fb},
        "limit": {"method": "background-free", "criterion": "N_signal >= 3",
                  "track_momentum_cut_MeV": p_cut_mev},
        "scope": {"in_band": in_band,
                  "interpretation": "theory/model variation diagnostic, not a statistical confidence band",
                  "idealized_partner_handoff": ["reconstruction", "detector response",
                                                "background", "statistics"]},
        "central_run": run,
        "band_sources": band_sources,
        "limitations_not_banded": lim,
        "known_features": ["Umu/Utau low-mass step where N->l pi closes (m_mu+m_pi=0.245 GeV) -- "
                           "real physics (REMAINING_WORK item 15)",
                           "island closes on real refined grid points (config_mass_grid 3.62-3.70) where "
                           "peak_N crosses 3: ~3.62 GeV (Ue/Umu), ~3.69 GeV (Utau) -- no extrapolation"],
        "generated": str(date.today()),
        "reproduce": {"env": f"HNL_TMP_DIR=tmp HNL_RUN_TAG=<run> HNL_P_CUT={p_cut_mev/1000:.3f} (hnl conda env)",
                      "steps": ["python -m analysis.run_sensitivity",
                                "run_variation_band.py + python -m analysis.combine_band  # FONLL band -> hnl_band.csv",
                                "python -m analysis.decay_model_band",
                                "python -m analysis.bc_nuisance",
                                "python -m analysis.channel_breakdown",
                                "python -m analysis.plot_money"],
                      "width_band_delta_table": "python -m analysis.width_band  (ROOT-enabled env)"},
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    # Defaults read the tracked published bundle (P>100 MeV cut), so a fresh clone
    # reproduces the figure with no tmp/ run tree. P>600 is legacy (pass
    # --central-csv .../analysis_exact_600/... --p-cut-mev 600 --decay-band ...
    # --bc-band ... to rebuild it).
    ap.add_argument("--run", default="central_newgrids_20260623")
    ap.add_argument("--central-csv", default=None)
    ap.add_argument("--fonll-band", default=str(BUNDLE / "hnl_band_fonll.csv"))
    ap.add_argument("--decay-band", default=str(BUNDLE / "decay_model_band.csv"))
    ap.add_argument("--bc-band", default=str(BUNDLE / "bc_nuisance_band.csv"))
    ap.add_argument("--kaon-band", default=str(BUNDLE / "kaon_desc_band.csv"))
    ap.add_argument("--breakdown", default=str(BUNDLE / "channel_breakdown_u2min.csv"))
    ap.add_argument("--l-int-fb", type=float, default=3000.0)
    ap.add_argument("--p-cut-mev", type=int, default=100)
    ap.add_argument("--out-dir", default=str(RUNS / "v1_bundle"))
    a = ap.parse_args(argv)

    central_csv = a.central_csv or str(PUBLISHED / "grendel_hnl_sensitivity.csv")
    cen = pd.read_csv(central_csv)
    have_fonll = Path(a.fonll_band).exists()
    have_kaon = Path(a.kaon_band).exists()
    fonll = pd.read_csv(a.fonll_band) if have_fonll else None
    dm = pd.read_csv(a.decay_band) if Path(a.decay_band).exists() else None
    bc = pd.read_csv(a.bc_band) if Path(a.bc_band).exists() else None
    kb = pd.read_csv(a.kaon_band) if have_kaon else None

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    for ax, fl in zip(axes, ["Ue", "Umu", "Utau"]):
        full = cen[cen.flavor == fl].sort_values("mass_GeV")
        sens = full[full.has_sensitivity == True]                            # noqa: E712
        cm = sens["mass_GeV"].to_numpy(float)
        c_min = sens["u2_min"].to_numpy(float)
        c_max = sens["u2_max"].to_numpy(float)
        # The island closes on REAL refined grid points (config_mass_grid
        # 3.62-3.70, where peak_N crosses 3): no synthetic pinch/extrapolation --
        # the contour and the bands both simply end on the last sensitive point.
        _plot_single_panel(ax, cen, fl, is_leftmost=(fl == "Ue"), band_df=None)

        # FONLL production band on both edges (hugs the dense central, follows nose).
        if have_fonll:
            bsub = fonll[fonll.flavor == fl].sort_values("mass_GeV")
            bm = bsub["mass_GeV"].to_numpy(float)
            labelled = False
            for col, c_edge in (("u2_min", c_min), ("u2_max", c_max)):
                m, lo, hi = _dex_densify(
                    bm, bsub[f"{col}_central"].to_numpy(float),
                    bsub[f"{col}_band_lo"].to_numpy(float),
                    bsub[f"{col}_band_hi"].to_numpy(float),
                    cm, c_edge)
                _ribbon(ax, m, lo, hi, "orange",
                        None if labelled else "FONLL production variation", zorder=4)
                labelled = True

        # Decay-model width band (upper edge) and Bc normalization band (lower edge).
        if dm is not None:
            dmf = dm[dm.flavor == fl].sort_values("mass_GeV")
            m, lo, hi = _dex_densify(
                dmf["mass_GeV"].to_numpy(float), dmf["u2_max"].to_numpy(float),
                dmf["u2_max_dm_lo"].to_numpy(float), dmf["u2_max_dm_hi"].to_numpy(float),
                cm, c_max)
            _ribbon(ax, m, lo, hi, "steelblue", "decay-model variation (upper)", zorder=4)
        if bc is not None:
            bcf = bc[bc.flavor == fl].sort_values("mass_GeV")
            m, lo, hi = _dex_densify(
                bcf["mass_GeV"].to_numpy(float), bcf["u2_min"].to_numpy(float),
                # Bc-up strengthens the limit (lower boundary); Bc-down weakens
                # it (upper boundary).  The historical column suffixes name the
                # nuisance direction, not the plotted vertical ordering.
                bcf["u2_min_bc_hi"].to_numpy(float), bcf["u2_min_bc_lo"].to_numpy(float),
                cm, c_min)
            _ribbon(ax, m, lo, hi, "teal", "Bc normalization variation (lower)", zorder=4)

        # Charged-kaon transport d_esc in [1, 3] m (nominal 1.5 m): the dominant
        # low-mass systematic (Ue/Umu below the K->l N threshold). Lower edge only.
        if kb is not None:
            kbf = kb[kb.flavor == fl].sort_values("mass_GeV")
            if len(kbf):
                m, lo, hi = _dex_densify(
                    kbf["mass_GeV"].to_numpy(float),
                    kbf["u2_min_desc1p5"].to_numpy(float),   # nominal d_esc = 1.5 m
                    kbf["u2_min_desc3p0"].to_numpy(float),   # d_esc = 3 m strengthens (lower)
                    kbf["u2_min_desc1p0"].to_numpy(float),   # d_esc = 1 m weakens (upper)
                    cm, c_min)
                _ribbon(ax, m, lo, hi, "purple",
                        "kaon transport d_esc [1,3] m (lower)", zorder=5)

        ax.set_title(f"HNL {LAB[fl]}", fontsize=13)
        if fl in ("Umu", "Utau"):
            ax.annotate(r"$N\to\ell\pi$ closes" + "\n→ weaker limit", xy=(0.25, 3e-6),
                        xytext=(0.33, 5e-5), fontsize=7, color="dimgray",
                        arrowprops=dict(arrowstyle="->", color="dimgray", lw=0.8))
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("GRENDEL HNL theory/model variation diagnostics (single-flavor, Majorana, 14 TeV, "
                 f"{a.l_int_fb:.0f} fb$^{{-1}}$, P>{a.p_cut_mev} MeV)", y=1.0, fontsize=13)
    fig.text(0.5, 0.005, _provenance(have_fonll, have_kaon), ha="center", va="bottom",
             fontsize=6.2, style="italic", wrap=True)
    plt.tight_layout(rect=[0, 0.075, 1, 0.97])

    out = Path(a.out_dir); out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "money_plot.png", dpi=130, bbox_inches="tight")
    fig.savefig(out / "money_plot.pdf", bbox_inches="tight")
    (out / "run_metadata.json").write_text(json.dumps(
        _metadata(a.run, a.l_int_fb, a.p_cut_mev, have_fonll, have_kaon), indent=2))
    for src, name in [(central_csv, "hnl_sensitivity_central.csv"),
                      (a.fonll_band, "hnl_band_fonll.csv"),
                      (a.decay_band, "decay_model_band.csv"),
                      (a.bc_band, "bc_nuisance_band.csv"),
                      (a.kaon_band, "kaon_desc_band.csv"),
                      (a.breakdown, "channel_breakdown_u2min.csv")]:
        if Path(src).exists():
            shutil.copy2(src, out / name)
    print(f"money plot + bundle -> {out}/  "
          f"(FONLL band: {'YES' if have_fonll else 'PENDING'}; "
          f"kaon band: {'YES' if have_kaon else 'PENDING'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

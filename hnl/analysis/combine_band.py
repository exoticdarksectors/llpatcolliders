#!/usr/bin/env python3
"""Combine per-variation HNL sensitivity curves into a theory-uncertainty band.

Reads the ``band_registry.json`` written by ``run_variation_band.py`` (one
``hnl_sensitivity.csv`` per coherent FONLL variation) and turns the spread of
each exclusion boundary into an uncertainty ribbon, per the agreed prescription:

* all combination is done in ``x = log10(U^2)`` (the boundary spans many
  decades, so min/max or std in linear U^2 would be meaningless);
* **scale** -> asymmetric envelope of the 7 coherent scale curves
  (max/min deviation of ``x`` from central), per boundary;
* **PDF** -> ``std(x, ddof=1)`` over the replica members (the NNPDF
  Monte-Carlo 1 sigma), per boundary;
* **mass** -> bottom and charm are independent: the per-quark maximum
  ``|x - x_central|`` are added in quadrature;
* **alpha_s** (optional, PDF4LHC) -> half the ``|x|`` spread between the
  as=0.119 and as=0.117 companion curves, folded symmetrically. The two
  companions are separate coherent grids (not registry replicas), passed via
  ``--alphas-lo``/``--alphas-hi``; omitting them reproduces the no-alpha_s band;
* the axes are combined in quadrature (uncorrelated) into an up/down
  ``sigma_x`` per boundary; the ribbon is ``10**(x_central +/- sigma_x)``.

Open-edge (``u2_*_open``) and no-sensitivity states are preserved: a boundary
is reported open if central or any contributing variation ran off the scan
edge, and masses where central has no sensitivity carry no band.  If any
contributing member creates or destroys an island/boundary relative to central,
the point is marked ``*_topology_changed`` and no numeric ribbon is reported;
topology changes are not silently reduced to a band from the surviving members.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HNL_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_REGISTRY = HNL_ROOT / "tmp" / "runs" / "band_registry.json"

BOUNDARIES = [("u2_min", "u2_min_open"), ("u2_max", "u2_max_open")]


def _load_curves(registry: dict) -> dict[str, pd.DataFrame]:
    curves = {}
    for v in registry["variations"]:
        df = pd.read_csv(v["sensitivity_csv"])
        df = df.set_index(["flavor", "mass_GeV"]).sort_index()
        curves[v["name"]] = df
    return curves


def _axis_members(registry: dict) -> dict[str, list[str]]:
    axes: dict[str, list[str]] = {"central": [], "scale": [], "pdf": [], "mass": []}
    for v in registry["variations"]:
        axes[v["axis"]].append(v["name"])
    return axes


def _read_curve(path) -> pd.DataFrame:
    """Load a standalone sensitivity CSV (e.g. an alpha_s companion) indexed
    like the registry curves."""
    return pd.read_csv(path).set_index(["flavor", "mass_GeV"]).sort_index()


def _boundary(curve: pd.DataFrame, key, mass_key, col, open_col):
    """Return (x=log10(boundary), is_open, has_sens) for one variation/point."""
    if key not in curve.index:
        return None, False, False
    row = curve.loc[key]
    if not bool(row.get("has_sensitivity", False)):
        return None, False, False
    is_open = bool(row.get(open_col, False))
    val = row.get(col, np.nan)
    if not np.isfinite(val) or val <= 0:
        return None, is_open, True
    return float(np.log10(val)), is_open, True


def combine_band(registry_path: Path, alphas_lo=None, alphas_hi=None) -> pd.DataFrame:
    registry = json.loads(registry_path.read_text())
    curves = _load_curves(registry)
    axes = _axis_members(registry)
    if not axes["central"]:
        raise ValueError("registry has no central variation")
    central_name = axes["central"][0]
    central = curves[central_name]
    registry_members = axes["scale"] + axes["pdf"] + axes["mass"]

    # alpha_s (PDF4LHC) is folded from two standalone companion grids, not
    # registry replicas; both must be present for the term to contribute.
    as_lo = _read_curve(alphas_lo) if alphas_lo else None
    as_hi = _read_curve(alphas_hi) if alphas_hi else None
    fold_alphas = as_lo is not None and as_hi is not None

    out_rows = []
    for (flavor, mass) in central.index:
        crow = central.loc[(flavor, mass)]
        rec = {"flavor": flavor, "mass_GeV": mass,
               "has_sensitivity": bool(crow.get("has_sensitivity", False))}
        if not rec["has_sensitivity"]:
            out_rows.append(rec)
            continue

        for col, open_col in BOUNDARIES:
            xc, c_open, _ = _boundary(central, (flavor, mass), mass, col, open_col)
            rec[f"{col}_central"] = crow.get(col, np.nan)
            rec[f"{col}_open"] = c_open
            member_states = [
                _boundary(curves[name], (flavor, mass), mass, col, open_col)
                for name in registry_members
            ]
            if fold_alphas:
                member_states.extend([
                    _boundary(as_lo, (flavor, mass), mass, col, open_col),
                    _boundary(as_hi, (flavor, mass), mass, col, open_col),
                ])
            rec[f"{col}_n_members_expected"] = len(member_states)
            rec[f"{col}_n_members_sensitive"] = sum(
                int(has_sens) for _, _, has_sens in member_states)
            rec[f"{col}_n_members_finite"] = sum(
                int(x is not None) for x, _, _ in member_states)
            rec[f"{col}_n_members_open"] = sum(
                int(is_open) for _, is_open, _ in member_states)
            topology_changed = any(
                (not has_sens)
                or (is_open != c_open)
                or ((x is not None) != (xc is not None))
                for x, is_open, has_sens in member_states
            )
            rec[f"{col}_topology_changed"] = topology_changed
            if xc is None:
                # central boundary open/undefined -> ribbon open on this side
                rec[f"{col}_band_lo"] = np.nan
                rec[f"{col}_band_hi"] = np.nan
                rec[f"{col}_open"] = True
                continue
            if topology_changed:
                rec[f"{col}_band_lo"] = np.nan
                rec[f"{col}_band_hi"] = np.nan
                rec[f"{col}_open"] = c_open or any(
                    is_open for _, is_open, _ in member_states)
                continue

            any_open = c_open

            # --- scale: asymmetric envelope of central + scale curves ---
            xs = [xc]
            for name in axes["scale"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xs.append(x)
            scale_up = max(xs) - xc
            scale_dn = xc - min(xs)

            # --- pdf: std over replica members in log10(U^2) ---
            xp = []
            for name in axes["pdf"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is not None:
                    xp.append(x)
            pdf_sigma = float(np.std(xp, ddof=1)) if len(xp) >= 2 else 0.0

            # --- mass: independent bottom/charm, per-quark max |dev|, in quadrature ---
            mdev = {"mb": 0.0, "mc": 0.0}
            for name in axes["mass"]:
                x, op, _ = _boundary(curves[name], (flavor, mass), mass, col, open_col)
                if op:
                    any_open = True
                if x is None:
                    continue
                quark = "mb" if name.startswith("mb") else "mc"
                mdev[quark] = max(mdev[quark], abs(x - xc))
            mass_dev = float(np.hypot(mdev["mb"], mdev["mc"]))

            # --- alpha_s: half the |x| spread of the two companion curves,
            #     folded symmetrically (PDF4LHC) ---
            alphas_dev = 0.0
            if fold_alphas:
                xlo, op_lo, _ = _boundary(as_lo, (flavor, mass), mass, col, open_col)
                xhi, op_hi, _ = _boundary(as_hi, (flavor, mass), mass, col, open_col)
                if op_lo or op_hi:
                    any_open = True
                if xlo is not None and xhi is not None:
                    alphas_dev = 0.5 * abs(xhi - xlo)

            sigma_up = float(np.sqrt(scale_up**2 + pdf_sigma**2 + mass_dev**2 + alphas_dev**2))
            sigma_dn = float(np.sqrt(scale_dn**2 + pdf_sigma**2 + mass_dev**2 + alphas_dev**2))

            rec[f"{col}_band_lo"] = 10.0 ** (xc - sigma_dn)
            rec[f"{col}_band_hi"] = 10.0 ** (xc + sigma_up)
            rec[f"{col}_open"] = any_open
            rec[f"{col}_scale_up_dex"] = scale_up
            rec[f"{col}_scale_dn_dex"] = scale_dn
            rec[f"{col}_pdf_sigma_dex"] = pdf_sigma
            rec[f"{col}_mass_dev_dex"] = mass_dev
            if fold_alphas:
                rec[f"{col}_alphas_dex"] = alphas_dev

        out_rows.append(rec)

    return pd.DataFrame(out_rows).sort_values(["flavor", "mass_GeV"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--registry", default=str(DEFAULT_REGISTRY),
                    help="band_registry.json from run_variation_band.py")
    ap.add_argument("--out", default=None,
                    help="output band CSV (default: <registry dir>/hnl_band.csv)")
    ap.add_argument("--alphas-lo", default=None,
                    help="as=0.117 companion sensitivity CSV (PDF4LHC alpha_s down); "
                         "folded only if --alphas-hi is also given")
    ap.add_argument("--alphas-hi", default=None,
                    help="as=0.119 companion sensitivity CSV (PDF4LHC alpha_s up)")
    args = ap.parse_args(argv)

    registry_path = Path(args.registry)
    if not registry_path.exists():
        print(f"registry not found: {registry_path}; run run_variation_band.py first")
        return 1

    df = combine_band(registry_path, args.alphas_lo, args.alphas_hi)
    out = Path(args.out) if args.out else registry_path.parent / "hnl_band.csv"
    df.to_csv(out, index=False)

    reg = json.loads(registry_path.read_text())
    axes = _axis_members(reg)
    n_sens = int(df["has_sensitivity"].sum())
    alphas_note = ", alpha_s" if (args.alphas_lo and args.alphas_hi) else ""
    print(f"combined {len(reg['variations'])} variations "
          f"(scale={len(axes['scale'])}, pdf={len(axes['pdf'])}, mass={len(axes['mass'])}{alphas_note})")
    print(f"band written: {out}  ({n_sens}/{len(df)} mass points with sensitivity)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

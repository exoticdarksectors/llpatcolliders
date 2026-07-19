#!/usr/bin/env python3
"""Read-only geometry-cache reuse audit for FONLL variation runs.

This study intentionally lives outside the production analysis package.  It does
not change any running campaign files.  It compares the central geometry cache
against an already-computed variation geometry cache and asks whether a row-wise
cache reuse would be safe.

The key diagnostic is the hit-mask mismatch.  If row ``i`` in a varied combined
CSV is assigned row ``i`` from the central geometry cache, missed/false hits
directly bias the acceptance and therefore the exclusion boundary.

Example:

    python geometry_reuse_study/geometry_reuse_audit.py \
      --central-run band_central_651bfaf4 \
      --variation-run band_scale_muR0p5_muF0p5_b0277fc5 \
      --channel-source-run central_newgrids_20260623 \
      --out geometry_reuse_study/scale_muR0p5_muF0p5_geometry_reuse.csv
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HNL_ROOT = Path(__file__).resolve().parents[1]
if str(HNL_ROOT) not in sys.path:
    sys.path.insert(0, str(HNL_ROOT))

from config_mass_grid import format_mass_for_filename  # noqa: E402
from production.combine_channels import CHANNELS       # noqa: E402


def _geom_path(run_tag: str, flavor: str, mass_label: str) -> Path:
    return (
        HNL_ROOT / "tmp" / "runs" / run_tag / "analysis" /
        "geometry_cache" / flavor / f"geom_{mass_label}.npz"
    )


def _channel_csv(run_tag: str, flavor: str, channel: str, mass_label: str) -> Path:
    return (
        HNL_ROOT / "tmp" / "runs" / run_tag / "llp_4vectors" /
        flavor / channel / f"mN_{mass_label}.csv"
    )


def _count_rows(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    # Headerless CSV.  Counting lines is faster and lower-memory than loading.
    with path.open("rb") as fh:
        return sum(1 for _ in fh)


def _channel_slices(source_run: str, flavor: str, mass_label: str) -> list[tuple[str, slice]]:
    out = []
    start = 0
    for ch in CHANNELS:
        n = _count_rows(_channel_csv(source_run, flavor, ch, mass_label))
        out.append((ch, slice(start, start + n)))
        start += n
    return out


def _safe_stats(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"median": np.nan, "p90": np.nan, "max": np.nan}
    return {
        "median": float(np.median(values)),
        "p90": float(np.quantile(values, 0.90)),
        "max": float(np.max(values)),
    }


def compare_one(
    central_run: str,
    variation_run: str,
    channel_source_run: str,
    flavor: str,
    mass: float,
) -> list[dict]:
    mass_label = format_mass_for_filename(mass)
    cp = _geom_path(central_run, flavor, mass_label)
    vp = _geom_path(variation_run, flavor, mass_label)
    if not cp.exists() or not vp.exists():
        return []

    c = np.load(cp)
    v = np.load(vp)
    chits = c["hits"].astype(bool)
    vhits = v["hits"].astype(bool)
    if len(chits) != len(vhits):
        return [{
            "flavor": flavor, "mass_GeV": mass, "channel": "ALL",
            "status": "length_mismatch", "central_rows": len(chits),
            "variation_rows": len(vhits),
        }]

    rows = []
    slices = _channel_slices(channel_source_run, flavor, mass_label)
    total_from_channels = sum(sl.stop - sl.start for _, sl in slices)
    if total_from_channels != len(chits):
        rows.append({
            "flavor": flavor, "mass_GeV": mass, "channel": "ALL",
            "status": "channel_count_mismatch",
            "central_rows": len(chits), "channel_rows": total_from_channels,
        })

    # Include ALL first, then channel-resolved rows.
    slices = [("ALL", slice(0, len(chits)))] + slices
    for ch, sl in slices:
        if sl.stop > len(chits):
            continue
        hc = chits[sl]
        hv = vhits[sl]
        n = len(hc)
        if n == 0:
            continue
        xor = hc ^ hv
        false_hits = hc & ~hv       # central says hit, variation says miss
        missed_hits = ~hc & hv      # central says miss, variation says hit
        common = hc & hv

        entry_diff = np.abs(c["entry_d"][sl][common] - v["entry_d"][sl][common])
        exit_diff = np.abs(c["exit_d"][sl][common] - v["exit_d"][sl][common])
        path_c = c["exit_d"][sl][common] - c["entry_d"][sl][common]
        path_v = v["exit_d"][sl][common] - v["entry_d"][sl][common]
        path_diff = np.abs(path_c - path_v)
        estats = _safe_stats(entry_diff)
        xstats = _safe_stats(exit_diff)
        pstats = _safe_stats(path_diff)

        vhit_n = int(hv.sum())
        rows.append({
            "flavor": flavor,
            "mass_GeV": mass,
            "mass_label": mass_label,
            "channel": ch,
            "status": "ok",
            "rows": n,
            "central_hits": int(hc.sum()),
            "variation_hits": vhit_n,
            "common_hits": int(common.sum()),
            "xor_hits": int(xor.sum()),
            "false_hits": int(false_hits.sum()),
            "missed_hits": int(missed_hits.sum()),
            "xor_per_row": float(xor.mean()),
            "xor_per_variation_hit": (
                float(xor.sum() / vhit_n) if vhit_n else np.nan
            ),
            "missed_per_variation_hit": (
                float(missed_hits.sum() / vhit_n) if vhit_n else np.nan
            ),
            "entry_absdiff_m_median": estats["median"],
            "entry_absdiff_m_p90": estats["p90"],
            "entry_absdiff_m_max": estats["max"],
            "exit_absdiff_m_median": xstats["median"],
            "exit_absdiff_m_p90": xstats["p90"],
            "exit_absdiff_m_max": xstats["max"],
            "path_absdiff_m_median": pstats["median"],
            "path_absdiff_m_p90": pstats["p90"],
            "path_absdiff_m_max": pstats["max"],
        })
    return rows


def default_masses(central_run: str, variation_run: str, flavors: list[str]) -> list[float]:
    labels = set()
    for run in (central_run, variation_run):
        for fl in flavors:
            d = HNL_ROOT / "tmp" / "runs" / run / "analysis" / "geometry_cache" / fl
            if d.exists():
                labels.update(p.stem.removeprefix("geom_") for p in d.glob("geom_*.npz"))
    masses = sorted(float(label.replace("p", ".")) for label in labels)
    return masses


def summarize(df: pd.DataFrame) -> dict:
    ok = df[(df["status"] == "ok") & (df["channel"] == "ALL")]
    by_channel = df[(df["status"] == "ok") & (df["channel"] != "ALL")]
    return {
        "n_points": int(len(ok)),
        "all_xor_per_variation_hit_median": float(ok["xor_per_variation_hit"].median()),
        "all_xor_per_variation_hit_max": float(ok["xor_per_variation_hit"].max()),
        "all_missed_per_variation_hit_median": float(ok["missed_per_variation_hit"].median()),
        "worst_points": ok.sort_values(
            "xor_per_variation_hit", ascending=False
        )[["flavor", "mass_GeV", "xor_per_variation_hit", "missed_per_variation_hit"]]
        .head(10)
        .to_dict(orient="records"),
        "by_channel_median_xor_per_variation_hit": (
            by_channel.groupby("channel")["xor_per_variation_hit"]
            .median()
            .sort_values(ascending=False)
            .to_dict()
        ),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--central-run", default="band_central_651bfaf4")
    ap.add_argument("--variation-run", required=True)
    ap.add_argument("--channel-source-run", default="central_newgrids_20260623")
    ap.add_argument("--flavor", nargs="+", default=["Ue", "Umu", "Utau"])
    ap.add_argument("--mass", nargs="+", type=float, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args(argv)

    masses = args.mass or default_masses(args.central_run, args.variation_run, args.flavor)
    records = []
    for fl in args.flavor:
        for mass in masses:
            records.extend(compare_one(
                args.central_run, args.variation_run, args.channel_source_run, fl, mass
            ))

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(records)
    df.to_csv(out, index=False)
    summary = summarize(df) if not df.empty else {}
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(f"wrote {out} ({len(df)} rows)")
    if summary:
        print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

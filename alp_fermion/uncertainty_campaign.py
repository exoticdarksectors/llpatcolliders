"""Definitions and provenance helpers for the exact BC10 uncertainty campaign."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


GRID_STEM = (
    "fonll_pp14tev_nnpdf40_nlo_as_01180_fonll_meson_"
    "dsdpTdy_pt0-50_y-3to3"
)
CBS_SCHEME_RELATIVE_AMPLITUDE = 0.20
BOUNDARIES = ("invf_min", "invf_max")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_tree(paths) -> str:
    """Hash file names and contents in stable order."""
    digest = hashlib.sha256()
    for path in sorted((Path(item) for item in paths), key=lambda item: str(item)):
        digest.update(path.name.encode())
        digest.update(b"\0")
        digest.update(sha256_file(path).encode())
        digest.update(b"\n")
    return digest.hexdigest()


def atomic_json(path: Path, payload) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def stable_seed(name: str) -> int:
    """Stable independent NumPy seed for one non-central production run."""
    digest = hashlib.sha256(f"GRENDEL-BC10-FONLL:{name}".encode()).digest()
    return int.from_bytes(digest[:4], "big") or 1


def discover_fonll_variations(grid_dir: Path) -> list[dict]:
    """Select the complete bottom-grid ensemble from its pinned manifest."""
    grid_dir = Path(grid_dir).expanduser().resolve()
    manifest_path = grid_dir / "variation_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    variations = []
    for entry in manifest["grids"]:
        if entry["quark"] != "bottom":
            continue
        tag = entry["variation_tag"]
        kind = entry["variation_kind"]
        if kind == "pdf" and int(entry["lhapdf_member"]) == 0:
            continue
        if kind not in {"central", "scale", "pdf", "mass"}:
            continue
        path = grid_dir / f"{GRID_STEM}_{tag}_bottom.dat"
        if not path.exists():
            raise FileNotFoundError(path)
        actual_sha = sha256_file(path)
        if actual_sha != entry["sha256"]:
            raise ValueError(f"FONLL checksum mismatch for {path}")
        axis = "mb" if kind == "mass" else kind
        variations.append({
            "name": tag,
            "axis": axis,
            "campaign_axis": "fonll",
            "grid_path": str(path),
            "grid_sha256": actual_sha,
            "muR": entry.get("muR"),
            "muF": entry.get("muF"),
            "pdf_member": entry.get("lhapdf_member"),
            "mb_GeV": entry.get("heavy_quark_mass_GeV"),
            "production_seed": 42 if kind == "central" else stable_seed(tag),
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": 1.0,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    order = {"central": 0, "scale": 1, "pdf": 2, "mb": 3}
    variations.sort(key=lambda item: (order[item["axis"]], item["name"]))
    counts = {axis: sum(v["axis"] == axis for v in variations) for axis in order}
    expected = {"central": 1, "scale": 6, "pdf": 100, "mb": 2}
    if counts != expected:
        raise ValueError(f"incomplete FONLL bottom ensemble: {counts} != {expected}")
    return variations


def auxiliary_variations(central_grid: Path) -> list[dict]:
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for flavor in ("u", "d", "s"):
        output.append({
            "name": f"gg_{flavor}",
            "axis": "decay_gg",
            "campaign_axis": "decay_gg",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": 42,
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": 1.0,
            "template_variant": f"gg_{flavor}",
            "production_mode": "central_vectors_exact_reuse",
        })
    for direction, scale in (
        ("down", 1.0 - CBS_SCHEME_RELATIVE_AMPLITUDE),
        ("up", 1.0 + CBS_SCHEME_RELATIVE_AMPLITUDE),
    ):
        output.append({
            "name": f"cbs_{direction}",
            "axis": "cbs",
            "campaign_axis": "cbs",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": 42,
            "reco_seed_offset": 0,
            "cbs_amplitude_scale": scale,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return output


def decay_structure_variation(central_grid: Path) -> dict:
    """Exact 2310 decay model as a one-sided structural comparison."""
    grid_path = str(Path(central_grid).resolve())
    return {
        "name": "decay_2310_structural",
        "axis": "decay_structure",
        "campaign_axis": "decay_structure",
        "grid_path": grid_path,
        "grid_sha256": sha256_file(Path(central_grid)),
        "production_seed": 42,
        "reco_seed_offset": 0,
        "cbs_amplitude_scale": 1.0,
        "template_variant": "decay_2310_structural",
        "decay_model": "2310_structural",
        "production_mode": "central_vectors_exact_reuse",
        "combination_role": "one_sided_structural_comparison_outside_halo",
    }


def numerical_control_variations(central_grid: Path) -> list[dict]:
    """Same-physics central repeats used only to measure numerical spread."""
    grid_path = str(Path(central_grid).resolve())
    grid_sha = sha256_file(Path(central_grid))
    output = []
    for index in (1, 2):
        name = f"central_repeat_{index}"
        output.append({
            "name": name,
            "axis": "numerical_control",
            "campaign_axis": "numerical_control",
            "grid_path": grid_path,
            "grid_sha256": grid_sha,
            "production_seed": stable_seed(f"{name}:production"),
            "reco_seed_offset": stable_seed(f"{name}:reconstruction"),
            "cbs_amplitude_scale": 1.0,
            "template_variant": "central",
            "production_mode": "fresh_600k",
        })
    return output


def all_variations(grid_dir: Path) -> list[dict]:
    fonll = discover_fonll_variations(grid_dir)
    central_grid = Path(fonll[0]["grid_path"])
    return [
        *fonll,
        *auxiliary_variations(central_grid),
        decay_structure_variation(central_grid),
        *numerical_control_variations(central_grid),
    ]


def _log_boundary(raw, mass, variation, boundary):
    selected = raw[
        (raw["mass_GeV"] == mass) & (raw["variation"] == variation)
    ]
    if len(selected) != 1:
        return None, False, False
    row = selected.iloc[0]
    sensitive = bool(row["has_sensitivity"])
    opened = bool(row.get(f"{boundary}_open", False))
    value = float(row.get(boundary, np.nan))
    if not sensitive or not np.isfinite(value) or value <= 0.0:
        return None, opened, sensitive
    return float(np.log10(value)), opened, sensitive


def _combine_campaign_band(raw: pd.DataFrame) -> pd.DataFrame:
    """Build a pointwise one-source-at-a-time variation envelope.

    Named alternatives use their extrema. The NNPDF ensemble uses its 16th and
    84th percentiles; its raw extrema are deliberately not used as a headline
    interval. No sources are combined in quadrature and this is not a
    confidence interval. Same-physics central repeats are excluded from the
    envelope and reported separately as numerical controls.
    """
    axes = {
        axis: list(raw.loc[raw["axis"] == axis, "variation"].unique())
        for axis in (
            "central", "scale", "pdf", "mb", "decay_gg", "cbs",
            "decay_structure", "numerical_control",
        )
    }
    expected_counts = {
        "central": 1, "scale": 6, "pdf": 100, "mb": 2,
        "decay_gg": 3, "cbs": 2, "decay_structure": 1,
        "numerical_control": 2,
    }
    counts = {axis: len(names) for axis, names in axes.items()}
    if counts != expected_counts:
        raise ValueError(f"incomplete contour ensemble: {counts} != {expected_counts}")
    halo_variations = [
        name for axis in ("scale", "pdf", "mb", "decay_gg", "cbs")
        for name in axes[axis]
    ]
    rows = []
    for mass in sorted(raw["mass_GeV"].unique()):
        central = raw[
            (raw["mass_GeV"] == mass) & (raw["variation"] == "central")
        ].iloc[0]
        structural_name = axes["decay_structure"][0]
        structural = raw[
            (raw["mass_GeV"] == mass)
            & (raw["variation"] == structural_name)
        ].iloc[0]
        central_sensitive = bool(central["has_sensitivity"])
        structural_sensitive = bool(structural["has_sensitivity"])
        numerical_rows = raw[
            (raw["mass_GeV"] == mass)
            & raw["variation"].isin(axes["numerical_control"])
        ].set_index("variation")
        numerical_sensitive = {
            name: bool(numerical_rows.loc[name, "has_sensitivity"])
            for name in axes["numerical_control"]
        }
        numerical_topology_differences = {
            name for name, sensitive in numerical_sensitive.items()
            if sensitive != central_sensitive
        }
        halo_rows = raw[
            (raw["mass_GeV"] == mass)
            & raw["variation"].isin(halo_variations)
        ].set_index("variation")
        halo_sensitive = {
            name: bool(halo_rows.loc[name, "has_sensitivity"])
            for name in halo_variations
        }
        halo_restores_sensitivity = {
            name for name, sensitive in halo_sensitive.items()
            if sensitive and not central_sensitive
        }
        halo_removes_sensitivity = {
            name for name, sensitive in halo_sensitive.items()
            if central_sensitive and not sensitive
        }
        halo_topology_differences = (
            halo_restores_sensitivity | halo_removes_sensitivity
        )
        any_halo_variation_sensitive = bool(
            raw.loc[
                (raw["mass_GeV"] == mass)
                & ~raw["axis"].isin(
                    ["numerical_control", "decay_structure"]
                ),
                "has_sensitivity",
            ].any()
        )
        record = {
            "mass_GeV": mass,
            "has_sensitivity": central_sensitive,
            # Backward-compatible name plus an explicit halo-only alias.
            "any_variation_sensitive": any_halo_variation_sensitive,
            "any_halo_variation_sensitive": any_halo_variation_sensitive,
            "decay_structure_variation": structural_name,
            "decay_structure_has_sensitivity": structural_sensitive,
            "decay_structure_restores_sensitivity": bool(
                structural_sensitive and not central_sensitive
            ),
            "decay_structure_removes_sensitivity": bool(
                central_sensitive and not structural_sensitive
            ),
            "decay_structure_included_in_halo": False,
            "numerical_control_n_sensitive": sum(numerical_sensitive.values()),
            "numerical_control_any_sensitive": any(numerical_sensitive.values()),
            "numerical_control_all_sensitive": all(numerical_sensitive.values()),
            "numerical_control_sensitive_variations": ";".join(
                name for name in axes["numerical_control"]
                if numerical_sensitive[name]
            ),
            "numerical_control_included_in_halo": False,
            "halo_restores_sensitivity": bool(halo_restores_sensitivity),
            "halo_restores_sensitivity_variations": ";".join(
                name for name in halo_variations
                if name in halo_restores_sensitivity
            ),
            "halo_removes_sensitivity": bool(halo_removes_sensitivity),
            "halo_removes_sensitivity_variations": ";".join(
                name for name in halo_variations
                if name in halo_removes_sensitivity
            ),
            "envelope_definition": "single_source_variation_envelope",
        }
        topology_differs = structural_sensitive != central_sensitive
        for boundary in BOUNDARIES:
            xc, central_open, _ = _log_boundary(
                raw, mass, "central", boundary
            )
            xs, structural_open, _ = _log_boundary(
                raw, mass, structural_name, boundary
            )
            record[f"{boundary}_central"] = float(central.get(boundary, np.nan))
            record[f"{boundary}_open"] = central_open
            record[f"{boundary}_decay_structure"] = (
                10.0 ** xs if xs is not None else np.nan
            )
            record[f"{boundary}_decay_structure_open"] = structural_open
            if central_sensitive and structural_sensitive:
                topology_differs = topology_differs or (
                    central_open != structural_open
                )
            for name in halo_variations:
                _, opened, sensitive = _log_boundary(
                    raw, mass, name, boundary
                )
                if (
                    central_sensitive and sensitive
                    and central_open != opened
                ):
                    halo_topology_differences.add(name)
            repeat_samples = []
            for name in axes["numerical_control"]:
                value, opened, sensitive = _log_boundary(
                    raw, mass, name, boundary
                )
                repeat_samples.append((name, value, opened, sensitive))
                record[f"{boundary}_{name}"] = (
                    10.0 ** value if value is not None else np.nan
                )
                if (
                    central_sensitive and sensitive
                    and central_open != opened
                ):
                    numerical_topology_differences.add(name)
            if xc is None:
                record[f"{boundary}_envelope_lo"] = np.nan
                record[f"{boundary}_envelope_hi"] = np.nan
                record[f"{boundary}_any_variation_open"] = central_open
                record[f"{boundary}_variation_missing"] = False
                record[f"{boundary}_envelope_lo_source"] = ""
                record[f"{boundary}_envelope_hi_source"] = ""
                record[f"{boundary}_physical_envelope_max_abs_dex"] = np.nan
                record[f"{boundary}_repeat_median_abs_dex"] = np.nan
                record[f"{boundary}_repeat_max_abs_dex"] = np.nan
                record[f"{boundary}_repeat_median_abs_fraction"] = np.nan
                record[f"{boundary}_repeat_max_abs_fraction"] = np.nan
                record[f"{boundary}_repeat_missing"] = False
                record[f"{boundary}_repeat_not_subdominant"] = False
                continue

            missing = False
            any_open = central_open

            def values(names):
                nonlocal missing, any_open
                result = []
                for name in names:
                    value, opened, sensitive = _log_boundary(
                        raw, mass, name, boundary
                    )
                    any_open = any_open or opened
                    if value is None:
                        missing = True
                    else:
                        result.append((name, value))
                return result

            source_intervals = {}
            for axis, label in (
                ("scale", "scale"),
                ("mb", "mb"),
                ("decay_gg", "gg"),
                ("cbs", "cbs"),
            ):
                candidates = [("central", xc), *values(axes[axis])]
                lo_name, lo = min(candidates, key=lambda item: item[1])
                hi_name, hi = max(candidates, key=lambda item: item[1])
                source_intervals[label] = (lo, hi)
                record.update({
                    f"{boundary}_{label}_lo": 10.0 ** lo,
                    f"{boundary}_{label}_hi": 10.0 ** hi,
                    f"{boundary}_{label}_lo_variation": lo_name,
                    f"{boundary}_{label}_hi_variation": hi_name,
                })

            pdf_pairs = values(axes["pdf"])
            pdf = np.asarray([value for _, value in pdf_pairs], dtype=float)
            if len(pdf):
                pdf_p16, pdf_p84 = np.quantile(pdf, [0.16, 0.84])
                pdf_std = float(np.std(pdf, ddof=1)) if len(pdf) >= 2 else np.nan
            else:
                pdf_p16 = pdf_p84 = xc
                pdf_std = np.nan
            pdf_lo = min(xc, float(pdf_p16))
            pdf_hi = max(xc, float(pdf_p84))
            source_intervals["pdf"] = (pdf_lo, pdf_hi)
            record.update({
                f"{boundary}_pdf_p16": 10.0 ** float(pdf_p16),
                f"{boundary}_pdf_p84": 10.0 ** float(pdf_p84),
                f"{boundary}_pdf_log10_std": pdf_std,
                f"{boundary}_pdf_n_finite": len(pdf),
            })

            lo_source, lo = min(
                ((source, interval[0]) for source, interval in source_intervals.items()),
                key=lambda item: item[1],
            )
            hi_source, hi = max(
                ((source, interval[1]) for source, interval in source_intervals.items()),
                key=lambda item: item[1],
            )
            record.update({
                f"{boundary}_envelope_lo": 10.0 ** lo,
                f"{boundary}_envelope_hi": 10.0 ** hi,
                f"{boundary}_any_variation_open": any_open,
                f"{boundary}_variation_missing": missing,
                f"{boundary}_envelope_lo_source": lo_source,
                f"{boundary}_envelope_hi_source": hi_source,
            })

            repeat_values = []
            repeat_missing = False
            for name, value, _, _ in repeat_samples:
                if value is None:
                    repeat_missing = True
                else:
                    repeat_values.append(value)
            repeat_abs_dex = np.abs(np.asarray(repeat_values) - xc)
            repeat_abs_fraction = np.abs(10.0 ** (
                np.asarray(repeat_values) - xc
            ) - 1.0)
            physical_max_abs_dex = max(xc - lo, hi - xc)
            repeat_max_abs_dex = (
                float(np.max(repeat_abs_dex)) if len(repeat_abs_dex) else np.nan
            )
            record.update({
                f"{boundary}_physical_envelope_max_abs_dex": physical_max_abs_dex,
                f"{boundary}_repeat_median_abs_dex": (
                    float(np.median(repeat_abs_dex))
                    if len(repeat_abs_dex) else np.nan
                ),
                f"{boundary}_repeat_max_abs_dex": repeat_max_abs_dex,
                f"{boundary}_repeat_median_abs_fraction": (
                    float(np.median(repeat_abs_fraction))
                    if len(repeat_abs_fraction) else np.nan
                ),
                f"{boundary}_repeat_max_abs_fraction": (
                    float(np.max(repeat_abs_fraction))
                    if len(repeat_abs_fraction) else np.nan
                ),
                f"{boundary}_repeat_missing": repeat_missing,
                f"{boundary}_repeat_not_subdominant": bool(
                    repeat_missing
                    or (
                        np.isfinite(repeat_max_abs_dex)
                        and repeat_max_abs_dex > 0.0
                        and repeat_max_abs_dex >= physical_max_abs_dex
                    )
                ),
            })
        record["decay_structure_topology_differs"] = bool(topology_differs)
        record["halo_topology_differs"] = bool(halo_topology_differences)
        record["halo_topology_difference_variations"] = ";".join(
            name for name in halo_variations
            if name in halo_topology_differences
        )
        record["numerical_control_topology_differs"] = bool(
            numerical_topology_differences
        )
        record["numerical_control_topology_difference_variations"] = ";".join(
            name for name in axes["numerical_control"]
            if name in numerical_topology_differences
        )
        rows.append(record)
    return pd.DataFrame(rows).sort_values("mass_GeV")


def combine_band(
    raw: pd.DataFrame,
    central_curve: pd.DataFrame | str | Path | None = None,
) -> pd.DataFrame:
    """Build the campaign band and optionally rebase it to a canonical contour.

    The expensive campaign measures one-source shifts relative to its independent
    nominal sample.  A higher-statistics central contour can therefore replace
    the nominal values without discarding those shifts.  Only the pointwise halo
    components are rebased; the exact 2310 structural contour and same-physics
    numerical repeats retain their independently simulated absolute values.
    """
    band = _combine_campaign_band(raw)
    if central_curve is None:
        return band

    reference = (
        pd.read_csv(central_curve)
        if isinstance(central_curve, (str, Path))
        else central_curve.copy()
    )
    if reference["mass_GeV"].duplicated().any():
        raise ValueError("canonical central curve contains duplicate masses")
    reference = reference.set_index("mass_GeV")
    missing = sorted(set(band["mass_GeV"]) - set(reference.index.astype(float)))
    if missing:
        raise ValueError(
            f"canonical central curve is missing {len(missing)} campaign masses; "
            f"first: {missing[0]}"
        )

    halo_value_columns = {
        boundary: [
            *(f"{boundary}_{source}_{edge}"
              for source in ("scale", "mb", "gg", "cbs")
              for edge in ("lo", "hi")),
            f"{boundary}_pdf_p16",
            f"{boundary}_pdf_p84",
            f"{boundary}_envelope_lo",
            f"{boundary}_envelope_hi",
        ]
        for boundary in BOUNDARIES
    }

    for index, campaign in band.iterrows():
        mass = float(campaign["mass_GeV"])
        canonical = reference.loc[mass]
        campaign_sensitive = bool(campaign["has_sensitivity"])
        canonical_sensitive = bool(canonical.get("has_sensitivity", False))
        band.at[index, "campaign_has_sensitivity"] = campaign_sensitive
        band.at[index, "has_sensitivity"] = canonical_sensitive
        band.at[index, "any_variation_sensitive"] = bool(
            campaign["any_variation_sensitive"] or canonical_sensitive
        )
        band.at[index, "any_halo_variation_sensitive"] = bool(
            campaign["any_halo_variation_sensitive"] or canonical_sensitive
        )

        topology_compatible = campaign_sensitive == canonical_sensitive
        for boundary in BOUNDARIES:
            campaign_value = float(campaign[f"{boundary}_central"])
            canonical_value = float(canonical.get(boundary, np.nan))
            campaign_open = bool(campaign[f"{boundary}_open"])
            canonical_open = bool(canonical.get(f"{boundary}_open", False))
            band.at[index, f"{boundary}_campaign_central"] = campaign_value
            band.at[index, f"{boundary}_central"] = canonical_value
            band.at[index, f"{boundary}_open"] = canonical_open

            compatible = (
                campaign_sensitive
                and canonical_sensitive
                and np.isfinite(campaign_value)
                and campaign_value > 0.0
                and np.isfinite(canonical_value)
                and canonical_value > 0.0
                and campaign_open == canonical_open
                and not campaign_open
            )
            topology_compatible = topology_compatible and (
                compatible or not campaign_sensitive
            )
            if compatible:
                factor = canonical_value / campaign_value
                for column in halo_value_columns[boundary]:
                    value = band.at[index, column]
                    if np.isfinite(value):
                        band.at[index, column] = float(value) * factor
            else:
                for column in halo_value_columns[boundary]:
                    band.at[index, column] = np.nan
                band.at[index, f"{boundary}_variation_missing"] = True
        band.at[index, "canonical_rebase_topology_compatible"] = bool(
            topology_compatible
        )

    band["envelope_reference"] = "canonical_high_statistics_central"
    return band


def structural_alternative_from_curves(
    central: pd.DataFrame | str | Path,
    structural: pd.DataFrame | str | Path,
    variation: str = "decay_2310_structural",
) -> pd.DataFrame:
    """Compare exact central/structural curves on one identical mass grid."""
    central = (
        pd.read_csv(central)
        if isinstance(central, (str, Path))
        else central.copy()
    )
    structural = (
        pd.read_csv(structural)
        if isinstance(structural, (str, Path))
        else structural.copy()
    )
    keep = [
        "mass_GeV", "has_sensitivity", "peak_N", "peak_invf",
        "invf_min", "invf_max", "invf_min_open", "invf_max_open",
    ]
    for label, frame in (("central", central), ("structural", structural)):
        missing = sorted(set(keep) - set(frame.columns))
        if missing:
            raise ValueError(f"{label} curve is missing columns: {missing}")
        if frame["mass_GeV"].duplicated().any():
            raise ValueError(f"{label} curve contains duplicate masses")
    central_masses = np.asarray(central["mass_GeV"], dtype=float)
    structural_masses = np.asarray(structural["mass_GeV"], dtype=float)
    if (
        len(central_masses) != len(structural_masses)
        or not np.array_equal(np.sort(central_masses), np.sort(structural_masses))
    ):
        raise ValueError("central/decay-structure mass grids differ")

    central_state = central[keep].rename(columns={
        column: f"central_{column}"
        for column in keep if column != "mass_GeV"
    })
    result = structural[keep].copy()
    result.insert(0, "variation", variation)
    result = result.merge(
        central_state, on="mass_GeV", validate="one_to_one"
    )
    result.insert(1, "axis", "decay_structure")
    result["restores_sensitivity"] = (
        result["has_sensitivity"] & ~result["central_has_sensitivity"]
    )
    result["removes_sensitivity"] = (
        ~result["has_sensitivity"] & result["central_has_sensitivity"]
    )
    open_differs = (
        result["has_sensitivity"]
        & result["central_has_sensitivity"]
        & (
            (result["invf_min_open"] != result["central_invf_min_open"])
            | (result["invf_max_open"] != result["central_invf_max_open"])
        )
    )
    result["topology_differs"] = (
        result["restores_sensitivity"]
        | result["removes_sensitivity"]
        | open_differs
    )
    result["comparison_role"] = "one_sided_structural_model_comparison"
    result["included_in_pointwise_halo"] = False
    return result.sort_values("mass_GeV")


def structural_alternative_table(raw: pd.DataFrame) -> pd.DataFrame:
    """Publish the campaign-grid structural contour with explicit topology state."""
    central = raw.loc[raw["axis"] == "central"].copy()
    structural = raw.loc[raw["axis"] == "decay_structure"].copy()
    if central["variation"].nunique() != 1 or len(central) != len(structural):
        raise ValueError("central/decay-structure contour counts are incomplete")
    if structural["variation"].nunique() != 1:
        raise ValueError("expected exactly one decay-structure variation")
    return structural_alternative_from_curves(
        central,
        structural,
        variation=str(structural["variation"].iloc[0]),
    )

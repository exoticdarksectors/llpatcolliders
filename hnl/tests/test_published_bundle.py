"""Integrity and topology checks for the tracked HNL diagnostic bundle."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"
BUNDLE = PUBLISHED / "bundle"


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def test_published_bundle_hashes_match_manifest():
    # The bundle is consistent with the published central: its central_curve hash
    # equals the live curve (the 0.305 GeV Ue/Umu band points were re-derived for
    # the 2026-07-18 kaon rerun), and every bundle file matches its recorded hash.
    manifest = json.loads((BUNDLE / "MANIFEST.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["central_curve"]["kaon_curve_publication_revision"] == "61aa73c"
    assert manifest["central_curve"]["sha256"] == _sha256(
        PUBLISHED / "grendel_hnl_sensitivity.csv")
    for name, metadata in manifest["files"].items():
        path = BUNDLE / name
        assert path.exists(), name
        assert metadata["sha256"] == _sha256(path), name

    fonll_manifest = json.loads((BUNDLE / "FONLL_MANIFEST.json").read_text())
    for metadata in fonll_manifest["files"].values():
        path = BUNDLE / metadata["file"]
        assert metadata["sha256"] == _sha256(path), metadata["file"]


def test_fonll_bundle_has_every_member_and_stable_topology():
    band = pd.read_csv(BUNDLE / "hnl_band_fonll.csv")
    raw = pd.read_csv(BUNDLE / "hnl_band_fonll_raw.csv")
    manifest = json.loads((BUNDLE / "FONLL_MANIFEST.json").read_text())
    assert manifest["schema_version"] == 2

    assert len(band) == 54
    assert len(raw) == 111 * 54
    assert raw["variation_name"].nunique() == 111
    assert set(raw.groupby("variation_name").size()) == {54}
    assert len(manifest["base_variations"]) == 111
    assert {entry["name"] for entry in manifest["base_variations"]} == set(
        raw["variation_name"])

    partial = manifest["partial_recomputes"]
    assert len(partial) == 1
    partial = partial[0]
    assert partial["affected_rows"] == 222
    assert partial["variation_count"] == 111
    replaced = (raw["mass_GeV"].eq(0.305)
                & raw["flavor"].isin(["Ue", "Umu"]))
    assert replaced.sum() == partial["affected_rows"]
    assert raw.loc[replaced, "run_tag"].str.endswith(
        partial["recompute_campaign_config_suffix"]).all()
    assert raw.loc[~replaced, "run_tag"].str.endswith(
        partial["base_campaign_config_suffix"]).all()
    assert set(raw.loc[replaced].groupby("variation_name").size()) == {2}

    masked = band[band["mass_GeV"].eq(0.305)
                  & band["flavor"].isin(["Ue", "Umu"])]
    assert len(masked) == 2
    assert masked[["u2_min_band_lo", "u2_min_band_hi"]].isna().all().all()
    for boundary in ("u2_min", "u2_max"):
        assert not band[f"{boundary}_topology_changed"].any()
        assert set(band[f"{boundary}_n_members_expected"]) == {110}
        assert set(band[f"{boundary}_n_members_sensitive"]) == {110}


def test_channel_fractions_and_endpoint_controls_are_consistent():
    from analysis.channel_breakdown import DEFAULT_MASSES

    channels = pd.read_csv(BUNDLE / "channel_breakdown_u2min.csv")
    assert min(DEFAULT_MASSES) == 0.305
    assert channels["mass_GeV"].min() == 0.305
    grouped = channels.groupby(["flavor", "mass_GeV", "boundary"])
    assert set(grouped["channel"].nunique()) == {8}
    assert np.allclose(grouped["frac"].sum(), 1.0, rtol=0.0, atol=1e-12)

    controls = json.loads((BUNDLE / "NUMERICAL_CONTROLS.json").read_text())
    assert controls["worker_count_independence"][
        "single_vs_two_worker_bitwise_equal"]
    assert controls["bc_endpoint_exact400"]["all_topologies_agree"]


def test_kaon_transport_band_is_anchored_to_published_central():
    band = pd.read_csv(BUNDLE / "kaon_desc_band.csv")
    central = pd.read_csv(PUBLISHED / "grendel_hnl_sensitivity.csv")
    merged = band.merge(
        central[["flavor", "mass_GeV", "u2_min"]],
        on=["flavor", "mass_GeV"], validate="one_to_one")

    assert np.array_equal(merged["u2_min_desc1p5"], merged["u2_min"])
    assert (merged["u2_min_desc1p0"] > merged["u2_min_desc1p5"]).all()
    assert (merged["u2_min_desc1p5"] > merged["u2_min_desc3p0"]).all()


def test_canonical_manifest_provenance():
    manifest = json.loads((PUBLISHED / "MANIFEST.json").read_text())
    # current publication is the 2026-07-18 kaon rerun
    assert manifest["producer_git_sha"] == "80668e6"
    # the earlier beta=p/E fix is preserved as prior_correction
    assert "fad3a8b" in manifest["prior_correction"]
    # The base FONLL campaign revision is preserved; the structured partial
    # recompute contract is checked above.
    assert manifest["diagnostic_bundle"]["campaign_git_sha"].startswith(
        "6482484")

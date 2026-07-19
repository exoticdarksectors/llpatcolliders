"""Output-path policy for the BC10 ALP package.

Generated artifacts default to ``alp_fermion/tmp`` so the source tree stays
small. Override the whole tree with ``ALP_TMP_DIR`` or select shared templates,
production vectors, analysis output, and geometry caches independently with
``ALP_TEMPLATE_DIR``, ``ALP_LLP_VECTORS_DIR``, ``ALP_ANALYSIS_DIR``, and
``ALP_GEOM_CACHE_DIR``. The uncertainty campaign uses these explicit paths for
large scratch runs; no repository symlinks are required. The shared GRENDEL
reconstruction and HNL FONLL sampler are imported from the sibling ``higgs/``
and ``hnl/`` packages (never copied).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ALP_ROOT = Path(__file__).resolve().parent
REPO_ROOT = ALP_ROOT.parent
HNL_ROOT = REPO_ROOT / "hnl"
HIGGS_DIR = REPO_ROOT / "higgs"

TMP_DIR = Path(os.environ.get("ALP_TMP_DIR", ALP_ROOT / "tmp")).expanduser()
LLP_VECTORS_DIR = Path(os.environ.get(
    "ALP_LLP_VECTORS_DIR", TMP_DIR / "llp_4vectors"
)).expanduser()                                  # a four-vector CSVs per mass
TEMPLATE_DIR = Path(os.environ.get(
    "ALP_TEMPLATE_DIR", TMP_DIR / "decay_templates"
)).expanduser()                                  # ALP rest-frame templates
ANALYSIS_DIR = Path(os.environ.get(
    "ALP_ANALYSIS_DIR", TMP_DIR / "analysis"
)).expanduser()                                  # island CSV + plots
GEOM_CACHE_DIR = Path(os.environ.get(
    "ALP_GEOM_CACHE_DIR", ANALYSIS_DIR / "geometry_cache"
)).expanduser()


def add_hnl_to_path():
    """Put the sibling ``hnl/`` package on sys.path so its production +
    analysis modules can be imported (they self-register ``higgs/``)."""
    if not HNL_ROOT.exists():
        raise FileNotFoundError(
            f"hnl/ not found at {HNL_ROOT}; the BC10 package reuses the HNL "
            "FONLL sampler + GRENDEL reconstruction. Check out PR #15's head."
        )
    p = str(HNL_ROOT)
    if p not in sys.path:
        sys.path.insert(0, p)

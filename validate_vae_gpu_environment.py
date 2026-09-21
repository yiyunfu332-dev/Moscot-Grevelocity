"""Fail-fast import/provenance check for the project's Python environment."""

from __future__ import annotations

import importlib
import importlib.metadata as metadata
import os
import sys
from pathlib import Path


CACHE_ROOT = Path("/tmp") / f"zebrafish_growth_{os.getuid()}"
for name in ("numba", "matplotlib", "xdg"):
    (CACHE_ROOT / name).mkdir(parents=True, exist_ok=True)
os.environ.setdefault("NUMBA_CACHE_DIR", str(CACHE_ROOT / "numba"))
os.environ.setdefault("MPLCONFIGDIR", str(CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_ROOT / "xdg"))

# mudata 0.3.8 accesses these private anndata 0.10.9 modules as top-level
# attributes. Importing and attaching them before mudata is deterministic and
# replaces the ad-hoc repair previously performed in interactive cells.
import anndata as ad

ad_core = importlib.import_module("anndata._core")
ad_file_backing = importlib.import_module("anndata._core.file_backing")
setattr(ad_core, "file_backing", ad_file_backing)
setattr(ad, "_core", ad_core)

import graphvelo
import mudata
import moscot
import scanpy
from moscot.problems.time import TemporalProblem


EXPECTED = {
    "anndata": "0.10.9",
    "scanpy": "1.11.5",
    "mudata": "0.3.8",
    "moscot": "0.5.1",
    "graphvelo": "0.1.9",
    "ott-jax": "0.6.0",
    "jax": "0.10.2",
    "numpy": "2.4.6",
    "scipy": "1.16.3",
    "pandas": "2.3.3",
}


rows = []
wrong = []
for package, expected in EXPECTED.items():
    distribution = metadata.distribution(package)
    actual = distribution.version
    location = str(distribution.locate_file(""))
    rows.append((package, expected, actual, location))
    if actual != expected:
        wrong.append(f"{package}: expected {expected}, found {actual}")

print("Python executable:", sys.executable)
print("Python version:", sys.version.replace("\n", " "))
for package, expected, actual, location in rows:
    print(f"{package:12s} expected={expected:8s} actual={actual:8s} {location}")
print("TemporalProblem:", TemporalProblem)
print("anndata private-module compatibility: OK")
print("Core imports: OK")
if wrong:
    raise RuntimeError("Version mismatch:\n" + "\n".join(wrong))

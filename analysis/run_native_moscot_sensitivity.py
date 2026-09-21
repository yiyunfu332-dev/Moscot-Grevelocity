"""One-factor-at-a-time sensitivity analysis for the current native moscot M1 model.

The script reuses the saved expanded-panel checkpoint and compares every refit with
the saved primary M1 couplings. It does not rescore genes or modify the checkpoint.
"""

from __future__ import annotations

import importlib
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_native_moscot_sensitivity")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_native_moscot_sensitivity")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg_native_moscot_sensitivity")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import jensenshannon
from scipy.stats import spearmanr

# Compatibility for the installed mudata/anndata combination.
ad_core = importlib.import_module("anndata._core")
ad_file_backing = importlib.import_module("anndata._core.file_backing")
setattr(ad_core, "file_backing", ad_file_backing)
setattr(ad, "_core", ad_core)

from moscot.problems.time import TemporalProblem


ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT = (
    ROOT
    / "blood_growth_graphvelo_moscot_v2/checkpoints/"
    / "zebrafish_growth_graphvelo_moscot_v2.h5ad"
)
RESULTS = ROOT / "blood_growth_graphvelo_moscot_v2/results"
COUPLINGS = RESULTS / "couplings_no_topk"

BASE = {"growth_scaling": 5.0, "epsilon": 1e-2, "tau_a": 0.9, "tau_b": 0.9}
RUNS = [
    {"run": "primary_refit", **BASE},
    *[
        {"run": f"growth_scaling_{value:g}", **BASE, "growth_scaling": value}
        for value in (1.0, 2.5, 10.0, 20.0)
    ],
    *[
        {"run": f"epsilon_{value:g}", **BASE, "epsilon": value}
        for value in (5e-3, 5e-2, 1e-1)
    ],
    *[
        {"run": f"tau_{value:g}", **BASE, "tau_a": value, "tau_b": value}
        for value in (0.8, 0.95, 1.0)
    ],
]


def as_csr(solution) -> sp.csr_matrix:
    matrix = solution.transport_matrix
    if sp.issparse(matrix):
        out = matrix.tocsr().astype(float)
    else:
        out = sp.csr_matrix(np.asarray(matrix, dtype=float))
    out.eliminate_zeros()
    assert np.isfinite(out.data).all() and (out.data >= 0).all()
    return out


def type_mass(matrix: sp.csr_matrix, source_codes, target_codes, n_types: int):
    coo = matrix.tocoo()
    aggregate = np.zeros((n_types, n_types), dtype=float)
    np.add.at(
        aggregate,
        (source_codes[coo.row], target_codes[coo.col]),
        coo.data,
    )
    total = aggregate.sum()
    assert total > 0
    joint = aggregate / total
    source = joint.sum(axis=1)
    conditional = np.divide(
        joint,
        source[:, None],
        out=np.zeros_like(joint),
        where=source[:, None] > 0,
    )
    return joint, source, conditional


adata = ad.read_h5ad(CHECKPOINT)
required_obs = {
    "hpf", "time_days", "official_fine_cell_type",
    "zf_proliferation_score", "zf_apoptosis_score",
    "moscot_posterior_growth_rate",
}
assert required_obs.issubset(adata.obs.columns)
assert "X_moscot" in adata.obsm

hpf = pd.to_numeric(adata.obs["hpf"]).to_numpy(float)
time_days = pd.to_numeric(adata.obs["time_days"]).to_numpy(float)
hpf_values = np.sort(np.unique(hpf))
hpf_pairs = list(zip(hpf_values[:-1], hpf_values[1:]))
day_pairs = [(a / 24.0, b / 24.0) for a, b in hpf_pairs]

cell_types = pd.Categorical(adata.obs["official_fine_cell_type"])
type_labels = cell_types.categories.astype(str).tolist()
type_codes = cell_types.codes
n_types = len(type_labels)
saved_posterior = pd.to_numeric(
    adata.obs["moscot_posterior_growth_rate"]
).to_numpy(float).reshape(-1)

reference = {}
for (t0, t1), (d0, d1) in zip(hpf_pairs, day_pairs):
    name = f"{t0:g}_to_{t1:g}_hpf"
    source_idx = np.flatnonzero(hpf == t0)
    target_idx = np.flatnonzero(hpf == t1)
    source_ids = pd.read_csv(COUPLINGS / f"{name}_source_cells.csv")["cell_id"].astype(str)
    target_ids = pd.read_csv(COUPLINGS / f"{name}_target_cells.csv")["cell_id"].astype(str)
    assert source_ids.tolist() == adata.obs_names[source_idx].astype(str).tolist()
    assert target_ids.tolist() == adata.obs_names[target_idx].astype(str).tolist()
    matrix = sp.load_npz(COUPLINGS / f"{name}_M1_raw.npz").tocsr()
    assert matrix.shape == (len(source_idx), len(target_idx))
    reference[(d0, d1)] = {
        "matrix": matrix,
        "source_idx": source_idx,
        "target_idx": target_idx,
        "row_mass": np.asarray(matrix.sum(axis=1)).ravel(),
        "type": type_mass(
            matrix,
            type_codes[source_idx],
            type_codes[target_idx],
            n_types,
        ),
    }

rows = []
for spec in RUNS:
    print("Fitting", spec, flush=True)
    problem = TemporalProblem(adata)
    problem.proliferation_key = "zf_proliferation_score"
    problem.apoptosis_key = "zf_apoptosis_score"
    problem = problem.prepare(
        time_key="time_days",
        joint_attr="X_moscot",
        policy="sequential",
        cost="sq_euclidean",
        marginal_kwargs={"scaling": spec["growth_scaling"]},
    )
    problem = problem.solve(
        epsilon=spec["epsilon"],
        tau_a=spec["tau_a"],
        tau_b=spec["tau_b"],
    )
    posterior_frame = problem.posterior_growth_rates
    assert posterior_frame.index.is_unique and posterior_frame.shape[1] == 1
    posterior = (
        posterior_frame.iloc[:, 0]
        .reindex(adata.obs_names)
        .to_numpy(dtype=float)
        .reshape(-1)
    )
    posterior_valid = np.isfinite(saved_posterior) & np.isfinite(posterior)
    posterior_rho = (
        spearmanr(
            saved_posterior[posterior_valid], posterior[posterior_valid]
        ).statistic
        if posterior_valid.sum() > 1 else np.nan
    )

    for (t0, t1), (d0, d1) in zip(hpf_pairs, day_pairs):
        ref = reference[(d0, d1)]
        solution = problem.solutions[(d0, d1)]
        matrix = as_csr(solution)
        assert matrix.shape == ref["matrix"].shape
        row_mass = np.asarray(matrix.sum(axis=1)).ravel()
        joint, source_mass, conditional = type_mass(
            matrix,
            type_codes[ref["source_idx"]],
            type_codes[ref["target_idx"]],
            n_types,
        )
        ref_joint, ref_source_mass, ref_conditional = ref["type"]
        valid = np.isfinite(row_mass) & np.isfinite(ref["row_mass"])
        row_rho = (
            spearmanr(ref["row_mass"][valid], row_mass[valid]).statistic
            if valid.sum() > 1 else np.nan
        )
        js_values = []
        js_weights = []
        for code in range(n_types):
            if ref_source_mass[code] > 0 and source_mass[code] > 0:
                js_values.append(
                    jensenshannon(
                        ref_conditional[code], conditional[code], base=2.0
                    ) ** 2
                )
                js_weights.append(ref_source_mass[code])
        rows.append(
            {
                **spec,
                "source_hpf": t0,
                "target_hpf": t1,
                "converged": bool(getattr(solution, "converged", True)),
                "cost": float(getattr(solution, "cost", np.nan)),
                "posterior_growth_spearman_vs_saved_primary": posterior_rho,
                "cell_source_mass_spearman_vs_saved_primary": row_rho,
                "type_joint_mass_l1_vs_saved_primary": float(
                    np.abs(joint - ref_joint).sum()
                ),
                "type_source_mass_max_abs_vs_saved_primary": float(
                    np.abs(source_mass - ref_source_mass).max()
                ),
                "type_conditional_js_vs_saved_primary": float(
                    np.average(js_values, weights=js_weights)
                ),
            }
        )

detail = pd.DataFrame(rows)
assert detail["converged"].all()
detail.to_csv(RESULTS / "native_moscot_parameter_sensitivity.csv", index=False)

summary = (
    detail.groupby(
        ["run", "growth_scaling", "epsilon", "tau_a", "tau_b"],
        observed=True,
    )
    .agg(
        all_pairs_converged=("converged", "all"),
        min_posterior_growth_spearman=(
            "posterior_growth_spearman_vs_saved_primary", "min"
        ),
        min_cell_source_mass_spearman=(
            "cell_source_mass_spearman_vs_saved_primary", "min"
        ),
        median_type_joint_mass_l1=("type_joint_mass_l1_vs_saved_primary", "median"),
        max_type_joint_mass_l1=("type_joint_mass_l1_vs_saved_primary", "max"),
        max_type_source_mass_abs=(
            "type_source_mass_max_abs_vs_saved_primary", "max"
        ),
        median_type_conditional_js=(
            "type_conditional_js_vs_saved_primary", "median"
        ),
        max_type_conditional_js=("type_conditional_js_vs_saved_primary", "max"),
    )
    .reset_index()
)
summary.to_csv(
    RESULTS / "native_moscot_parameter_sensitivity_summary.csv", index=False
)
print(summary.to_string(index=False))
print("Saved sensitivity results to", RESULTS)

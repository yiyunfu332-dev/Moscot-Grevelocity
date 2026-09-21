"""Direction-aware GraphVelo versus moscot validation.

Separates absolute directional agreement (cosine relative to zero) from performance
relative to a within-time, within-fish shuffled null. A negative cosine can therefore
never be labeled positive agreement merely because the shuffled null is more negative.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_graphvelo_directionality")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT = (
    ROOT
    / "blood_growth_graphvelo_moscot_v2/checkpoints/"
    / "zebrafish_growth_graphvelo_moscot_v2.h5ad"
)
RESULTS = ROOT / "blood_growth_graphvelo_moscot_v2/results"
COUPLINGS = RESULTS / "couplings_no_topk"
N_PERMUTATIONS = 2_000
N_BOOTSTRAPS = 5_000
SEED = 322


def cosine_rows(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(
        np.sum(left * right, axis=1),
        denominator,
        out=np.full(len(left), np.nan),
        where=denominator > 0,
    )


def equal_fish_mean(values: np.ndarray, fish: np.ndarray) -> float:
    frame = pd.DataFrame({"value": values, "fish": fish}).dropna()
    return float(frame.groupby("fish", observed=True)["value"].mean().mean())


rng = np.random.default_rng(SEED)
adata = ad.read_h5ad(CHECKPOINT, backed="r")
hpf = pd.to_numeric(adata.obs["hpf"]).to_numpy(float)
fish_all = adata.obs["fish_id"].astype(str).to_numpy()
x = np.asarray(adata.obsm["X_moscot"], dtype=float)
velocity = np.asarray(adata.obsm["gv_pca"], dtype=float)[:, : x.shape[1]]
hpf_values = np.sort(np.unique(hpf))

rows = []
for t0, t1 in zip(hpf_values[:-1], hpf_values[1:]):
    source_idx = np.flatnonzero(hpf == t0)
    target_idx = np.flatnonzero(hpf == t1)
    name = f"{t0:g}_to_{t1:g}_hpf"
    coupling = sp.load_npz(COUPLINGS / f"{name}_M1_cond.npz").tocsr()
    assert coupling.shape == (len(source_idx), len(target_idx))

    barycenter = np.asarray(coupling @ x[target_idx])
    velocity_moscot = (barycenter - x[source_idx]) / ((t1 - t0) / 24.0)
    velocity_graphvelo = velocity[source_idx]
    fish = fish_all[source_idx]
    cosine = cosine_rows(velocity_graphvelo, velocity_moscot)
    valid = np.isfinite(cosine)
    assert valid.any()
    cosine = cosine[valid]
    velocity_graphvelo = velocity_graphvelo[valid]
    velocity_moscot = velocity_moscot[valid]
    fish = fish[valid]

    fish_means = (
        pd.DataFrame({"cosine": cosine, "fish": fish})
        .groupby("fish", observed=True)["cosine"]
        .mean()
        .to_numpy()
    )
    observed = float(fish_means.mean())
    bootstrap = rng.choice(
        fish_means, size=(N_BOOTSTRAPS, len(fish_means)), replace=True
    ).mean(axis=1)

    null = np.empty(N_PERMUTATIONS, dtype=float)
    fish_indices = [np.flatnonzero(fish == label) for label in np.unique(fish)]
    for iteration in range(N_PERMUTATIONS):
        permuted = velocity_graphvelo.copy()
        for indices in fish_indices:
            permuted[indices] = velocity_graphvelo[rng.permutation(indices)]
        null_cosine = cosine_rows(permuted, velocity_moscot)
        null[iteration] = equal_fish_mean(null_cosine, fish)

    p_greater = float((1 + np.sum(null >= observed)) / (N_PERMUTATIONS + 1))
    p_less = float((1 + np.sum(null <= observed)) / (N_PERMUTATIONS + 1))
    rows.append(
        {
            "source_hpf": t0,
            "target_hpf": t1,
            "n_valid_cells": len(cosine),
            "n_fish": len(fish_means),
            "cell_weighted_mean_cosine": float(cosine.mean()),
            "fish_weighted_mean_cosine": observed,
            "fish_bootstrap_ci_low": float(np.quantile(bootstrap, 0.025)),
            "fish_bootstrap_ci_high": float(np.quantile(bootstrap, 0.975)),
            "within_fish_permuted_mean_cosine": float(null.mean()),
            "permutation_p_greater": p_greater,
            "permutation_p_less": p_less,
            "permutation_p_two_sided": min(1.0, 2.0 * min(p_greater, p_less)),
        }
    )

result = pd.DataFrame(rows)
for column in ["permutation_p_greater", "permutation_p_less", "permutation_p_two_sided"]:
    result[column.replace("_p_", "_q_")] = multipletests(
        result[column], method="fdr_bh"
    )[1]

positive = (
    result["fish_weighted_mean_cosine"].gt(0)
    & result["fish_bootstrap_ci_low"].gt(0)
    & result["permutation_q_greater"].lt(0.05)
)
negative = (
    result["fish_weighted_mean_cosine"].lt(0)
    & result["fish_bootstrap_ci_high"].lt(0)
    & result["permutation_q_less"].lt(0.05)
)
result["directional_evidence"] = np.select(
    [positive, negative],
    ["positive_concordance", "negative_discordance"],
    default="insufficient_or_mixed",
)
result["interpretation"] = np.select(
    [
        positive,
        negative,
        result["permutation_q_greater"].lt(0.05)
        & result["fish_weighted_mean_cosine"].le(0),
    ],
    [
        "positive direction, fish CI above zero, and better than within-fish shuffle",
        "negative direction, fish CI below zero, and worse than within-fish shuffle",
        "better than shuffled null but not positive absolute directional agreement",
    ],
    default="no robust direction-level evidence",
)

result.to_csv(RESULTS / "moscot_graphvelo_directional_evidence.csv", index=False)
print(result.to_string(index=False))
print("Saved", RESULTS / "moscot_graphvelo_directional_evidence.csv")
adata.file.close()

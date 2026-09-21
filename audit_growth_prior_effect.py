"""Quantify how much the proliferation/apoptosis prior changes native moscot."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def _tvd(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    if a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    return float(0.5 * np.abs(a / a.sum() - b / b.sum()).sum())


def _sparse_tvd(a: sparse.spmatrix, b: sparse.spmatrix) -> float:
    a = a.tocsr().astype(float)
    b = b.tocsr().astype(float)
    if a.sum() <= 0 or b.sum() <= 0:
        return np.nan
    difference = a / a.sum() - b / b.sum()
    return float(0.5 * np.abs(difference.data).sum())


def audit_growth_prior_effect(
    checkpoint_h5ad: str | Path,
    coupling_index_csv: str | Path,
    results_dir: str | Path,
    *,
    n_fish_bootstraps: int = 5000,
    dominant_tvd_threshold: float = 0.05,
    seed: int = 0,
) -> dict[str, pd.DataFrame]:
    """Measure M0→M1 changes and bootstrap fish-level descriptive effects.

    The bootstrap resamples per-fish summaries; it does not refit moscot and is
    therefore an uncertainty summary of represented fish, not full model-fit
    uncertainty.
    """

    results_dir = Path(results_dir)
    checkpoint = ad.read_h5ad(Path(checkpoint_h5ad), backed="r")
    coupling_index = pd.read_csv(coupling_index_csv)
    hpf = pd.to_numeric(checkpoint.obs["hpf"], errors="raise").to_numpy(float)
    fish_ids = checkpoint.obs["fish_id"].astype(str).to_numpy()
    rng = np.random.default_rng(seed)

    fish_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    pairs = coupling_index[["source_hpf", "target_hpf"]].drop_duplicates()
    for pair in pairs.itertuples(index=False):
        t0, t1 = float(pair.source_hpf), float(pair.target_hpf)
        records = coupling_index.loc[
            coupling_index["source_hpf"].eq(t0)
            & coupling_index["target_hpf"].eq(t1)
        ].set_index("matrix")
        required = {"M0_raw", "M1_raw", "M0_cond", "M1_cond"}
        if not required.issubset(records.index):
            raise KeyError(f"Missing coupling matrices for {t0:g}→{t1:g} hpf")
        m0_raw = sparse.load_npz(records.loc["M0_raw", "path"]).tocsr()
        m1_raw = sparse.load_npz(records.loc["M1_raw", "path"]).tocsr()
        m0_cond = sparse.load_npz(records.loc["M0_cond", "path"]).tocsr()
        m1_cond = sparse.load_npz(records.loc["M1_cond", "path"]).tocsr()
        source_global = np.flatnonzero(np.isclose(hpf, t0))
        if len(source_global) != m0_raw.shape[0]:
            raise ValueError(
                f"Source-cell order/shape mismatch for {t0:g}→{t1:g}: "
                f"{len(source_global)} metadata cells vs {m0_raw.shape[0]} rows"
            )
        source_fish = fish_ids[source_global]
        source_mass0 = np.asarray(m0_raw.sum(axis=1)).ravel()
        source_mass1 = np.asarray(m1_raw.sum(axis=1)).ravel()
        row_fate_tvd = 0.5 * np.asarray(
            np.abs(m0_cond - m1_cond).sum(axis=1)
        ).ravel()

        for fish_id in np.unique(source_fish):
            local = np.flatnonzero(source_fish == fish_id)
            fish_rows.append({
                "source_hpf": t0,
                "target_hpf": t1,
                "fish_id": fish_id,
                "n_source_cells": len(local),
                "source_mass_tvd_M0_M1": _tvd(
                    source_mass0[local], source_mass1[local]
                ),
                "mean_cell_conditional_fate_tvd_M0_M1": float(
                    np.mean(row_fate_tvd[local])
                ),
            })

        current_fish = pd.DataFrame(fish_rows).loc[
            lambda x: x["source_hpf"].eq(t0) & x["target_hpf"].eq(t1)
        ]
        fish_source_tvd = current_fish["source_mass_tvd_M0_M1"].to_numpy(float)
        fish_fate_tvd = current_fish[
            "mean_cell_conditional_fate_tvd_M0_M1"
        ].to_numpy(float)
        boot_index = rng.integers(
            0, len(current_fish), size=(n_fish_bootstraps, len(current_fish))
        )
        boot_source = fish_source_tvd[boot_index].mean(axis=1)
        boot_fate = fish_fate_tvd[boot_index].mean(axis=1)
        aggregate_source_tvd = _tvd(source_mass0, source_mass1)
        joint_tvd = _sparse_tvd(m0_raw, m1_raw)
        ci_low, ci_high = np.quantile(boot_source, [0.025, 0.975])
        summary_rows.append({
            "source_hpf": t0,
            "target_hpf": t1,
            "n_source_cells": len(source_global),
            "n_fish": len(current_fish),
            "aggregate_source_mass_tvd_M0_M1": aggregate_source_tvd,
            "joint_coupling_tvd_M0_M1": joint_tvd,
            "mean_fish_source_mass_tvd_M0_M1": float(fish_source_tvd.mean()),
            "fish_bootstrap_source_tvd_ci_low": float(ci_low),
            "fish_bootstrap_source_tvd_ci_high": float(ci_high),
            "mean_fish_conditional_fate_tvd_M0_M1": float(fish_fate_tvd.mean()),
            "fish_bootstrap_fate_tvd_ci_low": float(np.quantile(boot_fate, 0.025)),
            "fish_bootstrap_fate_tvd_ci_high": float(np.quantile(boot_fate, 0.975)),
            "dominant_growth_effect_threshold_tvd": dominant_tvd_threshold,
            "supports_dominant_growth_redistribution": bool(
                aggregate_source_tvd >= dominant_tvd_threshold
                and ci_low >= dominant_tvd_threshold
            ),
        })

    fish_effect = pd.DataFrame(fish_rows)
    summary = pd.DataFrame(summary_rows)
    summary["interpretation"] = np.where(
        summary["supports_dominant_growth_redistribution"],
        "growth prior produces >=5% source-mass redistribution",
        "no evidence that growth prior dominates transport",
    )
    fish_effect.to_csv(results_dir / "growth_prior_effect_by_fish.csv", index=False)
    summary.to_csv(results_dir / "growth_prior_effect_summary.csv", index=False)
    checkpoint.file.close()
    return {"by_fish": fish_effect, "summary": summary}


if __name__ == "__main__":
    project = Path(__file__).resolve().parent
    results = project / "blood_growth_graphvelo_moscot_v2" / "results"
    outputs = audit_growth_prior_effect(
        project / "blood_growth_graphvelo_moscot_v2" / "checkpoints"
        / "zebrafish_growth_graphvelo_moscot_v2.h5ad",
        results / "coupling_file_index.csv",
        results,
    )
    print(outputs["summary"].to_string(index=False))

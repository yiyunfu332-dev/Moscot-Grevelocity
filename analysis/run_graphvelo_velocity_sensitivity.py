"""Retrain GraphVelo across stricter kinetic-gene QC thresholds.

This tests whether the saved velocity direction and its moscot concordance depend on
the permissive gamma-R2 cutoff used in the primary run.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_graphvelo_sensitivity")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_graphvelo_sensitivity")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp/xdg_graphvelo_sensitivity")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from graphvelo.graph_velocity import GraphVelo


ROOT = Path(__file__).resolve().parent.parent
CHECKPOINT = (
    ROOT
    / "blood_growth_graphvelo_moscot_v2/checkpoints/"
    / "zebrafish_growth_graphvelo_moscot_v2.h5ad"
)
RESULTS = ROOT / "blood_growth_graphvelo_moscot_v2/results"
COUPLINGS = RESULTS / "couplings_no_topk"
QC_FILE = RESULTS / "velocity_gene_qc.csv"
THRESHOLDS = [0.01, 0.05, 0.10, 0.20]
N_PCS = 30


def row_cosine(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    denominator = np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1)
    return np.divide(
        np.sum(left * right, axis=1),
        denominator,
        out=np.full(len(left), np.nan),
        where=denominator > 0,
    )


adata = ad.read_h5ad(CHECKPOINT)
assert {"M_s", "velocity_S"}.issubset(adata.layers.keys())
assert {"X_pca", "X_moscot", "gv_pca"}.issubset(adata.obsm.keys())
assert "indices" in adata.uns["neighbors"]

qc = pd.read_csv(QC_FILE, index_col=0).reindex(adata.var_names)
assert qc.index.equals(adata.var_names) and qc.notna().any(axis=1).all()
base_velocity = np.asarray(adata.obsm["gv_pca"], dtype=float)[:, :N_PCS]
hpf = pd.to_numeric(adata.obs["hpf"]).to_numpy(float)
hpf_values = np.sort(np.unique(hpf))
hpf_pairs = list(zip(hpf_values[:-1], hpf_values[1:]))

moscot_velocity = {}
for t0, t1 in hpf_pairs:
    name = f"{t0:g}_to_{t1:g}_hpf"
    source_idx = np.flatnonzero(hpf == t0)
    target_idx = np.flatnonzero(hpf == t1)
    coupling = sp.load_npz(COUPLINGS / f"{name}_M1_cond.npz").tocsr()
    assert coupling.shape == (len(source_idx), len(target_idx))
    source_x = np.asarray(adata.obsm["X_moscot"])[source_idx, :N_PCS]
    target_x = np.asarray(adata.obsm["X_moscot"])[target_idx, :N_PCS]
    barycenter = np.asarray(coupling @ target_x)
    moscot_velocity[(t0, t1)] = (
        source_idx,
        (barycenter - source_x) / ((t1 - t0) / 24.0),
    )

summary_rows = []
interval_rows = []
gene_rows = []
for threshold in THRESHOLDS:
    valid = (
        qc["gamma"].gt(0)
        & qc["gamma_r2"].ge(threshold)
        & qc["velocity_finite_fraction"].ge(0.95)
        & qc["spliced_detected_cells"].ge(25)
        & qc["unspliced_detected_cells"].ge(25)
    )
    genes = qc.index[valid].astype(str).tolist()
    assert len(genes) >= 50
    for gene in genes:
        gene_rows.append({"gamma_r2_threshold": threshold, "gene": gene})

    if threshold == THRESHOLDS[0]:
        projected = base_velocity.copy()
        fit_cosine_median = np.nan
        source = "saved_primary"
    else:
        print(f"Training GraphVelo with gamma_r2 >= {threshold}: {len(genes)} genes", flush=True)
        graphvelo = GraphVelo(
            adata,
            gene_subset=genes,
            xkey="M_s",
            vkey="velocity_S",
        )
        graphvelo.train()
        projected = np.asarray(
            graphvelo.project_velocity(adata.obsm["X_pca"]), dtype=float
        )[:, :N_PCS]
        fit_cosine_median = float(
            np.nanmedian(
                row_cosine(
                    np.asarray(graphvelo.V),
                    np.asarray(graphvelo.project_velocity(graphvelo.X)),
                )
            )
        )
        source = "strict_qc_refit"

    stability = row_cosine(base_velocity, projected)
    summary_rows.append(
        {
            "gamma_r2_threshold": threshold,
            "n_velocity_genes": len(genes),
            "velocity_source": source,
            "graphvelo_fit_cosine_median": fit_cosine_median,
            "cell_cosine_vs_primary_median": float(np.nanmedian(stability)),
            "cell_cosine_vs_primary_q05": float(np.nanquantile(stability, 0.05)),
            "cell_fraction_positive_vs_primary": float(np.nanmean(stability > 0)),
            "cell_fraction_cosine_ge_0.5_vs_primary": float(np.nanmean(stability >= 0.5)),
        }
    )

    for (t0, t1), (source_idx, velocity_moscot) in moscot_velocity.items():
        cosine = row_cosine(projected[source_idx], velocity_moscot)
        interval_rows.append(
            {
                "gamma_r2_threshold": threshold,
                "n_velocity_genes": len(genes),
                "source_hpf": t0,
                "target_hpf": t1,
                "n_valid_cells": int(np.isfinite(cosine).sum()),
                "mean_moscot_graphvelo_cosine": float(np.nanmean(cosine)),
                "median_moscot_graphvelo_cosine": float(np.nanmedian(cosine)),
                "fraction_positive_moscot_graphvelo_cosine": float(
                    np.nanmean(cosine > 0)
                ),
            }
        )

summary = pd.DataFrame(summary_rows)
interval = pd.DataFrame(interval_rows)
genes = pd.DataFrame(gene_rows)
summary.to_csv(RESULTS / "graphvelo_gene_qc_sensitivity_summary.csv", index=False)
interval.to_csv(RESULTS / "graphvelo_gene_qc_sensitivity_by_interval.csv", index=False)
genes.to_csv(RESULTS / "graphvelo_gene_qc_sensitivity_gene_sets.csv", index=False)
print(summary.to_string(index=False))
print("Saved GraphVelo sensitivity results to", RESULTS)

"""Fish-aware pseudobulk validation of zebrafish cell-type markers.

The notebook's Scanpy rank_genes_groups table is useful for exploration, but its
cell-level p-values treat cells from the same fish as independent observations.
This module forms within-fish/time cell-type-versus-rest pseudobulk contrasts
and performs inference across fish instead.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t as student_t
from statsmodels.stats.multitest import multipletests


def _log2_cpm(counts: np.ndarray) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.float64).ravel()
    library_size = counts.sum()
    if library_size <= 0:
        raise ValueError("Encountered an empty pseudobulk library.")
    return np.log2(1.0 + counts * (1_000_000.0 / library_size))


def run_fish_aware_marker_qc(
    full_expression_h5ad: str | Path,
    checkpoint_h5ad: str | Path | ad.AnnData,
    results_dir: str | Path,
    *,
    min_target_cells: int = 5,
    min_rest_cells: int = 20,
    min_fish: int = 3,
    candidate_log2fc: float = 1.0,
    candidate_fdr: float = 0.05,
) -> dict[str, pd.DataFrame]:
    """Run paired pseudobulk marker QC with fish as the replicate.

    Each contrast compares one cell type with all other observed cell types in
    the same fish and hpf. If a fish contributes more than one hpf, its contrast
    is averaged before testing so that it still contributes one replicate.
    """

    full_expression_h5ad = Path(full_expression_h5ad)
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    full = ad.read_h5ad(full_expression_h5ad)
    owns_checkpoint = not isinstance(checkpoint_h5ad, ad.AnnData)
    checkpoint = (
        ad.read_h5ad(Path(checkpoint_h5ad), backed="r")
        if owns_checkpoint
        else checkpoint_h5ad
    )

    required_obs = {"fish_id", "hpf", "official_fine_cell_type"}
    missing_obs = required_obs.difference(checkpoint.obs.columns)
    if missing_obs:
        raise KeyError(f"Checkpoint is missing metadata: {sorted(missing_obs)}")
    if "counts" not in full.layers:
        raise KeyError("Full-expression AnnData has no raw 'counts' layer.")

    full_index = full.obs_names.get_indexer(checkpoint.obs_names)
    if np.any(full_index < 0):
        missing = checkpoint.obs_names[full_index < 0].tolist()[:5]
        raise KeyError(f"Checkpoint cells absent from full expression data: {missing}")

    counts = full.layers["counts"][full_index]
    counts = sparse.csr_matrix(counts)
    obs = checkpoint.obs[
        ["fish_id", "hpf", "official_fine_cell_type"]
    ].copy()
    obs["fish_id"] = obs["fish_id"].astype(str)
    obs["hpf"] = pd.to_numeric(obs["hpf"], errors="raise").astype(float)
    obs["official_fine_cell_type"] = obs["official_fine_cell_type"].astype(str)
    gene_names = full.var_names.astype(str).to_numpy()

    # One target-vs-rest contrast per fish/time/type.
    contrast_rows: list[dict[str, object]] = []
    contrast_vectors: dict[str, list[tuple[str, float, np.ndarray]]] = {}
    grouped_time = obs.groupby(["fish_id", "hpf"], observed=True).indices
    for (fish_id, hpf), time_indices_raw in grouped_time.items():
        time_indices = np.asarray(time_indices_raw, dtype=int)
        time_counts = np.asarray(counts[time_indices].sum(axis=0)).ravel()
        time_total_cells = len(time_indices)
        time_labels = obs.iloc[time_indices]["official_fine_cell_type"].to_numpy()
        for cell_type in np.unique(time_labels):
            target_indices = time_indices[time_labels == cell_type]
            n_target = len(target_indices)
            n_rest = time_total_cells - n_target
            if n_target < min_target_cells or n_rest < min_rest_cells:
                continue
            target_counts = np.asarray(counts[target_indices].sum(axis=0)).ravel()
            rest_counts = time_counts - target_counts
            contrast = _log2_cpm(target_counts) - _log2_cpm(rest_counts)
            contrast_vectors.setdefault(cell_type, []).append(
                (str(fish_id), float(hpf), contrast)
            )
            contrast_rows.append({
                "fish_id": str(fish_id),
                "hpf": float(hpf),
                "official_fine_cell_type": cell_type,
                "target_cells": n_target,
                "rest_cells": n_rest,
            })

    all_statistics: list[pd.DataFrame] = []
    summary_rows: list[dict[str, object]] = []
    for cell_type, entries in sorted(contrast_vectors.items()):
        # Average repeated time contrasts within fish, if they occur.
        fish_to_vectors: dict[str, list[np.ndarray]] = {}
        fish_to_hpf: dict[str, list[float]] = {}
        for fish_id, hpf, vector in entries:
            fish_to_vectors.setdefault(fish_id, []).append(vector)
            fish_to_hpf.setdefault(fish_id, []).append(hpf)
        fish_ids = sorted(fish_to_vectors)
        if len(fish_ids) < min_fish:
            summary_rows.append({
                "official_fine_cell_type": cell_type,
                "n_independent_fish": len(fish_ids),
                "n_fish_time_contrasts": len(entries),
                "status": "insufficient_fish",
                "significant_positive_markers": 0,
            })
            continue

        fish_matrix = np.vstack([
            np.mean(fish_to_vectors[fish_id], axis=0) for fish_id in fish_ids
        ])
        n_fish_used = fish_matrix.shape[0]
        mean_difference = fish_matrix.mean(axis=0)
        sd_difference = fish_matrix.std(axis=0, ddof=1)
        standard_error = sd_difference / np.sqrt(n_fish_used)
        t_statistic = np.divide(
            mean_difference,
            standard_error,
            out=np.zeros_like(mean_difference),
            where=standard_error > 0,
        )
        zero_se_nonzero = (standard_error == 0) & (mean_difference != 0)
        t_statistic[zero_se_nonzero] = np.sign(mean_difference[zero_se_nonzero]) * np.inf
        p_value = 2.0 * student_t.sf(np.abs(t_statistic), df=n_fish_used - 1)
        p_value[(standard_error == 0) & (mean_difference == 0)] = 1.0
        fdr = multipletests(p_value, method="fdr_bh")[1]

        stat = pd.DataFrame({
            "official_fine_cell_type": cell_type,
            "gene": gene_names,
            "n_independent_fish": n_fish_used,
            "mean_within_fish_log2cpm_difference": mean_difference,
            "sd_within_fish_log2cpm_difference": sd_difference,
            "t_statistic": t_statistic,
            "p_value": p_value,
            "fdr_bh": fdr,
        })
        all_statistics.append(stat)
        n_positive = int(
            ((fdr < candidate_fdr) & (mean_difference > candidate_log2fc)).sum()
        )
        summary_rows.append({
            "official_fine_cell_type": cell_type,
            "n_independent_fish": n_fish_used,
            "n_fish_time_contrasts": len(entries),
            "status": "tested",
            "significant_positive_markers": n_positive,
        })

    if not all_statistics:
        raise RuntimeError("No cell type had enough independent fish for marker QC.")

    statistics = pd.concat(all_statistics, ignore_index=True)
    candidates = statistics.loc[
        statistics["fdr_bh"].lt(candidate_fdr)
        & statistics["mean_within_fish_log2cpm_difference"].gt(candidate_log2fc)
    ].copy()
    candidates = candidates.sort_values(
        ["official_fine_cell_type", "mean_within_fish_log2cpm_difference"],
        ascending=[True, False],
    )
    top_candidates = candidates.groupby(
        "official_fine_cell_type", observed=True, group_keys=False
    ).head(20)
    summary = pd.DataFrame(summary_rows).sort_values("official_fine_cell_type")
    contrasts = pd.DataFrame(contrast_rows).sort_values(
        ["official_fine_cell_type", "hpf", "fish_id"]
    )

    statistics.to_csv(
        results_dir / "fish_aware_celltype_marker_statistics.csv.gz",
        index=False,
        compression="gzip",
    )
    candidates.to_csv(results_dir / "fish_aware_celltype_marker_candidates.csv", index=False)
    top_candidates.to_csv(
        results_dir / "fish_aware_celltype_marker_top20.csv", index=False
    )
    summary.to_csv(results_dir / "fish_aware_celltype_marker_summary.csv", index=False)
    contrasts.to_csv(
        results_dir / "fish_aware_celltype_marker_contrast_support.csv", index=False
    )

    if owns_checkpoint:
        checkpoint.file.close()
    return {
        "summary": summary,
        "top_candidates": top_candidates,
        "candidates": candidates,
        "contrast_support": contrasts,
    }


if __name__ == "__main__":
    project = Path(__file__).resolve().parent
    outputs = run_fish_aware_marker_qc(
        project / "zf_atlas_hematopoetic_endothelial_v4_release.h5ad",
        project
        / "blood_growth_graphvelo_moscot_v2"
        / "checkpoints"
        / "zebrafish_growth_graphvelo_moscot_v2.h5ad",
        project / "blood_growth_graphvelo_moscot_v2" / "results",
    )
    print(outputs["summary"].to_string(index=False))
    print("\nTop fish-aware marker candidates:")
    print(outputs["top_candidates"].to_string(index=False))

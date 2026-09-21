"""Audit time- and cell-type support for transport summaries."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def audit_transport_support(
    metadata: ad.AnnData | pd.DataFrame,
    type_transitions: pd.DataFrame,
    results_dir: str | Path,
    *,
    min_time_cells: int = 50,
    min_type_cells: int = 20,
    min_fish: int = 3,
) -> dict[str, pd.DataFrame]:
    """Annotate transport rows using explicit cell and fish support gates."""

    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    obs = metadata.obs.copy() if isinstance(metadata, ad.AnnData) else metadata.copy()
    required = {"hpf", "fish_id", "official_fine_cell_type"}
    missing = required.difference(obs.columns)
    if missing:
        raise KeyError(f"Missing support metadata: {sorted(missing)}")
    obs["hpf"] = pd.to_numeric(obs["hpf"], errors="raise").astype(float)
    obs["fish_id"] = obs["fish_id"].astype(str)
    obs["official_fine_cell_type"] = obs["official_fine_cell_type"].astype(str)

    time_support = (
        obs.groupby("hpf", observed=True)
        .agg(total_cells=("fish_id", "size"), total_fish=("fish_id", "nunique"))
        .reset_index()
        .sort_values("hpf")
    )
    time_support["passes_time_support"] = (
        time_support["total_cells"].ge(min_time_cells)
        & time_support["total_fish"].ge(min_fish)
    )

    type_support = (
        obs.groupby(["hpf", "official_fine_cell_type"], observed=True)
        .agg(type_cells=("fish_id", "size"), type_fish=("fish_id", "nunique"))
        .reset_index()
        .sort_values(["hpf", "official_fine_cell_type"])
    )
    type_support["passes_type_support"] = (
        type_support["type_cells"].ge(min_type_cells)
        & type_support["type_fish"].ge(min_fish)
    )

    annotated = type_transitions.copy()
    source_support = type_support.rename(columns={
        "hpf": "source_hpf",
        "official_fine_cell_type": "source",
        "type_cells": "source_type_cells",
        "type_fish": "source_type_fish",
        "passes_type_support": "source_type_supported",
    })
    target_support = type_support.rename(columns={
        "hpf": "target_hpf",
        "official_fine_cell_type": "target",
        "type_cells": "target_type_cells",
        "type_fish": "target_type_fish",
        "passes_type_support": "target_type_supported",
    })
    source_time = time_support.rename(columns={
        "hpf": "source_hpf",
        "total_cells": "source_time_cells",
        "total_fish": "source_time_fish",
        "passes_time_support": "source_time_supported",
    })
    target_time = time_support.rename(columns={
        "hpf": "target_hpf",
        "total_cells": "target_time_cells",
        "total_fish": "target_time_fish",
        "passes_time_support": "target_time_supported",
    })
    annotated = annotated.merge(source_support, on=["source_hpf", "source"], how="left")
    annotated = annotated.merge(target_support, on=["target_hpf", "target"], how="left")
    annotated = annotated.merge(source_time, on="source_hpf", how="left")
    annotated = annotated.merge(target_time, on="target_hpf", how="left")
    bool_columns = [
        "source_type_supported", "target_type_supported",
        "source_time_supported", "target_time_supported",
    ]
    for column in bool_columns:
        annotated[column] = annotated[column].fillna(False).astype(bool)
    annotated["supported_for_type_level_claim"] = annotated[bool_columns].all(axis=1)

    def reason(row: pd.Series) -> str:
        failures = []
        if not row["source_time_supported"]:
            failures.append("source_time")
        if not row["target_time_supported"]:
            failures.append("target_time")
        if not row["source_type_supported"]:
            failures.append("source_type")
        if not row["target_type_supported"]:
            failures.append("target_type")
        return "supported" if not failures else "+".join(failures)

    annotated["support_status"] = annotated.apply(reason, axis=1)
    supported = annotated.loc[annotated["supported_for_type_level_claim"]].copy()

    pair_support = (
        annotated.groupby(["source_hpf", "target_hpf"], observed=True)
        .agg(
            total_rows=("model", "size"),
            supported_rows=("supported_for_type_level_claim", "sum"),
            source_time_cells=("source_time_cells", "first"),
            source_time_fish=("source_time_fish", "first"),
            target_time_cells=("target_time_cells", "first"),
            target_time_fish=("target_time_fish", "first"),
        )
        .reset_index()
    )
    pair_support["supported_row_fraction"] = (
        pair_support["supported_rows"] / pair_support["total_rows"]
    )
    pair_support["has_any_supported_type_transition"] = pair_support[
        "supported_rows"
    ].gt(0)

    time_support.to_csv(results_dir / "transport_time_support.csv", index=False)
    type_support.to_csv(results_dir / "transport_cell_type_time_support.csv", index=False)
    pair_support.to_csv(results_dir / "transport_pair_support_summary.csv", index=False)
    annotated.to_csv(
        results_dir / "moscot_M0_M1_type_transitions_support_annotated.csv",
        index=False,
    )
    supported.to_csv(
        results_dir / "moscot_M0_M1_type_transitions_supported_only.csv",
        index=False,
    )
    return {
        "time_support": time_support,
        "type_support": type_support,
        "pair_support": pair_support,
        "annotated_transitions": annotated,
        "supported_transitions": supported,
    }


if __name__ == "__main__":
    project = Path(__file__).resolve().parent
    checkpoint_path = (
        project / "blood_growth_graphvelo_moscot_v2" / "checkpoints"
        / "zebrafish_growth_graphvelo_moscot_v2.h5ad"
    )
    results = project / "blood_growth_graphvelo_moscot_v2" / "results"
    adata = ad.read_h5ad(checkpoint_path, backed="r")
    transitions = pd.read_csv(results / "moscot_M0_M1_type_transitions.csv")
    outputs = audit_transport_support(adata, transitions, results)
    print(outputs["time_support"].to_string(index=False))
    print("\nPair support:")
    print(outputs["pair_support"].to_string(index=False))
    adata.file.close()

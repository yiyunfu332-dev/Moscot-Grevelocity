# Reproducing the MOSCOT analysis

## Setup and input data

Run all commands from the repository root. Create the environment from `environment.yml`, activate `zebrafish-growth-reproducible`, register its Jupyter kernel, and run `python validate_vae_gpu_environment.py`. The specification records the previously used versions; a fresh solver/install has not been tested in this update. The original environment mixed Conda and user-site packages, so an import check alone does not establish a clean rebuild. See `ENVIRONMENT_NOTES.md`.

Set `ZEBRAFISH_DATA_DIR` to the directory containing `zebrahub_full_velocity.h5ad`, `lineage.h5ad` (lowercase), and `zf_atlas_hematopoetic_endothelial_v4_release.h5ad`. If unset, the repository root is used. Outputs always go to the repository's `blood_growth_graphvelo_moscot_v2/`. Reference tables are resolved within the repository.

The cross-species workflow additionally requires the Popescu and GSE189161 inputs at the exact paths in `DATA_MANIFEST.csv`. `ZEBRAFISH_PROJECT_DIR` optionally changes its project root; normally leave it unset and launch from the repository root. The raw full-expression file respects `ZEBRAFISH_DATA_DIR`. Do not substitute a selected-feature velocity matrix for full expression.

## Run order

1. Execute `zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb` from a fresh kernel. This regenerates baseline results, couplings, aligned source/target IDs and the checkpoint.
2. Run `python analysis/run_native_moscot_sensitivity.py` to generate the complete one-factor sensitivity detail and summary, including tau variation.
3. Run `python analysis/run_graphvelo_velocity_sensitivity.py`.
4. Run `python analysis/reanalyze_graphvelo_directionality.py` to write the authoritative corrected fish-aware directional evidence table. The notebook's earlier exploratory direction summaries must not replace this corrected report.
5. Run `python transport_support_audit.py` and `python audit_growth_prior_effect.py` after the checkpoint/coupling index is available. These also run from the main notebook as applicable.
6. Run `python zebrafish_human_cross_species_validation.py`, or execute its notebook, after supplying the human data. It reads the checkpoint; it does not refit moscot/GraphVelo.

The optional notebook sensitivity block now writes `native_moscot_notebook_sensitivity.csv`. Previously it collided with the standalone script's more complete `native_moscot_parameter_sensitivity.csv`. Run the standalone script for the published 11-configuration audit.

## Snapshot versus new runs

`results_snapshot/2026-09-15/` is an immutable copy of saved source results, verified by the manifest. Working outputs remain gitignored. Do not rewrite the snapshot when rerunning. Source CSV paths are preserved as provenance and may point to the original machine. The newly generated working coupling index references the newly generated matrices and is the input to the growth audit.

The snapshot includes exploratory and corrected outputs for completeness. Prefer `moscot_graphvelo_directional_evidence.csv` over older cell-level/permutation tables, supported-only transitions for type claims, and fish-aware markers over cell-level candidate p-values. The archived raw source-mass numbers and fish bootstrap estimands differ; see `MOSCOT_SUMMARY.md`.

## Known limitations

No new full fit or fresh environment reconstruction was performed for the publication update. The notebook comments note that GraphVelo may need a SciPy sparse `.A` compatibility correction in some installations; the current environment is not proof that an unmodified package install will work. Do not silently change scientific settings to pass an import or data assertion. Store any required package patch/version in a future reproducibility update.

GraphVelo uses expression and PCA shared with transport. Its agreement is complementary computational evidence, not independent validation. The 12 hpf cohort remains too sparse for type-level claims. All fish are cross-sectional. Derived mass is not measured growth and no clonal ground truth is present.

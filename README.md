# Zebrafish hematopoiesis with GraphVelo and moscot

Growth-informed native moscot transport across nine developmental stages, with fish-aware audits and complementary GraphVelo direction analysis.

**Current finding:** the growth prior changes source mass most at 120→240 hpf (TVD 0.0565), while robust positive velocity concordance appears only at 16→19 and 48→72 hpf. The late interval has negative absolute direction agreement. These are model comparisons, not experimentally validated lineage transitions.

## Read the results

- [Complete MOSCOT work summary](docs/MOSCOT_SUMMARY.md): methods, quantitative results, corrections and interpretation limits.
- [Saved results and figures](results_snapshot/2026-09-15/README.md): 60 original result/figure files with SHA-256 provenance, including negative tests.
- [September repair report](docs/PROJECT_REPAIR_REPORT.md): historical execution and scientific corrections.
- [Reproduction guide](docs/REPRODUCIBILITY.md): input files, run order, portability changes and checks.

## Workflow

1. Exact cell-ID joins to official cell type, hpf and fish metadata.
2. Dynamo preprocessing, moments and stochastic RNA velocity; GraphVelo refinement.
3. Full-count expression scoring with 95 reviewed zebrafish proliferation and 46 pro-death genes.
4. Native M0 uniform-input-marginal and M1 growth-informed unbalanced transport, using identical solver settings.
5. Separate raw transport mass from conditional allocation; apply fish/type/time support gates.
6. Compare native M1 displacement with GraphVelo without changing the coupling; audit fish-level uncertainty and parameter sensitivity.
7. Evaluate human expression-program concordance separately.

The entry points are [the main notebook](zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb), [cross-species notebook](zebrafish_human_cross_species_validation.ipynb), and its [Python script](zebrafish_human_cross_species_validation.py). Standalone refit and direction audits are in `analysis/`. The frozen marker sources are in `growth_marker_reference/`.

## Data and execution status

The dated snapshot publishes the saved MOSCOT and cross-species result tables and figures. Multi-gigabyte raw expression data, the 267 MB checkpoint, and approximately 59 MB of coupling files remain outside Git; [the data manifest](docs/DATA_MANIFEST.csv) describes required inputs and generated dependencies. Recompute them with the main notebook before running checkpoint-based audits. Frozen snapshot path fields preserve the historical machine locations for provenance; they are not portable input paths.

The original local main notebook has a recorded 41-cell execution with zero saved errors. Published notebooks have outputs cleared after portability edits; this update does not claim a fresh scientific rerun. A clean rebuild of the pinned environment remains unverified. See [publication verification](docs/PUBLICATION_VERIFICATION.md).

The September 2 bilingual narrative is retained with an explicit historical notice. Its M2/M3 and smaller-panel descriptions are superseded. The current method is native M0/M1, with no Top-K pruning and no velocity reweighting.

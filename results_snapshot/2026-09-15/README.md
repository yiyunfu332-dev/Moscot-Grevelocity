# Frozen results underlying the repaired analysis

Published 24 September 2026; copied byte-for-byte from the local repaired-analysis outputs. The date labels the repaired result set, not a newly executed run.

- `moscot/`: every top-level MOSCOT result table, the full compressed fish-aware marker statistics, and four saved figures.
- `cross_species/results/`: all current cross-species result tables and the validation report.
- `cross_species/figures/`: all current cross-species figures.
- `manifest.csv`: source-relative filenames, published paths, byte sizes and SHA-256 checksums for all 60 files.

Start with `moscot/growth_prior_effect_summary.csv`, `moscot/moscot_graphvelo_directional_evidence.csv`, `moscot/transport_pair_support_summary.csv`, `moscot/native_moscot_parameter_sensitivity_summary.csv`, and `cross_species/results/validation_report.md`.

The legacy permutation/bootstrap tables and exploratory marker rankings remain for provenance. They are superseded for inferential claims by corrected directional evidence and fish-aware marker summaries. Source paths in CSVs are historical provenance; generated coupling matrices and input datasets are excluded from Git. Use the reproduction guide to regenerate working outputs.

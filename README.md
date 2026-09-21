# Zebrafish hematopoiesis with GraphVelo and moscot

This is the core growth and transport workflow from `~/dynamo`.

## Parts

- **Main analysis:** `zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb`. The source workspace saved a complete run of all 41 nonempty code cells without error outputs. The notebook imports the three helper modules at the repository root.
- **Sensitivity analyses:** `analysis/` contains native moscot refits, GraphVelo velocity QC sensitivity, and direction agreement analysis. Run these from the repository root after supplying the checkpoint and saved couplings.
- **Marker and cross-species code:** `build_expanded_growth_marker_panel.py` and `zebrafish_human_cross_species_validation.py`, with reviewed ortholog and pathway reference tables at the root and in `growth_marker_reference/`. The ZFIN ortholog dump and gene list are included for rebuilding the marker panel.
- **Setup and interpretation:** `environment.yml`, `validate_vae_gpu_environment.py`, and `docs/`.

## Data needed for a full rerun

The AnnData checkpoint, full expression matrices, cross-species datasets, transport couplings, and generated results are not included. The main notebook uses the original `~/dynamo` data paths; adjust its `PROJECT_DIR` and the input paths before running it elsewhere. The sensitivity scripts expect the checkpoint under `blood_growth_graphvelo_moscot_v2/checkpoints/` and couplings under that project's `results/` directory. The cross-species workflow expects data under `cross_species_validation/data/`.

The main notebook and environment were checked in the source workspace. The exploratory notebook with a saved error, unrelated dynamo analyses, repair scripts, SLURM launchers, and large generated outputs were deliberately left out. See `docs/PROJECT_REPAIR_REPORT.md` for results and limitations.

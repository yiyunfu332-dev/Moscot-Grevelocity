# Zebrafish growth/transport project repair report

Date: 2026-09-15

## Final execution status

The repaired main notebook was executed from the first cell to the last in a
fresh `vae_gpu` Jupyter kernel. All 41 non-empty code cells ran, no code cell was
skipped, and no error output was present. The executed notebook is:

`zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb`

## Problems resolved

### 1. Zebrafish genes were incorrectly treated as human markers

- Replaced the old 38-proliferation/12-apoptosis cross-species panel with the
  reviewed expanded zebrafish-to-human ortholog table.
- The primary model now uses 95 zebrafish proliferation genes mapping to 91
  unique human orthologs and 46 zebrafish apoptosis genes mapping to 38 unique
  human orthologs.
- Fixed exact zebrafish symbol resolution before case-insensitive fallback.
- Cross-species coverage is 100% for the required panels in zebrafish, Popescu,
  and GSE189161.

### 2. Old and current outputs were mixed

- Archived legacy 38/12 cross-species results.
- Archived old mapping and custom M2/reweighted-transport results.
- Archived the superseded cell-level marker output after introducing the
  explicitly named exploratory result.
- Active results now describe native moscot M0/M1 only.

### 3. Native moscot sensitivity had not been tested

- Refit all eight intervals across growth-scaling, epsilon, and tau settings.
- The primary refit reproduces the saved M1 result essentially exactly.
- Growth scaling 2.5–20 is reasonably stable, but epsilon and tau can strongly
  change the coupling. Conclusions must therefore be restricted to the tested
  neighborhood around epsilon 0.01 and tau 0.9.

### 4. GraphVelo agreement was overinterpreted

- Tightened velocity-gene QC and compared refits with the primary velocity.
- Moderate gamma-R2 thresholds 0.05 and 0.10 are stable; the extreme 0.20
  threshold destabilizes the lower tail and is not recommended.
- Replaced cell-level inference with fish-weighted means, fish bootstrap
  intervals, within-fish/time velocity permutations, and FDR correction.
- Robust positive direction concordance occurs only at 16→19 and 48→72 hpf.
  The other six intervals are labeled insufficient or mixed.

### 5. Marker p-values used cells as independent replicates

- Retained Scanpy Wilcoxon rankings only as exploratory visualization.
- Added paired fish-aware pseudobulk contrasts: each type is compared with the
  other cells from the same fish and hpf, followed by inference across fish.
- All eight modeled types have at least 7 independent fish in this analysis.
- The recovered top markers include expected erythroid, macrophage, microglia,
  and neutrophil genes, supporting the annotation QC.

### 6. Sparse times and cell types appeared equally reliable

- Added explicit support gates requiring at least 50 cells per endpoint time,
  at least 20 cells per source/target type, and at least three fish per type.
- The 12 hpf sample has only 19 cells; consequently 12→14 hpf has no supported
  type-level transition and is no longer plotted as a reliable fate heatmap.
- Later intervals retain only supported rows, reaching full type-level support
  for 120→240 hpf.

### 7. The growth prior's effect size was unclear

- Added M0→M1 total-variation-distance summaries and descriptive fish
  bootstrap intervals.
- From 12→120 hpf, source-mass TVD is 0.0004–0.0184: mostly modest changes.
- At 120→240 hpf, source-mass TVD is 0.0565 with a fish-bootstrap interval of
  about 0.0528–0.0610. This is the only interval meeting the operational 5%
  dominant-redistribution reporting rule.
- The notebook now reports the late effect without claiming that growth
  dominates all transport intervals.

### 8. The Python environment caused import and cache failures

- Pinned the verified package versions in `environment.yml`.
- Added `validate_vae_gpu_environment.py` and verified all core imports.
- Added writable Numba, Matplotlib, and XDG cache directories before imports.
- Kept the required `anndata 0.10.9`/`mudata 0.3.8` compatibility shim.
- Pinned both project notebooks to the `vae_gpu` Jupyter kernel.

### 9. Notebook ordering/state had not been verified

- Executed the complete repaired notebook in a clean kernel.
- Result: 41/41 non-empty code cells executed, 0 skipped, 0 error outputs.

## Cross-species evidence after repair

- Popescu primary broad-lineage proliferation concordance: Spearman rho 0.283,
  donor/fish bootstrap interval 0.077–0.418. Exact ortholog-label specificity
  does not pass the permutation test (p=0.0969).
- Popescu apoptosis concordance is weak: rho 0.083, bootstrap interval
  0.019–0.171.
- GSE189161 exploratory proliferation concordance is stronger (rho 0.734,
  permutation FDR about 0.004) but spans only three harmonized broad lineages.
- The appropriate grade remains **moderate program-level transcriptomic
  support**, not causal, lineage, transport, or velocity validation in humans.

## Limitations that code cannot remove

1. **Missing early cells:** the 12 hpf sample has only 19 modeled cells. No
   statistical method can recover biological replication that was not sampled.
2. **Rare/absent lineages:** a special blood subtype cannot be analyzed reliably
   unless it is present in the official annotation with enough cells and fish.
3. **No direct outcome measurements:** proliferation/apoptosis expression
   programs are proxies. EdU/BrdU, Ki-67/protein, TUNEL/caspase, cell counts, or
   perturbation data would be needed to validate actual division/death rates.
4. **No clonal ground truth:** transport mass is not a measured descendant
   count. Lineage tracing is needed to validate fate couplings.
5. **Method sensitivity:** epsilon and tau materially affect moscot. A single
   parameter setting must not be presented as uniquely correct.
6. **Partial velocity agreement:** GraphVelo supports only two of eight intervals
   under the corrected test and cannot be used as universal validation.
7. **Cross-species scope:** expression-program conservation does not establish
   conserved cell-to-cell couplings, RNA velocity, or absolute growth rates.

## Primary outputs

- `blood_growth_graphvelo_moscot_v2/results/fish_aware_celltype_marker_summary.csv`
- `blood_growth_graphvelo_moscot_v2/results/transport_pair_support_summary.csv`
- `blood_growth_graphvelo_moscot_v2/results/growth_prior_effect_summary.csv`
- `blood_growth_graphvelo_moscot_v2/results/native_moscot_parameter_sensitivity_summary.csv`
- `blood_growth_graphvelo_moscot_v2/results/graphvelo_gene_qc_sensitivity_summary.csv`
- `blood_growth_graphvelo_moscot_v2/results/moscot_graphvelo_directional_evidence.csv`
- `cross_species_validation/results/validation_report.md`
- `cross_species_validation/results/validation_scorecard.csv`

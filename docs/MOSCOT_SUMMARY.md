# MOSCOT work summary — 24 September 2026

## What was built

The current workflow models zebrafish hematopoietic development using native moscot, with a uniform-input-marginal baseline (M0) and proliferation/apoptosis-informed input marginals (M1). Both use the same unbalanced solver settings; M0 is not a balanced-OT control. GraphVelo supplies a separate directional comparison and does not modify either coupling. It shares expression data and PCA with transport, so it is complementary computational evidence, not independent experimental validation.

The analysis aligns cells by exact IDs with official annotations, preserves spliced/unspliced counts, reruns Dynamo preprocessing, computes moments and stochastic velocity, filters velocity genes, and fits GraphVelo. Growth scores use a separate full-count, log-normalized expression object so selected velocity features do not truncate the marker panel. The frozen ZFIN/Seurat/KEGG/Reactome panel contains 95 zebrafish proliferation genes and 46 pro-death genes, mapping to 91 and 38 unique human genes. Human symbols are annotations, never zebrafish expression indices. Anti-apoptotic pathway genes are excluded from the death score.

The final cohort is 5,426 cells, eight modeled types, 36 fish, nine times (12, 14, 16, 19, 24, 48, 72, 120, 240 hpf), and eight adjacent intervals. Four different fish represent each time; these are cross-sectional samples. The 5,441-cell input annotation loses 15 cells through modeling filters. The upstream full velocity atlas and 5 dpf Revelio experiments are distinct analyses, not additional independent replicates of this cohort.

## Results supported by saved tables

| Finding | Saved evidence | Interpretation |
|---|---|---|
| M0 and M1 solved all intervals | 16/16 convergence flags true | Numerical convergence, not biological accuracy |
| Largest growth effect at 120→240 hpf | Source-mass TVD 0.056526; mean within-fish TVD 0.056851, bootstrap interval 0.052802–0.061037 | Only interval passing the operational 0.05 effect-size rule |
| Earlier growth effects are small | Source-mass TVD 0.000424–0.018374 | No general claim that growth dominates transport |
| Conditional fate changes are smaller at 120→240 | Mean fish conditional-fate TVD 0.026367 | Source mass and conditional allocation answer different questions |
| Positive velocity concordance at 16→19 and 48→72 | Fish-weighted cosine 0.157612 and 0.200617; positive bootstrap lower bound and permutation FDR <0.05 | Complementary direction evidence at two intervals |
| Late interval is directionally unresolved | 120→240 cosine −0.091062; bootstrap interval −0.137811 to −0.033078 | Better than shuffled null does not mean positive absolute agreement |
| Early stage is under-sampled | 12 hpf has 19 cells; no supported 12→14 type transition | Do not use an early fate heatmap as a reliable biological result |
| Parameters matter | 11 configurations × 8 intervals in standalone sensitivity results | Growth scaling 2.5–20 is relatively stable; epsilon and tau can materially change transport |
| Marker QC uses biological replicates | Paired fish-aware pseudobulk; all modeled types represented in at least seven fish overall | Cell-level Wilcoxon rankings remain exploratory |

The bootstrap interval for source TVD describes a mean of within-fish normalized effects, not an interval for the aggregate source TVD. Neither bootstrap refits transport, so neither measures full model uncertainty. Support gates require at least 50 cells per endpoint time and at least 20 cells and three fish for each endpoint type. These are reporting screens, not proof of adequate power.

Primary settings: growth scaling 5, epsilon 0.01, tau_a=tau_b=0.9. The standalone sensitivity sweep is one factor at a time; joint parameter interactions remain untested. Moderate GraphVelo gamma-R2 thresholds 0.05 and 0.10 were stable relative to the saved primary fit; 0.20 destabilized the lower tail. Positive directional support and sufficient type support must both be checked before interpreting any particular transition.

## Cross-species work

The validation script reuses the zebrafish checkpoint and restores missing marker genes from full counts. Popescu E-MTAB-7407 contributes 113,063 annotated fetal-liver cells from 14 donors; GSE189161 contributes 58,041 CD34-enriched cells from 26 samples, with nine fetal-liver donors used for replication.

Popescu proliferation concordance is rho=0.283 (donor/fish bootstrap 0.077–0.418), but exact ortholog-label specificity fails the permutation test (p=0.0969). Apoptosis concordance is weak, rho=0.083 (0.019–0.171), with specificity p=0.05195. Exploratory GSE189161 proliferation concordance is rho=0.734 (FDR≈0.004), but covers only three harmonized broad lineages. Independent human HSC/MPP proliferation slope magnitude rankings do not replicate (rho=−0.036, FDR=0.7327); slope directions agree in 60.4% of genes (FDR=0.02929). Apoptosis slope ranks agree more modestly (rho=0.365, FDR=0.04819), with 78.9% same-direction slopes.

Overall: moderate expression-program evidence with important negative tests. These comparisons do not validate zebrafish fate couplings, velocity, absolute division/death rates, or causal mechanisms.

## Work superseded

Earlier notebooks explored hand-assigned/subtype annotations, smaller 38/12-era marker panels, custom M2 velocity reweighting, optional M3 priors, and Top-K transport. They are not the current primary model. September 2 narrative/slides must not be cited as the final method. The September 15 repair report and the tables summarized here supersede them. Historical exploratory notebooks contain saved failures; they are not the reproducible entry point.

## Execution and remaining limitations

The local primary notebook contains 41 nonempty code cells with saved execution counts and no saved error outputs. The September 15 report records a clean-kernel execution. The current cross-species notebook has no saved executions; its script and saved result tables provide the analysis evidence. This September 24 update does not claim a new full fit, a clean-environment rebuild, or new biological experiments.

Next scientific needs are refitted fish-level uncertainty, leakage-free predictive comparisons, cell-cycle and time-gap ablations, and external lineage/proliferation measurements. Gene expression is a growth proxy; transport mass is not a measured descendant count. Four fish per time and sparse early types limit inference.

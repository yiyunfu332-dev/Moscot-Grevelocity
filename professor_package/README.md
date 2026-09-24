# Zebrafish hematopoietic development: MOSCOT, growth programs, and GraphVelo

This project investigated how blood-cell populations change during zebrafish development by combining time-resolved single-cell expression, proliferation/apoptosis programs, optimal transport, and RNA-velocity analysis. This page records the completed computational work, its development, and its limitations. It includes successful runs, partial analyses, negative findings, and superseded methods.

**Project status:** preserved record of work completed through September 2026. The sharing package was assembled on 24 September; no scientific model was rerun for publication. Further experiments discussed during development have not been performed or established as results.

## Start here

1. [Six-page discussion summary](MOSCOT_summary_for_professor.pdf), [GitHub figure index](FIGURE_INDEX.md), and [complete illustrated figure report](MOSCOT_project_figures_for_professor.pdf): six new overview charts/pages drawn from saved results, followed by the historical/current figure appendix. Every page states its evidence status. The overview pages are newly drawn; the underlying fits are not new.
2. [Previously executed primary MOSCOT notebook](executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb): exact copy of the local saved run, **41/41 nonempty code cells executed and zero saved errors**. Original figures, tables, execution counts, kernel metadata and paths are preserved. Use this to review the past run; use the clean notebook at the repository root for a new run.
3. [Current scientific summary](../docs/MOSCOT_SUMMARY.md): corrected conclusions, including weak/negative results.
4. [Browse the figure gallery](gallery.html): download the package and open this file in a browser; GitHub does not render standalone HTML as a website.

## Project explanation

### Research question and motivation

The main question was: **how much does adding proliferation and apoptosis information change inferred developmental transport, and are the resulting directions consistent with RNA velocity?** A separate analysis asked whether the corresponding expression programs have support in human hematopoietic datasets.

Single-cell RNA sequencing measures cells at the time they are collected. Cells sampled at a later stage are not known descendants of the cells sampled earlier. A change in the observed population can reflect differentiation, proliferation, death, migration, or differences in sampling. The analysis therefore distinguishes three quantities:

| Quantity | Question | What this project can provide |
|---|---|---|
| Molecular state change | In which direction is a cell's expression changing? | A local RNA-velocity estimate from Dynamo/GraphVelo |
| Relative transported mass | How does the growth prior redistribute modeled source contributions? | An M0/M1 optimal-transport comparison |
| Conditional allocation | Given a source cell's transported mass, where is it allocated? | A normalized coupling between sampled states |

None of these is a directly measured parent–descendant relationship. The project used existing computational tools to build and audit a biological workflow; it did not develop a new OT solver or demonstrate a new universally accurate fate-prediction method.

### Data, cell identity, and biological replication

The project uses zebrafish developmental single-cell data from ZEBRAHUB, with spliced/unspliced expression for velocity estimation and official annotation, developmental time, and fish identity. A separate full-count expression object supplies genes missing from the selected-feature velocity object. The input files and sizes are recorded in the [data manifest](../docs/DATA_MANIFEST.csv).

The official annotation subset contains **5,441 cells**. The final model contains **5,426 cells**, after excluding the 15-cell thymocyte group under the minimum cell-type size filter. The eight modeled labels are erythroblast, HSPC, macrophage, immature macrophage, neutrophil, microglia, innate lymphoid cell, and non-specific hematopoietic. These names come from the official annotation rather than a manual interpretation of cluster numbers. See the [input label inventory](../results_snapshot/2026-09-15/moscot/official_fine_cell_type_counts.csv).

There are nine modeled stages—**12, 14, 16, 19, 24, 48, 72, 120, and 240 hours post fertilization (hpf)**—and eight adjacent time intervals. Each stage has four different fish, giving **36 fish** in total. This is a cross-sectional design: different animals contribute different time points. The analysis does not follow the same cells or fish longitudinally.

![Cell numbers and sampling support](overview/01_cohort_support.png)

The 12 hpf stage has only 19 modeled cells. Although a coupling can be fitted computationally, no 12→14 hpf type-level transition passes the reporting support gates. A numerical output is not automatically a supported biological result.

### How the analysis was constructed

**1. Audit inputs and align cells.** Velocity, annotation, and full-expression objects are joined using exact cell IDs. Official labels, fish IDs, and experimental time are preserved. Inherited embeddings are discarded before constructing the blood-cell analysis so that unrelated preprocessing does not silently determine the result.

**2. Estimate local expression dynamics.** Dynamo preprocessing, PCA, a neighbor graph, moments, and stochastic RNA velocity provide the initial dynamic representation. Velocity genes undergo numerical and expression-support checks. GraphVelo then refines the velocity representation, which is projected into the PCA coordinates shared with transport. These are estimates of local expression change, not calibrated speeds in biological hours.

**3. Construct species-appropriate growth scores.** A dedicated full-count, log-normalized expression object is used for scoring. Restricting scores to selected velocity features would omit many relevant genes. The reviewed panel contains **95 zebrafish proliferation genes** and **46 zebrafish pro-death genes**, corresponding to 91 and 38 unique human orthologs. Frozen ZFIN orthology and Seurat/Tirosh, KEGG, and Reactome evidence are documented in the [marker reference files](../growth_marker_reference/).

Human gene symbols are used for mapping and interpretation, not directly as zebrafish expression indices. Whole apoptosis pathways are not treated as a uniform death signature because they also contain survival/anti-apoptotic members. Even with these corrections, expression scores remain proxies; no direct division or death measurements calibrate them.

**4. Fit two native MOSCOT models.** Both models connect adjacent experimental stages in the same state space and use the same unbalanced-OT solver settings:

| Model | Input weighting | Role |
|---|---|---|
| **M0** | Uniform input marginals | Baseline without the gene-derived growth prior |
| **M1** | Proliferation/apoptosis-informed input marginals | Primary growth-informed model |

M0 is not a balanced-OT control. Primary settings are epsilon 0.01, tau_a=tau_b=0.9, and growth scaling 5. Experimental time is represented in days for the growth calculation and in hpf for reporting. The current model retains the solver's positive coupling entries without Top-K pruning. GraphVelo does not reweight the coupling.

**5. Separate mass from conditional allocation.** If the raw coupling is $P_{ij}$, its source-row mass is $m_i=\sum_j P_{ij}$. For a nonzero row, the conditional allocation is $C_{ij}=P_{ij}/m_i$. The first quantity describes modeled source contribution; the second describes its distribution among targets. Normalizing every row removes information about its total mass, so these quantities cannot be interpreted interchangeably. “Fate” in older figure titles refers to this inferred conditional allocation, not observed ancestry.

**6. Compare direction without modifying the model.** The conditional coupling gives each source cell a weighted target-state average. The displacement from its current state to that average is compared with GraphVelo using cosine similarity. Positive cosine indicates alignment, zero indicates orthogonality, and negative cosine indicates opposing components. GraphVelo and MOSCOT share expression data and PCA, so this is complementary computational evidence rather than independent experimental validation.

**7. Audit replication, support, and sensitivity.** Fish-level summaries, bootstrap intervals, within-fish/time velocity permutations, multiple-testing correction, gene-QC refits, and native MOSCOT sensitivity analyses test how much confidence to place in the outputs. The [executed primary notebook](executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb) records the main workflow; the [standalone analysis scripts](../analysis/) supply additional corrected audits.

### Main result: growth changes late source mass most

All 16 primary interval fits—eight M0 and eight M1—have true convergence flags in the [saved convergence table](../results_snapshot/2026-09-15/moscot/moscot_convergence_M0_M1.csv). This establishes numerical convergence, not biological accuracy.

The M0→M1 change is quantified with total-variation distance (TVD). For normalized distributions, TVD ranges from zero for identical distributions to one for disjoint distributions. A value of 0.0565 represents roughly 5.65% redistribution of normalized modeled mass; it is not a measurement of 5.65% population growth.

| Interval | Aggregate source-mass TVD, M0→M1 |
|---|---:|
| 12→14 hpf | 0.000428 |
| 14→16 hpf | 0.000424 |
| 16→19 hpf | 0.001205 |
| 19→24 hpf | 0.002147 |
| 24→48 hpf | 0.009821 |
| 48→72 hpf | 0.011103 |
| 72→120 hpf | 0.018374 |
| 120→240 hpf | **0.056526** |

![Growth-prior effects on source mass](overview/02_growth_effect.png)

At 120→240 hpf, the mean within-fish source-mass TVD is 0.056851, with a descriptive bootstrap interval of 0.052802–0.061037. This is a different estimand from the aggregate TVD in the table. The mean fish conditional-allocation TVD is smaller, 0.026367. The effect therefore concerns source weighting and conditional allocation to different degrees.

Only the late interval meets the operational 0.05 source-redistribution reporting rule. Earlier effects are mostly modest. The late interval is also the longest, and the growth prior depends on time duration. These results alone do not establish that late biological proliferation is uniquely strong. Source: [growth-effect table](../results_snapshot/2026-09-15/moscot/growth_prior_effect_summary.csv).

### Main result: velocity concordance is limited

Positive directional evidence requires all three conditions: positive fish-weighted cosine, a positive fish-bootstrap lower bound, and a within-fish permutation FDR below 0.05. Only **16→19 hpf** and **48→72 hpf** meet these criteria.

![Corrected directional concordance](overview/03_directional_agreement.png)

At 120→240 hpf, the fish-weighted cosine is **−0.091062**, with bootstrap interval **−0.137811 to −0.033078**. The observed result is better than a more-negative shuffled null, but that does not make its absolute direction positive. Earlier analyses that emphasized only improvement over the null could overstate agreement; the [corrected directional evidence table](../results_snapshot/2026-09-15/moscot/moscot_graphvelo_directional_evidence.csv) is the authoritative result.

The late disagreement could involve the difference between local velocity and a five-day displacement, curved/cyclic state-space paths, composition differences, or model misspecification. These are possible explanations, not mechanisms established by the completed analysis. This project did not resolve that disagreement experimentally.

### Robustness and statistical corrections

**Parameter sensitivity.** The standalone sweep tested 11 configurations across eight intervals. Growth scaling 2.5–20 was relatively stable compared with the primary fit, whereas epsilon and tau changes could substantially alter the coupling. This was a one-factor-at-a-time analysis, not a complete interaction study. See the [sensitivity overview](overview/04_parameter_sensitivity.png) and [saved parameter summary](../results_snapshot/2026-09-15/moscot/native_moscot_parameter_sensitivity_summary.csv).

**Velocity-gene sensitivity.** Moderate gamma-R2 thresholds 0.05 and 0.10 largely preserved the primary GraphVelo direction. At 0.20, the lower tail of cell-level agreement became negative. Stability relative to a primary fit is not proof that the primary fit is correct. See the [gene-QC overview](overview/05_velocity_qc.png).

**Fish-aware marker inference.** Cell-level Wilcoxon marker rankings remain exploratory. The corrected analysis compares a type against other cells from the same fish/time using pseudobulk contrasts, then performs inference across fish. All eight retained modeled types have at least seven independent fish across the dataset; the excluded thymocyte group has only one and is explicitly marked insufficient. This does not mean every type at every stage has seven fish. See the [marker support summary](../results_snapshot/2026-09-15/moscot/fish_aware_celltype_marker_summary.csv).

**Transport support gates.** Type-level claims require at least 50 cells per endpoint stage and at least 20 cells from at least three fish for each endpoint type. Unsupported outputs remain available for transparency but should not be read as reliable transitions. These thresholds are reporting screens, not a formal power calculation. The [supported-only transitions](../results_snapshot/2026-09-15/moscot/moscot_M0_M1_type_transitions_supported_only.csv) identify which rows pass.

**Uncertainty scope.** The saved fish bootstrap resamples fitted summaries; it does not refit preprocessing, velocity, and transport. It describes variation among represented fish conditional on the fitted model, rather than full model uncertainty. Cells from one fish are not independent biological replicates.

### Human comparison: expression-program support, with negative findings

The cross-species analysis reused the zebrafish checkpoint and compared mapped programs with Popescu E-MTAB-7407 fetal-liver data and GSE189161. Popescu contributed 113,063 annotated cells from 14 donors. GSE189161 included 58,041 CD34-enriched cells from 26 samples; the replication analysis used nine fetal-liver donors. The human analysis did not rerun MOSCOT or GraphVelo.

| Comparison | Saved result | Interpretation |
|---|---|---|
| Popescu proliferation profile | Spearman rho 0.283; bootstrap interval 0.077–0.418 | Modest program-level concordance |
| Popescu exact proliferation ortholog specificity | Permutation p=0.0969 | Did not pass the specificity test |
| Popescu apoptosis profile | rho 0.083; interval 0.019–0.171 | Weak concordance; specificity p=0.05195 |
| Exploratory GSE189161 proliferation | rho 0.734; FDR≈0.004 | Stronger association, but only three harmonized broad lineages |
| Exploratory GSE189161 apoptosis | rho −0.220; interval −0.353–0.124; FDR≈0.473 | No convincing positive concordance |
| Human HSC/MPP proliferation age-slope ranks across cohorts | rho −0.036; FDR=0.7327 | Gene-specific slope magnitudes did not replicate |

The saved analysis also found some agreement in the directions of human age-associated gene trends; direction agreement and reproduction of effect-size rankings are different tests. The overall interpretation is limited expression-program support with substantial qualifications, especially for apoptosis. It does not demonstrate conserved cell-to-cell transport, RNA velocity, absolute growth, or causal regulation in humans. See the [cross-species report](../results_snapshot/2026-09-15/cross_species/results/validation_report.md) and [profile comparison table](../results_snapshot/2026-09-15/cross_species/results/cross_species_lineage_profile_concordance.csv).

### How the project evolved—and why imperfect results are included

| Stage of work | What was explored | Status in this record |
|---|---|---|
| Whole-atlas GraphVelo and cluster-23 pilot | Velocity fields, potential, gene trends, and candidate trajectories | Exploratory; cluster labels and mixed blood/non-blood markers do not establish blood lineage |
| Early blood-cell notebooks | Annotation handling and GraphVelo analysis | Some partial runs and saved errors; retained as development history |
| Early combined transport/velocity models | Coupling reweighting and older transition summaries | Superseded; current GraphVelo comparison does not alter native transport |
| Smaller marker-panel and human analyses | Initial ortholog mapping and program comparisons | Superseded by the expanded reviewed panel |
| Repaired native MOSCOT workflow | Exact IDs, official labels, separate scoring, native M0/M1 | Current primary saved run |
| Corrected audits | Fish-level inference, support filters, direction/null distinction, parameter effects | Basis for the conclusions on this page |

The original primary run contains **41 executed nonempty code cells and no saved errors**. Earlier notebooks include a categorical-median aggregation error, an undefined `SUBTYPE_RULES` error, and a partially executed GraphVelo notebook. Their outputs are preserved because they document what was attempted. An image generated before a later error is a historical output, not proof that the full notebook completed successfully.

Saved notebook figures can precede standalone corrections. Current interpretation should follow the corrected tables and overview figures, even when an older title or caption sounds more conclusive. The [repair report](../docs/PROJECT_REPAIR_REPORT.md) describes the earlier corrections in detail.

### What the completed project contributes

The completed work provides an executed native MOSCOT workflow, a documented species-appropriate growth panel, a comparison of growth-informed and uniform-input transport, fish-aware reporting, complementary velocity checks, and a preserved record of sensitivity and negative results. It shows **where the growth prior affects this dataset and where supporting evidence is weak**.

It does not establish a causal lineage map, an absolute division/death model, a demonstrated predictive advantage for M1, or an experimentally validated mechanism for the late-interval disagreement. No new wet-lab or lineage-tracing experiments were completed as part of this record. These are boundaries of the completed work, not omitted successes or promised future results.

### Reading the figures and reproducing the work

For an initial discussion, read the [six-page summary](MOSCOT_summary_for_professor.pdf), then the [executed primary notebook](executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb). For the full history, use the [80-page report](MOSCOT_project_figures_for_professor.pdf) or [figure index](FIGURE_INDEX.md). Every displayed figure has a status label; original vector files remain available for close inspection.

The notebooks under `executed_notebooks/` are exact historical copies with saved outputs and original environment/path assumptions. The [root notebook](../zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb) is the clean, portable entry point for rerunning. Its cleared outputs should not be confused with the executed record. Large raw datasets, the H5AD checkpoint, and coupling binaries are excluded from Git and listed in the data manifest. A new environment installation and a new end-to-end run were not performed during publication. See the [reproduction guide](../docs/REPRODUCIBILITY.md) and [package verification](VERIFICATION.md).

## How to interpret the history

- **Current:** native M0/M1, expanded 95/46 zebrafish marker panels, and the corrected fish-aware audits. Saved notebook plots may precede standalone corrected audits; the new summary charts and current scientific summary govern interpretation.
- **Historical/superseded:** older combined/reweighted transport, smaller marker panels and older marker inference. Retained to show development of the project; not additional evidence that the current model is correct.
- **Exploratory:** whole-atlas GraphVelo and cluster-23 trajectories. Numeric cluster identity is not a validated lineage label; some plotted genes are non-blood markers.
- **Saved errors/partial runs:** preserved and disclosed below. A figure produced before a later error is still a historical output, but the notebook is not a clean completed pipeline.

## Saved notebook history

| Notebook | Executed / nonempty cells | Saved errors | Status |
|---|---:|---:|---|
| [zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb](executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb) | 41/41 | 0 | CURRENT SAVED RUN |
| [hematopoiesis_GraphVelo_moscot.ipynb](executed_notebooks/hematopoiesis_GraphVelo_moscot.ipynb) | 40/40 | 1 | HISTORICAL / SAVED ERROR |
| [hematopoiesis_graphvelo_2.ipynb](executed_notebooks/hematopoiesis_graphvelo_2.ipynb) | 27/31 | 0 | HISTORICAL / PARTIAL RUN |
| [hematopoiesis_graphvelo.ipynb](executed_notebooks/hematopoiesis_graphvelo.ipynb) | 14/43 | 1 | HISTORICAL / SAVED ERROR |
| [zebrafish_human_cross_species_validation.before_expanded_panel.ipynb](executed_notebooks/zebrafish_human_cross_species_validation.before_expanded_panel.ipynb) | 19/19 | 0 | SUPERSEDED HUMAN VALIDATION |

The full error messages and source hashes are in [the notebook execution manifest](notebook_execution_manifest.json). Execution counts describe saved state and do not independently prove a fresh-kernel run. The prior repair report documents the primary clean-kernel execution.

## Figure coverage and provenance

Includes every PNG/PDF/SVG/JPEG found in the current MOSCOT figures, current human validation, earlier official blood GraphVelo figures, relevant archive figures, cluster-23 trajectories and whole-atlas GraphVelo figure folder. Also extracts every saved PNG output from the five notebooks above. Original copies are byte-identical to local sources. The appendix renders every page of each unique PDF and retains unique PNG outputs; identical content and paired SVG/PNG display versions are not repeated in the appendix. All original formats and duplicates remain in the artifact manifest. Charts from unrelated Revelio projects, the separate non-OT project, and the integrated proposal are not uploaded here.

[Artifact manifest](artifact_manifest.csv): source, path, status, interpretation note, size and SHA-256. [Gallery manifest](gallery_manifest.json): displayed page order and original links. PDF previews are for reading; original PDFs/SVGs retain publication-quality vectors.

## Main takeaways for discussion

The growth-prior source-mass effect is largest at 120–240 hpf, but GraphVelo has negative absolute agreement there. Positive corrected directional evidence occurs only at 16–19 and 48–72 hpf. Early 12 hpf sampling is insufficient for type-level claims. Epsilon/tau materially influence transport. Human expression-program concordance does not validate lineage or absolute growth. These limitations are part of the project results, not omitted failures.

## Rebuilding the package

`python analysis/build_professor_package.py --source /path/to/original/dynamo --repository /path/to/this/repository` copies the preserved source files and redraws summaries from the dated result snapshot. Report generation requires matplotlib, numpy, Pillow and PyMuPDF; these are reporting dependencies, not changes to the scientific environment.

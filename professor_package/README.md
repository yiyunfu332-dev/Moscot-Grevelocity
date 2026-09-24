# MOSCOT project: professor sharing package

This package preserves previously executed notebooks and the complete discovered figure collection for the MOSCOT / blood-GraphVelo project, including unsuccessful, exploratory and superseded analyses. It was assembled on 24 September 2026; no scientific model was rerun.

## Start here

1. [Six-page discussion summary](MOSCOT_summary_for_professor.pdf), [GitHub figure index](FIGURE_INDEX.md), and [complete illustrated figure report](MOSCOT_project_figures_for_professor.pdf): six new overview charts/pages drawn from saved results, followed by the historical/current figure appendix. Every page states its evidence status. The overview pages are newly drawn; the underlying fits are not new.
2. [Previously executed primary MOSCOT notebook](executed_notebooks/zebrafish_hematopoiesis_growth_mapping_graphvelo_moscot.ipynb): exact copy of the local saved run, **41/41 nonempty code cells executed and zero saved errors**. Original figures, tables, execution counts, kernel metadata and paths are preserved. Use this to review the past run; use the clean notebook at the repository root for a new run.
3. [Current scientific summary](../docs/MOSCOT_SUMMARY.md): corrected conclusions, including weak/negative results.
4. [Browse the figure gallery](gallery.html): download the package and open this file in a browser; GitHub does not render standalone HTML as a website.

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

# Zebrafish-to-human cross-species validation report

## Scope

This analysis reused the saved zebrafish checkpoint and did not rerun GraphVelo,
moscot, or the earlier transport workflow. The full zebrafish count H5AD was used only
to recover marker genes absent from the checkpoint's 2,986-gene selected feature space.

## Cohorts

- Zebrafish: 5,426 saved cells, 36 fish.
- Popescu E-MTAB-7407: 113,063 annotated fetal-liver cells,
  14 donors, 7.86-17.00 PCW.
- GSE189161: 58,041 CD34-enriched cells from 26 samples;
  the replication analysis used 9 fetal-liver donors.

## Overall evidence grade

**Proliferation: moderate program-level transcriptomic support; apoptosis: moderate program-level transcriptomic support; not causal or dynamic validation**

## Primary zebrafish-Popescu lineage result

- Proliferation: Spearman rho=0.283,
  95% donor/fish bootstrap CI [0.077, 0.418],
  ortholog-label permutation p=0.0969.
- Apoptosis: Spearman rho=0.083,
  95% donor/fish bootstrap CI [0.019, 0.171],
  ortholog-label permutation p=0.05195.

## Independent human HSC/MPP replication

- Proliferation gene temporal slopes: rho=-0.036,
  rank-correlation FDR=0.7327, same-direction
  fraction=60.4%, direction FDR=0.02929.
- Apoptosis gene temporal slopes: rho=0.365,
  rank-correlation FDR=0.04819, same-direction
  fraction=78.9%, direction FDR=0.000472.

Slope-direction agreement and slope-magnitude rank correlation answer different
questions. Significant direction agreement means the two human cohorts tend to move
the same way with fetal age; a nonsignificant rank correlation means they do not
reliably reproduce the ordering of gene-specific effect sizes.

## Interpretation boundary

Agreement supports conservation of the mapped expression program. It does not validate
an absolute growth rate, a cell-to-cell transport coupling, or an RNA-velocity vector in
humans. Strong causal validation would require independent human lineage tracing,
perturbation, or measured proliferation/death outcomes.

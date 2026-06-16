# CZSH_FF pyCSH_Zn paper code package

This package provides a clean paper-code snapshot of the pyCSH_Zn workflow for
constructing, validating, screening, and selected quasi-static deformation
diagnostics of Zn-modified C-S-H atomistic models.

## Supported

- Brick-based C-S-H model construction based on pyCSH logic.
- Q1_Zn and Q2b_Zn motif construction.
- Single-Zn and multi-Zn model generation.
- CementFF4-Zn compatible typing and output.
- LAMMPS read/run0/static minimization workflow.
- Post-min Zn-O coordination diagnostics.
- Ensemble analysis and representative model selection.
- Selected x-direction quasi-static deformation diagnostics.

## Not Supported Or Not Claimed

- No finite-temperature MD.
- No `md_ready_candidate` label.
- No final elastic constants.
- No production mechanical properties.
- Overcoordinated models are minimum-valid, not ideal ZnO4 tetrahedral
  structures.

## Start Here

- `QUICK_START.md`: minimal command sequence.
- `REPRODUCIBILITY_COMMANDS.md`: lightweight verification commands.
- `pyCSH_Zn/FINAL_WORKFLOW.md`: full workflow reference.
- `pyCSH_Zn/SUPPLEMENTARY_WORKFLOW_TABLES.md`: supplement-ready tables.
- `checks/PAPER_CODE_AUDIT.md`: package audit summary.

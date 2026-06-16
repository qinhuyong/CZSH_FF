# pyCSH-Zn Static Candidate Workflow

This directory contains a conservative pyCSH-based workflow for Zn-modified
C-S-H static candidate generation, CementFF4/CementFF4-Zn-compatible LAMMPS
data output, validation, static minimization, ensemble screening, and selected
quasi-static mechanics smoke tests.

v1.8 completed the main generation, validation, ensemble, and selected mechanics
workflow. v1.9 is documentation and reproducibility cleanup only.

## Start Here

Use `FINAL_WORKFLOW.md` for the complete reproducible workflow, script index,
validation label table, coordination-quality table, output directory map, and
recommended paper-scale commands.

Use `MANUSCRIPT_METHODS_DRAFT.md` as the manuscript methods starting point.
Use `MANUSCRIPT_FIGURE_PLAN.md`, `SUPPLEMENTARY_WORKFLOW_TABLES.md`, and
`RELEASE_AUDIT_v1.9.md` for manuscript figures, supplement tables, and release
audit notes.

## Supported Workflows

- Pure C-S-H static generation.
- Q2b_Zn and Q1_Zn single-Zn static candidates.
- Q1_Zn motif diagnostics and screening.
- Single-Zn ensemble generation and representative selection.
- Single-structure multi-Zn candidates, including independent Q1_Zn + Q2b_Zn
  mixed-motif structures.
- Multi-Zn ensemble generation and representative selection.
- Selected multi-Zn x-direction quasi-static mechanics.

## Boundaries

This project does not claim:

- finite-temperature MD;
- MD-ready classification;
- final elastic constants;
- production mechanical properties;
- experimental uniqueness of a Zn local motif;
- support for the old `mixed_Q1_Q2b_Zn` site type.

CementFF4-Zn parameters, the 2.5 Angstrom Zn-O validation gate, and validation
semantics are not changed by v1.9.

## Default And Opt-In Mechanics

`examples/07_run_quasistatic_mechanics.py` remains the default mechanics smoke
test and its default targets are only:

- `pure_csh`
- `q2b_zn`

Q1_Zn and multi-Zn mechanics are opt-in:

- `examples/11_run_q1_quasistatic_mechanics.py`
- `examples/18_run_selected_multi_zn_mechanics.py`

## Minimal Smoke Test

Run from the repository root:

```bash
python pyCSH_Zn/examples/01_generate_pure_csh.py
python pyCSH_Zn/examples/02_generate_q2b_zn.py
python pyCSH_Zn/examples/08_generate_q1_zn.py
python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py
```

## Common Commands

Single-Zn ensemble:

```bash
python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 20 --seed-start 1000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation
python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 5 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots
```

Multi-Zn ensemble:

```bash
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 9 --seed-start 8400 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
```

Selected multi-Zn mechanics:

```bash
python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 4 --include-overcoordinated --write-plots
```

## Coordination Language

For multi-Zn models, post-min valid means minimum-valid: every Zn center has at
least four O neighbors within the unchanged 2.5 Angstrom Zn-O gate.

`ideal_fourfold` means all centers have coordination exactly 4.
`overcoordinated` means all centers have coordination >= 4 and at least one
center has coordination > 4. `undercoordinated_failed` means at least one center
has coordination < 4.

Overcoordinated models are minimum-valid candidates, not ideal ZnO4 tetrahedral
structures.

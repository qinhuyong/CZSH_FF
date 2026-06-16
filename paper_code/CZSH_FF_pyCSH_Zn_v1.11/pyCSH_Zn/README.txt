pyCSH-Zn workflow
=================

pyCSH-Zn is a focused pyCSH-based workflow for static Zn-modified C-S-H
candidate generation, CementFF4/CementFF4-Zn-compatible LAMMPS data output,
validation, static minimization, ensemble screening, and quasi-static mechanics
smoke tests.

Detailed reproducible commands are in FINAL_WORKFLOW.md.

Current scope
-------------

Supported:

- Pure C-S-H static generation.
- Q2b_Zn and Q1_Zn single-Zn static candidates.
- Single-Zn ensembles and representative selection.
- Single-structure multi-Zn candidates, including independent Q1_Zn + Q2b_Zn
  mixed-motif structures.
- Multi-Zn ensembles and selected multi-Zn quasi-static mechanics.

Not claimed:

- Finite-temperature MD.
- MD-ready classification.
- Final elastic constants.
- Production mechanical properties.
- Experimental uniqueness of any Zn local motif.
- The old mixed_Q1_Q2b_Zn site type.

Default mechanics
-----------------

examples/07_run_quasistatic_mechanics.py remains the default pure/Q2b
quasi-static mechanics runner. Its default active targets are only:

- pure_csh
- q2b_zn

Q1_Zn mechanics and multi-Zn mechanics are opt-in workflows:

- examples/11_run_q1_quasistatic_mechanics.py
- examples/18_run_selected_multi_zn_mechanics.py

Quick smoke test
----------------

Run from the repository root:

    python pyCSH_Zn/examples/01_generate_pure_csh.py
    python pyCSH_Zn/examples/02_generate_q2b_zn.py
    python pyCSH_Zn/examples/08_generate_q1_zn.py
    python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py

Recommended workflows
---------------------

Single-Zn ensemble:

    python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 20 --seed-start 1000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation
    python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 5 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots

Multi-Zn ensemble:

    python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 9 --seed-start 8400 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots

Selected multi-Zn mechanics:

    python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 4 --include-overcoordinated --write-plots

Main references
---------------

- FINAL_WORKFLOW.md: command recipes, script index, outputs, validation labels,
  coordination-quality definitions, and recommended paper-scale workflows.
- MANUSCRIPT_METHODS_DRAFT.md: manuscript-ready methods draft.
- MANUSCRIPT_FIGURE_PLAN.md: manuscript figure structure and data sources.
- SUPPLEMENTARY_WORKFLOW_TABLES.md: supplement-ready workflow tables.
- RELEASE_AUDIT_v1.9.md: v1.9 release audit and boundary checklist.
- v1.9-final-cleanup-and-manuscript-workflow-status.md: final version summary.

Coordination language
---------------------

For multi-Zn models, minimum-valid means every Zn center has at least four O
neighbors within the unchanged 2.5 Angstrom Zn-O validation gate.

- ideal_fourfold: every Zn center has coordination exactly 4.
- overcoordinated: all centers have coordination >= 4 and at least one center
  has coordination > 4.
- undercoordinated_failed: at least one Zn center has coordination < 4.

Overcoordinated models may be retained as minimum-valid candidates when clearly
labeled, but they are not ideal ZnO4 tetrahedral structures.

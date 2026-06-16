# Release Audit v1.9

This audit summarizes the frozen v1.9 documentation/reproducibility closure
release.

## Latest Release

- Commit: `6d187a7 docs: add final workflow and manuscript methods`
- Tag: `v1.9-final-cleanup-and-manuscript-workflow`

## Changed Files In v1.9

- `README.txt`
- `README_pyCSH_Zn_workflow.md`
- `FINAL_WORKFLOW.md`
- `MANUSCRIPT_METHODS_DRAFT.md`
- `v1.9-final-cleanup-and-manuscript-workflow-status.md`

## Scientific Boundaries Unchanged

- No scientific Python workflow logic was modified.
- No CementFF4-Zn force-field parameters were modified.
- The 2.5 Angstrom Zn-O coordination threshold was not relaxed.
- Single-Zn and multi-Zn validation semantics were not changed.
- The old `mixed_Q1_Q2b_Zn` site type remains unsupported.
- finite-temperature MD is not included.
- `md_ready_candidate` is not introduced.
- final elastic constants are not claimed.
- production mechanical properties are not claimed.

## Smoke-Test Commands

```bash
python -m py_compile pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py
python pyCSH_Zn/examples/01_generate_pure_csh.py
python pyCSH_Zn/examples/02_generate_q2b_zn.py
python pyCSH_Zn/examples/08_generate_q1_zn.py
python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py
```

## Default Mechanics Target Confirmation

The default `examples/07_run_quasistatic_mechanics.py` workflow remains limited
to:

- `pure_csh`
- `q2b_zn`

Q1_Zn mechanics and multi-Zn mechanics remain opt-in workflows.

## Ignored Output Confirmation

Generated outputs and runtime artifacts are not release files:

- `output_Y/`
- `__pycache__/`
- `log.lammps`
- LAMMPS dump files
- large intermediate outputs

These should remain ignored and uncommitted.

## Known Limitations

- The workflow constructs static Zn-C-S-H candidates and selected quasi-static
  deformation diagnostics; it is not a finite-temperature MD workflow.
- Quasi-static mechanics outputs are pipeline diagnostics and should not be
  treated as final elastic constants.
- Overcoordinated models are minimum-valid candidates, not ideal ZnO4
  tetrahedral structures.
- Candidate motifs are construction and screening motifs, not claims of unique
  experimentally resolved Zn environments.

## Appropriate Use

This release can be used for:

- reproducible static pyCSH-Zn candidate construction;
- CementFF4/CementFF4-Zn-compatible LAMMPS data generation;
- validation and post-minimization screening;
- single-Zn and multi-Zn ensemble generation;
- representative model selection;
- selected quasi-static deformation diagnostics;
- manuscript methods and supplement preparation.

This release should not be used as evidence for:

- production MD equilibration;
- final elastic constants;
- production mechanical properties;
- MD-ready Zn-C-S-H structures.

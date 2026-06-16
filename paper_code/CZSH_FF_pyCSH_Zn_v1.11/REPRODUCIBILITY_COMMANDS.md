# Reproducibility Commands

Run from this package root.

```bash
python -m py_compile pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py
python pyCSH_Zn/examples/01_generate_pure_csh.py
python pyCSH_Zn/examples/02_generate_q2b_zn.py
python pyCSH_Zn/examples/08_generate_q1_zn.py
python pyCSH_Zn/examples/04_build_lammps_inputs.py
python pyCSH_Zn/examples/06_run_static_relaxation.py
python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py
```

In a clean package checkout, `examples/07_run_quasistatic_mechanics.py` needs
the post-minimized pure C-S-H and Q2b_Zn references produced by
`04_build_lammps_inputs.py` and `06_run_static_relaxation.py`.

After running `examples/07_run_quasistatic_mechanics.py`, inspect:

```text
pyCSH_Zn/output_Y/workflow_v1/quasistatic_mechanics/mechanics_summary.json
```

The actual target list should be:

```text
pure_csh
q2b_zn
```

Q1_Zn mechanics and multi-Zn mechanics are opt-in workflows and should not
appear in the default `07` target list.

# Paper Code Audit

Package:

```text
paper_code/CZSH_FF_pyCSH_Zn_v1.11/
```

## Audit Summary

- Package folder created.
- Current stable `pyCSH_Zn/` source snapshot copied into the package.
- Ignored and generated outputs excluded.
- No Python workflow logic changed.
- No CementFF4-Zn force-field parameters changed.
- No validation semantics changed.
- No new Zn motif added.
- `examples/07_run_quasistatic_mechanics.py` default targets unchanged.
- Q1_Zn mechanics and multi-Zn mechanics remain opt-in.
- No final elastic constants claimed.
- No production mechanical properties claimed.

## Excluded From Package

- `output_Y/`
- `__pycache__/`
- `log.lammps`
- `*.lammpstrj`
- LAMMPS dump files
- large temporary/intermediate run outputs
- `.git/`
- local caches

## Required Checks

```bash
git status --short
git diff --stat
git diff --check
git status --ignored
find paper_code/CZSH_FF_pyCSH_Zn_v1.11 -type d -name "__pycache__"
find paper_code/CZSH_FF_pyCSH_Zn_v1.11 -name "log.lammps"
find paper_code/CZSH_FF_pyCSH_Zn_v1.11 -name "*.lammpstrj"
find paper_code/CZSH_FF_pyCSH_Zn_v1.11 -path "*/output_Y"
rg -n "md_ready_candidate|finite-temperature MD input|production elastic constants|final elastic constants|mixed_Q1_Q2b_Zn" paper_code/CZSH_FF_pyCSH_Zn_v1.11
```

Keyword hits are acceptable only when they are negative boundary statements,
unsupported-site statements, or audit text.

## Lightweight Package Verification

```bash
python -m py_compile paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/01_generate_pure_csh.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/02_generate_q2b_zn.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/08_generate_q1_zn.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/04_build_lammps_inputs.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/06_run_static_relaxation.py
python paper_code/CZSH_FF_pyCSH_Zn_v1.11/pyCSH_Zn/examples/07_run_quasistatic_mechanics.py
```

The `07` target list should remain:

```text
pure_csh
q2b_zn
```

## Use Boundary

This package is suitable as manuscript supplementary code for construction,
validation, screening, ensemble analysis, and selected quasi-static deformation
diagnostics. It is not a finite-temperature MD package, does not provide final
elastic constants, and does not provide production mechanical properties.

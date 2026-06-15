# pyCSH-Zn Static CementFF4 Workflow

This branch adds a conservative pyCSH-Zn workflow for generating CementFF4/CementFF4-Zn-compatible atomistic candidate structures.

## Supported

- Pure C-S-H generation with `enable_zinc = False`.
- Q2b_Zn candidate generation as the main Zn path.
- Q1_Zn candidate generation as a conservative static path.
- Charge-balanced ZnO2(OH)2-style substitutional candidate motif for Q2b_Zn.
- Fixed CementFF4 atom, bond, and angle type maps.
- LAMMPS data output with molecule IDs and `CS-Info`.
- JSON summaries for mapping, water, Zn, validation, and post-processing.
- LAMMPS input templates for read checks, run0 checks, static minimization, shell static relaxation, and quasi-static elastic deformation.
- LAMMPS static-relaxation runner for pure C-S-H first and Q2b_Zn second.

## Not Supported Yet

- Production finite-temperature CementFF4 core-shell MD.
- Claims that generated Zn structures are experimentally unique or chemically proven.
- Interlayer Zn or Ca-substitution-control Zn modes.
- mixed_Q1_Q2b_Zn.
- Automatic rescue of mixed_Q1_Q2b_Zn candidates.
- Silent lowering of W/Si when water placement is difficult.

## Main Files

- `forcefields/CementFF4_Zn_parameters.json`: authoritative CementFF4-Zn parameter database.
- `forcefields/build_cementff4_zn.py`: generates `in.CementFF4_Zn` and `cementff4_type_map.json`.
- `mod_zinc.py`: Q1_Zn and Q2b_Zn candidate motif construction and Zn summary.
- `mod_water.py`: water topology/contact helpers.
- `mod_write_Y.py`: CementFF4 fixed-map LAMMPS data writer.
- `validate_cementff_data.py`: static data validator.
- `lammps_templates/build_inputs.py`: LAMMPS input template generator.
- `postprocess/analyze_structure.py`: normalized RDF, coordination, angle, and contact summaries.
- `forcefields/validate_forcefield.py`: pair_coeff syntax and pair coverage audit.

## Quick Start

```bash
python examples/01_generate_pure_csh.py
python examples/02_generate_q2b_zn.py
python examples/03_validate_outputs.py
python examples/04_build_lammps_inputs.py
python examples/05_postprocess_q2b_zn.py
python examples/06_run_static_relaxation.py
python examples/07_run_quasistatic_mechanics.py
python examples/10_screen_q1_motifs.py
```

Outputs are written under:

```text
output_Y/workflow_v1/
```

`examples/06_run_static_relaxation.py` requires a LAMMPS executable available as `lmp` or via `LAMMPS_EXE`.

## Fixed Type Maps

Atom types:

1. Ca
2. Si
3. O_core
4. O_shell / O(S)
5. Ow
6. Oh
7. Hw
8. Hoh
9. Zn
10. Al optional
11. Cl optional

Bond types:

1. O_core-O_shell
2. Ow-Hw
3. Oh-Hoh

Angle types:

1. Hw-Ow-Hw
2. O-Si-O / Oh-Si-O / Oh-Si-Oh
3. Si-Oh-H
4. O-Zn-O / Oh-Zn-O / Oh-Zn-Oh
5. Zn-Oh-H

## Classification Policy

The validator uses static candidate classifications:

- `valid_static_candidate`
- `valid_q1_zn_candidate`
- `valid_q2b_zn_candidate`
- `needs_static_relaxation`
- `failed_charge`
- `failed_charge_assignment`
- `failed_topology`
- `failed_water_contacts`
- `failed_zinc_geometry`
- `failed_csinfo`
- `experimental_md_only`

`valid_q1_zn_candidate` and `valid_q2b_zn_candidate` both mean static CementFF4-Zn candidates only. Neither is an MD-ready label.

An MD-ready classification is intentionally not used. Finite-temperature core-shell MD must be validated separately before that label is introduced.

## v1.1 Static Relaxation

The v1.1 static-relaxation workflow builds and tests:

- `in.read_check`
- `in.run0`
- `in.minimize_static`
- `in.elastic_x_plus`
- `in.elastic_x_minus`

It runs pure C-S-H first, then Q2b_Zn. LAMMPS `write_data` does not preserve the custom `CS-Info` section, so the runner reattaches original `CS-Info` entries by atom ID before post-minimization validation. No finite-temperature MD is run.

The x-strain templates use scale factors `1.0+strain` and `1.0-strain`. They are quasi-static elastic input-validation smoke tests only, not final elastic constants.

## v1.2 Quasi-Static Mechanics

`examples/07_run_quasistatic_mechanics.py` starts from the v1.1 post-minimized structures and runs x-direction strain cases at +/-0.001, +/-0.002, and +/-0.003 for pure C-S-H and Q2b_Zn. It writes mechanics CSV/JSON summaries and simple SVG plots.

Q1_Zn mechanics remains opt-in and should only be run after the Q1 static-relaxation path validates.

This is a controlled quasi-static mechanics pipeline validation only. It is not a final elastic constants workflow, not a production mechanical-property calculation, and not finite-temperature MD.

## v1.3.2 Q1 Motif Screening

`examples/10_screen_q1_motifs.py` screens multiple topology-valid Q1_Zn motifs from the same generated base structure. Each candidate is independently minimized and validated; failures are kept in the summary rather than hidden.

The screening writes CSV/JSON rankings under `output_Y/workflow_v1/q1_motif_screening/`. Q1_Zn remains outside the default mechanics workflow unless a screened post-minimized candidate validates as `valid_q1_zn_candidate`.

## v1.3.3 Q1 Selected Static Candidate

Default Q1_Zn generation now uses `PYCSH_ZN_Q1_SELECTION_MODE=ranked_static`. The generator enumerates topology-valid Q1 candidates and selects a deterministic static candidate using pre-min motif geometry, hydroxylation safety, tetrahedral angle deviation, O-O separation, and secondary-shell crowding diagnostics.

This default path does not run post-min screening and does not hard-code a site atom ID. Post-min validation is still required before any opt-in Q1 mechanics workflow.

## Opt-In Q1 Mechanics Smoke Test

`examples/11_run_q1_quasistatic_mechanics.py` runs the x-direction quasi-static smoke-test strain set for Q1_Zn only. It refuses to start unless the Q1 reference is the post-minimized `valid_q1_zn_candidate` from `examples/09_run_q1_static_relaxation.py`.

Outputs are written under `output_Y/workflow_v1/mechanics_q1_zn/`. This is not a final elastic-constants workflow and not a production mechanical-property calculation. `examples/07_run_quasistatic_mechanics.py` still defaults to pure C-S-H and Q2b_Zn only.

## CS-Info Policy

`CS-Info` contains entries for all atoms. Bonded `O_core`/`O_shell` pairs share the same CSID. Non-core-shell atoms have singleton CSIDs. The validator checks both complete CS-Info coverage and bonded core-shell pair consistency.

## Force-Field Audit

`examples/04_build_lammps_inputs.py` generates `forcefield_validation_report.json` next to `in.CementFF4_Zn`. The report checks static `pair_coeff` syntax and pair coverage. The generated `lammps_inputs/in.read_check` should be run with the target LAMMPS executable before production calculations.

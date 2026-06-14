# pyCSH-Zn Static CementFF4 Workflow

This branch adds a conservative pyCSH-Zn workflow for generating CementFF4/CementFF4-Zn-compatible atomistic candidate structures.

## Supported

- Pure C-S-H generation with `enable_zinc = False`.
- Q2b_Zn candidate generation as the main Zn path.
- Q1_Zn and mixed_Q1_Q2b_Zn remain prototype/debug modes.
- Charge-balanced ZnO2(OH)2-style substitutional candidate motif for Q2b_Zn.
- Fixed CementFF4 atom, bond, and angle type maps.
- LAMMPS data output with molecule IDs and `CS-Info`.
- JSON summaries for mapping, water, Zn, validation, and post-processing.
- LAMMPS input templates for read checks, static minimization, shell static relaxation, and quasi-static elastic deformation.
- Experimental short-MD template, clearly marked as not validated.

## Not Supported Yet

- Production finite-temperature CementFF4 core-shell MD.
- Claims that generated Zn structures are experimentally unique or chemically proven.
- Interlayer Zn or Ca-substitution-control Zn modes.
- Automatic rescue of Q1_Zn or mixed_Zn candidates.
- Silent lowering of W/Si when water placement is difficult.

## Main Files

- `forcefields/CementFF4_Zn_parameters.json`: authoritative CementFF4-Zn parameter database.
- `forcefields/build_cementff4_zn.py`: generates `in.CementFF4_Zn` and `cementff4_type_map.json`.
- `mod_zinc.py`: Q2b_Zn candidate motif construction and Zn summary.
- `mod_water.py`: water topology/contact helpers.
- `mod_write_Y.py`: CementFF4 fixed-map LAMMPS data writer.
- `validate_cementff_data.py`: static data validator.
- `lammps_templates/build_inputs.py`: LAMMPS input template generator.
- `postprocess/analyze_structure.py`: RDF, coordination, angle, and contact summaries.

## Quick Start

```bash
python examples/01_generate_pure_csh.py
python examples/02_generate_q2b_zn.py
python examples/03_validate_outputs.py
python examples/04_build_lammps_inputs.py
python examples/05_postprocess_q2b_zn.py
```

Outputs are written under:

```text
output_Y/workflow_v1/
```

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
- `valid_q2b_zn_candidate`
- `needs_static_relaxation`
- `failed_charge`
- `failed_topology`
- `failed_water_contacts`
- `failed_zinc_geometry`
- `failed_csinfo`
- `experimental_md_only`

`md_ready_candidate` is intentionally not used. Finite-temperature core-shell MD must be validated separately before that label is introduced.

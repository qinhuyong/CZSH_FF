# pyCSH-Zn v1.0-static Release Acceptance Report

Commit accepted for freeze: `76a22d1`

Release scope:

- Static pyCSH-Zn structure generation.
- Pure C-S-H generation.
- Q2b_Zn generation as the main Zn path.
- CementFF4/CementFF4-Zn LAMMPS data output.
- CementFF4-Zn force-field include generation.
- Static validation, force-field syntax audit, and normalized RDF post-processing.

Out of scope for this release:

- Finite-temperature MD.
- MD-ready classification.
- Changes to Zn motif generation.
- Q1_Zn or mixed_Zn rescue logic.

## Acceptance Commands

Run from `pyCSH_Zn/`:

```bash
python examples/01_generate_pure_csh.py
python examples/02_generate_q2b_zn.py
python examples/03_validate_outputs.py
python examples/04_build_lammps_inputs.py
python examples/05_postprocess_q2b_zn.py
```

All commands completed successfully.

## Pure C-S-H Validation

Output:

```text
output_Y/workflow_v1/pure_csh/pure_csh_cementff.data
```

Result:

```text
classification: valid_static_candidate
atoms: 418
bonds: 164
angles: 227
charge_assignment.n_bad: 0
total_charge: -1.8207657603852567e-14
csinfo.n_entries: 418
csinfo.n_bad_pairs: 0
water.n_bad_water: 0
zinc.n_zinc: 0
```

## Q2b_Zn Validation

Output:

```text
output_Y/workflow_v1/q2b_zn/q2b_zn_cementff_zn.data
```

Result:

```text
classification: valid_q2b_zn_candidate
atoms: 418
bonds: 164
angles: 229
charge_assignment.n_bad: 0
total_charge: -1.8207657603852567e-14
csinfo.n_entries: 418
csinfo.n_bad_pairs: 0
water.n_bad_water: 0
zinc.n_zinc: 1
zinc.coordination_2p3: 4
zinc.coordination_2p5: 4
```

Nearest Zn-O distances from validator:

```text
1.587104644675014 A  O_core
1.6367823778554067 A O_core
1.9500001851828117 A Oh
1.9500004419050767 A Oh
```

## Force-Field Audit

Output:

```text
output_Y/workflow_v1/q2b_zn/forcefield_validation_report.json
```

Result:

```text
ok: true
n_pair_coeff_lines: 45
expected_pair_coeff_lines: 45
missing_pairs: []
duplicates: []
```

Optional Al/Cl pair coefficients are reported separately as optional missing
pairs because Al and Cl are optional types and are not used by the v1.0-static
default workflow.

The generated read-check input is:

```text
output_Y/workflow_v1/q2b_zn/lammps_inputs/in.read_check
```

## Post-Processing Acceptance

Output:

```text
output_Y/workflow_v1/q2b_zn/postprocess/structure_analysis.json
```

Result:

```text
rdf_definition.normalized: true
rdf_definition.box_volume: 4681.991606199814
zinc_coordination_2p3: [4]
zinc_coordination_2p5: [4]
```

Generated normalized RDF files:

```text
rdf_Zn_O.csv
rdf_Zn_Si.csv
rdf_Zn_Ca.csv
rdf_Si_O.csv
rdf_Ca_O.csv
```

## Freeze Statement

`valid_q2b_zn_candidate` is a static candidate classification only. It is not
an MD-ready label. An MD-ready classification is intentionally not used by
the supported workflow.

This release is frozen as `pyCSH-Zn v1.0-static`.

pyCSH-Zn v1.0-static
====================

pyCSH-Zn is a focused extension of pyCSH for generating Zn-modified
C-S-H static candidate structures and exporting CementFF4/CementFF4-Zn
compatible LAMMPS files.

Scope
-----

This release supports:

1. Pure C-S-H generation from the pyCSH brick workflow.
2. Q2b_Zn candidate generation as the main Zn path.
3. Charge-balanced ZnO2(OH)2-style Q2b_Zn candidate environments.
4. CementFF4/CementFF4-Zn LAMMPS data output.
5. CementFF4-Zn force-field include generation from
   forcefields/CementFF4_Zn_parameters.json.
6. Static validation of type maps, charge assignment, CS-Info, water topology,
   and Zn coordination.
7. Normalized RDF output, Zn-O coordination, selected angle summaries, and
   water/contact summaries.

This release does not reproduce the original CementFF4 finite-temperature MD
workflow and does not reproduce the Morales-Melgares simulation protocol.
Finite-temperature core-shell MD is not a v1.0-static acceptance criterion.

Do not interpret valid_q2b_zn_candidate as MD-ready. It means only that the
generated structure is a static CementFF4-Zn candidate with valid charge
assignment, topology, CS-Info, water topology, and Zn coordination checks.

Quick start
-----------

Run from this directory:

    python examples/01_generate_pure_csh.py
    python examples/02_generate_q2b_zn.py
    python examples/03_validate_outputs.py
    python examples/04_build_lammps_inputs.py
    python examples/05_postprocess_q2b_zn.py

For the v1.1-static-relaxation workflow, run after the first four commands:

    python examples/06_run_static_relaxation.py

For the v1.2-quasistatic-mechanics workflow, run after v1.1:

    python examples/07_run_quasistatic_mechanics.py

Expected validation result:

    pure_csh_cementff.data  -> valid_static_candidate
    q2b_zn_cementff_zn.data -> valid_q2b_zn_candidate

Outputs are written under:

    output_Y/workflow_v1/

CementFF4-Zn charges
--------------------

forcefields/CementFF4_Zn_parameters.json stores CementFF4 SI Table S1 charges:

    O_core      +0.84819
    O_shell     -2.84819
    Ow          -1.1128
    Hw          +0.5564
    Oh          -1.4
    Hoh/H       +0.4
    Zn          +2.0

validate_cementff_data.py checks every atom charge against this table. A
nonzero total charge is not the only charge gate; per-type charge assignment
must also pass.

CS-Info policy
--------------

The LAMMPS data file contains a CS-Info entry for every atom.

- O_core/O_shell bonded pairs share the same CSID.
- Non-core-shell atoms receive singleton CSIDs.
- The validator checks both CS-Info coverage and bonded core-shell pair
  consistency.

Main files
----------

- mod_zinc.py
  Q2b_Zn site selection, substitution, hydroxylation, charge balance, and
  zinc_summary.json generation.

- mod_write_Y.py
  CementFF4/CementFF4-Zn LAMMPS data writer with fixed type maps, molecule
  IDs, water topology, and CS-Info.

- validate_cementff_data.py
  Static validator with per-atom charge assignment checks and orthogonal or
  triclinic minimum-image distances.

- forcefields/CementFF4_Zn_parameters.json
  CementFF4-Zn parameter database.

- forcefields/build_cementff4_zn.py
  Generates in.CementFF4_Zn, cementff4_type_map.json, and
  forcefield_validation_report.json.

- forcefields/validate_forcefield.py
  Audits pair_coeff syntax and pair coverage. The generated LAMMPS
  lammps_inputs/in.read_check file should also be run with the target LAMMPS
  executable to catch runtime pair_coeff syntax errors.

- lammps_templates/build_inputs.py
  Generates read-check, run0, static minimization, static shell relaxation,
  and quasi-static elastic templates. It does not generate finite-temperature
  MD inputs.

- postprocess/analyze_structure.py
  Produces normalized RDF CSV files using PBC, shell volume, number density,
  and box volume; also writes coordination, angle, and contact summaries.

Validation classes
------------------

- valid_static_candidate
- valid_q2b_zn_candidate
- needs_static_relaxation
- failed_charge
- failed_charge_assignment
- failed_topology
- failed_water_contacts
- failed_zinc_geometry
- failed_csinfo
- experimental_md_only

The label md_ready_candidate is intentionally not used.

Post-processing outputs
-----------------------

examples/05_postprocess_q2b_zn.py writes:

- structure_analysis.json
- rdf_Zn_O.csv
- rdf_Zn_Si.csv
- rdf_Zn_Ca.csv
- rdf_Si_O.csv
- rdf_Ca_O.csv

These RDF files are normalized g(r), not raw distance histograms.

v1.1-static-relaxation
----------------------

examples/06_run_static_relaxation.py runs LAMMPS for pure C-S-H first and then
Q2b_Zn. For each target it tests:

- in.read_check
- in.run0
- in.minimize_static
- in.elastic_x_plus
- in.elastic_x_minus

LAMMPS write_data output does not preserve the custom CS-Info section, so the
runner reattaches the original CS-Info by atom ID before post-minimization
validation. The validator semantics are unchanged.

The x-direction strain templates are small-strain smoke tests for quasi-static
input validation only. They do not claim final elastic constants or final
mechanical-property results.

The runner writes:

    output_Y/workflow_v1/static_relaxation_report.json

Finite-temperature MD is not run.

v1.2-quasistatic-mechanics
--------------------------

examples/07_run_quasistatic_mechanics.py starts from the v1.1 post-minimized
pure C-S-H and Q2b_Zn structures. It generates x-direction +/-0.001, +/-0.002,
and +/-0.003 strain cases, runs LAMMPS run0/deform/minimize/run0, reattaches
CS-Info, validates each deformed minimized data file, and writes CSV/JSON
summaries plus simple SVG plots.

This is quasi-static mechanics pipeline validation only. It is not a final
elastic-constant calculation and not a production mechanical-property result.
Finite-temperature MD is not generated or run.

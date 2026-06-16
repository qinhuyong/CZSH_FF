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

## v1.4 Zn-C-S-H Ensemble Generator

`examples/12_generate_zn_csh_ensemble.py` generates constrained random Zn-C-S-H static candidate ensembles from independent seeds.

Example:

```text
python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 20 --seed-start 1000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation
```

Outputs are written under `output_Y/workflow_v1/zn_csh_ensemble/`:

- `ensemble_manifest.json`
- `ensemble_summary.csv`
- `ensemble_summary.json`
- `accepted_models.csv`
- `rejected_models.csv`
- `models/model_000001/`

Accepted models are post-minimized structures classified as `valid_q1_zn_candidate` or `valid_q2b_zn_candidate`. Rejected models include generation failures, failed initial validation, and failed post-min validation. The ensemble run continues after individual model failures.

Supported modes are `q2b_only`, `q1_only`, and `q1_q2b_mixture`. The `q1_q2b_mixture` mode is an ensemble-level mixture: each independent model is assigned either Q1_Zn or Q2b_Zn according to `--q1-fraction`. It is not the old `mixed_Q1_Q2b_Zn` site type and does not create a single local mixed motif.

v1.4 supports one Zn motif per generated structure. Multiple Zn motifs within the same structure are not yet supported because shared O/H topology, Zn-Zn separation, and coupled charge-balance constraints require additional screening.

This is static candidate ensemble generation, not finite-temperature MD, not final elastic constants, and not production mechanical-property calculation.

## v1.5 Ensemble Analysis And Selection

After generating a v1.4 ensemble, run:

```text
python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 5 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots
```

The analysis writes:

- `output_Y/workflow_v1/zn_csh_ensemble_analysis/ensemble_analysis_summary.json`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/ensemble_analysis_summary.csv`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/motif_survival_summary.csv`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/failure_reason_summary.csv`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/mechanics_ready_models.csv`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/representative_models.json`
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/plots/*.svg`

`mechanics_ready_models.csv` is a filtered input manifest for a future opt-in v1.6 batch mechanics workflow. It does not mean production mechanics, final elastic constants, or final mechanical properties have been calculated.

`representative_models.json` records selected model IDs, seeds, motif types, model directories, post-min data paths, validation labels, selection scores, and the reason each model was selected.

v1.5 only performs ensemble statistics, failure reason summaries, and representative selection. It does not add finite-temperature MD, single-structure mixed Q1+Q2b, single-structure multi-Zn, final elastic constants, or production mechanical-property workflows.

## v1.6-alpha Single-Structure Multi-Zn

`examples/15_generate_multi_zn_structure.py` generates a minimal single-structure multi-Zn alpha candidate.

Supported modes:

- `multi_q2b`: two independent Q2b_Zn motifs in one C-S-H structure.
- `multi_q1`: two independent Q1_Zn motifs in one C-S-H structure.
- `q1_q2b_single_structure_mixture`: one Q1_Zn motif plus one Q2b_Zn motif in the same C-S-H structure.

Example commands:

```text
python pyCSH_Zn/examples/15_generate_multi_zn_structure.py --mode multi_q2b --n-q2b 2 --seed 6100 --run-static-relaxation
python pyCSH_Zn/examples/15_generate_multi_zn_structure.py --mode multi_q1 --n-q1 2 --seed 6200 --run-static-relaxation
python pyCSH_Zn/examples/15_generate_multi_zn_structure.py --mode q1_q2b_single_structure_mixture --n-q1 1 --n-q2b 1 --seed 6300 --run-static-relaxation
```

Concepts are distinct:

- Ensemble-level `q1_q2b_mixture`: different models in the ensemble are Q1_Zn or Q2b_Zn.
- Single-structure Q1+Q2b mixed motif: one structure contains independent Q1_Zn and Q2b_Zn motifs.
- `mixed_Q1_Q2b_Zn` site type: old prototype site type; still unsupported.

v1.6-alpha requires independent motif records, no duplicate Si site, PBC Zn-Zn separation screening, no repeated hydroxylated O core-shell pair, charge balance, CS-Info validation, water topology validation, and per-center Zn-O coordination diagnostics.

This alpha stage does not add finite-temperature MD, final elastic constants, `md_ready_candidate`, production mechanical-property calculation, or Q1_Zn to the default `examples/07_run_quasistatic_mechanics.py` targets.

## v1.6-beta Multi-Zn Site Pair Screening

`examples/16_screen_multi_zn_combinations.py` screens multiple single-structure
multi-Zn site combinations for `multi_q2b`, `multi_q1`, and
`q1_q2b_single_structure_mixture`. Each combination starts independently from
the same pure C-S-H parent structure; failed candidates are recorded and do not
stop the screen.

Post-min valid multi-Zn candidates are minimum-valid candidates: every Zn center
has at least four O neighbors within the unchanged 2.5 Angstrom validation
threshold. This is distinct from an ideal ZnO4 fourfold result.

Coordination quality labels are:

- `ideal_fourfold`: every Zn center has coordination exactly 4.
- `overcoordinated`: at least one Zn center has coordination greater than 4 and no center is undercoordinated.
- `undercoordinated_failed`: at least one Zn center has coordination less than 4, so the structure is `failed_multi_zn_candidate`.
- `minimum_valid`: all centers satisfy coordination >= 4 when no more specific label applies.

The best candidate files only promote post-min valid candidates with the matching
`valid_multi_*` validation label. Undercoordinated cases such as 3;5 or 5;3 are
kept as failed candidates and are not best candidates.

Current v1.6-beta screening identifies an ideal-fourfold `multi_q2b` best
candidate, while the best `multi_q1` and single-structure Q1+Q2b mixed-motif
candidates are overcoordinated minimum-valid candidates. They should not be
described as ideal fourfold ZnO4 motifs.

## v1.7 Multi-Zn Ensemble Generator

`examples/17_generate_multi_zn_ensemble.py` performs constrained multi-Zn
ensemble generation and analysis using the v1.6-beta single-structure screening
logic.

Example commands:

```text
python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 20 --seed-start 1000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation
python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 5 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots
python pyCSH_Zn/examples/16_screen_multi_zn_combinations.py --mode q1_q2b_single_structure_mixture --n-q1 1 --n-q2b 1 --seed 7300 --max-combinations 10 --run-static-relaxation
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 9 --seed-start 8400 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
```

Supported v1.7 modes:

- `multi_q2b_ensemble`
- `multi_q1_ensemble`
- `q1_q2b_single_structure_mixed_ensemble`
- `mixed_multi_zn_ensemble`

Outputs are written under `output_Y/workflow_v1/multi_zn_ensemble/`, including
`multi_zn_ensemble_summary.csv`, accepted/rejected model CSV files,
`representative_multi_zn_models.json`, `mechanics_ready_multi_zn_models.csv`,
survival and coordination-quality summaries, failure-reason summaries, per-model
diagnostics, and optional SVG plots.

`mechanics_ready_multi_zn_models.csv` is only a downstream manifest for future
opt-in mechanics. v1.7 does not run batch mechanics, finite-temperature MD,
production mechanical-property calculations, or final elastic-constant
workflows.

The same coordination-quality language from v1.6-beta applies: `ideal_fourfold`
means all Zn centers have coordination exactly 4, while `overcoordinated`
candidates are only minimum-valid and must not be described as ideal ZnO4
tetrahedral structures.

## v1.8 Selected Multi-Zn Batch Mechanics

`examples/18_run_selected_multi_zn_mechanics.py` reads the v1.7
`mechanics_ready_multi_zn_models.csv` manifest and runs x-direction
quasi-static mechanics for selected post-min-valid multi-Zn models.

Example:

```text
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 6 --seed-start 8500 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 4 --include-overcoordinated --write-plots
```

The runner refuses `failed_multi_zn_candidate`, `undercoordinated_failed`, and
`postmin_valid = false` models. Each strain case starts independently from the
same post-min reference structure; sequential strain accumulation is not used.

Outputs are written under
`output_Y/workflow_v1/selected_multi_zn_mechanics/`, including batch summary
CSV/JSON files, accepted/failed mechanics model CSV files, stress-strain
summaries by model and mode, failure-reason summaries, per-model strain case
directories, and optional SVG plots.

v1.8 is a selected multi-Zn quasi-static mechanics pipeline. It is not
finite-temperature MD, not final elastic constants, and not a production
mechanical-property workflow. `ideal_fourfold` and `overcoordinated` references
are recorded separately in mechanics summaries; overcoordinated models are
minimum-valid candidates, not ideal ZnO4 tetrahedra.

## CS-Info Policy

`CS-Info` contains entries for all atoms. Bonded `O_core`/`O_shell` pairs share the same CSID. Non-core-shell atoms have singleton CSIDs. The validator checks both complete CS-Info coverage and bonded core-shell pair consistency.

## Force-Field Audit

`examples/04_build_lammps_inputs.py` generates `forcefield_validation_report.json` next to `in.CementFF4_Zn`. The report checks static `pair_coeff` syntax and pair coverage. The generated `lammps_inputs/in.read_check` should be run with the target LAMMPS executable before production calculations.

# Supplementary Workflow Tables

These tables are written in a supplement-ready style and summarize scripts,
labels, coordination-quality definitions, output columns, and reproducibility
commands. They describe the workflow as of v1.9 and do not introduce new
scientific functionality.

## Table S1. Script Index

| Script | Main input | Main output | Default or opt-in | Runs LAMMPS | Runtime |
| --- | --- | --- | --- | --- | --- |
| `01_generate_pure_csh.py` | pyCSH defaults | pure C-S-H data and summaries | default | no | short |
| `02_generate_q2b_zn.py` | generated C-S-H parent | Q2b_Zn data and summaries | default | no | short |
| `06_run_static_relaxation.py` | pure/Q2b data | read/run0/minimize reports | default | yes | medium |
| `07_run_quasistatic_mechanics.py` | post-min pure/Q2b references | pure/Q2b mechanics summaries | default | yes | medium |
| `08_generate_q1_zn.py` | generated C-S-H parent | Q1_Zn data and summaries | opt-in | no | short |
| `09_run_q1_static_relaxation.py` | Q1_Zn data | Q1 post-min validation | opt-in | yes | medium |
| `10_screen_q1_motifs.py` | Q1 candidate sites | Q1 motif ranking and diagnostics | opt-in | yes | long |
| `11_run_q1_quasistatic_mechanics.py` | valid Q1 post-min reference | Q1 mechanics summaries | opt-in | yes | medium |
| `12_generate_zn_csh_ensemble.py` | independent seeds | single-Zn ensemble | opt-in | optional | medium/long |
| `13_analyze_zn_csh_ensemble.py` | single-Zn ensemble summaries | statistics and representatives | opt-in | no | short |
| `15_generate_multi_zn_structure.py` | pure parent seed | one multi-Zn structure | opt-in | optional | medium |
| `16_screen_multi_zn_combinations.py` | pure parent seed | multi-Zn combination screening | opt-in | yes | long |
| `17_generate_multi_zn_ensemble.py` | independent seeds | multi-Zn ensemble | opt-in | yes | long |
| `18_run_selected_multi_zn_mechanics.py` | mechanics-ready multi-Zn CSV | selected mechanics summaries | opt-in | yes | medium/long |

## Table S2. Validation Labels

| Label | Interpretation |
| --- | --- |
| `valid_static_candidate` | Pure C-S-H candidate passes static validation. |
| `valid_q1_zn_candidate` | Single Q1_Zn candidate passes charge, topology, CS-Info, water, and Zn coordination gates. |
| `valid_q2b_zn_candidate` | Single Q2b_Zn candidate passes charge, topology, CS-Info, water, and Zn coordination gates. |
| `valid_multi_q1_zn_candidate` | Multi-Zn structure with Q1_Zn centers passes all per-center Zn gates. |
| `valid_multi_q2b_zn_candidate` | Multi-Zn structure with Q2b_Zn centers passes all per-center Zn gates. |
| `valid_multi_q1_q2b_zn_candidate` | Single structure containing independent Q1_Zn and Q2b_Zn motifs passes all per-center Zn gates. |
| `failed_multi_zn_candidate` | One or more multi-Zn validation gates failed. |
| `needs_static_relaxation` | Initial candidate is not post-min valid but is close enough to justify static relaxation. |

No validation label in this workflow means MD-ready.

## Table S3. Coordination Quality

| Quality | Definition | Use |
| --- | --- | --- |
| `ideal_fourfold` | All Zn centers have coordination exactly 4 within the 2.5 Angstrom Zn-O gate. | Preferred representative when available. |
| `overcoordinated` | All Zn centers have coordination >= 4 and at least one center has coordination > 4. | Minimum-valid, but not ideal ZnO4. |
| `undercoordinated_failed` | At least one Zn center has coordination < 4. | Rejected for mechanics. |
| `minimum-valid` | All Zn centers have coordination >= 4. Includes ideal_fourfold and overcoordinated. | Eligibility criterion, not a geometry claim. |

## Table S4. Accepted/Rejected Model Table Template

| Column | Meaning |
| --- | --- |
| `model_id` | Unique model identifier. |
| `seed` | Random seed used for the independent parent structure. |
| `requested_mode` | Requested ensemble mode. |
| `internal_mode` | Actual motif class used for the model. |
| `n_q1_actual` | Number of Q1_Zn motifs in the accepted structure. |
| `n_q2b_actual` | Number of Q2b_Zn motifs in the accepted structure. |
| `n_Zn_total` | Total Zn centers in the structure. |
| `postmin_validation_label` | Post-min validation classification. |
| `coordination_quality` | ideal_fourfold, overcoordinated, or undercoordinated_failed. |
| `per_center_coordination` | Semicolon-separated Zn-O coordination values. |
| `accepted` | Whether the model passed the post-min acceptance gates. |
| `mechanics_ready` | Whether the model is eligible for selected quasi-static mechanics. |
| `failure_reason` | Rejection reason, if any. |

## Table S5. Mechanics Output Columns

| Column | Meaning |
| --- | --- |
| `strain` | Prescribed engineering strain. |
| `actual_strain` | `final_lx / initial_lx - 1`. |
| `initial_lx` | x box length before deformation. |
| `final_lx` | x box length after deformation and minimization. |
| `stress_xx_bar` | LAMMPS Pxx in bar. |
| `stress_xx_GPa` | `stress_xx_bar / 10000`. |
| `pressure_bar` | LAMMPS pressure in bar. |
| `pressure_GPa` | `pressure_bar / 10000`. |
| `energy_initial` | Initial potential energy from thermo output. |
| `energy_final` | Final potential energy from thermo output. |
| `validation_label_after_strain` | Validation label for the minimized strained structure. |
| `per_center_coordination_after_strain` | Per-Zn coordination after strained minimization. |
| `coordination_quality_after_strain` | Coordination quality after strained minimization. |
| `case_ok` | Whether the strain case passed the workflow gates. |
| `failure_reason` | Case failure reason, if any. |

## Table S6. Reproducibility Commands

| Purpose | Command |
| --- | --- |
| Minimal smoke test | `python pyCSH_Zn/examples/01_generate_pure_csh.py` |
| Minimal Q2b_Zn | `python pyCSH_Zn/examples/02_generate_q2b_zn.py` |
| Minimal Q1_Zn | `python pyCSH_Zn/examples/08_generate_q1_zn.py` |
| Default pure/Q2b mechanics | `python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py` |
| Single-Zn ensemble | `python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 20 --seed-start 1000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation` |
| Single-Zn analysis | `python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 5 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots` |
| Multi-Zn ensemble smoke test | `python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 3 --seed-start 9000 --n-q1 1 --n-q2b 1 --max-combinations-per-model 5 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots` |
| Selected mechanics | `python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 4 --include-overcoordinated --write-plots` |

## Table S7. Release Boundary

| Boundary | Status |
| --- | --- |
| CementFF4-Zn parameters modified | No |
| 2.5 Angstrom Zn-O gate relaxed | No |
| Validation semantics changed | No |
| finite-temperature MD included | No |
| MD-ready label introduced | No |
| final elastic constants claimed | No |
| production mechanical properties claimed | No |
| Q1 added to default `07` mechanics | No |

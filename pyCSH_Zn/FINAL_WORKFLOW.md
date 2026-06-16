# pyCSH-Zn Final Workflow

This document is the reproducible workflow map for pyCSH-Zn through v1.9. The
workflow builds static C-S-H and Zn-modified C-S-H candidates, validates
CementFF4/CementFF4-Zn-compatible data files, runs static minimization, screens
single-Zn and multi-Zn ensembles, and runs selected quasi-static mechanics smoke
tests.

It does not claim finite-temperature MD, final elastic constants, production
mechanical properties, or MD-ready structures.

## A. Minimal Smoke-Test Workflow

Run from the repository root:

```bash
python pyCSH_Zn/examples/01_generate_pure_csh.py
python pyCSH_Zn/examples/02_generate_q2b_zn.py
python pyCSH_Zn/examples/08_generate_q1_zn.py
python pyCSH_Zn/examples/07_run_quasistatic_mechanics.py
```

`examples/07_run_quasistatic_mechanics.py` defaults to `pure_csh` and `q2b_zn`
only. Q1_Zn and multi-Zn mechanics are opt-in.

## B. Single-Zn Workflow

```bash
python pyCSH_Zn/examples/01_generate_pure_csh.py
python pyCSH_Zn/examples/02_generate_q2b_zn.py
python pyCSH_Zn/examples/08_generate_q1_zn.py
python pyCSH_Zn/examples/09_run_q1_static_relaxation.py
python pyCSH_Zn/examples/10_screen_q1_motifs.py
```

Optional Q1 mechanics:

```bash
python pyCSH_Zn/examples/11_run_q1_quasistatic_mechanics.py
```

## C. Multi-Zn Workflow

Single-structure multi-Zn generation and screening:

```bash
python pyCSH_Zn/examples/15_generate_multi_zn_structure.py --mode q1_q2b_single_structure_mixture --n-q1 1 --n-q2b 1 --seed 6300 --run-static-relaxation
python pyCSH_Zn/examples/16_screen_multi_zn_combinations.py --mode q1_q2b_single_structure_mixture --n-q1 1 --n-q2b 1 --seed 7300 --max-combinations 10 --run-static-relaxation
```

Multi-Zn ensemble:

```bash
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 9 --seed-start 8400 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
```

## D. Selected Mechanics Workflow

```bash
python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 4 --include-overcoordinated --write-plots
```

Every strain case starts independently from the same post-minimized reference
structure. Sequential strain accumulation is not used.

## E. Recommended Paper-Scale Workflow

Smoke test:

```bash
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode mixed_multi_zn_ensemble --n-models 3 --seed-start 9000 --n-q1 1 --n-q2b 1 --max-combinations-per-model 5 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
```

Paper-scale single-Zn ensemble, start with 20 to 50 models:

```bash
python pyCSH_Zn/examples/12_generate_zn_csh_ensemble.py --n-models 50 --seed-start 10000 --mode q1_q2b_mixture --q1-fraction 0.5 --target-zn-si 0.05 --run-static-relaxation
python pyCSH_Zn/examples/13_analyze_zn_csh_ensemble.py --ensemble-dir pyCSH_Zn/output_Y/workflow_v1/zn_csh_ensemble --top-n 10 --select-for-mechanics --prefer-balanced-q1-q2b --write-plots
```

Paper-scale multi-Zn ensemble, start with about 20 models per class before
scaling:

```bash
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode multi_q2b_ensemble --n-models 20 --seed-start 11000 --n-q2b 2 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode multi_q1_ensemble --n-models 20 --seed-start 12000 --n-q1 2 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
python pyCSH_Zn/examples/17_generate_multi_zn_ensemble.py --mode q1_q2b_single_structure_mixed_ensemble --n-models 20 --seed-start 13000 --n-q1 1 --n-q2b 1 --max-combinations-per-model 10 --min-zn-zn-distance 5.0 --run-static-relaxation --prefer-ideal-fourfold --write-plots
```

Selected mechanics should be limited to representative or mechanics-ready
models:

```bash
python pyCSH_Zn/examples/18_run_selected_multi_zn_mechanics.py --models-csv pyCSH_Zn/output_Y/workflow_v1/multi_zn_ensemble/mechanics_ready_multi_zn_models.csv --max-models 6 --include-overcoordinated --write-plots
```

## F. Output Directory Map

- `output_Y/workflow_v1/pure_csh/`: pure C-S-H generation outputs.
- `output_Y/workflow_v1/q2b_zn/`: Q2b_Zn single-Zn outputs.
- `output_Y/workflow_v1/q1_zn/`: Q1_Zn single-Zn outputs.
- `output_Y/workflow_v1/static_relaxation_report.json`: pure/Q2b static relaxation report.
- `output_Y/workflow_v1/quasistatic_mechanics/`: default pure/Q2b mechanics.
- `output_Y/workflow_v1/mechanics_q1_zn/`: opt-in Q1 mechanics.
- `output_Y/workflow_v1/q1_motif_screening/`: Q1 motif screening.
- `output_Y/workflow_v1/zn_csh_ensemble/`: single-Zn ensemble.
- `output_Y/workflow_v1/zn_csh_ensemble_analysis/`: single-Zn ensemble analysis.
- `output_Y/workflow_v1/multi_zn_structure/`: single multi-Zn structure.
- `output_Y/workflow_v1/multi_zn_screening/`: multi-Zn site-combination screening.
- `output_Y/workflow_v1/multi_zn_ensemble/`: multi-Zn ensemble.
- `output_Y/workflow_v1/selected_multi_zn_mechanics/`: selected multi-Zn mechanics.

## G. Script Index

| Script | Input | Output | Default or opt-in | Runs LAMMPS | Runtime |
| --- | --- | --- | --- | --- | --- |
| `01_generate_pure_csh.py` | seed/config defaults | pure C-S-H data | default | no | short |
| `02_generate_q2b_zn.py` | generated pyCSH parent | Q2b_Zn data | default | no | short |
| `06_run_static_relaxation.py` | pure/Q2b data | read/run0/minimize reports | default | yes | medium |
| `07_run_quasistatic_mechanics.py` | post-min pure/Q2b | pure/Q2b mechanics summaries | default | yes | medium |
| `08_generate_q1_zn.py` | generated pyCSH parent | Q1_Zn data | opt-in | no | short |
| `09_run_q1_static_relaxation.py` | Q1_Zn data | Q1 post-min validation | opt-in | yes | medium |
| `10_screen_q1_motifs.py` | Q1 candidate pool | ranked Q1 candidates | opt-in | yes | long |
| `11_run_q1_quasistatic_mechanics.py` | valid Q1 post-min data | Q1 mechanics summaries | opt-in | yes | medium |
| `12_generate_zn_csh_ensemble.py` | independent seeds | single-Zn ensemble | opt-in | optional | medium/long |
| `13_analyze_zn_csh_ensemble.py` | single-Zn ensemble CSV/JSON | analysis and representatives | opt-in | no | short |
| `15_generate_multi_zn_structure.py` | pure parent seed | one multi-Zn structure | opt-in | optional | medium |
| `16_screen_multi_zn_combinations.py` | pure parent seed | multi-Zn combination ranking | opt-in | yes | long |
| `17_generate_multi_zn_ensemble.py` | independent seeds | multi-Zn ensemble | opt-in | yes | long |
| `18_run_selected_multi_zn_mechanics.py` | mechanics-ready multi-Zn CSV | selected mechanics summaries | opt-in | yes | medium/long |

## H. Validation Labels

| Label | Meaning |
| --- | --- |
| `valid_static_candidate` | Pure C-S-H static candidate passes charge, topology, CS-Info, and water checks. |
| `valid_q1_zn_candidate` | Single Q1_Zn static candidate passes static validation. |
| `valid_q2b_zn_candidate` | Single Q2b_Zn static candidate passes static validation. |
| `valid_multi_q1_zn_candidate` | Multi-Zn structure with Q1_Zn centers passes all per-center gates. |
| `valid_multi_q2b_zn_candidate` | Multi-Zn structure with Q2b_Zn centers passes all per-center gates. |
| `valid_multi_q1_q2b_zn_candidate` | Single structure with independent Q1_Zn and Q2b_Zn motifs passes all per-center gates. |
| `failed_multi_zn_candidate` | At least one multi-Zn gate fails, often a Zn center below four O neighbors. |
| `needs_static_relaxation` | Initial single-Zn candidate is close enough to justify static relaxation but is not post-min valid yet. |

None of these labels means MD-ready.

## I. Coordination Quality

| Quality | Meaning |
| --- | --- |
| `ideal_fourfold` | All Zn centers have coordination exactly 4. |
| `overcoordinated` | All Zn centers have coordination >= 4 and at least one center has coordination > 4. |
| `undercoordinated_failed` | At least one Zn center has coordination < 4. |
| `minimum-valid` | All Zn centers have coordination >= 4; includes ideal_fourfold and overcoordinated. |

Overcoordinated models are minimum-valid candidates, not ideal ZnO4 tetrahedral
structures.

## J. What This Workflow Does Not Claim

- It does not run finite-temperature MD.
- It does not reproduce a production CementFF4 MD protocol.
- It does not report final elastic constants.
- It does not report production mechanical properties.
- It does not use an `md_ready_candidate` label.
- It does not enable the old `mixed_Q1_Q2b_Zn` site type.
- It does not modify CementFF4-Zn parameters or relax the 2.5 Angstrom Zn-O gate.

"""v2.6 Q2b_Zn water sanitizer audit and short equilibration gate."""

from __future__ import print_function

import json
import math
import os
import re
import subprocess

from audit_v25_water_q2b import (
    AUDIT as V25_AUDIT,
    LABELS,
    LMP_EXE,
    OUT,
    Q2B_DATA,
    Q2B_FF,
    Q2B_SUMMARY,
    parse_data,
    parse_lammps_log,
    top_force_atoms,
    water_audit,
    zinc_geometry,
)


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
AUDIT = os.path.join(OUT, "v26_water_gate")
PURE_DATA = os.path.join(OUT, "example_pure_csh_cementff1.data")


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def write_input(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines))
        f.write("\n")


def run_lammps(input_path, cwd):
    log_path = os.path.splitext(input_path)[0] + ".log"
    if os.path.exists(log_path):
        os.remove(log_path)
    if not os.path.exists(LMP_EXE):
        return {"code": None, "skipped": "LAMMPS executable not found", "log": log_path}
    proc = subprocess.run(
        [LMP_EXE, "-in", input_path, "-log", log_path],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=180,
    )
    return {"code": proc.returncode, "log": log_path, "stdout_tail": proc.stdout[-2000:]}


def common_header(data_file, ff_file):
    return [
        "clear",
        "units metal",
        "dimension 3",
        "atom_style full",
        "boundary p p p",
        "box tilt large",
        "read_data {}".format(data_file),
        "include {}".format(ff_file),
        "neighbor 2.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        "thermo 100",
        "thermo_style custom step pe ebond eangle evdwl ecoul elong temp fnorm fmax press",
    ]


def make_inputs():
    data = os.path.abspath(Q2B_DATA)
    ff = os.path.abspath(Q2B_FF)
    inputs = {}

    stage0 = os.path.join(AUDIT, "q2b_stage0_run0.in")
    write_input(
        stage0,
        common_header(data, ff)
        + [
            "dump fdump all custom 1 q2b_stage0_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "write_data q2b_stage0.data nocoeff",
        ],
    )
    inputs["stage0_run0"] = stage0

    stage1 = os.path.join(AUDIT, "q2b_stage1_water_only.in")
    write_input(
        stage1,
        common_header(data, ff)
        + [
            "group mobile_water type 5 7",
            "group frozen subtract all mobile_water",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "compute twater mobile_water temp",
            "velocity mobile_water create 50.0 246813 dist gaussian mom yes rot yes",
            "fix lang mobile_water langevin 50.0 100.0 0.02 97531",
            "fix nvew mobile_water nve",
            "fix_modify lang temp twater",
            "dump fdump all custom 100 q2b_stage1_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data q2b_stage1_water.data nocoeff",
        ],
    )
    inputs["stage1_water_only"] = stage1

    stage2 = os.path.join(AUDIT, "q2b_stage2_water_znoh.in")
    write_input(
        stage2,
        common_header(os.path.join(AUDIT, "q2b_stage1_water.data"), ff)
        + [
            "group mobile type 5 7 6 8 9",
            "group frozen subtract all mobile",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "compute tmob mobile temp",
            "velocity mobile create 50.0 13579 dist gaussian mom yes rot yes",
            "fix lang mobile langevin 50.0 100.0 0.02 86420",
            "fix nvew mobile nve",
            "fix_modify lang temp tmob",
            "dump fdump all custom 100 q2b_stage2_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data q2b_stage2_water_znoh.data nocoeff",
        ],
    )
    inputs["stage2_water_znoh"] = stage2

    stage3 = os.path.join(AUDIT, "q2b_stage3_full_short.in")
    write_input(
        stage3,
        common_header(os.path.join(AUDIT, "q2b_stage2_water_znoh.data"), ff)
        + [
            "compute tall all temp",
            "velocity all create 50.0 97531 dist gaussian mom yes rot yes",
            "fix lang all langevin 50.0 100.0 0.02 12345",
            "fix nvea all nve",
            "fix_modify lang temp tall",
            "dump fdump all custom 100 q2b_stage3_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 20000",
            "write_data q2b_stage3_full_short.data nocoeff",
        ],
    )
    inputs["stage3_full_short"] = stage3

    freeze_water = os.path.join(AUDIT, "q2b_diag_freeze_water.in")
    write_input(
        freeze_water,
        common_header(data, ff)
        + [
            "group water_atoms type 5 7",
            "fix freeze water_atoms setforce 0.0 0.0 0.0",
            "dump fdump all custom 100 q2b_diag_freeze_water_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "min_style cg",
            "min_modify dmax 0.01 line quadratic",
            "minimize 1e-6 1e-8 300 3000",
            "write_data q2b_diag_freeze_water.data nocoeff",
        ],
    )
    inputs["diag_freeze_water"] = freeze_water
    return inputs


def run_gate(inputs):
    results = {}
    sequence = [
        "stage0_run0",
        "stage1_water_only",
        "stage2_water_znoh",
        "stage3_full_short",
        "diag_freeze_water",
    ]
    prereq = {
        "stage2_water_znoh": os.path.join(AUDIT, "q2b_stage1_water.data"),
        "stage3_full_short": os.path.join(AUDIT, "q2b_stage2_water_znoh.data"),
    }
    dump_names = {
        "stage0_run0": "q2b_stage0_forces.lammpstrj",
        "stage1_water_only": "q2b_stage1_forces.lammpstrj",
        "stage2_water_znoh": "q2b_stage2_forces.lammpstrj",
        "stage3_full_short": "q2b_stage3_forces.lammpstrj",
        "diag_freeze_water": "q2b_diag_freeze_water_forces.lammpstrj",
    }
    base = parse_data(Q2B_DATA)
    for name in sequence:
        if name in prereq and not os.path.exists(prereq[name]):
            results[name] = {"skipped": "missing prerequisite {}".format(prereq[name])}
            continue
        result = run_lammps(inputs[name], AUDIT)
        result["parsed"] = parse_lammps_log(result.get("log", ""))
        dump_path = os.path.join(AUDIT, dump_names[name])
        result["top_forces"] = top_force_atoms(dump_path, base, 20)
        results[name] = result
    return results


def summarize_data_file(path):
    if not os.path.exists(path):
        return None
    data = parse_data(path)
    wa = water_audit(data, highlighted=[58, 106, 268, 324, 380])
    return {
        "data_file": path,
        "water_topology": {
            "n_water": wa["n_water"],
            "n_tip4p_compatible": wa["n_tip4p_compatible"],
            "n_bad_water": wa["n_bad_water"],
            "water_summary": wa["all_water_summary"],
            "highlighted": wa["highlighted"],
        },
        "zinc_geometry": zinc_geometry(data),
    }


def classify(summary, protocols):
    sanitizer = summary.get("zinc_summary", {}).get("water_sanitizer", {})
    if summary["q2b_initial_water_audit"]["n_bad_water"] != 0:
        return "failed_tip4p_topology"
    if sanitizer.get("rejected"):
        return "failed_water_sanitization"
    stage3 = protocols.get("stage3_full_short", {})
    parsed = stage3.get("parsed", {})
    if stage3.get("code") != 0 or parsed.get("errors") or parsed.get("lost"):
        return "failed_short_equilibration"
    final = parsed.get("final") or {}
    if final.get("fmax", 1.0e9) > 50.0:
        return "failed_short_equilibration"
    final_geom = summary.get("post_stage_summaries", {}).get("stage3_full_short", {}).get("zinc_geometry", {})
    if final_geom.get("coordination_by_cutoff", {}).get("2.3") != 4:
        return "debug_only_bad_zinc_geometry"
    return "md_ready_candidate"


def main():
    ensure_dir(AUDIT)
    q2b_summary = json.load(open(Q2B_SUMMARY)) if os.path.exists(Q2B_SUMMARY) else {}
    pure_audit = water_audit(parse_data(PURE_DATA), highlighted=[58, 106, 268, 324, 380])
    q2b_audit = water_audit(parse_data(Q2B_DATA), highlighted=[58, 106, 268, 324, 380])
    inputs = make_inputs()
    protocols = run_gate(inputs)
    stage_files = {
        "initial": Q2B_DATA,
        "stage0_run0": os.path.join(AUDIT, "q2b_stage0.data"),
        "stage1_water_only": os.path.join(AUDIT, "q2b_stage1_water.data"),
        "stage2_water_znoh": os.path.join(AUDIT, "q2b_stage2_water_znoh.data"),
        "stage3_full_short": os.path.join(AUDIT, "q2b_stage3_full_short.data"),
        "diag_freeze_water": os.path.join(AUDIT, "q2b_diag_freeze_water.data"),
    }
    post = {name: summarize_data_file(path) for name, path in stage_files.items() if summarize_data_file(path)}
    report = {
        "inputs": {
            "pure_data": PURE_DATA,
            "q2b_data": Q2B_DATA,
            "q2b_forcefield": Q2B_FF,
            "q2b_summary": Q2B_SUMMARY,
        },
        "zinc_summary": q2b_summary,
        "pure_water_audit": {
            "n_water": pure_audit["n_water"],
            "n_tip4p_compatible": pure_audit["n_tip4p_compatible"],
            "n_bad_water": pure_audit["n_bad_water"],
            "water_summary": pure_audit["all_water_summary"],
        },
        "q2b_initial_water_audit": {
            "n_water": q2b_audit["n_water"],
            "n_tip4p_compatible": q2b_audit["n_tip4p_compatible"],
            "n_bad_water": q2b_audit["n_bad_water"],
            "water_summary": q2b_audit["all_water_summary"],
            "highlighted": q2b_audit["highlighted"],
        },
        "protocol_inputs": inputs,
        "protocol_results": protocols,
        "post_stage_summaries": post,
    }
    report["final_classification"] = classify(report, protocols)
    out = os.path.join(AUDIT, "v26_water_gate_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print("Wrote {}".format(out))
    print("Final classification: {}".format(report["final_classification"]))


if __name__ == "__main__":
    main()

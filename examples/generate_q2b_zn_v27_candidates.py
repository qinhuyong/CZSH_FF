"""Generate water-screened Q2b_Zn candidates for v2.7.

This script keeps the Zn motif unchanged.  It moves water screening upstream:
candidate structures with hard water contacts are repaired or rejected before
LAMMPS data output.
"""

from __future__ import print_function

import json
import os
import random
import shutil
import subprocess
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from mod_construct_brick_Y import get_all_bricks, pieces
from mod_construct_supercell_Y import get_angles, get_full_coordinates, resize_crystal
from mod_sample import fill_water, sample_Ca_Si_ratio
from mod_water import screen_and_repair_waters
from mod_write_Y import write_output
from mod_zinc import apply_zinc_modification, finalize_zinc_summary

from audit_v25_water_q2b import parse_data, parse_lammps_log, top_force_atoms, water_audit, zinc_geometry


OUT = os.path.join(ROOT, "output_Y")
AUDIT = os.path.join(OUT, "v27_q2b_water_resampling")
LMP_EXE = r"C:\Program Files\LAMMPS 64-bit 4Feb2025-MSMPI\bin\lmp.exe"


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def build_supercell(size):
    unitcell = np.array(
        [
            [6.7352, 0.0, 0.0],
            [-4.071295, 6.209521, 0.0],
            [0.7037701, -6.2095578, 13.9936836],
        ]
    )
    supercell = np.zeros((3, 3))
    for i in range(3):
        supercell[i, :] = unitcell[i, :] * size[i]
    return unitcell, supercell


def generate_candidate(seed, index, sorted_bricks):
    np.random.seed(seed)
    random.seed(seed + 10)
    size = (2, 2, 2)
    ca_si_ratio = 1.7
    w_si_ratio = 0.2
    widths = [0.06, 0.08, 0.08]
    unitcell, supercell = build_supercell(size)
    n_brick = size[0] * size[1] * size[2]
    crystal, n_ca, n_si, r_sioh, r_caoh, mcl, n_water, _ = sample_Ca_Si_ratio(
        sorted_bricks, ca_si_ratio, w_si_ratio, n_brick, widths
    )
    water_in_crystal = fill_water(crystal, n_water)
    crystal_rs, water_in_crystal_rs = resize_crystal(crystal, water_in_crystal, size)
    entries_crystal, entries_bonds, crystal_dict, water_dict = get_full_coordinates(
        crystal_rs, water_in_crystal_rs, size, pieces, False, [0]
    )
    entries_angle = get_angles(crystal_dict, water_dict, size)
    entries_crystal, crystal_dict, zinc_summary = apply_zinc_modification(
        entries_crystal,
        crystal_dict,
        supercell,
        0.03,
        "Q2b_Zn",
        seed,
        n_ca / n_si,
        "hydroxylate_two_oxygens",
        entries_bonds,
        entries_angle,
        False,
        True,
        1.95,
    )
    zinc_summary = finalize_zinc_summary(
        entries_crystal,
        entries_bonds,
        entries_angle,
        supercell,
        zinc_summary,
        "hydroxylate_two_oxygens",
        False,
    )
    entries_crystal, water_screen = screen_and_repair_waters(entries_crystal, entries_bonds, supercell)
    zn_geom = zinc_summary.get("pre_minimization_geometry", {})
    selected = zinc_summary.get("selected_sites", [{}])[0]
    candidate = {
        "index": index,
        "seed": seed,
        "N_water": n_water,
        "W_Si_ratio_actual": float(n_water) / float(n_si),
        "N_Si": n_si,
        "N_Ca": n_ca,
        "Zn_site_atom_id": selected.get("atom_id"),
        "Zn_site_piece": selected.get("piece"),
        "Zn_site_cell": selected.get("cell"),
        "water_screen": water_screen,
        "zinc_summary": zinc_summary,
        "classification": classify_candidate(water_screen, zinc_summary),
        "score": score_candidate(water_screen, zinc_summary),
    }
    return {
        "candidate": candidate,
        "entries_crystal": entries_crystal,
        "entries_bonds": entries_bonds,
        "entries_angle": entries_angle,
        "size": size,
        "crystal_rs": crystal_rs,
        "water_in_crystal_rs": water_in_crystal_rs,
        "supercell": supercell,
        "unitcell": unitcell,
        "n_ca": n_ca,
        "n_si": n_si,
        "r_sioh": r_sioh,
        "r_caoh": r_caoh,
        "mcl": mcl,
        "zinc_summary": zinc_summary,
    }


def classify_candidate(water_screen, zinc_summary):
    if water_screen.get("rejected"):
        return "failed_water_resampling"
    zn_class = zinc_summary.get("output_classification")
    if zn_class == "debug_only_bad_zinc_geometry":
        return "debug_only_bad_zinc_geometry"
    if water_screen.get("warning_contacts"):
        return "needs_short_equilibration_test"
    return "water_screened_candidate"


def score_candidate(water_screen, zinc_summary):
    rejected = len(water_screen.get("rejected", []))
    hard = len(water_screen.get("hard_fail_before", []))
    warn = sum(len(x.get("violations", [])) for x in water_screen.get("warning_contacts", []))
    repaired = len(water_screen.get("repaired_by_rotation", [])) + len(water_screen.get("repaired_by_translation", []))
    zn_bad = 1 if zinc_summary.get("output_classification") == "debug_only_bad_zinc_geometry" else 0
    return [rejected, zn_bad, warn, hard, repaired]


def write_best_candidate(best):
    prefix = "v27_q2b_best"
    best["zinc_summary"]["v27_water_screen"] = best["candidate"]["water_screen"]
    best["zinc_summary"]["output_classification"] = best["candidate"]["classification"]
    best["zinc_summary"].setdefault("classification_reasons", []).append(
        "v2.7 water-screened Q2b candidate; short equilibration gate still required"
    )
    write_output(
        0,
        best["entries_crystal"],
        best["entries_bonds"],
        best["entries_angle"],
        best["size"],
        best["crystal_rs"],
        best["water_in_crystal_rs"],
        best["supercell"],
        best["n_ca"],
        best["n_si"],
        best["r_sioh"],
        best["r_caoh"],
        best["mcl"],
        False,
        False,
        False,
        False,
        prefix,
        best["unitcell"],
        False,
        False,
        True,
        None,
        False,
        [0],
        False,
        True,
        best["zinc_summary"],
        True,
    )
    return {
        "data": os.path.join(OUT, prefix + "_cementff1_zn.data"),
        "ff": os.path.join(OUT, prefix + "_in.CementFF4_Zn_1"),
        "summary": os.path.join(OUT, prefix + "_zinc_summary_1.json"),
        "water_summary": os.path.join(OUT, prefix + "_cementff_water_summary_1.json"),
    }


def write_input(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines))
        f.write("\n")


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


def make_gate_inputs(files):
    data = os.path.abspath(files["data"])
    ff = os.path.abspath(files["ff"])
    inputs = {}
    stage0 = os.path.join(AUDIT, "v27_stage0_run0.in")
    write_input(
        stage0,
        common_header(data, ff)
        + [
            "dump fdump all custom 1 v27_stage0_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "write_data v27_stage0.data nocoeff",
        ],
    )
    inputs["stage0_run0"] = stage0
    stage1 = os.path.join(AUDIT, "v27_stage1_water_only.in")
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
            "dump fdump all custom 100 v27_stage1_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data v27_stage1_water.data nocoeff",
        ],
    )
    inputs["stage1_water_only"] = stage1
    stage2 = os.path.join(AUDIT, "v27_stage2_water_znoh.in")
    write_input(
        stage2,
        common_header(os.path.join(AUDIT, "v27_stage1_water.data"), ff)
        + [
            "group mobile type 5 7 6 8 9",
            "group frozen subtract all mobile",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "compute tmob mobile temp",
            "velocity mobile create 50.0 13579 dist gaussian mom yes rot yes",
            "fix lang mobile langevin 50.0 100.0 0.02 86420",
            "fix nvew mobile nve",
            "fix_modify lang temp tmob",
            "dump fdump all custom 100 v27_stage2_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data v27_stage2_water_znoh.data nocoeff",
        ],
    )
    inputs["stage2_water_znoh"] = stage2
    stage3 = os.path.join(AUDIT, "v27_stage3_full_short.in")
    write_input(
        stage3,
        common_header(os.path.join(AUDIT, "v27_stage2_water_znoh.data"), ff)
        + [
            "compute tall all temp",
            "velocity all create 50.0 97531 dist gaussian mom yes rot yes",
            "fix lang all langevin 50.0 100.0 0.02 12345",
            "fix nvea all nve",
            "fix_modify lang temp tall",
            "dump fdump all custom 100 v27_stage3_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data v27_stage3_full_short.data nocoeff",
        ],
    )
    inputs["stage3_full_short"] = stage3
    return inputs


def run_gate(inputs, base_data):
    results = {}
    sequence = ["stage0_run0", "stage1_water_only", "stage2_water_znoh", "stage3_full_short"]
    prereq = {
        "stage2_water_znoh": os.path.join(AUDIT, "v27_stage1_water.data"),
        "stage3_full_short": os.path.join(AUDIT, "v27_stage2_water_znoh.data"),
    }
    dumps = {
        "stage0_run0": "v27_stage0_forces.lammpstrj",
        "stage1_water_only": "v27_stage1_forces.lammpstrj",
        "stage2_water_znoh": "v27_stage2_forces.lammpstrj",
        "stage3_full_short": "v27_stage3_forces.lammpstrj",
    }
    base = parse_data(base_data)
    for name in sequence:
        if name in prereq and not os.path.exists(prereq[name]):
            results[name] = {"skipped": "missing prerequisite {}".format(prereq[name])}
            continue
        result = run_lammps(inputs[name], AUDIT)
        result["parsed"] = parse_lammps_log(result.get("log", ""))
        result["top_forces"] = top_force_atoms(os.path.join(AUDIT, dumps[name]), base, 20)
        results[name] = result
    return results


def summarize_stage(path):
    if not os.path.exists(path):
        return None
    data = parse_data(path)
    wa = water_audit(data, highlighted=[])
    return {
        "water": {
            "n_water": wa["n_water"],
            "n_tip4p_compatible": wa["n_tip4p_compatible"],
            "n_bad_water": wa["n_bad_water"],
            "summary": wa["all_water_summary"],
        },
        "zinc_geometry": zinc_geometry(data),
    }


def classify_gate(best, results, stage_summaries):
    if best["candidate"]["classification"] == "failed_water_resampling":
        return "failed_water_resampling"
    stage1 = results.get("stage1_water_only", {}).get("parsed", {})
    stage3_result = results.get("stage3_full_short", {})
    stage3 = stage3_result.get("parsed", {})
    if stage1.get("lost") or stage1.get("errors"):
        return "failed_short_equilibration"
    if stage3_result.get("code") != 0 or stage3.get("lost") or stage3.get("errors"):
        return "failed_short_equilibration"
    geom = stage_summaries.get("stage3_full_short", {}).get("zinc_geometry", {})
    if geom.get("coordination_by_cutoff", {}).get("2.3") != 4:
        return "debug_only_bad_zinc_geometry"
    return "md_ready_candidate"


def main():
    ensure_dir(AUDIT)
    attempts = 12
    base_seed = 23137
    _, sorted_bricks = get_all_bricks(pieces)
    candidates = []
    full = []
    for i in range(attempts):
        seed = base_seed + 101 * i
        try:
            generated = generate_candidate(seed, i, sorted_bricks)
            full.append(generated)
            candidates.append(generated["candidate"])
            print(
                "candidate {} seed {} classification {} score {}".format(
                    i, seed, generated["candidate"]["classification"], generated["candidate"]["score"]
                ),
                flush=True,
            )
        except Exception as exc:
            candidates.append({"index": i, "seed": seed, "classification": "failed", "error": str(exc), "score": [999]})
            print("candidate {} seed {} failed: {}".format(i, seed, exc), flush=True)

    viable = [item for item in full if item["candidate"]["classification"] != "failed_water_resampling"]
    if viable:
        best = sorted(viable, key=lambda x: x["candidate"]["score"])[0]
    else:
        best = sorted(full, key=lambda x: x["candidate"]["score"])[0] if full else None

    files = None
    gate_results = {}
    stage_summaries = {}
    final_classification = "failed_water_resampling"
    if best is not None:
        files = write_best_candidate(best)
        inputs = make_gate_inputs(files)
        gate_results = run_gate(inputs, files["data"])
        stage_paths = {
            "initial": files["data"],
            "stage0_run0": os.path.join(AUDIT, "v27_stage0.data"),
            "stage1_water_only": os.path.join(AUDIT, "v27_stage1_water.data"),
            "stage2_water_znoh": os.path.join(AUDIT, "v27_stage2_water_znoh.data"),
            "stage3_full_short": os.path.join(AUDIT, "v27_stage3_full_short.data"),
        }
        stage_summaries = {
            name: summarize_stage(path)
            for name, path in stage_paths.items()
            if summarize_stage(path) is not None
        }
        final_classification = classify_gate(best, gate_results, stage_summaries)

    report = {
        "workflow": {
            "water_origin": "Water sites are selected from each Brick.elegible_water list in mod_sample.fill_water(); coordinates and Ow-Hw-Hw topology are expanded in mod_construct_supercell_Y.get_coordinates_brick().",
            "target_water": "sample_Ca_Si_ratio sets N_water = round(N_Si * W_Si_ratio); v2.7 does not silently lower W/Si.",
            "attempts": attempts,
        },
        "candidates": candidates,
        "best_candidate": best["candidate"] if best is not None else None,
        "best_files": files,
        "gate_results": gate_results,
        "stage_summaries": stage_summaries,
        "final_classification": final_classification,
    }
    out = os.path.join(AUDIT, "v27_q2b_water_resampling_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print("Wrote {}".format(out))
    print("Final classification: {}".format(final_classification))


if __name__ == "__main__":
    main()

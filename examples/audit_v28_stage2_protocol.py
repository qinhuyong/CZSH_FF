"""v2.8 missing-angle diagnosis and Stage-2 protocol repair for best Q2b_Zn."""

from __future__ import print_function

import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from examples.audit_v25_water_q2b import parse_data, parse_lammps_log, top_force_atoms, water_audit, zinc_geometry


OUT = os.path.join(ROOT, "output_Y")
AUDIT = os.path.join(OUT, "v28_stage2_protocol")
BEST_DATA = os.path.join(OUT, "v27_q2b_best_cementff1_zn.data")
BEST_FF = os.path.join(OUT, "v27_q2b_best_in.CementFF4_Zn_1")
BEST_SUMMARY = os.path.join(OUT, "v27_q2b_best_zinc_summary_1.json")
V27_STAGE1 = os.path.join(OUT, "v27_q2b_water_resampling", "v27_stage1_water.data")
V27_STAGE2 = os.path.join(OUT, "v27_q2b_water_resampling", "v27_stage2_water_znoh.data")
V27_STAGE3_LOG = os.path.join(OUT, "v27_q2b_water_resampling", "v27_stage3_full_short.log")
LMP_EXE = r"C:\Program Files\LAMMPS 64-bit 4Feb2025-MSMPI\bin\lmp.exe"

ANGLE_LABELS = {
    1: "Hw-Ow-Hw",
    2: "O-Si-O / Oh-Si-O / Oh-Si-Oh",
    3: "Si-Oh-H",
    4: "O-Zn-O / Oh-Zn-O / Oh-Zn-Oh",
    5: "Zn-Oh-H",
}


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def atom_roles(data, zinc_summary, atom_id):
    atom = data["atoms"].get(atom_id)
    if atom is None:
        return ["missing"]
    roles = []
    t = atom["type"]
    if t in (5, 7):
        roles.append("water")
    if t == 9:
        roles.append("Zn motif")
    if t in (1, 2, 3, 4):
        roles.append("framework")
    if t == 3:
        roles.append("O core")
    if t == 4:
        roles.append("O shell")
    if t in (6, 8):
        roles.append("hydroxyl")
    if t == 1:
        roles.append("Ca framework")
    if t == 2:
        roles.append("Si framework")
    converted = set()
    for rec in zinc_summary.get("hydroxylation_records", []):
        for oxy in rec.get("hydroxylated_oxygens", []):
            converted.add(oxy.get("oxygen_atom_id"))
            converted.add(oxy.get("reused_shell_as_H_atom_id"))
    if atom_id in converted:
        roles.append("converted Zn hydroxyl")
    return roles


def bonds_angles_for(data, atom_ids):
    atom_ids = set(atom_ids)
    bonds = [b for b in data["bonds"] if atom_ids & {b["a1"], b["a2"]}]
    angles = [a for a in data["angles"] if atom_ids & {a["a1"], a["a2"], a["a3"]}]
    return bonds, angles


def pair_distances(data, ids):
    pairs = {}
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            pairs["{}-{}".format(ids[i], ids[j])] = distance(data, ids[i], ids[j])
    return pairs


def distance(data, id1, id2):
    from examples.audit_v25_water_q2b import dist

    if id1 not in data["atoms"] or id2 not in data["atoms"]:
        return None
    return dist(data, id1, id2)


def find_exact_angle(data, ids):
    target = tuple(ids)
    found = []
    for angle in data["angles"]:
        triplet = (angle["a1"], angle["a2"], angle["a3"])
        if triplet == target or triplet == tuple(reversed(target)):
            found.append(angle)
    return found


def mobile_ids_from_summary(data, summary, include_water=True, include_converted=True):
    ids = set()
    if include_water:
        ids.update(atom_id for atom_id, atom in data["atoms"].items() if atom["type"] in (5, 7))
    if include_converted:
        for site in summary.get("selected_sites", []):
            ids.add(site.get("atom_id"))
        for rec in summary.get("hydroxylation_records", []):
            ids.add(rec.get("zn_atom_id"))
            for oxy in rec.get("hydroxylated_oxygens", []):
                ids.add(oxy.get("oxygen_atom_id"))
                ids.add(oxy.get("reused_shell_as_H_atom_id"))
    return sorted(x for x in ids if x is not None)


def group_membership(ids, mobile_ids):
    mobile = set(mobile_ids)
    return {str(atom_id): ("mobile" if atom_id in mobile else "frozen") for atom_id in ids}


def topology_split_report(data, mobile_ids):
    mobile = set(mobile_ids)
    split_bonds = []
    split_angles = []
    for bond in data["bonds"]:
        states = [(bond["a1"] in mobile), (bond["a2"] in mobile)]
        if any(states) and not all(states):
            split_bonds.append(bond)
    for angle in data["angles"]:
        states = [(angle["a1"] in mobile), (angle["a2"] in mobile), (angle["a3"] in mobile)]
        if any(states) and not all(states):
            split_angles.append(angle)
    return {"split_bonds": split_bonds, "split_angles": split_angles}


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
        "comm_modify cutoff 14.0",
        "thermo 25",
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


def make_protocols(summary, stage1_data):
    data = parse_data(stage1_data)
    water_ids = mobile_ids_from_summary(data, summary, include_water=True, include_converted=False)
    motif_ids = mobile_ids_from_summary(data, summary, include_water=False, include_converted=True)
    precise_ids = sorted(set(water_ids + motif_ids))
    inputs = {}

    def group_line(ids):
        return "group mobile id " + " ".join(str(x) for x in ids)

    specs = {
        "A_precise_small_dt": {
            "ids": precise_ids,
            "data": stage1_data,
            "steps": 5000,
            "dt": "0.000025",
            "temp": "50.0",
            "out": "v28_A_precise_small_dt.data",
        },
        "B_water_only_1ps": {
            "ids": water_ids,
            "data": stage1_data,
            "steps": 10000,
            "dt": "0.0001",
            "temp": "100.0",
            "out": "v28_B_water_only_1ps.data",
        },
        "C_motif_only": {
            "ids": motif_ids,
            "data": stage1_data,
            "steps": 5000,
            "dt": "0.000025",
            "temp": "25.0",
            "out": "v28_C_motif_only.data",
        },
    }
    for name, spec in specs.items():
        path = os.path.join(AUDIT, name + ".in")
        lines = common_header(os.path.abspath(spec["data"]), os.path.abspath(BEST_FF))
        lines += [
            group_line(spec["ids"]),
            "group frozen subtract all mobile",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "compute tmob mobile temp",
            "velocity mobile create 1.0 13579 dist gaussian mom yes rot yes",
            "fix lang mobile langevin 1.0 {} 0.05 86420".format(spec["temp"]),
            "fix nvew mobile nve",
            "fix_modify lang temp tmob",
            "dump fdump all custom 25 {} id mol type q x y z fx fy fz".format(name + "_forces.lammpstrj"),
            "timestep {}".format(spec["dt"]),
            "run {}".format(spec["steps"]),
            "write_data {} nocoeff".format(spec["out"]),
        ]
        write_input(path, lines)
        inputs[name] = {"input": path, "ids": spec["ids"], "output": os.path.join(AUDIT, spec["out"])}

    # Stage-3 repair: start from successful precise Stage 2 and avoid immediately heating all core/shell atoms.
    stage3 = os.path.join(AUDIT, "E_stage3_water_motif_only.in")
    lines = common_header(os.path.join(AUDIT, "v28_A_precise_small_dt.data"), os.path.abspath(BEST_FF))
    lines += [
        group_line(precise_ids),
        "group frozen subtract all mobile",
        "fix freeze frozen setforce 0.0 0.0 0.0",
        "compute tmob mobile temp",
        "velocity mobile create 10.0 97531 dist gaussian mom yes rot yes",
        "fix lang mobile langevin 10.0 50.0 0.05 12345",
        "fix nvew mobile nve",
        "fix_modify lang temp tmob",
        "dump fdump all custom 5 E_stage3_water_motif_only_forces.lammpstrj id mol type q x y z fx fy fz",
        "timestep 0.000025",
        "run 5000",
        "write_data v28_E_stage3_water_motif_only.data nocoeff",
    ]
    write_input(stage3, lines)
    inputs["E_stage3_water_motif_only"] = {"input": stage3, "ids": precise_ids, "output": os.path.join(AUDIT, "v28_E_stage3_water_motif_only.data")}
    return inputs, {"water_ids": water_ids, "motif_ids": motif_ids, "precise_ids": precise_ids}


def run_protocols(inputs):
    results = {}
    for name, info in inputs.items():
        if name == "E_stage3_water_motif_only" and not os.path.exists(os.path.join(AUDIT, "v28_A_precise_small_dt.data")):
            results[name] = {"skipped": "A_precise_small_dt did not produce data"}
            continue
        result = run_lammps(info["input"], AUDIT)
        result["parsed"] = parse_lammps_log(result.get("log", ""))
        dump = os.path.join(AUDIT, name + "_forces.lammpstrj")
        result["top_forces"] = top_force_atoms(dump, parse_data(V27_STAGE1), 20)
        results[name] = result
    return results


def summarize_output(path):
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


def classify(results):
    a = results.get("A_precise_small_dt", {}).get("parsed", {})
    c = results.get("C_motif_only", {}).get("parsed", {})
    e = results.get("E_stage3_water_motif_only", {}).get("parsed", {})
    if a.get("errors") or a.get("lost"):
        return "needs_timestep_protocol_fix"
    if c.get("errors") or c.get("lost"):
        return "failed_converted_hydroxyl_topology"
    if e.get("errors") or e.get("lost"):
        return "needs_short_equilibration_test"
    return "needs_short_equilibration_test"


def main():
    ensure_dir(AUDIT)
    summary = json.load(open(BEST_SUMMARY))
    best = parse_data(BEST_DATA)
    stage1 = parse_data(V27_STAGE1)
    target_ids = [340, 339, 346]
    exact = find_exact_angle(best, target_ids)
    bonds, angles = bonds_angles_for(best, target_ids)
    old_mobile = [atom_id for atom_id, atom in stage1["atoms"].items() if atom["type"] in (5, 7, 6, 8, 9)]
    inputs, groups = make_protocols(summary, V27_STAGE1)
    precise_mobile = groups["precise_ids"]
    report = {
        "target_angle_query": {
            "requested_atoms": target_ids,
            "exact_angle_records": exact,
            "interpretation": "No exact angle 340-339-346 exists in the current v27 best data; this failure signature is stale for the current best candidate.",
            "atoms": {
                str(atom_id): {
                    "atom": best["atoms"].get(atom_id),
                    "roles": atom_roles(best, summary, atom_id),
                    "old_stage2_membership": "mobile" if atom_id in old_mobile else "frozen",
                    "precise_stage2_membership": "mobile" if atom_id in precise_mobile else "frozen",
                }
                for atom_id in target_ids
            },
            "pair_distances_initial": pair_distances(best, target_ids),
            "pair_distances_stage1": pair_distances(stage1, target_ids),
            "bonds_involving_atoms": bonds,
            "angles_involving_atoms": [
                dict(angle, label=ANGLE_LABELS.get(angle["type"], "unknown")) for angle in angles
            ],
        },
        "old_stage2_group": {
            "definition": "group mobile type 5 7 6 8 9",
            "mobile_count": len(old_mobile),
            "target_membership": group_membership(target_ids, old_mobile),
            "split_topology": {
                "n_split_bonds": len(topology_split_report(stage1, old_mobile)["split_bonds"]),
                "n_split_angles": len(topology_split_report(stage1, old_mobile)["split_angles"]),
                "examples_split_bonds": topology_split_report(stage1, old_mobile)["split_bonds"][:10],
                "examples_split_angles": topology_split_report(stage1, old_mobile)["split_angles"][:10],
            },
        },
        "precise_stage2_group": {
            "definition": "water atoms plus Zn and converted Oh/Hoh from zinc_summary only",
            "water_ids": groups["water_ids"],
            "motif_ids": groups["motif_ids"],
            "mobile_count": len(precise_mobile),
            "target_membership": group_membership(target_ids, precise_mobile),
            "split_topology": {
                "n_split_bonds": len(topology_split_report(stage1, precise_mobile)["split_bonds"]),
                "n_split_angles": len(topology_split_report(stage1, precise_mobile)["split_angles"]),
                "examples_split_bonds": topology_split_report(stage1, precise_mobile)["split_bonds"][:10],
                "examples_split_angles": topology_split_report(stage1, precise_mobile)["split_angles"][:10],
            },
        },
        "protocol_inputs": inputs,
    }
    results = run_protocols(inputs)
    report["protocol_results"] = results
    report["output_summaries"] = {
        name: summarize_output(info["output"])
        for name, info in inputs.items()
        if summarize_output(info["output"]) is not None
    }
    report["current_v27_stage3_failure"] = {
        "log": V27_STAGE3_LOG,
        "note": "Current v27 report fails at Stage 3 with missing core-shell bond atoms 300 301, not at Stage 2 with angle 340-339-346.",
    }
    report["final_classification"] = classify(results)
    out = os.path.join(AUDIT, "v28_stage2_protocol_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print("Wrote {}".format(out))
    print("Final classification: {}".format(report["final_classification"]))


if __name__ == "__main__":
    main()

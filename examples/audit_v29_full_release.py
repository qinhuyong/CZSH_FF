"""v2.9 full-system release diagnosis for best Q2b_Zn candidate."""

from __future__ import print_function

import json
import math
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from examples.audit_v25_water_q2b import parse_data, parse_lammps_log, top_force_atoms, water_audit, zinc_geometry, dist


OUT = os.path.join(ROOT, "output_Y")
AUDIT = os.path.join(OUT, "v29_full_release")
BEST_DATA = os.path.join(OUT, "v27_q2b_best_cementff1_zn.data")
BEST_FF = os.path.join(OUT, "v27_q2b_best_in.CementFF4_Zn_1")
BEST_SUMMARY = os.path.join(OUT, "v27_q2b_best_zinc_summary_1.json")
STAGE2_DATA = os.path.join(OUT, "v28_stage2_protocol", "v28_A_precise_small_dt.data")
PURE_DATA = os.path.join(OUT, "example_pure_csh_cementff1.data")
PURE_FF = os.path.join(OUT, "v22_audit", "generated_CementFF4_noZn_forcefield_only.in")
LMP_EXE = r"C:\Program Files\LAMMPS 64-bit 4Feb2025-MSMPI\bin\lmp.exe"


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def run_lammps(input_path, cwd, timeout=180):
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
        timeout=timeout,
    )
    return {"code": proc.returncode, "log": log_path, "stdout_tail": proc.stdout[-2000:]}


def write_input(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines))
        f.write("\n")


def common_header(data_file, ff_file, thermo=25):
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
        "thermo {}".format(thermo),
        "thermo_style custom step pe ebond eangle evdwl ecoul elong temp fnorm fmax press",
    ]


def atom_info(data, atom_id, summary=None):
    atom = data["atoms"].get(atom_id)
    if atom is None:
        return None
    roles = []
    if atom["type"] == 1:
        roles.append("Ca framework")
    elif atom["type"] == 2:
        roles.append("Si framework")
    elif atom["type"] == 3:
        roles.append("O core/framework")
    elif atom["type"] == 4:
        roles.append("O shell/core-shell")
    elif atom["type"] in (5, 7):
        roles.append("water")
    elif atom["type"] in (6, 8):
        roles.append("hydroxyl")
    elif atom["type"] == 9:
        roles.append("Zn motif")
    converted = set()
    if summary:
        for rec in summary.get("hydroxylation_records", []):
            for oxy in rec.get("hydroxylated_oxygens", []):
                converted.add(oxy.get("oxygen_atom_id"))
                converted.add(oxy.get("reused_shell_as_H_atom_id"))
    if atom_id in converted:
        roles.append("converted Zn hydroxyl")
    return {"atom": atom, "roles": roles}


def bonds_angles_for(data, atom_ids):
    atom_ids = set(atom_ids)
    bonds = [b for b in data["bonds"] if atom_ids & {b["a1"], b["a2"]}]
    angles = [a for a in data["angles"] if atom_ids & {a["a1"], a["a2"], a["a3"]}]
    return bonds, angles


def mobile_membership(atom_id):
    return {
        "stage1_water_only": "mobile" if atom_id in WATER_IDS else "frozen",
        "stage2_precise": "mobile" if atom_id in PRECISE_IDS else "frozen",
        "stage3_full": "mobile",
    }


def ids_from_summary(data, summary):
    water = sorted(atom_id for atom_id, atom in data["atoms"].items() if atom["type"] in (5, 7))
    motif = set()
    for site in summary.get("selected_sites", []):
        motif.add(site.get("atom_id"))
    for rec in summary.get("hydroxylation_records", []):
        motif.add(rec.get("zn_atom_id"))
        for oxy in rec.get("hydroxylated_oxygens", []):
            motif.add(oxy.get("oxygen_atom_id"))
            motif.add(oxy.get("reused_shell_as_H_atom_id"))
    motif = sorted(x for x in motif if x is not None)
    return water, motif, sorted(set(water + motif))


def shell_distance_summary(data):
    vals = []
    for bond in data["bonds"]:
        if bond["type"] != 1:
            continue
        if bond["a1"] in data["atoms"] and bond["a2"] in data["atoms"]:
            vals.append(dist(data, bond["a1"], bond["a2"]))
    return {
        "count": len(vals),
        "min": min(vals) if vals else None,
        "mean": sum(vals) / len(vals) if vals else None,
        "max": max(vals) if vals else None,
    }


def summarize_data(path):
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
        "core_shell": shell_distance_summary(data),
        "zinc_geometry": zinc_geometry(data),
    }


def last_frame_displacements(dump_path, data, atom_ids):
    if not os.path.exists(dump_path):
        return {}
    frames = []
    with open(dump_path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(float(lines[i + 1].strip()))
            n = int(lines[i + 3].strip())
            cols = lines[i + 8].split()[2:]
            atoms = {}
            for row in lines[i + 9:i + 9 + n]:
                p = row.split()
                atom = {cols[j]: float(p[j]) if cols[j] not in ("id", "mol", "type") else int(p[j]) for j in range(len(cols))}
                if atom["id"] in atom_ids:
                    atoms[atom["id"]] = atom
            frames.append({"step": step, "atoms": atoms})
            i += 9 + n
        else:
            i += 1
    if not frames:
        return {}
    out = {"last_step": frames[-1]["step"], "atoms": {}}
    for atom_id in atom_ids:
        if atom_id not in frames[-1]["atoms"]:
            continue
        last = frames[-1]["atoms"][atom_id]
        prev = None
        for frame in reversed(frames[:-1]):
            if atom_id in frame["atoms"]:
                prev = frame["atoms"][atom_id]
                break
        initial = data["atoms"].get(atom_id)
        disp_from_initial = None
        disp_last_step = None
        if initial:
            disp_from_initial = math.sqrt((last["x"] - initial["x"]) ** 2 + (last["y"] - initial["y"]) ** 2 + (last["z"] - initial["z"]) ** 2)
        if prev:
            disp_last_step = math.sqrt((last["x"] - prev["x"]) ** 2 + (last["y"] - prev["y"]) ** 2 + (last["z"] - prev["z"]) ** 2)
        out["atoms"][str(atom_id)] = {
            "last": last,
            "disp_from_initial_raw": disp_from_initial,
            "disp_from_previous_dump_raw": disp_last_step,
        }
    return out


def make_protocols():
    protocols = {}
    # failing reference with dump every step
    path = os.path.join(AUDIT, "Q_ref_full_release.in")
    lines = common_header(os.path.abspath(STAGE2_DATA), os.path.abspath(BEST_FF), thermo=1)
    lines += [
        "compute tall all temp",
        "velocity all create 50.0 97531 dist gaussian mom yes rot yes",
        "fix lang all langevin 50.0 100.0 0.02 12345",
        "fix nvea all nve",
        "fix_modify lang temp tall",
        "dump fdump all custom 1 Q_ref_full_release.lammpstrj id mol type q x y z fx fy fz",
        "timestep 0.0001",
        "run 200",
        "write_data Q_ref_full_release.data nocoeff",
    ]
    write_input(path, lines)
    protocols["Q_ref_full_release"] = {"input": path, "data": os.path.join(AUDIT, "Q_ref_full_release.data"), "dump": os.path.join(AUDIT, "Q_ref_full_release.lammpstrj")}

    # pure C-S-H same protocol
    path = os.path.join(AUDIT, "P_pure_same_full_release.in")
    lines = common_header(os.path.abspath(PURE_DATA), os.path.abspath(PURE_FF), thermo=1)
    lines += [
        "compute tall all temp",
        "velocity all create 50.0 97531 dist gaussian mom yes rot yes",
        "fix lang all langevin 50.0 100.0 0.02 12345",
        "fix nvea all nve",
        "fix_modify lang temp tall",
        "dump fdump all custom 5 P_pure_same_full_release.lammpstrj id mol type q x y z fx fy fz",
        "timestep 0.0001",
        "run 200",
        "write_data P_pure_same_full_release.data nocoeff",
    ]
    write_input(path, lines)
    protocols["P_pure_same_full_release"] = {"input": path, "data": os.path.join(AUDIT, "P_pure_same_full_release.data"), "dump": os.path.join(AUDIT, "P_pure_same_full_release.lammpstrj")}

    specs = {
        "A_lower_timestep": {"dt": "0.000025", "steps": "5000", "temp0": "10.0", "temp1": "50.0"},
        "B_low_temp": {"dt": "0.000025", "steps": "5000", "temp0": "1.0", "temp1": "10.0"},
        "D_damped": {"dt": "0.000025", "steps": "5000", "temp0": "1.0", "temp1": "25.0"},
    }
    for name, spec in specs.items():
        path = os.path.join(AUDIT, name + ".in")
        lines = common_header(os.path.abspath(STAGE2_DATA), os.path.abspath(BEST_FF), thermo=100)
        lines += [
            "compute tall all temp",
            "velocity all create {} 13579 dist gaussian mom yes rot yes".format(spec["temp0"]),
            "fix lang all langevin {} {} 0.10 86420".format(spec["temp0"], spec["temp1"]),
            "fix nvea all nve",
            "fix_modify lang temp tall",
            "dump fdump all custom 25 {} id mol type q x y z fx fy fz".format(name + ".lammpstrj"),
            "timestep {}".format(spec["dt"]),
            "run {}".format(spec["steps"]),
            "write_data {} nocoeff".format(name + ".data"),
        ]
        write_input(path, lines)
        protocols[name] = {"input": path, "data": os.path.join(AUDIT, name + ".data"), "dump": os.path.join(AUDIT, name + ".lammpstrj")}

    # staged release: water + motif from Stage2, not full core-shell release
    path = os.path.join(AUDIT, "C_staged_release_no_framework.in")
    lines = common_header(os.path.abspath(STAGE2_DATA), os.path.abspath(BEST_FF), thermo=100)
    lines += [
        "group mobile id " + " ".join(str(x) for x in PRECISE_IDS),
        "group frozen subtract all mobile",
        "fix freeze frozen setforce 0.0 0.0 0.0",
        "compute tmob mobile temp",
        "velocity mobile create 25.0 24680 dist gaussian mom yes rot yes",
        "fix lang mobile langevin 25.0 100.0 0.10 86420",
        "fix nvea mobile nve",
        "fix_modify lang temp tmob",
        "dump fdump all custom 25 C_staged_release_no_framework.lammpstrj id mol type q x y z fx fy fz",
        "timestep 0.000025",
        "run 5000",
        "write_data C_staged_release_no_framework.data nocoeff",
    ]
    write_input(path, lines)
    protocols["C_staged_release_no_framework"] = {"input": path, "data": os.path.join(AUDIT, "C_staged_release_no_framework.data"), "dump": os.path.join(AUDIT, "C_staged_release_no_framework.lammpstrj")}

    # conservative minimization after Stage2
    path = os.path.join(AUDIT, "E_limited_minimize.in")
    lines = common_header(os.path.abspath(STAGE2_DATA), os.path.abspath(BEST_FF), thermo=10)
    lines += [
        "dump fdump all custom 10 E_limited_minimize.lammpstrj id mol type q x y z fx fy fz",
        "min_style cg",
        "min_modify dmax 0.001 line quadratic",
        "minimize 1e-6 1e-8 200 2000",
        "write_data E_limited_minimize.data nocoeff",
    ]
    write_input(path, lines)
    protocols["E_limited_minimize"] = {"input": path, "data": os.path.join(AUDIT, "E_limited_minimize.data"), "dump": os.path.join(AUDIT, "E_limited_minimize.lammpstrj")}
    return protocols


def run_protocols(protocols):
    results = {}
    for name, info in protocols.items():
        result = run_lammps(info["input"], AUDIT, timeout=180)
        result["parsed"] = parse_lammps_log(result.get("log", ""))
        base = parse_data(STAGE2_DATA if not name.startswith("P_") else PURE_DATA)
        result["top_forces"] = top_force_atoms(info["dump"], base, 20)
        result["summary"] = summarize_data(info["data"])
        if name == "Q_ref_full_release":
            result["tracked_displacements"] = last_frame_displacements(info["dump"], parse_data(STAGE2_DATA), [300, 301, 95, 96, 97, 98, 99])
        results[name] = result
    return results


def classify(results):
    pure = results.get("P_pure_same_full_release", {}).get("parsed", {})
    q_ref = results.get("Q_ref_full_release", {}).get("parsed", {})
    stable = []
    for name in ["A_lower_timestep", "B_low_temp", "C_staged_release_no_framework", "D_damped", "E_limited_minimize"]:
        parsed = results.get(name, {}).get("parsed", {})
        if results.get(name, {}).get("code") == 0 and not parsed.get("errors") and not parsed.get("lost"):
            stable.append(name)
    if pure.get("errors") or pure.get("lost"):
        return "failed_fullsystem_protocol"
    if q_ref.get("errors") and stable:
        return "needs_extended_short_equilibration"
    if stable:
        return "needs_extended_short_equilibration"
    return "needs_fullsystem_protocol_fix"


def main():
    ensure_dir(AUDIT)
    summary = json.load(open(BEST_SUMMARY))
    stage2 = parse_data(STAGE2_DATA)
    best = parse_data(BEST_DATA)
    global WATER_IDS, MOTIF_IDS, PRECISE_IDS
    WATER_IDS = sorted(atom_id for atom_id, atom in stage2["atoms"].items() if atom["type"] in (5, 7))
    motif = set()
    for site in summary.get("selected_sites", []):
        motif.add(site.get("atom_id"))
    for rec in summary.get("hydroxylation_records", []):
        motif.add(rec.get("zn_atom_id"))
        for oxy in rec.get("hydroxylated_oxygens", []):
            motif.add(oxy.get("oxygen_atom_id"))
            motif.add(oxy.get("reused_shell_as_H_atom_id"))
    MOTIF_IDS = sorted(x for x in motif if x is not None)
    PRECISE_IDS = sorted(set(WATER_IDS + MOTIF_IDS))

    bonds, angles = bonds_angles_for(stage2, [300, 301])
    bond_300_301 = [b for b in stage2["bonds"] if {b["a1"], b["a2"]} == {300, 301}]
    protocols = make_protocols()
    results = run_protocols(protocols)
    report = {
        "bond_300_301_identity": {
            "atoms": {
                "300": atom_info(stage2, 300, summary),
                "301": atom_info(stage2, 301, summary),
            },
            "bond_records": bond_300_301,
            "bond_label": "O core-shell" if bond_300_301 and bond_300_301[0]["type"] == 1 else "unknown",
            "bonds_involving_atoms": bonds,
            "angles_involving_atoms": angles,
            "distance_to_Zn": {
                "300": dist(stage2, 300, summary["selected_sites"][0]["atom_id"]),
                "301": dist(stage2, 301, summary["selected_sites"][0]["atom_id"]),
            },
            "bond_length_before_stage3": dist(stage2, 300, 301),
            "membership": {
                "300": mobile_membership(300),
                "301": mobile_membership(301),
            },
        },
        "groups": {
            "water_ids": WATER_IDS,
            "motif_ids": MOTIF_IDS,
            "precise_ids": PRECISE_IDS,
        },
        "protocol_inputs": protocols,
        "protocol_results": results,
        "stage2_input_summary": summarize_data(STAGE2_DATA),
    }
    report["final_classification"] = classify(results)
    out = os.path.join(AUDIT, "v29_full_release_report.json")
    with open(out, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print("Wrote {}".format(out))
    print("Final classification: {}".format(report["final_classification"]))


if __name__ == "__main__":
    WATER_IDS = []
    MOTIF_IDS = []
    PRECISE_IDS = []
    main()

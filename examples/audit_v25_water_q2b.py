"""v2.5 water/TIP4P audit and short Q2b_Zn relaxation diagnostics.

This script intentionally does not generate new Zn motifs.  It reads the
existing pure C-S-H and Q2b_Zn CementFF4 data files, audits water topology and
contacts, then runs short diagnostic LAMMPS protocols focused on water.
"""

from __future__ import print_function

import json
import math
import os
import re
import shutil
import subprocess


ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "output_Y")
AUDIT = os.path.join(OUT, "v25_water_audit")

PURE_DATA = os.path.join(OUT, "example_pure_csh_cementff1.data")
Q2B_DATA = os.path.join(OUT, "example_q2b_zn_cementff1_zn.data")
Q2B_FF = os.path.join(OUT, "example_q2b_zn_in.CementFF4_Zn_1")
Q2B_SUMMARY = os.path.join(OUT, "example_q2b_zn_zinc_summary_1.json")
V24_FORCE_DUMP = os.path.join(OUT, "v24_audit", "q2b_force_cg_forces.lammpstrj")

LMP_EXE = r"C:\Program Files\LAMMPS 64-bit 4Feb2025-MSMPI\bin\lmp.exe"

LABELS = {
    1: "Ca",
    2: "Si",
    3: "O_core",
    4: "O_shell",
    5: "Ow",
    6: "Oh",
    7: "Hw",
    8: "Hoh",
    9: "Zn",
}


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


def parse_data(path):
    data = {
        "path": path,
        "counts": {},
        "box": {},
        "masses": {},
        "atoms": {},
        "bonds": [],
        "angles": [],
    }
    section = None
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            m = re.match(r"^(\d+)\s+(atoms|bonds|angles|atom types|bond types|angle types)$", line)
            if m:
                data["counts"][m.group(2)] = int(m.group(1))
                continue
            parts = line.split()
            if len(parts) >= 4 and parts[2:4] == ["xlo", "xhi"]:
                data["box"]["xlo"], data["box"]["xhi"] = float(parts[0]), float(parts[1])
                continue
            if len(parts) >= 4 and parts[2:4] == ["ylo", "yhi"]:
                data["box"]["ylo"], data["box"]["yhi"] = float(parts[0]), float(parts[1])
                continue
            if len(parts) >= 4 and parts[2:4] == ["zlo", "zhi"]:
                data["box"]["zlo"], data["box"]["zhi"] = float(parts[0]), float(parts[1])
                continue
            if len(parts) >= 6 and parts[3:6] == ["xy", "xz", "yz"]:
                data["box"]["xy"], data["box"]["xz"], data["box"]["yz"] = (
                    float(parts[0]),
                    float(parts[1]),
                    float(parts[2]),
                )
                continue
            if line.startswith("Masses"):
                section = "Masses"
                continue
            if line.startswith("Atoms"):
                section = "Atoms"
                continue
            if line.startswith("Bonds"):
                section = "Bonds"
                continue
            if line.startswith("Angles"):
                section = "Angles"
                continue
            if line.startswith("Velocities"):
                section = "Velocities"
                continue
            if section == "Masses":
                p = line.split("#")[0].split()
                if len(p) >= 2:
                    data["masses"][int(p[0])] = float(p[1])
            elif section == "Atoms":
                p = line.split("#")[0].split()
                if len(p) >= 7:
                    atom_id = int(p[0])
                    atom_type = int(p[2])
                    data["atoms"][atom_id] = {
                        "id": atom_id,
                        "mol": int(p[1]),
                        "type": atom_type,
                        "label": LABELS.get(atom_type, "type{}".format(atom_type)),
                        "q": float(p[3]),
                        "x": float(p[4]),
                        "y": float(p[5]),
                        "z": float(p[6]),
                    }
            elif section == "Bonds":
                p = line.split("#")[0].split()
                if len(p) >= 4:
                    data["bonds"].append(
                        {"id": int(p[0]), "type": int(p[1]), "a1": int(p[2]), "a2": int(p[3])}
                    )
            elif section == "Angles":
                p = line.split("#")[0].split()
                if len(p) >= 5:
                    data["angles"].append(
                        {
                            "id": int(p[0]),
                            "type": int(p[1]),
                            "a1": int(p[2]),
                            "a2": int(p[3]),
                            "a3": int(p[4]),
                        }
                    )
    for k in ("xy", "xz", "yz"):
        data["box"].setdefault(k, 0.0)
    return data


def box_matrix(data):
    b = data["box"]
    lx = b["xhi"] - b["xlo"]
    ly = b["yhi"] - b["ylo"]
    lz = b["zhi"] - b["zlo"]
    return (
        (lx, b["xy"], b["xz"]),
        (0.0, ly, b["yz"]),
        (0.0, 0.0, lz),
    )


def inv_upper(m):
    a, b, c = m[0]
    _, d, e = m[1]
    _, _, f = m[2]
    return (
        (1.0 / a, -b / (a * d), (b * e - c * d) / (a * d * f)),
        (0.0, 1.0 / d, -e / (d * f)),
        (0.0, 0.0, 1.0 / f),
    )


def mat_vec(m, v):
    return (
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    )


def pbc_delta(data, a, b):
    m = box_matrix(data)
    inv = inv_upper(m)
    dr = (a["x"] - b["x"], a["y"] - b["y"], a["z"] - b["z"])
    frac = list(mat_vec(inv, dr))
    frac = [x - round(x) for x in frac]
    return mat_vec(m, frac)


def dist(data, id1, id2):
    d = pbc_delta(data, data["atoms"][id1], data["atoms"][id2])
    return math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])


def angle_deg(data, i, j, k):
    v1 = pbc_delta(data, data["atoms"][i], data["atoms"][j])
    v2 = pbc_delta(data, data["atoms"][k], data["atoms"][j])
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0.0 or n2 == 0.0:
        return None
    cosang = sum(v1[n] * v2[n] for n in range(3)) / (n1 * n2)
    cosang = max(-1.0, min(1.0, cosang))
    return math.degrees(math.acos(cosang))


def build_bond_map(data):
    by_atom = {}
    for bond in data["bonds"]:
        by_atom.setdefault(bond["a1"], []).append(bond)
        by_atom.setdefault(bond["a2"], []).append(bond)
    return by_atom


def build_angle_map(data):
    by_atom = {}
    for angle in data["angles"]:
        for atom_id in (angle["a1"], angle["a2"], angle["a3"]):
            by_atom.setdefault(atom_id, []).append(angle)
    return by_atom


def water_molecules(data):
    bond_map = build_bond_map(data)
    angles = data["angles"]
    waters = []
    for ow_id, atom in sorted(data["atoms"].items()):
        if atom["type"] != 5:
            continue
        hw_bonds = []
        for bond in bond_map.get(ow_id, []):
            other = bond["a2"] if bond["a1"] == ow_id else bond["a1"]
            if data["atoms"][other]["type"] == 7:
                hw_bonds.append((other, bond))
        hw_ids = sorted(x[0] for x in hw_bonds)
        water_angles = [
            a for a in angles
            if a["a2"] == ow_id and set((a["a1"], a["a3"])) == set(hw_ids)
        ]
        oh = [dist(data, ow_id, h) for h in hw_ids]
        hh = dist(data, hw_ids[0], hw_ids[1]) if len(hw_ids) == 2 else None
        hoh = angle_deg(data, hw_ids[0], ow_id, hw_ids[1]) if len(hw_ids) == 2 else None
        waters.append(
            {
                "Ow": ow_id,
                "Hw": hw_ids,
                "molecule_id": atom["mol"],
                "bond_types": [b["type"] for _, b in sorted(hw_bonds)],
                "angle_ids": [a["id"] for a in water_angles],
                "angle_types": [a["type"] for a in water_angles],
                "OH_lengths": oh,
                "HH_distance": hh,
                "HOH_angle": hoh,
                "has_two_H": len(hw_ids) == 2,
                "has_tip4p_bonds": len(hw_bonds) == 2 and all(b["type"] == 2 for _, b in hw_bonds),
                "has_tip4p_angle": len(water_angles) == 1 and water_angles[0]["type"] == 1,
                "consecutive_ids": len(hw_ids) == 2 and hw_ids == [ow_id + 1, ow_id + 2],
                "tip4p_compatible": (
                    len(hw_ids) == 2
                    and len(hw_bonds) == 2
                    and all(b["type"] == 2 for _, b in hw_bonds)
                    and len(water_angles) == 1
                    and water_angles[0]["type"] == 1
                ),
            }
        )
    return waters


def min_distance_to_types(data, atom_id, target_types, exclude=None):
    exclude = exclude or set()
    best = None
    for other_id, other in data["atoms"].items():
        if other_id == atom_id or other_id in exclude or other["type"] not in target_types:
            continue
        d = dist(data, atom_id, other_id)
        item = {"distance": d, "atom_id": other_id, "type": other["type"], "label": other["label"]}
        if best is None or d < best["distance"]:
            best = item
    return best


def water_contacts(data, waters):
    by_ow = {}
    for w in waters:
        water_ids = set([w["Ow"]] + w["Hw"])
        ow = w["Ow"]
        hw_contacts = []
        for h in w["Hw"]:
            hw_contacts.append(
                {
                    "H": h,
                    "H_O_nonbonded": min_distance_to_types(data, h, {3, 4, 5, 6}, exclude=water_ids),
                    "H_H_nonbonded": min_distance_to_types(data, h, {7, 8}, exclude=water_ids),
                    "H_Ca": min_distance_to_types(data, h, {1}, exclude=water_ids),
                    "H_Si": min_distance_to_types(data, h, {2}, exclude=water_ids),
                    "H_Zn": min_distance_to_types(data, h, {9}, exclude=water_ids),
                }
            )
        by_ow[ow] = {
            "Ow": ow,
            "Hw": w["Hw"],
            "Ow_Ca": min_distance_to_types(data, ow, {1}, exclude=water_ids),
            "Ow_Si": min_distance_to_types(data, ow, {2}, exclude=water_ids),
            "Ow_O": min_distance_to_types(data, ow, {3, 4, 5, 6}, exclude=water_ids),
            "Ow_Zn": min_distance_to_types(data, ow, {9}, exclude=water_ids),
            "Hw_contacts": hw_contacts,
        }
    return by_ow


def water_audit(data, highlighted=None):
    highlighted = highlighted or []
    waters = water_molecules(data)
    contacts = water_contacts(data, waters)
    bad = [
        w for w in waters
        if not w["has_two_H"] or not w["has_tip4p_bonds"] or not w["has_tip4p_angle"]
    ]
    oh_lengths = [x for w in waters for x in w["OH_lengths"]]
    hh = [w["HH_distance"] for w in waters if w["HH_distance"] is not None]
    hoh = [w["HOH_angle"] for w in waters if w["HOH_angle"] is not None]
    by_ow = {w["Ow"]: w for w in waters}
    return {
        "file": data["path"],
        "n_water": len(waters),
        "n_tip4p_compatible": sum(1 for w in waters if w["tip4p_compatible"]),
        "n_bad_water": len(bad),
        "bad_water": bad,
        "all_water_summary": {
            "OH_min": min(oh_lengths) if oh_lengths else None,
            "OH_max": max(oh_lengths) if oh_lengths else None,
            "HH_min": min(hh) if hh else None,
            "HH_max": max(hh) if hh else None,
            "HOH_min": min(hoh) if hoh else None,
            "HOH_max": max(hoh) if hoh else None,
        },
        "waters": waters,
        "highlighted": {
            str(ow): {"topology": by_ow.get(ow), "contacts": contacts.get(ow)}
            for ow in highlighted
        },
        "contacts": contacts,
    }


def parse_lammps_log(path):
    if not os.path.exists(path):
        return {"exists": False}
    columns = None
    rows = []
    stops = []
    errors = []
    lost = False
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("ERROR"):
                errors.append(line)
            if "Lost atoms" in line:
                lost = True
            if "Stopping criterion =" in line:
                stops.append(line.split("=", 1)[1].strip())
            if re.match(r"^(Step|Time)\s+", line):
                columns = line.split()
                continue
            if columns:
                p = line.split()
                if len(p) == len(columns):
                    try:
                        rows.append({columns[i].lower(): float(p[i]) for i in range(len(columns))})
                    except ValueError:
                        pass
    return {
        "exists": True,
        "errors": errors,
        "lost": lost,
        "nrows": len(rows),
        "initial": rows[0] if rows else None,
        "final": rows[-1] if rows else None,
        "stops": stops,
    }


def parse_last_dump_frame(path):
    if not os.path.exists(path):
        return None
    last = None
    with open(path, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            timestep = int(float(lines[i + 1].strip()))
            n = int(lines[i + 3].strip())
            header_i = i + 8
            columns = lines[header_i].strip().split()[2:]
            atoms = []
            for row in lines[header_i + 1: header_i + 1 + n]:
                p = row.split()
                atoms.append({columns[j]: float(p[j]) if columns[j] not in ("id", "mol", "type") else int(p[j]) for j in range(len(columns))})
            last = {"timestep": timestep, "atoms": atoms}
            i = header_i + 1 + n
        else:
            i += 1
    return last


def dump_atom_dict(frame):
    if frame is None:
        return {}
    return {
        atom["id"]: atom
        for atom in frame["atoms"]
    }


def frame_distance(data, frame_atoms, id1, id2):
    if id1 not in frame_atoms or id2 not in frame_atoms:
        return dist(data, id1, id2)
    a0 = data["atoms"][id1].copy()
    b0 = data["atoms"][id2].copy()
    a = frame_atoms[id1]
    b = frame_atoms[id2]
    a0.update({"x": a["x"], "y": a["y"], "z": a["z"]})
    b0.update({"x": b["x"], "y": b["y"], "z": b["z"]})
    tmp = {"atoms": {id1: a0, id2: b0}, "box": data["box"]}
    return dist(tmp, id1, id2)


def top_force_atoms(dump_path, data, n=20):
    frame = parse_last_dump_frame(dump_path)
    if frame is None:
        return []
    frame_atoms = dump_atom_dict(frame)
    out = []
    for atom in frame["atoms"]:
        fmag = math.sqrt(atom["fx"] ** 2 + atom["fy"] ** 2 + atom["fz"] ** 2)
        atom_id = atom["id"]
        base = data["atoms"].get(atom_id, {})
        entry = {
            "id": atom_id,
            "type": atom["type"],
            "label": LABELS.get(atom["type"], str(atom["type"])),
            "mol": atom.get("mol"),
            "force": fmag,
            "fx": atom["fx"],
            "fy": atom["fy"],
            "fz": atom["fz"],
            "distance_to_Zn": None,
        }
        zn_ids = [i for i, a in data["atoms"].items() if a["type"] == 9]
        if zn_ids and atom_id in data["atoms"]:
            entry["distance_to_Zn"] = frame_distance(data, frame_atoms, atom_id, zn_ids[0])
        neighbors = []
        if atom_id in data["atoms"]:
            for other_id in data["atoms"]:
                if other_id == atom_id:
                    continue
                neighbors.append((frame_distance(data, frame_atoms, atom_id, other_id), other_id))
            neighbors.sort()
            entry["nearest_neighbors"] = [
                {
                    "distance": d,
                    "atom_id": oid,
                    "type": data["atoms"][oid]["type"],
                    "label": data["atoms"][oid]["label"],
                }
                for d, oid in neighbors[:8]
            ]
        out.append(entry)
    out.sort(key=lambda x: x["force"], reverse=True)
    return out[:n]


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
        timeout=120,
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
        "thermo 50",
        "thermo_style custom step pe ebond eangle evdwl ecoul elong temp fnorm fmax press",
    ]


def make_protocol_inputs():
    q2b_data = os.path.abspath(Q2B_DATA)
    q2b_ff = os.path.abspath(Q2B_FF)
    inputs = {}

    a = os.path.join(AUDIT, "q2b_protocol_A_water_only_min.in")
    write_input(
        a,
        common_header(q2b_data, q2b_ff)
        + [
            "group mobile_water type 5 7",
            "group frozen subtract all mobile_water",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "dump fdump all custom 20 q2b_protocol_A_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "min_style cg",
            "min_modify dmax 0.01 line quadratic",
            "minimize 1e-6 1e-8 200 2000",
            "write_data q2b_protocol_A_water_only_min.data nocoeff",
        ],
    )
    inputs["A_water_only_min"] = a

    b = os.path.join(AUDIT, "q2b_protocol_B_water_only_langevin.in")
    write_input(
        b,
        common_header(q2b_data, q2b_ff)
        + [
            "group mobile_water type 5 7",
            "group frozen subtract all mobile_water",
            "fix freeze frozen setforce 0.0 0.0 0.0",
            "compute twater mobile_water temp",
            "velocity mobile_water create 50.0 246813 dist gaussian mom yes rot yes",
            "fix lang mobile_water langevin 50.0 100.0 0.02 97531",
            "fix nvew mobile_water nve",
            "fix_modify lang temp twater",
            "dump fdump all custom 100 q2b_protocol_B_forces.lammpstrj id mol type q x y z fx fy fz",
            "timestep 0.0001",
            "run 10000",
            "write_data q2b_protocol_B_water_relaxed.data nocoeff",
        ],
    )
    inputs["B_water_only_langevin"] = b

    c = os.path.join(AUDIT, "q2b_protocol_C_full_after_water.in")
    write_input(
        c,
        common_header(os.path.join(AUDIT, "q2b_protocol_B_water_relaxed.data"), q2b_ff)
        + [
            "dump fdump all custom 20 q2b_protocol_C_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "min_style cg",
            "min_modify dmax 0.01 line quadratic",
            "minimize 1e-6 1e-8 300 3000",
            "write_data q2b_protocol_C_full_min.data nocoeff",
        ],
    )
    inputs["C_full_after_water"] = c

    freeze = os.path.join(AUDIT, "q2b_diag_freeze_water_min.in")
    write_input(
        freeze,
        common_header(q2b_data, q2b_ff)
        + [
            "group mobile_water type 5 7",
            "fix freeze mobile_water setforce 0.0 0.0 0.0",
            "dump fdump all custom 20 q2b_diag_freeze_water_forces.lammpstrj id mol type q x y z fx fy fz",
            "run 0",
            "min_style cg",
            "min_modify dmax 0.01 line quadratic",
            "minimize 1e-6 1e-8 300 3000",
            "write_data q2b_diag_freeze_water_min.data nocoeff",
        ],
    )
    inputs["diag_freeze_water"] = freeze
    return inputs


def run_protocols(inputs):
    results = {}
    for name, path in inputs.items():
        if name == "C_full_after_water":
            prereq = os.path.join(AUDIT, "q2b_protocol_B_water_relaxed.data")
            if not os.path.exists(prereq):
                results[name] = {"skipped": "Protocol B did not produce water-relaxed data"}
                continue
        result = run_lammps(path, AUDIT)
        result["parsed"] = parse_lammps_log(result.get("log", ""))
        dump_name = {
            "A_water_only_min": "q2b_protocol_A_forces.lammpstrj",
            "B_water_only_langevin": "q2b_protocol_B_forces.lammpstrj",
            "C_full_after_water": "q2b_protocol_C_forces.lammpstrj",
            "diag_freeze_water": "q2b_diag_freeze_water_forces.lammpstrj",
        }.get(name)
        if dump_name:
            result["top_forces"] = top_force_atoms(os.path.join(AUDIT, dump_name), parse_data(Q2B_DATA), 20)
        results[name] = result
    return results


def zinc_geometry(data):
    zn_ids = [i for i, a in data["atoms"].items() if a["type"] == 9]
    if not zn_ids:
        return {}
    zn = zn_ids[0]
    # Treat type 4 as the shell of an O(S) core, not as an independent
    # chemically coordinating oxygen.  Counting both core and shell produces
    # artificial duplicate Zn-O contacts and near-zero O-Zn-O angles.
    oxygen_ids = [i for i, a in data["atoms"].items() if a["type"] in (3, 5, 6)]
    distances = sorted((dist(data, zn, i), i, data["atoms"][i]["type"], data["atoms"][i]["label"]) for i in oxygen_ids)
    near4 = distances[:4]
    near_by_cutoff = {
        str(c): sum(1 for d, _, _, _ in distances if d <= c)
        for c in (2.1, 2.3, 2.5)
    }
    angles = []
    for i in range(len(near4)):
        for j in range(i + 1, len(near4)):
            angles.append(angle_deg(data, near4[i][1], zn, near4[j][1]))
    oh_bonds = []
    for bond in data["bonds"]:
        t1 = data["atoms"][bond["a1"]]["type"]
        t2 = data["atoms"][bond["a2"]]["type"]
        if set((t1, t2)) == set((6, 8)):
            oh_bonds.append({"atoms": [bond["a1"], bond["a2"]], "type": bond["type"], "distance": dist(data, bond["a1"], bond["a2"])})
    return {
        "Zn_id": zn,
        "coordination_by_cutoff": near_by_cutoff,
        "nearest_Zn_O": [
            {"distance": d, "atom_id": i, "type": t, "label": lab}
            for d, i, t, lab in near4
        ],
        "O_Zn_O_angle_min": min(angles) if angles else None,
        "O_Zn_O_angle_max": max(angles) if angles else None,
        "OH_bond_min": min(x["distance"] for x in oh_bonds) if oh_bonds else None,
        "OH_bond_max": max(x["distance"] for x in oh_bonds) if oh_bonds else None,
    }


def summarize_protocol_data(results):
    out = {}
    candidates = {
        "initial": Q2B_DATA,
        "A_water_only_min": os.path.join(AUDIT, "q2b_protocol_A_water_only_min.data"),
        "B_water_only_langevin": os.path.join(AUDIT, "q2b_protocol_B_water_relaxed.data"),
        "C_full_after_water": os.path.join(AUDIT, "q2b_protocol_C_full_min.data"),
        "diag_freeze_water": os.path.join(AUDIT, "q2b_diag_freeze_water_min.data"),
    }
    for name, path in candidates.items():
        if os.path.exists(path):
            data = parse_data(path)
            wa = water_audit(data, highlighted=[268, 324])
            out[name] = {
                "data_file": path,
                "water_summary": wa["all_water_summary"],
                "water_bad_count": wa["n_bad_water"],
                "zinc_geometry": zinc_geometry(data),
            }
    return out


def final_classification(protocols, summaries):
    b = protocols.get("B_water_only_langevin", {})
    c = protocols.get("C_full_after_water", {})
    b_ok = b.get("code") == 0 and not b.get("parsed", {}).get("errors")
    c_ok = c.get("code") == 0 and not c.get("parsed", {}).get("errors")
    c_final = c.get("parsed", {}).get("final") or {}
    fmax = c_final.get("fmax")
    q2b_geom = summaries.get("C_full_after_water", {}).get("zinc_geometry", {})
    cn23 = q2b_geom.get("coordination_by_cutoff", {}).get("2.3")
    if not b_ok:
        return "needs_water_protocol_fix"
    if not c_ok:
        return "failed_water_relaxation"
    if fmax is None or fmax > 10.0:
        return "failed_water_relaxation"
    if cn23 != 4:
        return "debug_only_bad_zinc_geometry"
    return "md_ready_candidate"


def main():
    ensure_dir(AUDIT)
    q2b = parse_data(Q2B_DATA)
    pure = parse_data(PURE_DATA)
    audits = {
        "pure": water_audit(pure, highlighted=[]),
        "q2b": water_audit(q2b, highlighted=[268, 324]),
    }
    top_v24 = top_force_atoms(V24_FORCE_DUMP, q2b, 20)
    inputs = make_protocol_inputs()
    protocols = run_protocols(inputs)
    summaries = summarize_protocol_data(protocols)
    report = {
        "inputs": {
            "pure_data": PURE_DATA,
            "q2b_data": Q2B_DATA,
            "q2b_forcefield": Q2B_FF,
            "q2b_summary": Q2B_SUMMARY,
        },
        "water_topology_audit": audits,
        "v24_top_force_atoms": top_v24,
        "protocol_inputs": inputs,
        "protocol_results": protocols,
        "post_protocol_summaries": summaries,
        "tip4p_shake_diagnosis": {
            "pair_style": "lj/cut/tip4p/long 5 7 2 1 0.1546 10 10",
            "Ow_type": 5,
            "Hw_type": 7,
            "Ow_Hw_bond_type": 2,
            "Hw_Ow_Hw_angle_type": 1,
            "shake_line": "fix 1 water shake 1e-4 150 0 b 2 a 1",
            "topology_based_conclusion": "TIP4P topology is valid if every Ow has exactly two type-7 H bonded by type-2 bonds and one type-1 H-O-H angle.",
        },
    }
    report["final_classification"] = final_classification(protocols, summaries)
    out_json = os.path.join(AUDIT, "v25_water_audit_report.json")
    with open(out_json, "w") as f:
        json.dump(report, f, indent=2, sort_keys=True)
    print("Wrote {}".format(out_json))
    print("Final classification: {}".format(report["final_classification"]))


if __name__ == "__main__":
    main()

from __future__ import print_function

import argparse
import csv
import json
import math
import os

from validate_cementff_data import parse_data, nearest, distance


def rdf(data, type_a, type_b, rmax=8.0, dr=0.1):
    bins = int(rmax / dr)
    hist = [0] * bins
    atoms_a = [a for a in data["atoms"].values() if a["type"] in type_a]
    atoms_b = [a for a in data["atoms"].values() if a["type"] in type_b]
    for a in atoms_a:
        for b in atoms_b:
            if a["id"] == b["id"]:
                continue
            r = distance(a, b)
            idx = int(r / dr)
            if 0 <= idx < bins:
                hist[idx] += 1
    return [{"r": (i + 0.5) * dr, "count": hist[i]} for i in range(bins)]


def angle(a, b, c):
    v1 = (a["x"] - b["x"], a["y"] - b["y"], a["z"] - b["z"])
    v2 = (c["x"] - b["x"], c["y"] - b["y"], c["z"] - b["z"])
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return None
    cosv = max(-1.0, min(1.0, sum(v1[i] * v2[i] for i in range(3)) / (n1 * n2)))
    return math.degrees(math.acos(cosv))


def zinc_angles(data):
    zn_centered = []
    zn_oh_h = []
    for ang in data["angles"]:
        ids = [ang["a1"], ang["a2"], ang["a3"]]
        if any(i not in data["atoms"] for i in ids):
            continue
        atoms = [data["atoms"][i] for i in ids]
        val = angle(atoms[0], atoms[1], atoms[2])
        if val is None:
            continue
        if atoms[1]["type"] == 9:
            zn_centered.append({"angle_id": ang["id"], "type": ang["type"], "angle": val, "atoms": ids})
        if ang["type"] == 5:
            zn_oh_h.append({"angle_id": ang["id"], "type": ang["type"], "angle": val, "atoms": ids})
    return {"O_Zn_O": zn_centered, "Zn_Oh_H": zn_oh_h}


def summarize(values):
    if not values:
        return {"count": 0, "min": None, "mean": None, "max": None}
    vals = sorted(values)
    return {"count": len(vals), "min": vals[0], "mean": sum(vals) / len(vals), "max": vals[-1]}


def water_contacts(data):
    rows = []
    for atom in data["atoms"].values():
        if atom["type"] not in (5, 7):
            continue
        contacts = nearest(data, atom["id"], types={1, 2, 3, 4, 6, 8, 9})[:3]
        rows.append({"atom_id": atom["id"], "label": atom["label"], "nearest": contacts})
    return rows


def analyze(data_file, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    data = parse_data(data_file)
    zn_sites = [a["id"] for a in data["atoms"].values() if a["type"] == 9]
    zn_nn = []
    for zid in zn_sites:
        zn_nn.append({"Zn": zid, "nearest_oxygen": nearest(data, zid, types={3, 5, 6})[:8]})
    angle_records = zinc_angles(data)
    rdf_specs = {
        "Zn_O": ({9}, {3, 5, 6}),
        "Zn_Si": ({9}, {2}),
        "Zn_Ca": ({9}, {1}),
        "Si_O": ({2}, {3, 6}),
        "Ca_O": ({1}, {3, 5, 6}),
    }
    rdf_out = {name: rdf(data, *spec) for name, spec in rdf_specs.items()}
    summary = {
        "data_file": data_file,
        "n_atoms": len(data["atoms"]),
        "zinc_nearest_neighbors": zn_nn,
        "zinc_coordination_2p3": [sum(1 for x in rec["nearest_oxygen"] if x["distance"] <= 2.3) for rec in zn_nn],
        "zinc_coordination_2p5": [sum(1 for x in rec["nearest_oxygen"] if x["distance"] <= 2.5) for rec in zn_nn],
        "angle_summary": {
            "O_Zn_O": summarize([x["angle"] for x in angle_records["O_Zn_O"]]),
            "Zn_Oh_H": summarize([x["angle"] for x in angle_records["Zn_Oh_H"]]),
        },
        "water_contacts": water_contacts(data),
        "rdf_files": {},
    }
    for name, rows in rdf_out.items():
        path = os.path.join(out_dir, "rdf_{}.csv".format(name))
        with open(path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["r", "count"])
            writer.writeheader()
            writer.writerows(rows)
        summary["rdf_files"][name] = path
    out_json = os.path.join(out_dir, "structure_analysis.json")
    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = analyze(args.data_file, args.out)
    print(json.dumps({"out": os.path.join(args.out, "structure_analysis.json"), "n_atoms": result["n_atoms"]}, indent=2))


if __name__ == "__main__":
    main()

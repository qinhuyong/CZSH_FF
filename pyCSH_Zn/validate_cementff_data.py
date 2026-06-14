from __future__ import print_function

import argparse
import json
import math
import os
import re


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
    10: "Al",
    11: "Cl",
}
VALID_ATOM_TYPES = set(LABELS)
VALID_BOND_TYPES = {1, 2, 3}
VALID_ANGLE_TYPES = {1, 2, 3, 4, 5}


def parse_data(path):
    data = {"path": path, "counts": {}, "box": {}, "masses": {}, "atoms": {}, "bonds": [], "angles": [], "csinfo": {}}
    section = None
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
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
            header = line.split("#")[0].strip()
            if header in ("Masses", "Bonds", "Angles", "CS-Info") or header.startswith("Atoms"):
                section = "Atoms" if header.startswith("Atoms") else header
                continue
            p = line.split("#")[0].split()
            if section == "Masses" and len(p) >= 2:
                data["masses"][int(p[0])] = float(p[1])
            elif section == "Atoms" and len(p) >= 7:
                aid = int(p[0])
                typ = int(p[2])
                data["atoms"][aid] = {
                    "id": aid,
                    "mol": int(p[1]),
                    "type": typ,
                    "label": LABELS.get(typ, "type{}".format(typ)),
                    "q": float(p[3]),
                    "x": float(p[4]),
                    "y": float(p[5]),
                    "z": float(p[6]),
                }
            elif section == "Bonds" and len(p) >= 4:
                data["bonds"].append({"id": int(p[0]), "type": int(p[1]), "a1": int(p[2]), "a2": int(p[3])})
            elif section == "Angles" and len(p) >= 5:
                data["angles"].append({"id": int(p[0]), "type": int(p[1]), "a1": int(p[2]), "a2": int(p[3]), "a3": int(p[4])})
            elif section == "CS-Info" and len(p) >= 2:
                data["csinfo"][int(p[0])] = int(p[1])
    return data


def distance(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def nearest(data, atom_id, types=None, exclude=None):
    exclude = set(exclude or [])
    atom = data["atoms"][atom_id]
    out = []
    for other_id, other in data["atoms"].items():
        if other_id == atom_id or other_id in exclude:
            continue
        if types is not None and other["type"] not in types:
            continue
        out.append({"atom_id": other_id, "type": other["type"], "label": other["label"], "distance": distance(atom, other)})
    return sorted(out, key=lambda x: x["distance"])


def audit_csinfo(data):
    bad = []
    pairs = []
    for b in data["bonds"]:
        if b["a1"] in data["atoms"] and b["a2"] in data["atoms"]:
            t = {data["atoms"][b["a1"]]["type"], data["atoms"][b["a2"]]["type"]}
            if t == {3, 4}:
                same = data["csinfo"].get(b["a1"]) == data["csinfo"].get(b["a2"])
                rec = {"bond_id": b["id"], "a1": b["a1"], "a2": b["a2"], "same_csid": same, "distance": distance(data["atoms"][b["a1"]], data["atoms"][b["a2"]])}
                pairs.append(rec)
                if not same:
                    bad.append(rec)
    return {"n_pairs": len(pairs), "n_csinfo": len(data["csinfo"]), "bad_pairs": bad, "pairs": pairs}


def audit_water(data):
    waters = []
    bad = []
    for aid, atom in data["atoms"].items():
        if atom["type"] != 5:
            continue
        bonded_h = []
        for b in data["bonds"]:
            if b["type"] == 2 and aid in (b["a1"], b["a2"]):
                other = b["a2"] if b["a1"] == aid else b["a1"]
                if other in data["atoms"] and data["atoms"][other]["type"] == 7:
                    bonded_h.append(other)
        rec = {"Ow": aid, "Hw": sorted(bonded_h), "molecule_id": atom["mol"], "tip4p_compatible": len(bonded_h) == 2}
        waters.append(rec)
        if len(bonded_h) != 2:
            bad.append(rec)
    return {"n_water": len(waters), "n_bad_water": len(bad), "waters": waters, "bad_waters": bad}


def audit_zinc(data):
    zns = [aid for aid, a in data["atoms"].items() if a["type"] == 9]
    records = []
    for zid in zns:
        zn_o = nearest(data, zid, types={3, 5, 6})[:8]
        records.append({
            "Zn": zid,
            "coordination_2p3": sum(1 for x in zn_o if x["distance"] <= 2.3),
            "coordination_2p5": sum(1 for x in zn_o if x["distance"] <= 2.5),
            "nearest_oxygen": zn_o[:4],
        })
    return {"n_zinc": len(zns), "zinc_sites": records}


def validate(path):
    data = parse_data(path)
    total_charge = sum(a["q"] for a in data["atoms"].values())
    atom_type_bad = [a for a in data["atoms"].values() if a["type"] not in VALID_ATOM_TYPES]
    bond_type_bad = [b for b in data["bonds"] if b["type"] not in VALID_BOND_TYPES]
    angle_type_bad = [a for a in data["angles"] if a["type"] not in VALID_ANGLE_TYPES]
    cs = audit_csinfo(data)
    water = audit_water(data)
    zinc = audit_zinc(data)
    reasons = []
    classification = "valid_static_candidate"
    if abs(total_charge) > 1.0e-5:
        classification = "failed_charge"
        reasons.append("charge residual {:.6g}".format(total_charge))
    elif atom_type_bad or bond_type_bad or angle_type_bad:
        classification = "failed_topology"
        reasons.append("invalid atom/bond/angle type")
    elif cs["bad_pairs"] or cs["n_csinfo"] != len(data["atoms"]):
        classification = "failed_csinfo"
        reasons.append("CS-Info does not match core-shell pairs")
    elif water["n_bad_water"]:
        classification = "failed_water_contacts"
        reasons.append("invalid TIP4P water topology")
    elif zinc["n_zinc"]:
        if any(site["coordination_2p5"] < 4 for site in zinc["zinc_sites"]):
            classification = "failed_zinc_geometry"
            reasons.append("Zn has fewer than 4 O neighbors within 2.5 A")
        else:
            classification = "valid_q2b_zn_candidate"
    return {
        "data_file": path,
        "classification": classification,
        "reasons": reasons,
        "counts": {
            "atoms": len(data["atoms"]),
            "bonds": len(data["bonds"]),
            "angles": len(data["angles"]),
            "atom_types": sorted(set(a["type"] for a in data["atoms"].values())),
        },
        "total_charge": total_charge,
        "csinfo": {"n_pairs": cs["n_pairs"], "n_entries": cs["n_csinfo"], "n_bad_pairs": len(cs["bad_pairs"])},
        "water": {"n_water": water["n_water"], "n_bad_water": water["n_bad_water"]},
        "zinc": zinc,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("data_file")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()
    result = validate(args.data_file)
    out = args.out or os.path.splitext(args.data_file)[0] + "_validation.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")
    print(json.dumps({"classification": result["classification"], "out": out}, indent=2))


if __name__ == "__main__":
    main()

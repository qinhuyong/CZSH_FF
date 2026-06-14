import json
import math
import os
import copy

import numpy as np


ZN_SPECIE = 14
ZN_CHARGE = 2.0
CHARGE_TOLERANCE = 1.0e-6
O_CORE_CHARGE = 0.848
O_SHELL_CHARGE = -2.848
OH_CHARGE = -1.4
HOH_CHARGE = 0.4
OH_BOND_DISTANCE = 1.0
ZN_O_CUTOFF = 2.3
H_OVERLAP_CUTOFF = 0.75
H_MIN_DISTANCES = {
    "H_H": 1.2,
    "H_O_nonbonded": 1.2,
    "H_Si": 1.6,
    "H_Ca": 1.6,
    "H_Zn": 1.4,
}
CONVERTED_O_MIN_DISTANCES = {
    "O_Ca": 1.6,
    "O_Si": 1.45,
    "O_O_nonbonded": 1.2,
    "O_Zn": 1.7,
}

SUPPORTED_CHARGE_BALANCE_MODES = {
    "hydroxylate_two_oxygens",
    "fail_if_not_neutral",
    "allow_unbalanced_for_debug",
    "add_interlayer_Ca",
}

CEMENTFF4_TYPE_MAP = {
    1: {"lammps_type": 1, "label": "Ca", "source": "Ca1"},
    9: {"lammps_type": 1, "label": "Ca", "source": "Ca2"},
    2: {"lammps_type": 2, "label": "Si", "source": "Si1"},
    10: {"lammps_type": 2, "label": "Si", "source": "Si2"},
    3: {"lammps_type": 3, "label": "O", "source": "O core"},
    11: {"lammps_type": 3, "label": "O", "source": "OCa core"},
    4: {"lammps_type": 4, "label": "O(S)", "source": "O shell"},
    5: {"lammps_type": 5, "label": "Ow", "source": "water oxygen"},
    6: {"lammps_type": 6, "label": "Oh", "source": "hydroxide oxygen"},
    7: {"lammps_type": 7, "label": "Hw", "source": "water hydrogen"},
    8: {"lammps_type": 8, "label": "Hoh", "source": "hydroxide hydrogen"},
    14: {"lammps_type": 9, "label": "Zn", "source": "Zn substituted Si site"},
}

CEMENTFF4_ANGLE_MAP = {
    1: "Hw-Ow-Hw",
    2: "O-Si-O / Oh-Si-O / Oh-Si-Oh",
    3: "Si-Oh-H",
    4: "O-Zn-O / Oh-Zn-O / Oh-Zn-Oh",
    5: "Zn-Oh-H",
}

Q1_PIECES = {"<L", "<R", ">L", ">R", "<Lo", "<Ro", ">Lo", ">Ro"}
Q2B_PIECES = {"SU", "SD", "SUo", "SDo"}
SUPPORTED_SITE_TYPES = {"Q1_Zn", "Q2b_Zn", "mixed_Q1_Q2b_Zn"}
UNSUPPORTED_SITE_TYPES = {"interlayer_Zn", "Ca_substitution_control"}


def validate_zinc_site_type(site_type):
    if site_type in UNSUPPORTED_SITE_TYPES:
        raise NotImplementedError(site_type + " is not implemented in v0")
    if site_type not in SUPPORTED_SITE_TYPES:
        raise ValueError(
            "Unknown Zn_site_type {!r}. Expected one of {}".format(
                site_type, sorted(SUPPORTED_SITE_TYPES | UNSUPPORTED_SITE_TYPES)
            )
        )


def validate_charge_balance_mode(mode):
    if mode not in SUPPORTED_CHARGE_BALANCE_MODES:
        raise ValueError(
            "Unknown Zn_charge_balance_mode {!r}. Expected one of {}".format(
                mode, sorted(SUPPORTED_CHARGE_BALANCE_MODES)
            )
        )
    if mode == "add_interlayer_Ca":
        raise NotImplementedError("add_interlayer_Ca is not safely implemented in v2")


def inspect_zinc_candidates(crystal_dict):
    """Return Q1-like and Q2b-like Si centers from the expanded crystal."""
    q1_sites = []
    q2b_sites = []

    for cell, brick_dict in crystal_dict.items():
        for piece_name, piece_entries in brick_dict.items():
            if piece_name not in Q1_PIECES and piece_name not in Q2B_PIECES:
                continue

            si_entries = [entry for entry in piece_entries if entry[1] in (2, 10)]
            if not si_entries:
                raise ValueError(
                    "Cannot identify a Si center for Zn candidate piece "
                    "{!r} at cell {}".format(piece_name, cell)
                )

            for si_entry in si_entries:
                site = {
                    "atom_id": int(si_entry[0]),
                    "cell": [int(cell[0]), int(cell[1]), int(cell[2])],
                    "piece": piece_name,
                    "motif": "Q1_Zn" if piece_name in Q1_PIECES else "Q2b_Zn",
                    "coord": [float(si_entry[3]), float(si_entry[4]), float(si_entry[5])],
                    "original_specie": int(si_entry[1]),
                    "original_charge": float(si_entry[2]),
                }
                if piece_name in Q1_PIECES:
                    q1_sites.append(site)
                else:
                    q2b_sites.append(site)

    if not q1_sites and not q2b_sites:
        raise ValueError("No Q1-like or Q2b-like silicate sites were identified for Zn placement")

    return {"Q1_Zn": q1_sites, "Q2b_Zn": q2b_sites}


def count_species(entries):
    counts = {}
    for entry in entries:
        counts[int(entry[1])] = counts.get(int(entry[1]), 0) + 1
    return counts


def coords_by_atom_id(entries):
    return {int(entry[0]): np.array(entry[3:], dtype=float) for entry in entries}


def type_by_atom_id(entries):
    return {int(entry[0]): int(entry[1]) for entry in entries}


def bonded_atom_ids(entries_bonds, atom_id, bond_types=None):
    bonded = []
    for bond in entries_bonds:
        if bond_types is not None and int(bond[1]) not in bond_types:
            continue
        if int(bond[2]) == int(atom_id):
            bonded.append(int(bond[3]))
        elif int(bond[3]) == int(atom_id):
            bonded.append(int(bond[2]))
    return bonded


def minimum_periodic_distance(coord, selected_coords, supercell):
    if not selected_coords:
        return math.inf
    inv_supercell = np.linalg.inv(supercell)
    coord = np.array(coord, dtype=float)
    min_distance = math.inf
    for other in selected_coords:
        delta = coord - np.array(other, dtype=float)
        frac = np.dot(delta, inv_supercell)
        frac -= np.rint(frac)
        delta_pbc = np.dot(frac, supercell)
        min_distance = min(min_distance, float(np.linalg.norm(delta_pbc)))
    return min_distance


def minimum_distance_to_species(coord, entries, species, supercell):
    coords = [entry[3:] for entry in entries if int(entry[1]) in species]
    return minimum_periodic_distance(coord, coords, supercell)


def periodic_distance(coord1, coord2, supercell):
    return minimum_periodic_distance(coord1, [coord2], supercell)


def angle_degrees(v1, v2):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return None
    cosang = float(np.dot(v1, v2) / denom)
    cosang = max(-1.0, min(1.0, cosang))
    return float(np.degrees(np.arccos(cosang)))


def vector_pbc(origin, other, supercell):
    inv_supercell = np.linalg.inv(supercell)
    delta = np.array(other, dtype=float) - np.array(origin, dtype=float)
    frac = np.dot(delta, inv_supercell)
    frac -= np.rint(frac)
    return np.dot(frac, supercell)


def total_charge(entries_crystal):
    return float(sum(float(entry[2]) for entry in entries_crystal))


def next_atom_id(entries_crystal):
    return max(int(entry[0]) for entry in entries_crystal) + 1


def find_entry(entries, atom_id):
    for entry in entries:
        if int(entry[0]) == int(atom_id):
            return entry
    return None


def find_shell_for_core(entries_crystal, entries_bonds, core_id):
    for bond in entries_bonds:
        if int(bond[1]) == 3 and int(bond[2]) == int(core_id):
            shell = find_entry(entries_crystal, int(bond[3]))
            if shell is not None and int(shell[1]) == 4:
                return shell, bond
    return None, None


def count_nearby_si_neighbors(atom_id, entries_crystal, supercell, cutoff=2.2, exclude_si_ids=None):
    coords = coords_by_atom_id(entries_crystal)
    atom_types = type_by_atom_id(entries_crystal)
    o_coord = coords[int(atom_id)]
    count = 0
    neighbors = []
    exclude_si_ids = set(exclude_si_ids or [])
    for other_id, specie in atom_types.items():
        if other_id in exclude_si_ids:
            continue
        if specie not in (2, 10):
            continue
        distance = periodic_distance(o_coord, coords[other_id], supercell)
        if distance <= cutoff:
            count += 1
            neighbors.append({"atom_id": int(other_id), "distance": float(distance)})
    return count, neighbors


def classify_oxygen_for_hydroxylation(atom_id, entries_crystal, entries_bonds, supercell, exclude_si_ids=None):
    entry = find_entry(entries_crystal, atom_id)
    if entry is None:
        return {"oxygen_class": "missing", "safe_for_default_hydroxylation": False}

    specie = int(entry[1])
    if specie == 5:
        return {"oxygen_class": "water oxygen", "safe_for_default_hydroxylation": False}
    if specie == 6:
        return {"oxygen_class": "hydroxyl oxygen already", "safe_for_default_hydroxylation": False}
    if specie not in (3, 11):
        return {"oxygen_class": "not core/shell oxygen", "safe_for_default_hydroxylation": False}

    shell, shell_bond = find_shell_for_core(entries_crystal, entries_bonds, atom_id)
    if shell is None or shell_bond is None:
        return {"oxygen_class": "core oxygen without O(S) shell", "safe_for_default_hydroxylation": False}

    si_count, si_neighbors = count_nearby_si_neighbors(atom_id, entries_crystal, supercell, exclude_si_ids=exclude_si_ids)
    if si_count == 0:
        oxygen_class = "terminal/non-bridging O"
        safe = True
    else:
        oxygen_class = "bridging Zn/Si-O-Si oxygen"
        safe = False

    return {
        "oxygen_class": oxygen_class,
        "safe_for_default_hydroxylation": safe,
        "nearby_si_count": int(si_count),
        "nearby_si_neighbors": si_neighbors,
    }


def oxygen_candidates_for_site(site, entries_crystal, entries_bonds, supercell, allow_hydroxylate_bridging_oxygen=False):
    coords = coords_by_atom_id(entries_crystal)
    atom_types = type_by_atom_id(entries_crystal)
    zn_coord = np.array(site["coord"], dtype=float)
    candidates = []
    for atom_id, specie in atom_types.items():
        if specie not in (3, 11):
            continue
        shell, shell_bond = find_shell_for_core(entries_crystal, entries_bonds, atom_id)
        if shell is None:
            continue
        distance = periodic_distance(zn_coord, coords[atom_id], supercell)
        if distance <= ZN_O_CUTOFF:
            classification = classify_oxygen_for_hydroxylation(
                atom_id,
                entries_crystal,
                entries_bonds,
                supercell,
                exclude_si_ids={int(site["atom_id"])},
            )
            safe = bool(classification["safe_for_default_hydroxylation"])
            if classification["oxygen_class"] == "bridging Si-O-Si oxygen" and allow_hydroxylate_bridging_oxygen:
                safe = True
            candidates.append(
                {
                    "atom_id": atom_id,
                    "shell_id": int(shell[0]),
                    "shell_bond_id": int(shell_bond[0]),
                    "distance": distance,
                    "specie": int(specie),
                    "oxygen_class": classification["oxygen_class"],
                    "safe_for_default_hydroxylation": safe,
                    "nearby_si_count": classification.get("nearby_si_count"),
                    "nearby_si_neighbors": classification.get("nearby_si_neighbors", []),
                }
            )
    candidates.sort(key=lambda item: (not item["safe_for_default_hydroxylation"], item["distance"], item["atom_id"]))
    return candidates


def select_oxygens_for_hydroxylation(
    site,
    entries_crystal,
    entries_bonds,
    supercell,
    n_oxygen=2,
    allow_hydroxylate_bridging_oxygen=False,
):
    candidates = oxygen_candidates_for_site(
        site, entries_crystal, entries_bonds, supercell, allow_hydroxylate_bridging_oxygen
    )
    safe_candidates = [candidate for candidate in candidates if candidate["safe_for_default_hydroxylation"]]
    if len(safe_candidates) < n_oxygen:
        raise ValueError(
            "Selected Zn site atom_id={} has only {} safe terminal/non-bridging O core/shell "
            "candidates within {:.2f} A; need {}. Candidate classes: {}".format(
                site["atom_id"],
                len(safe_candidates),
                ZN_O_CUTOFF,
                n_oxygen,
                [
                    {
                        "atom_id": item["atom_id"],
                        "class": item["oxygen_class"],
                        "distance": item["distance"],
                    }
                    for item in candidates
                ],
            )
        )
    return safe_candidates[:n_oxygen]


def h_direction_away_from_zn(o_coord, zn_coord):
    direction = np.array(o_coord, dtype=float) - np.array(zn_coord, dtype=float)
    norm = np.linalg.norm(direction)
    if norm < 1.0e-8:
        return np.array([0.0, 0.0, 1.0])
    return direction / norm


def min_h_overlap(h_coord, entries_crystal, ignore_ids, supercell):
    distances = []
    for entry in entries_crystal:
        atom_id = int(entry[0])
        if atom_id in ignore_ids:
            continue
        distances.append(periodic_distance(h_coord, entry[3:], supercell))
    return None if not distances else float(min(distances))


def distance_to_type_set(h_coord, entries_crystal, type_set, supercell, ignore_ids=None):
    ignore_ids = set(ignore_ids or [])
    best = {"distance": None, "atom_id": None, "specie": None}
    for entry in entries_crystal:
        atom_id = int(entry[0])
        if atom_id in ignore_ids or int(entry[1]) not in type_set:
            continue
        distance = periodic_distance(h_coord, entry[3:], supercell)
        if best["distance"] is None or distance < best["distance"]:
            best = {"distance": float(distance), "atom_id": atom_id, "specie": int(entry[1])}
    return best


def h_contact_metrics(h_coord, entries_crystal, supercell, bonded_o_id, h_id, zn_id):
    ignore_o = {int(bonded_o_id), int(h_id)}
    return {
        "H_H": distance_to_type_set(h_coord, entries_crystal, {7, 8}, supercell, {int(h_id)}),
        "H_O_nonbonded": distance_to_type_set(h_coord, entries_crystal, {3, 4, 5, 6, 11, 12}, supercell, ignore_o),
        "H_Si": distance_to_type_set(h_coord, entries_crystal, {2, 10}, supercell, {int(h_id)}),
        "H_Ca": distance_to_type_set(h_coord, entries_crystal, {1, 9}, supercell, {int(h_id)}),
        "H_Zn": distance_to_type_set(h_coord, entries_crystal, {ZN_SPECIE}, supercell, {int(h_id)}),
    }


def h_contacts_are_safe(metrics):
    for key, cutoff in H_MIN_DISTANCES.items():
        distance = metrics.get(key, {}).get("distance")
        if distance is not None and distance < cutoff:
            return False
    return True


def oxygen_contact_metrics(o_coord, entries_crystal, supercell, o_id, h_id, zn_id):
    return {
        "O_Ca": distance_to_type_set(o_coord, entries_crystal, {1, 9}, supercell, {int(o_id), int(h_id)}),
        "O_Si": distance_to_type_set(o_coord, entries_crystal, {2, 10}, supercell, {int(o_id), int(h_id), int(zn_id)}),
        "O_O_nonbonded": distance_to_type_set(
            o_coord,
            entries_crystal,
            {3, 4, 5, 6, 11, 12},
            supercell,
            {int(o_id), int(h_id)},
        ),
        "O_Zn": distance_to_type_set(o_coord, entries_crystal, {ZN_SPECIE}, supercell, {int(o_id), int(h_id)}),
    }


def oxygen_contacts_are_safe(metrics):
    for key, cutoff in CONVERTED_O_MIN_DISTANCES.items():
        distance = metrics.get(key, {}).get("distance")
        if distance is not None and distance < cutoff:
            return False
    return True


def orthonormal_basis(axis):
    axis = np.array(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(axis, trial))) > 0.85:
        trial = np.array([0.0, 1.0, 0.0])
    u = trial - np.dot(trial, axis) * axis
    u = u / np.linalg.norm(u)
    v = np.cross(axis, u)
    return u, v


def trial_h_directions(base_direction):
    base = np.array(base_direction, dtype=float)
    base = base / np.linalg.norm(base)
    directions = [base]
    u, v = orthonormal_basis(base)
    for polar_deg in (30.0, 50.0, 70.0, 90.0, 110.0):
        polar = math.radians(polar_deg)
        for azimuth_deg in range(0, 360, 30):
            azimuth = math.radians(azimuth_deg)
            direction = (
                math.cos(polar) * base
                + math.sin(polar) * (math.cos(azimuth) * u + math.sin(azimuth) * v)
            )
            directions.append(direction / np.linalg.norm(direction))
    return directions


def choose_h_position(o_coord, zn_coord, entries_crystal, supercell, bonded_o_id, h_id, zn_id):
    base = h_direction_away_from_zn(o_coord, zn_coord)
    best = None
    for direction in trial_h_directions(base):
        h_coord = np.array(o_coord, dtype=float) + OH_BOND_DISTANCE * direction
        metrics = h_contact_metrics(h_coord, entries_crystal, supercell, bonded_o_id, h_id, zn_id)
        score = min(
            metrics[key]["distance"] if metrics[key]["distance"] is not None else 99.0
            for key in H_MIN_DISTANCES
        )
        record = {
            "coord": h_coord,
            "direction": direction,
            "metrics": metrics,
            "score": float(score),
            "safe": h_contacts_are_safe(metrics),
        }
        if best is None or record["score"] > best["score"]:
            best = record
        if record["safe"]:
            return record
    return best


def set_entry_coord(entry, coord):
    entry[3] = float(coord[0])
    entry[4] = float(coord[1])
    entry[5] = float(coord[2])


def precondition_converted_oxygen(o_entry, site, supercell, target_Zn_O_distance):
    zn_coord = np.array(site["coord"], dtype=float)
    o_coord = np.array(o_entry[3:], dtype=float)
    vector = vector_pbc(zn_coord, o_coord, supercell)
    distance = float(np.linalg.norm(vector))
    if distance < 1.0e-8:
        return {"moved": False, "old_distance": distance, "new_distance": distance}
    if distance >= target_Zn_O_distance:
        return {"moved": False, "old_distance": distance, "new_distance": distance}
    new_coord = zn_coord + vector / distance * float(target_Zn_O_distance)
    set_entry_coord(o_entry, new_coord)
    return {
        "moved": True,
        "old_distance": distance,
        "new_distance": float(target_Zn_O_distance),
        "old_coord": [float(x) for x in o_coord],
        "new_coord": [float(x) for x in new_coord],
    }


def hydroxylate_two_oxygens(
    entries_crystal,
    entries_bonds,
    entries_angle,
    selected_sites,
    supercell,
    allow_hydroxylate_bridging_oxygen=False,
    precondition_zinc_geometry=True,
    target_Zn_O_distance=1.95,
):
    records = []
    for site in selected_sites:
        selected_oxygens = select_oxygens_for_hydroxylation(
            site, entries_crystal, entries_bonds, supercell, 2, allow_hydroxylate_bridging_oxygen
        )
        site_records = []
        for oxygen in selected_oxygens:
            o_entry = find_entry(entries_crystal, oxygen["atom_id"])
            h_entry = find_entry(entries_crystal, oxygen["shell_id"])
            shell_bond = None
            for bond in entries_bonds:
                if int(bond[0]) == oxygen["shell_bond_id"]:
                    shell_bond = bond
                    break
            if o_entry is None or h_entry is None or shell_bond is None:
                raise ValueError("Internal error while hydroxylating O atom {}".format(oxygen["atom_id"]))

            precondition_record = {"moved": False}
            if precondition_zinc_geometry:
                precondition_record = precondition_converted_oxygen(
                    o_entry, site, supercell, target_Zn_O_distance
                )
            o_coord = np.array(o_entry[3:], dtype=float)
            o_contact_metrics = oxygen_contact_metrics(
                o_coord,
                entries_crystal,
                supercell,
                int(o_entry[0]),
                int(h_entry[0]),
                int(site["atom_id"]),
            )
            if not oxygen_contacts_are_safe(o_contact_metrics):
                if precondition_record.get("moved") and "old_coord" in precondition_record:
                    set_entry_coord(o_entry, precondition_record["old_coord"])
                raise ValueError(
                    "Preconditioning converted O atom {} for Zn atom {} creates unsafe contacts: {}".format(
                        o_entry[0], site["atom_id"], o_contact_metrics
                    )
                )
            placement = choose_h_position(
                o_coord,
                site["coord"],
                entries_crystal,
                supercell,
                int(o_entry[0]),
                int(h_entry[0]),
                int(site["atom_id"]),
            )
            if placement is None or not placement["safe"]:
                raise ValueError(
                    "Could not place Hoh safely for Zn atom {} and O atom {}. Best contact metrics: {}".format(
                        site["atom_id"], o_entry[0], None if placement is None else placement["metrics"]
                    )
                )
            h_coord = placement["coord"]
            direction = placement["direction"]

            o_entry[1] = 6
            o_entry[2] = OH_CHARGE
            h_entry[1] = 8
            h_entry[2] = HOH_CHARGE
            h_entry[3] = float(h_coord[0])
            h_entry[4] = float(h_coord[1])
            h_entry[5] = float(h_coord[2])

            shell_bond[1] = 1
            shell_bond[2] = int(o_entry[0])
            shell_bond[3] = int(h_entry[0])

            angle_id = max([int(angle[0]) for angle in entries_angle] or [0]) + 1
            entries_angle.append([angle_id, 5, int(site["atom_id"]), int(o_entry[0]), int(h_entry[0])])

            overlap = min_h_overlap(h_coord, entries_crystal, {int(o_entry[0]), int(h_entry[0]), int(site["atom_id"])}, supercell)
            site_records.append(
                {
                    "oxygen_atom_id": int(o_entry[0]),
                    "reused_shell_as_H_atom_id": int(h_entry[0]),
                    "modified_bond_id": int(shell_bond[0]),
                    "added_angle_id": int(angle_id),
                    "original_oxygen_specie": int(oxygen["specie"]),
                    "original_shell_atom_id": int(oxygen["shell_id"]),
                    "oxygen_class": oxygen["oxygen_class"],
                    "nearby_si_count_before_hydroxylation": oxygen.get("nearby_si_count"),
                    "nearby_si_neighbors_before_hydroxylation": oxygen.get("nearby_si_neighbors", []),
                    "Zn_O_distance_before_hydroxylation": float(oxygen["distance"]),
                    "Zn_O_preconditioning": precondition_record,
                    "converted_O_contact_metrics": o_contact_metrics,
                    "H_coord": [float(h_coord[0]), float(h_coord[1]), float(h_coord[2])],
                    "H_placement_vector": [float(direction[0]), float(direction[1]), float(direction[2])],
                    "O_H_distance": float(periodic_distance(o_coord, h_coord, supercell)),
                    "H_contact_metrics": placement["metrics"],
                    "min_H_overlap_distance": overlap,
                    "H_placement": "deterministic_overlap_checked_search",
                }
            )
        records.append({"zn_atom_id": int(site["atom_id"]), "hydroxylated_oxygens": site_records})
    return records


def apply_charge_balance(
    entries_crystal,
    entries_bonds,
    entries_angle,
    selected_sites,
    supercell,
    mode,
    allow_hydroxylate_bridging_oxygen=False,
    precondition_zinc_geometry=True,
    target_Zn_O_distance=1.95,
):
    validate_charge_balance_mode(mode)
    if mode == "hydroxylate_two_oxygens":
        return hydroxylate_two_oxygens(
            entries_crystal,
            entries_bonds,
            entries_angle,
            selected_sites,
            supercell,
            allow_hydroxylate_bridging_oxygen,
            precondition_zinc_geometry,
            target_Zn_O_distance,
        )
    if mode in ("fail_if_not_neutral", "allow_unbalanced_for_debug"):
        return []
    raise NotImplementedError(mode + " is not implemented in v2")


def select_zinc_sites(candidates, n_zinc, site_type, seed, supercell, min_zn_zn_distance=3.0, site_filter=None):
    validate_zinc_site_type(site_type)
    if n_zinc <= 0:
        return []

    if site_type == "Q1_Zn":
        pool = list(candidates["Q1_Zn"])
    elif site_type == "Q2b_Zn":
        pool = list(candidates["Q2b_Zn"])
    else:
        pool = list(candidates["Q1_Zn"]) + list(candidates["Q2b_Zn"])

    if not pool:
        raise ValueError("No candidate sites are available for " + site_type)
    if n_zinc > len(pool):
        raise ValueError(
            "Requested {} Zn sites, but only {} eligible {} sites were found".format(
                n_zinc, len(pool), site_type
            )
        )

    rng = np.random.default_rng(seed)
    order = rng.permutation(len(pool))
    selected = []
    selected_coords = []
    skipped_for_distance = 0

    for idx in order:
        site = dict(pool[int(idx)])
        if site_filter is not None and not site_filter(site):
            continue
        distance = minimum_periodic_distance(site["coord"], selected_coords, supercell)
        if distance < min_zn_zn_distance:
            skipped_for_distance += 1
            continue
        selected.append(site)
        selected_coords.append(site["coord"])
        if len(selected) == n_zinc:
            break

    if len(selected) != n_zinc:
        raise ValueError(
            "Could only select {} of {} Zn sites without Zn-Zn distances below {:.2f} A "
            "({} candidates skipped). Lower Zn_Si_ratio or increase the supercell.".format(
                len(selected), n_zinc, min_zn_zn_distance, skipped_for_distance
            )
        )

    return selected


def apply_zinc_sites(entries_crystal, crystal_dict, selected_sites):
    selected_ids = {site["atom_id"] for site in selected_sites}
    touched = set()

    for entry in entries_crystal:
        if int(entry[0]) in selected_ids:
            if int(entry[1]) in (1, 9):
                raise ValueError("Refusing to place Zn on Ca atom id {}".format(entry[0]))
            entry[1] = ZN_SPECIE
            entry[2] = ZN_CHARGE
            touched.add(int(entry[0]))

    for brick_dict in crystal_dict.values():
        for piece_entries in brick_dict.values():
            for entry in piece_entries:
                if int(entry[0]) in selected_ids:
                    entry[1] = ZN_SPECIE
                    entry[2] = ZN_CHARGE

    missing = selected_ids.difference(touched)
    if missing:
        raise ValueError("Selected Zn atom ids were not found in entries_crystal: {}".format(sorted(missing)))

    return entries_crystal, crystal_dict


def remap_zinc_angles(entries_crystal, entries_angle):
    """Convert Zn-centered old Si angle types to CementFF4 Zn angle types."""
    atom_types = type_by_atom_id(entries_crystal)
    selected_zn = {atom_id for atom_id, specie in atom_types.items() if specie == ZN_SPECIE}
    remapped_ozno = 0
    remapped_znohh = 0
    stale_zn_angles = []

    for angle in entries_angle:
        atom1, atom2, atom3 = int(angle[2]), int(angle[3]), int(angle[4])
        if atom2 in selected_zn and int(angle[1]) == 2:
            angle[1] = 4
            remapped_ozno += 1
        elif atom1 in selected_zn and int(angle[1]) == 3:
            angle[1] = 5
            remapped_znohh += 1
        elif atom2 in selected_zn and int(angle[1]) == 4:
            continue
        elif atom1 in selected_zn and int(angle[1]) == 5:
            continue
        elif atom1 in selected_zn or atom2 in selected_zn or atom3 in selected_zn:
            stale_zn_angles.append(list(angle))

    return {
        "remapped_O_Zn_O_angles": remapped_ozno,
        "remapped_Zn_Oh_H_angles": remapped_znohh,
        "stale_zn_angles": stale_zn_angles,
    }


def validate_no_zinc_bonds(entries_crystal, entries_bonds):
    atom_types = type_by_atom_id(entries_crystal)
    zn_ids = {atom_id for atom_id, specie in atom_types.items() if specie == ZN_SPECIE}
    zinc_bonds = []
    for bond in entries_bonds:
        if int(bond[2]) in zn_ids or int(bond[3]) in zn_ids:
            zinc_bonds.append(list(bond))
    if zinc_bonds:
        raise ValueError(
            "CementFF4 v1 does not define Zn-O bonds; refusing Zn-bonded topology: {}".format(zinc_bonds)
        )
    return zinc_bonds


def geometry_metrics(entries_crystal, entries_angle, supercell):
    coords = coords_by_atom_id(entries_crystal)
    atom_types = type_by_atom_id(entries_crystal)
    zn_ids = [atom_id for atom_id, specie in atom_types.items() if specie == ZN_SPECIE]
    o_species = {3, 4, 5, 6, 11, 12}
    silicate_o_species = {3, 6, 11, 12}
    ca_species = {1, 9}
    si_species = {2, 10}

    zn_o_distances = []
    zn_coordination = []
    for zn_id in zn_ids:
        distances = [
            periodic_distance(coords[zn_id], coord, supercell)
            for atom_id, coord in coords.items()
            if atom_types[atom_id] in silicate_o_species
        ]
        neighbors = [d for d in distances if d <= 2.3]
        zn_coordination.append(len(neighbors))
        zn_o_distances.extend(neighbors)

    zn_zn_distances = []
    for i, zn_i in enumerate(zn_ids):
        for zn_j in zn_ids[i + 1:]:
            zn_zn_distances.append(periodic_distance(coords[zn_i], coords[zn_j], supercell))

    si_o_distances = []
    for si_id, si_coord in coords.items():
        if atom_types[si_id] not in si_species:
            continue
        for o_id, o_coord in coords.items():
            if atom_types[o_id] in silicate_o_species:
                distance = periodic_distance(si_coord, o_coord, supercell)
                if distance <= 2.2:
                    si_o_distances.append(distance)

    ca_o_distances = []
    for ca_id, ca_coord in coords.items():
        if atom_types[ca_id] not in ca_species:
            continue
        for o_id, o_coord in coords.items():
            if atom_types[o_id] in o_species:
                distance = periodic_distance(ca_coord, o_coord, supercell)
                if distance <= 3.2:
                    ca_o_distances.append(distance)

    o_zn_o_angles = []
    for angle in entries_angle:
        if int(angle[1]) != 4:
            continue
        a1, center, a3 = int(angle[2]), int(angle[3]), int(angle[4])
        if center not in zn_ids:
            continue
        v1 = vector_pbc(coords[center], coords[a1], supercell)
        v2 = vector_pbc(coords[center], coords[a3], supercell)
        angle_value = angle_degrees(v1, v2)
        if angle_value is not None:
            o_zn_o_angles.append(angle_value)

    def stats(values):
        if not values:
            return {"min": None, "mean": None, "max": None}
        return {
            "min": float(min(values)),
            "mean": float(sum(values) / len(values)),
            "max": float(max(values)),
        }

    return {
        "Zn_O_cutoff_A": ZN_O_CUTOFF,
        "Zn_O_coordination_numbers_cutoff_2p3A": [int(x) for x in zn_coordination],
        "Zn_O_distance_A": stats(zn_o_distances),
        "O_Zn_O_angle_deg": stats(o_zn_o_angles),
        "minimum_Zn_Zn_distance": None if not zn_zn_distances else float(min(zn_zn_distances)),
        "minimum_Si_O_distance": None if not si_o_distances else float(min(si_o_distances)),
        "minimum_Ca_O_distance": None if not ca_o_distances else float(min(ca_o_distances)),
    }


def h_overlap_metrics(hydroxylation_records):
    overlaps = []
    for record in hydroxylation_records or []:
        for oxy in record["hydroxylated_oxygens"]:
            if oxy["min_H_overlap_distance"] is not None:
                overlaps.append(float(oxy["min_H_overlap_distance"]))
    if not overlaps:
        return {"minimum_added_H_overlap_distance": None, "has_severe_H_overlap": False}
    return {
        "minimum_added_H_overlap_distance": float(min(overlaps)),
        "has_severe_H_overlap": bool(min(overlaps) < H_OVERLAP_CUTOFF),
    }


def hydroxylation_topology_audit(entries_crystal, entries_bonds, entries_angle, hydroxylation_records):
    atom_types = type_by_atom_id(entries_crystal)
    audit = []
    bad_records = []
    for record in hydroxylation_records or []:
        for oxy in record["hydroxylated_oxygens"]:
            o_id = int(oxy["oxygen_atom_id"])
            h_id = int(oxy["reused_shell_as_H_atom_id"])
            bonds = [list(bond) for bond in entries_bonds if int(bond[2]) in (o_id, h_id) or int(bond[3]) in (o_id, h_id)]
            shell_bonds = [
                bond for bond in bonds
                if int(bond[1]) == 3 and {atom_types.get(int(bond[2])), atom_types.get(int(bond[3]))} == {3, 4}
            ]
            oh_bonds = [bond for bond in bonds if int(bond[1]) == 1 and {int(bond[2]), int(bond[3])} == {o_id, h_id}]
            angles = [list(angle) for angle in entries_angle if o_id in [int(angle[2]), int(angle[3]), int(angle[4])] or h_id in [int(angle[2]), int(angle[3]), int(angle[4])]]
            item = {
                "oxygen_atom_id": o_id,
                "hoh_atom_id": h_id,
                "oxygen_type": atom_types.get(o_id),
                "hoh_type": atom_types.get(h_id),
                "has_correct_oh_hoh_bond": bool(oh_bonds),
                "remaining_core_shell_bonds": shell_bonds,
                "bonds_involving_pair": bonds,
                "angles_involving_pair": angles,
            }
            if item["oxygen_type"] != 6 or item["hoh_type"] != 8 or not item["has_correct_oh_hoh_bond"] or shell_bonds:
                bad_records.append(item)
            audit.append(item)
    return {"records": audit, "bad_records": bad_records}


def classify_zinc_output(summary, allow_unbalanced_for_debug=False):
    classification = "needs_minimization"
    reasons = []

    charge_residual = float(summary.get("charge_residual_final", summary.get("total_charge_residual", 0.0)))
    if abs(charge_residual) > CHARGE_TOLERANCE:
        classification = "debug_only_unbalanced"
        reasons.append("non-neutral charge residual")

    geometry = summary.get("pre_minimization_geometry", summary.get("geometry_validation", {}))
    zn_o = geometry.get("Zn_O_distance_A", {})
    if zn_o.get("min") is None:
        classification = "debug_only_bad_zinc_geometry"
        reasons.append("Zn has no O neighbors within 2.3 A")
    elif zn_o["min"] < 1.7 or zn_o["max"] > 2.3:
        classification = "debug_only_bad_zinc_geometry"
        reasons.append("pre-minimization Zn-O distances outside v2 sanity window 1.7-2.3 A")

    h_overlap = geometry.get("added_H_overlap", {})
    if h_overlap.get("has_severe_H_overlap"):
        classification = "debug_only_bad_hydrogen_geometry"
        reasons.append("newly added H has severe overlap")
    if h_overlap.get("has_v21_H_contact_violation"):
        classification = "debug_only_bad_hydrogen_geometry"
        reasons.append("newly added H violates v2.1 minimum contact cutoffs")

    topology = summary.get("topology_validation", {})
    if topology.get("stale_zn_angles"):
        if classification == "md_ready_candidate":
            classification = "debug_only_missing_topology"
        reasons.append("unmapped Zn angle entries remain")
    if topology.get("zinc_bonds"):
        if classification == "md_ready_candidate":
            classification = "debug_only_missing_topology"
        reasons.append("Zn bonded topology has no defined CementFF4 bond type")
    if topology.get("hydroxylation_topology_audit", {}).get("bad_records"):
        classification = "debug_only_bad_topology"
        reasons.append("converted O(S)->Oh-Hoh pair still has shell/topology remnants")

    if classification == "needs_minimization":
        reasons.append("charge/topology mapped; minimization required before MD-ready classification")

    summary["output_classification"] = classification
    summary["classification_reasons"] = reasons
    return summary


def finalize_zinc_summary(
    entries_crystal,
    entries_bonds,
    entries_angle,
    supercell,
    summary,
    charge_balance_mode="fail_if_not_neutral",
    allow_unbalanced_for_debug=False,
):
    if summary is None:
        return None

    summary["total_charge_after_hydroxylation"] = total_charge(entries_crystal)
    summary["charge_residual_final"] = total_charge(entries_crystal)
    topology = remap_zinc_angles(entries_crystal, entries_angle)
    topology["zinc_bonds"] = validate_no_zinc_bonds(entries_crystal, entries_bonds)
    topology["hydroxylation_topology_audit"] = hydroxylation_topology_audit(
        entries_crystal, entries_bonds, entries_angle, summary.get("hydroxylation_records", [])
    )
    summary["topology_validation"] = topology
    summary["pre_minimization_geometry"] = geometry_metrics(entries_crystal, entries_angle, supercell)
    summary["pre_minimization_geometry"]["added_H_overlap"] = h_overlap_metrics(
        summary.get("hydroxylation_records", [])
    )
    h_violations = []
    for record in summary.get("hydroxylation_records", []):
        for oxy in record["hydroxylated_oxygens"]:
            for key, cutoff in H_MIN_DISTANCES.items():
                distance = oxy.get("H_contact_metrics", {}).get(key, {}).get("distance")
                if distance is not None and distance < cutoff:
                    h_violations.append(
                        {
                            "oxygen_atom_id": oxy["oxygen_atom_id"],
                            "hoh_atom_id": oxy["reused_shell_as_H_atom_id"],
                            "metric": key,
                            "distance": float(distance),
                            "cutoff": float(cutoff),
                        }
                    )
    summary["pre_minimization_geometry"]["added_H_overlap"]["v21_H_contact_violations"] = h_violations
    summary["pre_minimization_geometry"]["added_H_overlap"]["has_v21_H_contact_violation"] = bool(h_violations)
    summary["Zn_charge_balance_mode"] = charge_balance_mode
    summary["allow_unbalanced_for_debug"] = bool(allow_unbalanced_for_debug)
    summary["cementff4_type_mapping"] = CEMENTFF4_TYPE_MAP
    summary["cementff4_angle_mapping"] = CEMENTFF4_ANGLE_MAP
    summary["status_note"] = (
        "v2 uses a charge-balanced ZnO2(OH)2 substitutional motif, but minimization is required before MD."
    )
    summary = classify_zinc_output(summary, allow_unbalanced_for_debug)
    return summary


def build_zinc_summary(
    entries_crystal,
    selected_sites,
    candidates,
    target_zinc_si_ratio,
    ca_si_ratio,
    supercell,
    site_type,
    seed,
):
    counts = count_species(entries_crystal)
    n_si = counts.get(2, 0) + counts.get(10, 0)
    n_zn = counts.get(ZN_SPECIE, 0)
    n_ca = counts.get(1, 0) + counts.get(9, 0)
    n_si_original = n_si + n_zn
    min_zn_zn = minimum_periodic_distance(
        selected_sites[0]["coord"] if selected_sites else [0.0, 0.0, 0.0],
        [site["coord"] for site in selected_sites[1:]] if len(selected_sites) > 1 else [],
        supercell,
    )
    min_zn_o = math.inf
    for site in selected_sites:
        min_zn_o = min(
            min_zn_o,
            minimum_distance_to_species(site["coord"], entries_crystal, {3, 4, 5, 6, 11, 12}, supercell),
        )

    total_charge = float(sum(float(entry[2]) for entry in entries_crystal))
    denominator = n_si + n_zn

    return {
        "enable_zinc": True,
        "Zn_site_type": site_type,
        "Zn_seed": int(seed),
        "target_Zn_Si_ratio": float(target_zinc_si_ratio),
        "actual_Zn_Si_ratio": float(n_zn / n_si) if n_si else None,
        "actual_Zn_Si_original_ratio": float(n_zn / n_si_original) if n_si_original else None,
        "Ca_Si_ratio": float(ca_si_ratio),
        "Ca_over_Si_plus_Zn_ratio": float(n_ca / denominator) if denominator else None,
        "N_Si_original": int(n_si_original),
        "N_Si_final": int(n_si),
        "N_Si": int(n_si),
        "N_Zn": int(n_zn),
        "N_Ca": int(n_ca),
        "N_Q1_Zn": int(sum(1 for site in selected_sites if site["motif"] == "Q1_Zn")),
        "N_Q2b_Zn": int(sum(1 for site in selected_sites if site["motif"] == "Q2b_Zn")),
        "N_Q1_candidates": int(len(candidates["Q1_Zn"])),
        "N_Q2b_candidates": int(len(candidates["Q2b_Zn"])),
        "selected_sites": [
            {
                "atom_id": site["atom_id"],
                "motif": site["motif"],
                "cell": site["cell"],
                "piece": site["piece"],
                "coord": site["coord"],
                "original_specie": site["original_specie"],
            }
            for site in selected_sites
        ],
        "minimum_Zn_Zn_distance": None if math.isinf(min_zn_zn) else float(min_zn_zn),
        "minimum_Zn_O_distance": None if math.isinf(min_zn_o) else float(min_zn_o),
        "total_charge_before_zinc": None,
        "total_charge_after_zinc_before_hydroxylation": total_charge,
        "total_charge_after_hydroxylation": total_charge,
        "total_charge_residual": total_charge,
        "charge_residual_final": total_charge,
        "N_Os_converted_to_Oh": 0,
        "N_H_added_for_Zn_OH": 0,
        "hydroxylation_records": [],
        "Ca_Si_original": float(n_ca / n_si_original) if n_si_original else None,
        "Ca_Si_final": float(n_ca / n_si) if n_si else None,
        "Zn_Si_original": float(n_zn / n_si_original) if n_si_original else None,
        "Zn_Si_final": float(n_zn / n_si) if n_si else None,
        "notes": [
            "v2 creates a charge-balanced ZnO2(OH)2 substitutional candidate.",
            "v2 does not use guest_ions/substitute and does not randomly replace Ca-layer atoms.",
            "Zn parameters are taken from CementFF4 supplementary information.",
        ],
    }


def apply_zinc_modification(
    entries_crystal,
    crystal_dict,
    supercell,
    Zn_Si_ratio,
    Zn_site_type,
    Zn_seed,
    ca_si_ratio,
    charge_balance_mode="hydroxylate_two_oxygens",
    entries_bonds=None,
    entries_angle=None,
    allow_hydroxylate_bridging_oxygen=False,
    precondition_zinc_geometry=True,
    target_Zn_O_distance=1.95,
):
    validate_zinc_site_type(Zn_site_type)
    validate_charge_balance_mode(charge_balance_mode)
    if charge_balance_mode == "hydroxylate_two_oxygens" and (entries_bonds is None or entries_angle is None):
        raise ValueError("hydroxylate_two_oxygens requires entries_bonds and entries_angle")
    candidates = inspect_zinc_candidates(crystal_dict)
    counts = count_species(entries_crystal)
    n_si_initial = counts.get(2, 0) + counts.get(10, 0)
    if n_si_initial <= 0:
        raise ValueError("Cannot place Zn because the generated structure contains no Si atoms")

    n_zinc = int(round(float(Zn_Si_ratio) * n_si_initial))
    if float(Zn_Si_ratio) > 0.0 and n_zinc == 0:
        n_zinc = 1

    site_filter = None
    if charge_balance_mode == "hydroxylate_two_oxygens":
        def site_filter(site):
            try:
                trial_entries = copy.deepcopy(entries_crystal)
                trial_bonds = copy.deepcopy(entries_bonds)
                trial_angles = copy.deepcopy(entries_angle)
                for entry in trial_entries:
                    if int(entry[0]) == int(site["atom_id"]):
                        entry[1] = ZN_SPECIE
                        entry[2] = ZN_CHARGE
                        break
                hydroxylate_two_oxygens(
                    trial_entries,
                    trial_bonds,
                    trial_angles,
                    [site],
                    supercell,
                    allow_hydroxylate_bridging_oxygen,
                    precondition_zinc_geometry,
                    target_Zn_O_distance,
                )
                return True
            except ValueError:
                return False

    selected_sites = select_zinc_sites(candidates, n_zinc, Zn_site_type, Zn_seed, supercell, site_filter=site_filter)
    charge_before_zinc = total_charge(entries_crystal)
    entries_crystal, crystal_dict = apply_zinc_sites(entries_crystal, crystal_dict, selected_sites)
    charge_after_zinc = total_charge(entries_crystal)
    hydroxylation_records = apply_charge_balance(
        entries_crystal,
        entries_bonds,
        entries_angle,
        selected_sites,
        supercell,
        charge_balance_mode,
        allow_hydroxylate_bridging_oxygen,
        precondition_zinc_geometry,
        target_Zn_O_distance,
    )
    summary = build_zinc_summary(
        entries_crystal,
        selected_sites,
        candidates,
        Zn_Si_ratio,
        ca_si_ratio,
        supercell,
        Zn_site_type,
        Zn_seed,
    )
    summary["total_charge_before_zinc"] = charge_before_zinc
    summary["total_charge_after_zinc_before_hydroxylation"] = charge_after_zinc
    summary["total_charge_after_hydroxylation"] = total_charge(entries_crystal)
    summary["total_charge_residual"] = total_charge(entries_crystal)
    summary["charge_residual_final"] = total_charge(entries_crystal)
    summary["N_Os_converted_to_Oh"] = sum(len(record["hydroxylated_oxygens"]) for record in hydroxylation_records)
    summary["N_H_added_for_Zn_OH"] = summary["N_Os_converted_to_Oh"]
    summary["hydroxylation_records"] = hydroxylation_records
    summary["allow_hydroxylate_bridging_oxygen"] = bool(allow_hydroxylate_bridging_oxygen)
    summary["precondition_zinc_geometry"] = bool(precondition_zinc_geometry)
    summary["target_Zn_O_distance"] = float(target_Zn_O_distance)
    return entries_crystal, crystal_dict, summary


def write_zinc_summary(path, summary):
    with open(path, "w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

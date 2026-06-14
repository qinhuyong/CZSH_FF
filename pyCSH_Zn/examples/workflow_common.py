from __future__ import print_function

import os
import random
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from mod_construct_brick_Y import get_all_bricks, pieces
from mod_construct_supercell_Y import check_move_water_hydrogens, get_angles, get_full_coordinates, resize_crystal
from mod_sample import fill_water, sample_Ca_Si_ratio
from mod_write_Y import get_lammps_input_cementff, write_cementff4_mapping_json, write_cementff4_zinc_input
from mod_zinc import apply_zinc_modification, finalize_zinc_summary, write_zinc_summary


UNITCELL = np.array(
    [
        [6.7352, 0.0, 0.0],
        [-4.071295, 6.209521, 0.0],
        [0.7037701, -6.2095578, 13.9936836],
    ]
)


def generate_structure(output_dir, prefix, enable_zinc=False, zn_ratio=0.03, site_type="Q2b_Zn", seed=23743):
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    np.random.seed(seed)
    random.seed(seed + 10)
    size = (2, 2, 2)
    supercell = np.zeros((3, 3))
    for i in range(3):
        supercell[i, :] = UNITCELL[i, :] * size[i]
    bricks, sorted_bricks = get_all_bricks(pieces)
    crystal, n_ca, n_si, r_sioh, r_caoh, mcl, n_water, _ = sample_Ca_Si_ratio(
        sorted_bricks, 1.7, 0.2, size[0] * size[1] * size[2], [0.06, 0.08, 0.08]
    )
    water_in_crystal = fill_water(crystal, n_water)
    crystal_rs, water_in_crystal_rs = resize_crystal(crystal, water_in_crystal, size)
    entries, bonds, crystal_dict, water_dict = get_full_coordinates(
        crystal_rs, water_in_crystal_rs, size, pieces, False, [0]
    )
    angles = get_angles(crystal_dict, water_dict, size)
    zinc_summary = None
    if enable_zinc:
        entries, crystal_dict, zinc_summary = apply_zinc_modification(
            entries, crystal_dict, supercell, zn_ratio, site_type, seed, n_ca / n_si,
            "hydroxylate_two_oxygens", bonds, angles, False, True, 1.95
        )
        zinc_summary = finalize_zinc_summary(
            entries, bonds, angles, supercell, zinc_summary, "hydroxylate_two_oxygens", False
        )
    entries, _, _ = check_move_water_hydrogens(entries)
    data_name = prefix + ("_cementff_zn.data" if enable_zinc else "_cementff.data")
    data_path = os.path.join(output_dir, data_name)
    water_summary = get_lammps_input_cementff(data_path, entries, bonds, angles, supercell, zinc_summary)
    mapping_path = os.path.join(output_dir, "cementff_mapping_summary.json")
    write_cementff4_mapping_json(mapping_path, enable_zinc)
    water_path = os.path.join(output_dir, "water_summary.json")
    with open(water_path, "w") as f:
        import json
        json.dump(water_summary, f, indent=2, sort_keys=True)
        f.write("\n")
    result = {
        "data_file": data_path,
        "mapping_summary": mapping_path,
        "water_summary": water_path,
        "zinc_summary": None,
        "forcefield": None,
    }
    if enable_zinc:
        zinc_path = os.path.join(output_dir, "zinc_summary.json")
        write_zinc_summary(zinc_path, zinc_summary)
        ff_path = os.path.join(output_dir, "in.CementFF4_Zn")
        write_cementff4_zinc_input(ff_path)
        result["zinc_summary"] = zinc_path
        result["forcefield"] = ff_path
    return result

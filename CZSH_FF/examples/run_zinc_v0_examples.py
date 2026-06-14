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
from mod_write_Y import write_output
from mod_zinc import apply_zinc_modification, finalize_zinc_summary


def run_case(name, enable_zinc=False, zinc_ratio=0.0, site_type="mixed_Q1_Q2b_Zn"):
    seed = 23137
    np.random.seed(seed)
    random.seed(seed + 10)

    size = (2, 2, 2)
    ca_si_ratio = 1.7
    w_si_ratio = 0.2
    widths = [0.06, 0.08, 0.08]
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

    bricks, sorted_bricks = get_all_bricks(pieces)
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

    zinc_summary = None
    if enable_zinc:
        entries_crystal, crystal_dict, zinc_summary = apply_zinc_modification(
            entries_crystal,
            crystal_dict,
            supercell,
            zinc_ratio,
            site_type,
            seed,
            n_ca / n_si,
            "hydroxylate_two_oxygens",
            entries_bonds,
            entries_angle,
            False,
            True,
            1.95,
        )

    if zinc_summary is not None:
        zinc_summary = finalize_zinc_summary(
            entries_crystal,
            entries_bonds,
            entries_angle,
            supercell,
            zinc_summary,
            "hydroxylate_two_oxygens",
            False,
        )
        print("{} classification: {}".format(name, zinc_summary["output_classification"]))
    entries_crystal, _, _ = check_move_water_hydrogens(entries_crystal)
    write_output(
        0,
        entries_crystal,
        entries_bonds,
        entries_angle,
        size,
        crystal_rs,
        water_in_crystal_rs,
        supercell,
        n_ca,
        n_si,
        r_sioh,
        r_caoh,
        mcl,
        False,
        False,
        False,
        False,
        name,
        unitcell,
        False,
        False,
        True,
        None,
        False,
        [0],
        False,
        True,
        zinc_summary,
        True,
    )
    return zinc_summary


if __name__ == "__main__":
    os.makedirs("output_Y", exist_ok=True)
    cases = [
        ("example_pure_csh", False, 0.0, "mixed_Q1_Q2b_Zn"),
        ("example_q1_zn", True, 0.03, "Q1_Zn"),
        ("example_q2b_zn", True, 0.03, "Q2b_Zn"),
        ("example_mixed_zn", True, 0.03, "mixed_Q1_Q2b_Zn"),
    ]
    for case in cases:
        summary = run_case(*case)
        if summary is None:
            print(case[0] + ": pure C-S-H")
        else:
            print(
                "{}: N_Zn={} actual_Zn_Si={:.6f} Q1={} Q2b={}".format(
                    case[0],
                    summary["N_Zn"],
                    summary["actual_Zn_Si_ratio"],
                    summary["N_Q1_Zn"],
                    summary["N_Q2b_Zn"],
                )
            )

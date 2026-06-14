from __future__ import print_function

import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from forcefields.build_cementff4_zn import build as build_forcefield
from lammps_templates.build_inputs import build as build_inputs


if __name__ == "__main__":
    q2b_dir = os.path.join("output_Y", "workflow_v1", "q2b_zn")
    ff_result = build_forcefield(q2b_dir)
    inputs = build_inputs(
        os.path.join(q2b_dir, "q2b_zn_cementff_zn.data"),
        ff_result["forcefield"],
        os.path.join(q2b_dir, "lammps_inputs"),
    )
    print(json.dumps({"forcefield": ff_result, "inputs": inputs}, indent=2, sort_keys=True))

from __future__ import print_function

import argparse
import json
import os


def write(path, lines):
    with open(path, "w") as f:
        f.write("\n".join(lines))
        f.write("\n")


def header(data_file, ff_file, run_dir):
    data_ref = os.path.relpath(os.path.abspath(data_file), os.path.abspath(run_dir))
    ff_ref = os.path.relpath(os.path.abspath(ff_file), os.path.abspath(run_dir))
    return [
        "clear",
        "units metal",
        "dimension 3",
        "atom_style full",
        "boundary p p p",
        "box tilt large",
        "fix csinfo all property/atom i_CSID",
        "read_data {}".format(data_ref) + " fix csinfo NULL CS-Info",
        "include {}".format(ff_ref),
        "neighbor 2.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        "comm_modify vel yes cutoff 14.0",
        "thermo 100",
        "thermo_style custom step pe ebond eangle evdwl ecoul elong fnorm fmax press",
    ]


def build(data_file, ff_file, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    outputs = {}
    outputs["read_check"] = os.path.join(out_dir, "in.read_check")
    write(outputs["read_check"], header(data_file, ff_file, out_dir) + ["run 0"])

    outputs["minimize_static"] = os.path.join(out_dir, "in.minimize_static")
    write(outputs["minimize_static"], header(data_file, ff_file, out_dir) + [
        "min_style cg",
        "min_modify dmax 0.01 line quadratic",
        "minimize 1e-6 1e-8 1000 10000",
        "write_data minimized_static.data nocoeff",
    ])

    outputs["static_relax_shell"] = os.path.join(out_dir, "in.static_relax_shell")
    write(outputs["static_relax_shell"], header(data_file, ff_file, out_dir) + [
        "group mobile type 4",
        "group fixed subtract all mobile",
        "fix freeze fixed setforce 0.0 0.0 0.0",
        "min_style fire",
        "min_modify dmax 0.001",
        "minimize 1e-8 1e-10 500 5000",
        "unfix freeze",
        "write_data shell_relaxed_static.data nocoeff",
    ])

    outputs["elastic_quasistatic"] = os.path.join(out_dir, "in.elastic_quasistatic")
    write(outputs["elastic_quasistatic"], header(data_file, ff_file, out_dir) + [
        "# Quasi-static deformation template; edit variable strain before use.",
        "variable strain equal 0.001",
        "change_box all x scale ${strain} remap",
        "min_style cg",
        "min_modify dmax 0.002",
        "minimize 1e-6 1e-8 1000 10000",
        "write_data elastic_step.data nocoeff",
    ])

    outputs["short_md_experimental"] = os.path.join(out_dir, "in.short_md_experimental")
    write(outputs["short_md_experimental"], header(data_file, ff_file, out_dir) + [
        "# EXPERIMENTAL: finite-temperature core-shell MD protocol is not validated.",
        "# Use only after separately validating the core-shell dynamics protocol.",
        "compute CStemp all temp/cs cores shells",
        "thermo_modify temp CStemp",
        "timestep 0.000028",
        "velocity all create 10 134 dist gaussian mom yes rot no bias yes temp CStemp",
        "velocity all scale 10 temp CStemp",
        "fix thermostat all temp/berendsen 10 10 0.028",
        "fix_modify thermostat temp CStemp",
        "fix nve_all all nve",
        "run 1000",
    ])
    manifest = os.path.join(out_dir, "lammps_input_manifest.json")
    with open(manifest, "w") as f:
        json.dump(outputs, f, indent=2, sort_keys=True)
        f.write("\n")
    outputs["manifest"] = manifest
    return outputs


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True)
    parser.add_argument("--ff", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.data, args.ff, args.out), indent=2))


if __name__ == "__main__":
    main()

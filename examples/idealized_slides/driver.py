# Set up numerical experiment. Create one directory for each simulation and place
# a params.yaml file within it.
# Generate run_directory.csv, a summary csv listing all simulations.

import itertools
import os
from distutils.dir_util import copy_tree
from pathlib import Path

import numpy as np
import yaml

run_list = []
run_id = 1

AMR = False


def make_sim(run_id, kr, m0, phi, depth, src2method, order, riemann_method, params):

    sim_name = f"IDEAL_{str(run_id).zfill(3)}_K{int(np.log10(kr)*-1*100)}_m{int(m0*100)}_P{int(phi)}_depth{depth}_src{src2method[0]}_O{order}_R{riemann_method}"

    # Create unique simulation directory within results
    dst = f"results/{sim_name}"

    run_list.append((run_id, sim_name))

    print(f"Creating {sim_name}")

    if not Path(dst).exists():
        Path(dst).mkdir(parents=True)

    copy_tree("template", dst)

    os.chdir(dst)

    # Create yaml file containing parameters for this simulation.
    # these parameters will be read in by a template setrun.py to
    # run the simulation.
    params2 = {
        "run_id": run_id,
        "kr": kr,
        "phi": phi,
        "m0": m0,
        "depth": depth,
        "src2method": src2method[0],
        "alphamethod": src2method[1],
        "order": order,
        "riemann_method": riemann_method,
        "amr": AMR,
    }

    with open("params.yml", "w") as file:
        yaml.dump(params | params2, file, sort_keys=False)
    os.chdir("../..")


# params for all simulations
with open("params.yml", "r") as f:
    params = yaml.safe_load(f)

# Material parameter values
src2methods = [(-1, 0), (0, 0), (2, 1)]
orders = [1, 2]
riemann_methods = [0, 1]
krs = [
    1e-10,
    1e-9,
    1e-8,
]
m0s = [0.62]
phis = [
    35,
]
depths = [1, 5, 10]


krs = [1e-9]
m0s = [0.62]
phis = [35]
depths = [5]

for kr, phi, m0, depth, src2method, order, riemann_method in itertools.product(
    krs, phis, m0s, depths, src2methods, orders, riemann_methods
):
    make_sim(
        run_id,
        kr,
        m0,
        phi,
        depth,
        src2method,
        order,
        riemann_method,
        params=params,
    )
    run_id += 1

# write out a csv file that lists the simulations considered. This is used to
# manage job submission as a batch job in slurm.
with open("run_directory.csv", "w") as f:
    f.write("run_id,run_no,run_name\n")
    for run_id, run_name in run_list:
        f.write(f"{run_id},{run_id},{run_name}\n")

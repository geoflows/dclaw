# Set up numerical experiment. Create one directory for each simulation and place
# a params.yaml file within it.
# Generate run_directory.csv, a summary csv listing all simulations.

import os
from distutils.dir_util import copy_tree
from pathlib import Path

import numpy as np
import yaml

run_list = []
run_id = 1

AMR = False


def make_sim(run_id, kr, m0, phi, depth, src2method, order, params):

    sim_name = f"IDEAL_{str(run_id).zfill(3)}_K{int(np.log10(kr)*-1*100)}_m{int(m0*100)}_P{int(phi)}_depth{depth}_src{src2method[0]}_O{order}"

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
        "amr": AMR,
    }

    with open("params.yml", "w") as file:
        yaml.dump(params | params2, file, sort_keys=False)
    os.chdir("../..")


params = dict(
    mc=0.64,
    # geometry
    slope=40,
    length_scale=100,
    depth=5,
    xlower=0,
    xupper=500,
    ylower=-100,
    yupper=100,
    xls=50,
    yls=0,
    rls=10,
    # for interpolation of initial input data
    nxpoints=201,
    nypoints=501,
    # initial grids
    num_cells_x=30,
    num_cells_y=20,
)

# Material parameter values
# src2methods = [(-1, 0), (0,0), (2,1)]
src2methods = [
    (2, 1),
]

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

for kr in krs:
    for phi in phis:
        for m0 in m0s:
            for depth in depths:
                for src2method in src2methods:
                    make_sim(run_id, kr, m0, phi, depth, src2method, order=1, params=params)
                    run_id += 1

# write out a csv file that lists the simulations considered. This is used to
# manage job submission as a batch job in slurm.
with open("run_directory.csv", "w") as f:
    f.write("run_id,run_no,run_name\n")
    for run_id, run_name in run_list:
        f.write(f"{run_id},{run_id},{run_name}\n")

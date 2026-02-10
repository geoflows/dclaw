# Run all simulations, one at a time.

import os
import time

import pandas as pd

df = pd.read_csv("run_directory.csv")

for ind, row in df.iterrows():
    print(row.run_id, row.run_no, row.run_name, time.ctime())
    os.chdir(f"results/{row.run_name}")

    if not os.path.exists("conservation.csv"):

        print("input")
        os.system("make input > log_input.txt")
        print("data")
        os.system("make .data > log_data.txt")
        print("output")
        os.system("make .output > log_output.txt")
        print("plot")
        os.system("make .plots > log_plots.txt")
        print("post")
        os.system("python pyvista_plotting.py > log_pyvista.txt")
        print("post")
        os.system("make postprocess > log_post.txt")

        # os.system("rm -r _plots")
        # os.system("rm -r _frames")
        # os.system("rm _output/*fort*")

    os.chdir("../..")

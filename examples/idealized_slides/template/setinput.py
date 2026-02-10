import numpy as np
import yaml
from clawpack.clawutil.data import ClawData
from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import Rearth  # radius of earth
from pylab import *
from scipy.interpolate import interp1d

with open("params.yml", "r") as file:
    params = yaml.safe_load(file)

# problem geometry from yaml
depth = params["depth"]
length_scale = params["length_scale"]

slope = params["slope"]
xlower = params["xlower"]
ylower = params["ylower"]
xupper = params["xupper"]
yupper = params["yupper"]
xls = params["xls"]
yls = params["yls"]
rls = params["rls"]
nxpoints = params["nxpoints"]
nypoints = params["nypoints"]


def make_plots():

    basal = topotools.Topography("basal_topo.tt3", 3)
    basal.plot()
    title("Basal topo")
    fname = "basal_topo.png"
    savefig(fname)
    print("Created ", fname)

    eta = topotools.Topography("surface_topo.tt3", 3)
    eta.plot()
    title("Surface topo eta")
    fname = "surface_topo.png"
    savefig(fname)
    print("Created ", fname)

    h = eta.Z - basal.Z
    figure()
    pcolormesh(eta.X, eta.Y, h, cmap="Reds")
    axis("equal")
    colorbar()
    title("Landslide depth")
    fname = "landslide_depth.png"
    savefig(fname)
    print("Created ", fname)


def basal(x, y):
    """
    Cartesian: x,y in meters
    """
    b = length_scale
    a = b * np.tan(np.deg2rad(slope))
    z = a * np.exp(-x / b)
    return z


def eta(x, y):
    """
    Cartesian: x,y in meters
    """

    r = np.sqrt((x - xls) ** 2 + (y - yls) ** 2)
    is_ls = r <= rls

    z = basal(x, y)
    z[is_ls] += depth
    return z


def maketopo():
    """
    Output topography file for the entire domain
    """
    outfile = "basal_topo.tt3"
    topotools.topo3writer(
        outfile, basal, xlower, xupper, ylower, yupper, nxpoints, nypoints
    )


def make_surface():
    """
    Output surface topography file for the entire domain
    (Could be for smaller region)
    """
    outfile = "surface_topo.tt3"
    topotools.topo3writer(
        outfile, eta, xlower, xupper, ylower, yupper, nxpoints, nypoints
    )


if __name__ == "__main__":
    maketopo()
    make_surface()
    make_plots()

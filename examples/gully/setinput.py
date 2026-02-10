import numpy as np
from clawpack.clawutil.data import ClawData
from clawpack.geoclaw import topotools
from clawpack.geoclaw.data import Rearth  # radius of earth
from clawpack.visclaw import colormaps
from pylab import *
from scipy.interpolate import interp1d

rotate = False

# Create a gully that exits onto a fan
# landslide located on a side slop of the gully
dx = 25


if rotate:
    x0 = 0  # right side of domain
    x1 = 3000  # location of gully-fan intersection
    x2 = 6000  # right side of fan

    y0 = 0  # bottom of domain
    y1 = 1000  # center of channel
    y2 = 2000  # top of domain
    yc1 = 1025  # left and right side of channels
    yc2 = 9075

else:
    y0 = 0  # right side of domain
    y1 = 3000  # location of gully-fan intersection
    y2 = 6000  # right side of fan

    x0 = 0  # bottom of domain
    x1 = 1000  # center of channel
    x2 = 2000  # top of domain
    xc1 = 1025  # left and right side of channels
    xc2 = 9075

nx = int((x2 - x0) / dx)
ny = int((y2 - y0) / dx)

alpha_g = 45  # gully side slope angle
alpha_c = 15  # gully channel slope angle
alpha_f = 5  # fan slope angle

# landslide
if rotate:
    xl1, xl2 = 2100, 2200  # landslide x extent
    yl1, yl2 = 700, 800  # landslide y extent
else:
    yl1, yl2 = 2100, 2200  # landslide x extent
    xl1, xl2 = 700, 800  # landslide y extent

depth = 4  # landslide depth


def make_plots():

    basal = topotools.Topography("basal_topo.tt3", 3)
    basal.plot(long_lat=False)
    title("Basal topo")
    fname = "basal_topo.png"
    savefig(fname)
    print("Created ", fname)

    eta = topotools.Topography("surface_topo.tt3", 3)
    eta.plot(long_lat=False)
    title("Surface topo eta")
    fname = "surface_topo.png"
    savefig(fname)
    print("Created ", fname)

    h = eta.Z - basal.Z
    figure()
    pcolormesh(eta.X, eta.Y, h, cmap='bone')
    xlim(x0, x2)
    ylim(y0, y2)
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
    z = np.zeros_like(x)

    # slope down at the fan slope to the right of x1
    if rotate:
        fan_area = x > x1
        z[fan_area] -= np.sin(np.deg2rad(alpha_f)) * (x[fan_area] - x1)

        # slope up at the channel slope to the left of x1
        gully_area = x <= x1
        z[gully_area] -= np.sin(np.deg2rad(alpha_c)) * (x[gully_area] - x1)

        # add the channel side slopes
        # only add them outside the channl
        side_slopes = (x <= x1) * ((y < yc1) | (y > yc2))
        z[gully_area] += np.sin(np.deg2rad(alpha_g)) * np.abs(y[gully_area] - y1)
    else:
        fan_area = y > y1
        z[fan_area] -= np.sin(np.deg2rad(alpha_f)) * (y[fan_area] - y1)

        # slope up at the channel slope to the left of y1
        gully_area = y <= y1
        z[gully_area] -= np.sin(np.deg2rad(alpha_c)) * (y[gully_area] - y1)

        # add the channel side slopes
        # only add them outside the channl
        side_slopes = (y <= y1) * ((x < xc1) | (x > xc2))
        z[gully_area] += np.sin(np.deg2rad(alpha_g)) * np.abs(x[gully_area] - x1)

    return z


def thickness(x, y):
    z = np.zeros_like(x)

    landslide_location = (x > xl1) * (x <= xl2) * (y > yl1) * (y <= yl2)
    z[landslide_location] = depth

    dx = np.diff(x, axis=1)[0, 0]
    dy = np.diff(y, axis=0)[0, 0]
    volume = np.sum(z) * dx * dy

    print(f"volume: {volume} m**2")

    return z


def eta(x, y):
    """
    Cartesian: x,y in meters
    """
    z = basal(x, y) + thickness(x, y)
    return z


def maketopo():
    """
    Output topography file for the entire domain
    """

    for filename, function in [
        ("basal_topo.tt3", basal),
        ("surface_topo.tt3", eta),
        ("thickness.tt3", thickness),
    ]:
        topography = topotools.Topography(topo_func=function)
        topography.x = np.linspace(x0, x2, nx * 2 + 1)
        topography.y = np.linspace(y0, y2, ny * 2 + 1)
        topography.write(filename, topo_type=3)


if __name__ == "__main__":
    maketopo()
    make_plots()

"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import

import numpy as np
from clawpack.geoclaw.topotools import Topography

# from pyclaw.data import Data
# probdata = Data('setprob.data')

a = 1.0
sigma = 0.5
h0 = 0.1
grav = 9.81
omega = np.sqrt(2.0 * grav * h0) / a


def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 200
    nypoints = 200
    xupper = 2.0e0
    yupper = 2.0e0
    xlower = -2.0e0
    ylower = -2.0e0
    outfile = "bowl.topotype2"

    topography = Topography(topo_func=topo)
    topography.x = np.linspace(xlower, xupper, nxpoints)
    topography.y = np.linspace(ylower, yupper, nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")


def topo(x, y):
    """
    Parabolic bowl
    """
    z = h0 * (x**2 + y**2) / a**2 - h0
    return z


if __name__ == "__main__":
    maketopo()

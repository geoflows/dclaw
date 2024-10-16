
from pylab import *
from clawpack.geoclaw import topotools
from clawpack.visclaw import colormaps
from clawpack.geoclaw.data import Rearth  # radius of earth
from clawpack.clawutil.data import ClawData
from scipy.interpolate import interp1d
import numpy as np


# Create a gully that exits onto a fan
# landslide located on a side slop of the gully

x0 = 0 # right side of domain
x1 = 3000 # location of gully-fan intersection
x2 = 6000 # right side of fan

y0 = 0 # bottom of domain
y1 = 1000 # center of channel
y2 = 2000 # top of domain
yc1 = 1025 # left and right side of channels
yc2 = 9075

alpha_g = 45 # gully side slope angle
alpha_c = 15 # gully channel slope angle
alpha_f = 5 # fan slope angle

# landslide
xl1, xl2 = 2100, 2150 # landslide x extent
yl1, yl2 = 700, 800 # landslide y extent
depth = 4 # landslide depth

def make_plots():

    basal = topotools.Topography('basal_topo.tt3',3)
    basal.plot()
    title('Basal topo')
    fname = 'basal_topo.png'
    savefig(fname)
    print('Created ',fname)

    eta = topotools.Topography('surface_topo.tt3',3)
    eta.plot()
    title('Surface topo eta')
    fname = 'surface_topo.png'
    savefig(fname)
    print('Created ',fname)

    h = eta.Z - basal.Z
    figure()
    pcolormesh(eta.X,eta.Y,h,cmap=colormaps.white_red)
    axis('equal')
    colorbar()
    title('Landslide depth')
    fname = 'landslide_depth.png'
    savefig(fname)
    print('Created ',fname)


def basal(x,y):
    """
    Cartesian: x,y in meters
    """
    z = np.zeros_like(x)

    # slope down at the fan slope to the right of x1
    fan_area = x>x1
    z[fan_area] -= np.sin(np.deg2rad(alpha_f))*(x[fan_area]-x1)

    # slope up at the channel slope to the left of x1
    gully_area = x<=x1
    z[gully_area] -= np.sin(np.deg2rad(alpha_c))*(x[gully_area]-x1)

    # add the channel side slopes
    # only add them outside the channl
    side_slopes = (x<=x1)*((y<yc1)|(y>yc2))
    z[gully_area] += np.sin(np.deg2rad(alpha_g))*np.abs(y[gully_area]-y1)

    return z


def eta(x,y):
    """
    Cartesian: x,y in meters
    """
    z = basal(x,y)
    landslide_location = (x>xl1)*(x<=xl2)*(y>yl1)*(y<yl2)
    z[landslide_location] += depth

    return z


def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 501
    nypoints = 501
    xlower= x0
    xupper= x2
    ylower= y0
    yupper= y2
    outfile= "basal_topo.tt3"
    topotools.topo3writer(outfile,basal,xlower,xupper,ylower,yupper,nxpoints,nypoints)
    outfile= "surface_topo.tt3"
    topotools.topo3writer(outfile,eta,xlower,xupper,ylower,yupper,nxpoints,nypoints)

if __name__=='__main__':
    maketopo()
    make_plots()

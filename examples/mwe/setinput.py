
import numpy as np
from scipy.interpolate import interp1d
import os

from clawpack.geoclaw import topotools, marching_front, dtopotools
from clawpack.visclaw import plottools


inputdir = 'input_files'
os.system('mkdir -p %s' % inputdir)
print('Input files will be put in directory %s' % inputdir)


### Create topography as a simple slope

# Define piecewise linear function (unequally spaced):
xland = np.array([0, 100, 200, 300, 400])
zland = np.array([100, 50, 25, 12, 0])

# Set dx, ymin, and ymax
dx = 10
ymin = 0
ymax = 50

# Interpolate to equally spaced grid for topofile:
xo = np.arange(xland.min(), xland.max(), dx)

yo = np.array([ymin, ymax])
zfunc = interp1d(xland,zland,fill_value="extrapolate")
zo = zfunc(xo)

# Convert to 2d arrays:
Xo,Yo = np.meshgrid(xo,yo)
Zo = np.vstack((zo,zo))

# ### Save as a topofile:
topo = topotools.Topography()
topo.set_xyZ(xo,yo,Zo)

topofile = '%s/topo.tt3' % inputdir
topo.write(topofile, topo_type=3, Z_format="%11.3e")
print('Created ', topofile)

### Create initial value for q1, depth

# Define piecewise linear function (unequally spaced):
xdepth = np.array([100, 120])
zdepth = np.array([5, 5])

# Set dx, ymin, and ymax
dx = 2
ymin = 20
ymax = 30

# Interpolate to equally spaced grid for topofile:
xo = np.arange(xdepth.min(), xdepth.max(), dx)

yo = np.array([ymin, ymax])
zfunc = interp1d(xdepth,zdepth,fill_value="extrapolate")
zo = zfunc(xo)

# Convert to 2d arrays:
Xo,Yo = np.meshgrid(xo,yo)
Zo = np.vstack((zo,zo))

# ### Save as a topofile:
topo = topotools.Topography()
topo.set_xyZ(xo,yo,Zo)

topofile = '%s/q1.tt3' % inputdir
topo.write(topofile, topo_type=3, Z_format="%11.3e")
print('Created ', topofile)


### Create initial value for aux5, depth of entrainable material

# Define piecewise linear function (unequally spaced):
xdepth = np.array([25, 50])
zdepth = np.array([5, 5])

# Set dx, ymin, and ymax
dx = 2
ymin = 10
ymax = 40

# Interpolate to equally spaced grid for topofile:
xo = np.arange(xdepth.min(), xdepth.max(), dx)

yo = np.array([ymin, ymax])
zfunc = interp1d(xdepth,zdepth,fill_value="extrapolate")
zo = zfunc(xo)

# Convert to 2d arrays:
Xo,Yo = np.meshgrid(xo,yo)
Zo = np.vstack((zo,zo))

# ### Save as a topofile:
topo = topotools.Topography()
topo.set_xyZ(xo,yo,Zo)

topofile = '%s/aux5.tt3' % inputdir
topo.write(topofile, topo_type=3, Z_format="%11.3e")
print('Created ', topofile)

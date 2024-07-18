from pylab import *
import glob, os
import pyvista as pv
from clawpack.geoclaw import topotools
from clawpack.visclaw import animation_tools
from clawpack.geoclaw import fgout_tools
from clawpack.visclaw import colormaps

tan = [0.8,0.5,0.2]
light_blue = [0.0,1.0,1.0]
cmap_massfrac = colormaps.make_colormap({0:light_blue, 0.1:light_blue,
                                         1:tan})

def plot_topo():
    topo = topotools.Topography('basal_topo.tt3')
    z = array([0.])
    x = topo.x
    y = -topo.y
    X,Y,Z = meshgrid(x, y, z, indexing='ij')
    topoxyz = pv.StructuredGrid(X,Y,Z)

    B = flipud(topo.Z)
    B = fliplr(B)
    topoxyz.point_data['B'] = B.flatten(order='C')

    warpfactor = 5  # amplification of elevations

    topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

    # plot topo alone:
    p = pv.Plotter()
    #p.add_mesh(topoxyz,cmap='gist_earth',clim=(-10,100))
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    p.add_title('Topography')
    p.show(window_size=(1500,1500))

def plot_topo_fgout():
    fgno = 1  # which fgout grid

    outdir = '_output1'
    format = 'binary32'  # format of fgout grid output
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format, qmap='dclaw')

    fgout = fgout_grid.read_frame(1)

    x = fgout.x
    y = fgout.y
    z = array([0.])
    X,Y,Z = meshgrid(x, y, z, indexing='ij')
    topoxyz = pv.StructuredGrid(X,Y,Z)

    B = fgout.B
    topoxyz.point_data['B'] = B.flatten(order='C')
    warpfactor = 5  # amplification of elevations
    topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

    h = fgout.h
    topoxyz.point_data['h'] = h.flatten(order='C')
    warpfactor = 5  # amplification of elevations
    hwarp = topoxyz.warp_by_scalar('h', factor=warpfactor)

    eta = fgout.eta
    topoxyz.point_data['eta'] = B.flatten(order='C')
    warpfactor = 5  # amplification of elevations
    etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)

    p = pv.Plotter()
    #p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    p.add_mesh(topowarp,color='g')
    #p.add_mesh(hwarp,color='r')
    #p.add_mesh(etawarp,color='b')
    p.add_title('Topography fgout.B')
    p.show(window_size=(1500,1500))


def plot_fgout(save_fig=False):

    global etamesh

    fgno = 1  # which fgout grid
    outdir = '_output'
    format = 'binary32'  # format of fgout grid output
    warpfactor = 5  # amplification of elevations

    fgout_frames = glob.glob(os.path.join(outdir, \
                              'fgout%s.t*' % str(fgno).zfill(4)))
    nfgout = len(fgout_frames)
    fgframes = range(1, nfgout+1)
    print('Found %i fgout frames' % nfgout)

    # Instantiate object for reading fgout frames:
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format, qmap='dclaw')

    fgout = fgout_grid.read_frame(1)

    x = fgout.x
    y = fgout.y
    z = array([0.])
    X,Y,Z = meshgrid(x, y, z, indexing='ij')
    topoxyz = pv.StructuredGrid(X,Y,Z)

    B = fgout.B
    topoxyz.point_data['B'] = B.flatten(order='C')
    topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

    p = pv.Plotter()
    #p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    p.add_mesh(topowarp,color='lightgreen')


    def set_frame(tsec):
        global etamesh
        fgframeno = int(round(tsec)) + 1
        fgout = fgout_grid.read_frame(fgframeno)
        tsec = fgout.t
        print('Frame %i, t = %.1f seconds' % (fgframeno, fgout.t))

        # replace land surface by nan in eta so it only shows water:
        # (big tolerance for landslide, would normally be smaller)
        eta = where(fgout.h>1, fgout.eta, nan)

        topoxyz.point_data['eta'] = eta.flatten(order='C')
        etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
        p.remove_actor(etamesh)
        #etamesh = p.add_mesh(etawarp,color='c')  # dirt and water both cyan

        # color the mesh based on mass fraction, using cmap_massfrac
        h = fgout.h
        massfrac = divide(fgout.hm, h, where=h>0, out=nan*ones(h.shape))
        etamesh = p.add_mesh(etawarp,scalars=massfrac,
                             colormap=cmap_massfrac, clim=(0,0.6))

        p.add_title('Frame %i at time %.1f seconds' % (fgframeno,tsec))
        if save_fig:
            # to capture each frame after moving slider:
            fname = 'PyVistaFrame%s.png' % str(fgframeno).zfill(4)
            p.screenshot(fname)
            print('Saving ', fname)

    p.add_slider_widget(set_frame, [0,100], value=0,
                        pointa=(0.4,0.8), pointb=(0.9,0.8),
                        title='Time (sec)')

    p.camera_position = 'xz'
    p.camera.azimuth = 45
    p.camera.elevation = 20
    p.show(window_size=(1500,1500))

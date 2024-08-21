from pylab import *
import glob, os
import pyvista as pv
from clawpack.geoclaw import topotools
from clawpack.visclaw import animation_tools
from clawpack.geoclaw import fgout_tools
from clawpack.visclaw import colormaps

warpfactor = 5  # vertical amplification of elevations in 3D plots

outdir  = '_output'
fgno = 1 # fgout grid number

if 1:
    # determine how many fgout frames are in outdir:
    fgout_frames = glob.glob(os.path.join(outdir, \
                              'fgout%s.t*' % str(fgno).zfill(4)))
    nfgout = len(fgout_frames)
    fgframenos = array(range(1, nfgout+1))
    print('Found %i fgout frames' % nfgout)
else:
    # specify which frames to use:
    fgframenos = array(range(1,50,2))


# where to save a png file for each frame, for constructing an animation:
framedir = '_frames'
os.system('mkdir -p %s' % framedir)
os.system('rm %s/*' % framedir)  # remove frames from previous version

window_size = (1200,1200)
        
tan = [0.8,0.5,0.2]
light_blue = [0.0,1.0,1.0]
cmap_massfrac = colormaps.make_colormap({0:light_blue, 0.1:light_blue,
                                         1:tan})

def plot_topo():
    """
    plot topo from topofile
    """
    topo = topotools.Topography('basal_topo.tt3')
    z = array([0.])
    x = topo.x
    y = -topo.y
    X,Y,Z = meshgrid(x, y, z, indexing='ij')
    topoxyz = pv.StructuredGrid(X,Y,Z)

    B = flipud(topo.Z)
    B = fliplr(B)
    topoxyz.point_data['B'] = B.flatten(order='C')

    topowarp = topoxyz.warp_by_scalar('B', factor=warpfactor)

    p = pv.Plotter()
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    p.add_title('Topography')
    p.show(window_size=window_size)

def plot_topo_fgout():
    """
    plot topo from an fgout file
    """
    fgno = 1  # which fgout grid

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
    p.window_size = window_size
    p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    #p.add_mesh(topowarp,color='g')
    p.add_title('Topography fgout.B')
    p.show()


def plot_fgout(make_animation=False):
    """
    If make_animation == False, plot one fgout frame with a slider bar to
    change frame number. Also prints the camera position after changing frame,
    which you can copy and paste into this function in order to set the
    camera position for making an animation.
    
    If make_animation == True, make an mp4 animation of all frames
    specified by fgframenos set above.
    """
    
    global etamesh

    fgno = 1  # which fgout grid
    warpfactor = 5  # amplification of elevations

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

    p = pv.Plotter(off_screen=make_animation, lighting='three lights')
    p.window_size = window_size
    #p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))
    p.add_mesh(topowarp,color='lightgreen')
    etamesh = None


    def set_frameno(fgframeno):
        global etamesh
        fgframeno = int(round(fgframeno))
        fgout = fgout_grid.read_frame(fgframeno)
        tsec = fgout.t
        print('Frame %i, t = %.1f seconds' % (fgframeno, fgout.t))

        # replace land surface by nan in eta so it only shows water:
        # (big tolerance for landslide, would normally be smaller)
        eta = where(fgout.h>1, fgout.eta, nan)

        topoxyz.point_data['eta'] = eta.flatten(order='C')
        etawarp = topoxyz.warp_by_scalar('eta', factor=warpfactor)
        if etamesh:
            p.remove_actor(etamesh)

        # color the mesh based on mass fraction, using cmap_massfrac
        h = fgout.h
        massfrac = divide(fgout.hm, h, where=h>0, out=nan*ones(h.shape))
        etamesh = p.add_mesh(etawarp,scalars=massfrac,
                             colormap=cmap_massfrac, clim=(0,0.6))

        p.add_title('Frame %i at time %.1f seconds' % (fgframeno,tsec))
        
        if not make_animation:
            # print camera position so that this can be copied and pasted
            # into this script after adjusting (and then sliding frameno)
            print('p.camera_position = ', p.camera_position)


    # initial camera position:
    p.camera_position =  [
         (9034.796206317897, -2484.6171987620633, 3703.6128928929047),
         (1500.0, 1500.0, 950.89111328125),
         (-0.29984440650740857, 0.08916816605502331, 0.949811755059182)]
 
 
    if not make_animation:
        fgfr1 = fgframenos[0]
        fgfr2 = fgframenos[-1]
        p.add_slider_widget(set_frameno, [fgfr1,fgfr2],
                        value=fgfr1, title='Frame',
                        pointa=(0.4,0.85), pointb=(0.9,0.85), color='blue',
                        slider_width=0.02, tube_width=0.005)

        p.show()
        
    else:

        # make a png file for each frame:
        for fgframeno in fgframenos:
            set_frameno(fgframeno)
            fname_png = '%s/PyVistaFrame%s.png' \
                        % (framedir, str(fgframeno).zfill(4))
            p.screenshot(fname_png)
            print('Created ',fname_png)
            
        p.close()

        # combine png files into mp4 animation:
        anim = animation_tools.make_anim(framedir,
               fname_pattern='PyVistaFrame*.png')
        fname_mp4 = 'radial_slide_animation.mp4'
        animation_tools.make_mp4(anim, fname_mp4)
        print('Created ',fname_mp4)

if __name__=='__main__':
    
    plot_fgout(make_animation=False)

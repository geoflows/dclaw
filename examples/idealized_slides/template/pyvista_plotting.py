import glob
import os

import cmocean
import pyvista as pv
import yaml
from clawpack.geoclaw import fgout_tools, topotools
from clawpack.visclaw import animation_tools
from matplotlib.colors import LinearSegmentedColormap
from pylab import *

with open("params.yml", "r") as file:
    params = yaml.safe_load(file)

hmax = params["depth"]

make_animation = True
fps = 10

camera_position = None
# xyz of camera position
# xyz of where you are looking
# unit vector of up

# position 1 looking up fiord at the landslide
fname_mp4 = "idealized_runout.mp4"
camera_position = [
    (948.0862737449888, -381.7006596750718, 463.6465382450876),
    (223.62423088340623, 38.927734486851804, -21.986370544155236),
    (-0.4229897516013297, 0.27008747172574316, 0.8649464883199661),
]


warpfactor = 1  # vertical amplification of elevations in 3D plots

outdir = "_output"
fgno = 1  # fgout grid number

if 1:
    # determine how many fgout frames are in outdir:
    fgout_frames = glob.glob(os.path.join(outdir, "fgout%s.t*" % str(fgno).zfill(4)))
    nfgout = len(fgout_frames)
    fgframenos = array(range(1, nfgout + 1))
    print("Found %i fgout frames" % nfgout)
else:
    # specify which frames to use:
    fgframenos = array(range(1, 600, 1))
    fgframenos = array(range(1, 10, 1))


# where to save a png file for each frame, for constructing an animation:
framedir = "_frames"
os.system("mkdir -p %s" % framedir)
os.system("rm %s/*" % framedir)  # remove frames from previous version

window_size = (1200, 1200)
tan = [140 / 255, 81 / 255, 10 / 255]
light_blue = [171 / 255, 217 / 255, 233 / 255]
nodes = [0.0, 0.1, 0.6, 1.0]
colors = [light_blue, light_blue, tan, tan]
cmap_massfrac = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))
cmap_depth = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))


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
    global bmesh
    fgno = 1  # which fgout grid

    # Instantiate object for reading fgout frames:
    fgout_grid = fgout_tools.FGoutGrid(fgno, outdir, format, qmap="dclaw")
    fgout_grid.read_fgout_grids_data()
    fgout = fgout_grid.read_frame(1)

    x = fgout.x
    y = fgout.y
    z = array([0.0])
    X, Y, Z = meshgrid(x, y, z, indexing="ij")
    topoxyz = pv.StructuredGrid(X, Y, Z)

    p = pv.Plotter(off_screen=make_animation, lighting="three lights")
    p.window_size = window_size
    # p.add_mesh(topowarp,cmap='gist_earth',clim=(-20,400))

    B = fgout.B
    topoxyz.point_data["B"] = B.flatten(order="F")
    topowarp = topoxyz.warp_by_scalar("B", factor=warpfactor)
    bmesh = p.add_mesh(topowarp, color="#dfc27d")

    etamesh = None

    def set_frameno(fgframeno):
        global etamesh
        global bmesh

        fgframeno = int(round(fgframeno))
        fgout = fgout_grid.read_frame(fgframeno)
        tsec = fgout.t
        print("Frame %i, t = %.1f seconds" % (fgframeno, fgout.t))

        if bmesh:
            p.remove_actor(bmesh)

        B = fgout.B
        topoxyz.point_data["B"] = B.flatten(order="F")
        topowarp = topoxyz.warp_by_scalar("B", factor=warpfactor)
        bmesh = p.add_mesh(topowarp, color="#dfc27d")

        # replace land surface by nan in eta so it only shows water:
        # (big tolerance for landslide, would normally be smaller)
        eta = where(fgout.h > 0.0001, fgout.eta, nan)

        topoxyz.point_data["eta"] = eta.flatten(order="F")
        etawarp = topoxyz.warp_by_scalar("eta", factor=warpfactor)
        if etamesh:
            p.remove_actor(etamesh)

        # color the mesh based on mass fraction, using cmap_massfrac
        h = fgout.h.flatten(order="F")
        hm = fgout.hm.flatten(order="F")
        massfrac = divide(hm, h, where=h > 0, out=nan * ones(h.shape))
        etamesh = p.add_mesh(
            etawarp,
            scalars=h,
            colormap=cmocean.cm.turbid,
            clim=(0, 10),
            scalar_bar_args={"title": "Depth, m"},
        )

        p.add_title("Time: %.1f seconds" % (tsec))

        if not make_animation:
            # print camera position so that this can be copied and pasted
            # into this script after adjusting (and then sliding frameno)
            print("p.camera_position = ", p.camera_position)

    # initial camera position:
    p.camera_position = camera_position

    if not make_animation:
        fgfr1 = fgframenos[0]
        fgfr2 = fgframenos[-1]
        p.add_slider_widget(
            set_frameno,
            [fgfr1, fgfr2],
            value=fgfr1,
            title="Frame",
            pointa=(0.4, 0.85),
            pointb=(0.9, 0.85),
            color="blue",
            slider_width=0.02,
            tube_width=0.005,
        )

        p.show()

    else:

        # make a png file for each frame:
        for fgframeno in fgframenos:
            set_frameno(fgframeno)
            fname_png = "%s/PyVistaFrame%s.png" % (framedir, str(fgframeno).zfill(4))
            p.screenshot(fname_png)
            print("Created ", fname_png)

        p.close()
        # combine png files into mp4 animation:
        anim = animation_tools.make_anim(framedir, fname_pattern="PyVistaFrame*.png")
        animation_tools.make_mp4(anim, fname_mp4, fps=fps)
        print("Created ", fname_mp4)


if __name__ == "__main__":

    plot_fgout(make_animation=make_animation)

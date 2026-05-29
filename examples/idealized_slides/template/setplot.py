"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import os
import sys

import clawpack.dclaw.plot as dplot
import cmocean
import matplotlib as mpl
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
import yaml
from clawpack.visclaw import colormaps, geoplot, gridtools

with open("params.yml", "r") as file:
    params = yaml.safe_load(file)


xlower = params["xlower"]
ylower = params["ylower"]
xupper = params["xupper"]
yupper = params["yupper"]

xlimits = [xlower, xupper]
ylimits = [ylower, yupper]

outdir2 = None

bouss = False
# setplot uses python indexing
if bouss:
    i_eta = 10 - 1
    i_hm = 6 - 1
    i_pb = 7 - 1
else:
    i_eta = 8 - 1
    i_hm = 4 - 1
    i_pb = 5 - 1


# --------------------------
def setplot(plotdata=None):
    # --------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    import clawpack.dclaw.plot as dplot
    from numpy import linspace, mod
    from pylab import gca, ticklabel_format, title, xticks

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData

        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = "binary"

    def timeformat(t):
        hours = int(t / 3600.0)
        tmin = mod(t, 3600.0)
        min = int(tmin / 60.0)
        sec = int(mod(tmin, 60.0))
        timestr = "%s:%s:%s" % (hours, str(min).zfill(2), str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        t = current_data.t
        timestr = timeformat(t)
        title("t = %s" % timestr)

    def aa(current_data):
        gca().set_aspect(1.0)
        title_hours(current_data)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)

    def aa_notime(current_data):
        gca().set_aspect(1.0)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)
        title("")

    # -----------------------------------------
    # Figure for state variables
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="Computational domain", figno=0)
    plotfigure.kwargs = {"figsize": (8, 7), "dpi": 600}
    plotfigure.show = True

    # Panel 1: Hillshade and Depth
    plotaxes = plotfigure.new_plotaxes("depth")
    plotaxes.title = "Depth"
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(221)"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = aa
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.depth
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "Depth (m)"
    plotitem.imshow_cmap = "Purples"
    plotitem.imshow_norm = mpl.colors.LogNorm(vmin=0.001, vmax=4, clip=True)

    plotitem.patchedges_show = True

    # Panel 2 : Hillshade and Velocity
    plotaxes = plotfigure.new_plotaxes("velocity")
    plotaxes.title = ""
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(222)"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = aa_notime
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.velocity_magnitude
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "Velocity (m/s)"
    plotitem.imshow_cmap = "Greens"
    plotitem.imshow_norm = mpl.colors.Normalize(vmin=0, vmax=5)

    # Panel 2 : Hillshade and M
    plotaxes = plotfigure.new_plotaxes("solidfrac")
    plotaxes.title = ""
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(223)"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = aa_notime
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.solid_frac
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "Solid fraction (-)"
    plotitem.imshow_cmap = "pink_r"
    plotitem.imshow_norm = mpl.colors.Normalize(vmin=0, vmax=1)

    # Panel 4: Hillshade and Pressure
    plotaxes = plotfigure.new_plotaxes("Pressure")
    plotaxes.title = "Pressure Ratio"
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(224)"
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ylimits
    plotaxes.afteraxes = aa
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.basal_pressure_over_hydrostatic
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "P_b/Hydrostatic (-)"
    plotitem.imshow_cmap = "BrBG"
    plotitem.imshow_norm = mpl.colors.TwoSlopeNorm(vcenter=1, vmin=0, vmax=2.5)
    plotitem.patchedges_show = True

    # -----------------------
    # Figure for rheology
    # -----------------------

    plotfigure = plotdata.new_plotfigure(name="apparent_fc", figno=1)
    plotfigure.show = True

    def aa_apparent_fc(current_data):
        gca().set_ylim([0, 1])
        gca().set_xscale('log')
        gca().set_xlim([1e-7, 1])

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.title = "Apparent FC"
    plotaxes.grid = True
    plotaxes.afteraxes = aa_apparent_fc
    plotaxes.xlabel = "Modified Inertial Number [-]"
    plotaxes.ylabel = "Apparent friction coefficient [-]"

    def shear_rate_fc(current_data):

        phi = 35
        mu = np.tan(np.deg2rad(phi))

        shear = dplot.shear(current_data)
        depth = dplot.depth(current_data)

        N = dplot.N(current_data)
        N[N==0] = 1e-10
        sigma_e_over_lithostatic = dplot.sigma_e_over_lithostatic(current_data)

        mu_eff = mu*sigma_e_over_lithostatic

        return N, mu_eff, depth
        # ADD color mapping for plots

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")
    plotitem.map_2d_to_1d = shear_rate_fc

    plotitem.plotstyle = "o"
    plotitem.map_color = True
    plotitem.plot_cmap = "viridis"
    plotitem.plot_norm = mpl.colors.LogNorm(vmin=0.1, vmax=2)
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "right",
        "orientation": "vertical",
    }
    plotitem.colorbar_label = "Depth (m)"

    plotitem.amr_data_show = [True, False, False]
    plotitem.kwargs = {"markersize": 1}
    plotitem.show = True  # show on plot?



    # -------------------------------------
    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        import os

        from clawpack.visclaw import plot_timing_stats

        try:
            timing_plotdir = plotdata.plotdir + "/_timing_figures"
            os.system("mkdir -p %s" % timing_plotdir)
            units = {"comptime": "minutes", "simtime": "minutes", "cell": "millions"}
            plot_timing_stats.make_plots(
                outdir=plotdata.outdir,
                make_pngs=True,
                plotdir=timing_plotdir,
                units=units,
            )
            os.system("cp %s/timing.* %s" % (plotdata.outdir, timing_plotdir))
        except:
            print("*** Error making timing plots")

    otherfigure = plotdata.new_otherfigure(
        name="timing", fname="_timing_figures/timing.html"
    )
    otherfigure.makefig = make_timing_plots

    # -----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True  # print figures
    plotdata.print_format = "png"  # file format
    plotdata.print_framenos = "all"  # list of frames to print
    plotdata.print_gaugenos = "all"  # list of gauges to print
    plotdata.print_fignos = "all"  # list of figures to print
    plotdata.html = True  # create html files of plots?
    plotdata.html_homelink = "../README.html"  # pointer for top of index
    plotdata.latex = True  # create latex file of plots?
    plotdata.latex_figsperline = 2  # layout of plots
    plotdata.latex_framesperline = 1  # layout of plots
    plotdata.latex_makepdf = False  # also run pdflatex?
    plotdata.parallel = True  # make multiple frame png's at once
    plotdata.mp4_movie = True

    return plotdata

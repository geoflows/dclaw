"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from clawpack.visclaw import colormaps, geoplot, gridtools
from setinput import xl1, xl2, yl1, yl2

xlimits = [2000, 4000]
ylimits = [500, 1500]
xlimits = [xl1 - 200, xl2 + 200]
ylimits = [yl1 - 200, yl2 + 200]
import os
import sys


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
    plotfigure.kwargs = {"figsize": (8, 7)}
    plotfigure.show = True

    for i in [0, 1, 2]:
        amr_data = [0, 0, 0]
        amr_data[i] = 1

        plotind = i * 3

        # Panel 1: Hillshade and Depth
        plotaxes = plotfigure.new_plotaxes(f"depth_{i}")
        plotaxes.title = "Depth"
        plotaxes.scaled = True
        plotaxes.axescmd = f"subplot(3,3,{plotind+1})"
        plotaxes.xlimits = xlimits
        plotaxes.ylimits = ylimits
        plotaxes.afteraxes = aa
        plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
        plotitem.show = True
        plotitem.plot_var = dplot.eta
        plotitem.add_colorbar = False
        plotitem.amr_data_show = amr_data

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
        plotitem.amr_data_show = amr_data

        # Panel 2 : Hillshade and Velocity
        plotaxes = plotfigure.new_plotaxes(f"velocity_{i}")
        plotaxes.title = ""
        plotaxes.scaled = True
        plotaxes.axescmd = f"subplot(3,3,{plotind+2})"
        plotaxes.xlimits = xlimits
        plotaxes.ylimits = ylimits
        plotaxes.afteraxes = aa_notime
        plotitem.amr_data_show = amr_data

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
        plotitem.imshow_norm = mpl.colors.LogNorm(vmin=0.1, vmax=10, clip=True)
        plotitem.amr_data_show = amr_data

        # Panel 2 : Hillshade and M
        plotaxes = plotfigure.new_plotaxes(f"solidfrac_{i}")
        plotaxes.title = ""
        plotaxes.scaled = True
        plotaxes.axescmd = f"subplot(3,3,{plotind+3})"
        plotaxes.xlimits = xlimits
        plotaxes.ylimits = ylimits
        plotaxes.afteraxes = aa_notime
        plotitem.amr_data_show = amr_data

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
        plotitem.amr_data_show = amr_data

    # eventually add segregation and entrained depth

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

    return plotdata

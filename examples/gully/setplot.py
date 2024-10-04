
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import numpy as np
import matplotlib.pyplot as plt
import cmocean
import matplotlib as mpl

from clawpack.visclaw import geoplot
from clawpack.visclaw import colormaps
from clawpack.visclaw import gridtools


import os,sys


outdir2 = None
#outdir2 = '_output_SL50_order2_trans0'

#--------------------------
def setplot(plotdata=None):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    import clawpack.dclaw.plot as dplot
    from numpy import linspace

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = 'binary'



    def timeformat(t):
        from numpy import mod
        hours = int(t/3600.)
        tmin = mod(t,3600.)
        min = int(tmin/60.)
        sec = int(mod(tmin,60.))
        timestr = '%s:%s:%s' % (hours,str(min).zfill(2),str(sec).zfill(2))
        return timestr

    def title_hours(current_data):
        from pylab import title
        t = current_data.t
        timestr = timeformat(t)
        title('t = %s' % timestr)

    def aa(current_data):
        from pylab import ticklabel_format, xticks, gca, cos, pi, savefig
        gca().set_aspect(1.)
        title_hours(current_data)
        ticklabel_format(useOffset=False)
        xticks(rotation=20)


    #-----------------------------------------
    # Figure for surface and velocity
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Computational domain', figno=0)
    plotfigure.kwargs = {'figsize':(8,7)}
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('depth')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(211)"
    plotaxes.xlimits = [2000, 4000]
    plotaxes.ylimits = [500, 1500]
    plotaxes.afteraxes = aa

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False

    # Surface
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.depth
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "Surface (m)"
    plotitem.imshow_cmap = cmocean.cm.matter
    plotitem.imshow_cmin = 0
    plotitem.imshow_cmax = 2

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('velocity')
    plotaxes.title = 'Velocity'
    plotaxes.scaled = True
    plotaxes.axescmd = "subplot(212)"
    plotaxes.xlimits = [2000, 4000]
    plotaxes.ylimits = [500, 1500]
    plotaxes.afteraxes = aa

    # Hillshade
    plotitem = plotaxes.new_plotitem(plot_type="2d_hillshade")
    plotitem.show = True
    plotitem.plot_var = dplot.eta
    plotitem.add_colorbar = False

    # Surface
    plotitem = plotaxes.new_plotitem(plot_type="2d_imshow")
    plotitem.plot_var = dplot.velocity_magnitude
    plotitem.add_colorbar = True
    plotitem.colorbar_kwargs = {
        "shrink": 0.5,
        "location": "bottom",
        "orientation": "horizontal",
    }
    plotitem.colorbar_label = "Velocity (m/s)"
    plotitem.imshow_cmap = cmocean.cm.speed
    plotitem.imshow_cmin = 0
    plotitem.imshow_cmax = 10
    #-------------------------------------
    # Plots of timing (CPU and wall time):

    def make_timing_plots(plotdata):
        import os
        from clawpack.visclaw import plot_timing_stats
        try:
            timing_plotdir = plotdata.plotdir + '/_timing_figures'
            os.system('mkdir -p %s' % timing_plotdir)
            units = {'comptime':'minutes', 'simtime':'minutes', 'cell':'millions'}
            plot_timing_stats.make_plots(outdir=plotdata.outdir, make_pngs=True,
                                          plotdir=timing_plotdir, units=units)
            os.system('cp %s/timing.* %s' % (plotdata.outdir, timing_plotdir))
        except:
            print('*** Error making timing plots')

    otherfigure = plotdata.new_otherfigure(name='timing',
                    fname='_timing_figures/timing.html')
    otherfigure.makefig = make_timing_plots


    #-----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'        # list of frames to print
    plotdata.print_gaugenos = 'all'          # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata

"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

from __future__ import absolute_import

import numpy as np

a = 1.0
sigma = 0.5
h0 = 0.1
grav = 9.81
omega = np.sqrt(2.0 * grav * h0) / a


# --------------------------
def setplot(plotdata=None):
    # --------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    from clawpack.visclaw import geoplot

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData

        plotdata = ClawPlotData()

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user["drytol"] = 1.0e-3

    plotdata.beforeframe = set_drytol

    def surface(current_data):
        """
        Return a masked array containing the surface elevation only in wet cells.
        Surface is eta = h+topo, assumed to be output as 4th column of fort.q
        files.
        """
        drytol = current_data.user.get("drytol", 1e-3)
        q = current_data.q
        h = q[0, ...]
        eta = q[3, ...]  # DIG --  NO! should still be -1 component

        water = np.ma.masked_where(h <= drytol, eta)

        try:
            # Use mask covering coarse regions if it's set:
            m = current_data.mask_coarse
            water = np.ma.masked_where(m, water)
        except:
            pass

        return water

    # -----------------------------------------
    # Figure for pcolor plot
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="pcolor", figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes("pcolor")
    plotaxes.title = "Surface"
    plotaxes.scaled = True

    # Water
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.surface
    # plotitem.plot_var = surface
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.1
    plotitem.pcolor_cmax = 0.1
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type="2d_pcolor")
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0, 0, 0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [-2, 2]
    plotaxes.ylimits = [-2, 2]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type="2d_contour")
    plotitem.plot_var = geoplot.topo

    plotitem.contour_levels = np.linspace(-0.1, 0.5, 20)
    plotitem.amr_contour_colors = ["k"]  # color on each level
    plotitem.kwargs = {"linestyles": "solid"}
    plotitem.amr_contour_show = [1]
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # -----------------------------------------
    # Figure for cross section
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="cross-section", figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2, 2]
    plotaxes.ylimits = [-0.15, 0.3]
    plotaxes.title = "Cross section at y=0"

    def plot_topo_xsec(current_data):
        from pylab import cos, legend, nan, plot, sin, where

        t = current_data.t

        x = np.linspace(-2, 2, 201)
        y = 0.0
        B = h0 * (x**2 + y**2) / a**2 - h0
        eta1 = (
            sigma
            * h0
            / a**2
            * (2.0 * x * cos(omega * t) + 2.0 * y * sin(omega * t) - sigma)
        )
        etatrue = where(eta1 > B, eta1, nan)
        plot(x, etatrue, "r", label="true solution", linewidth=2)
        plot(x, B, "g", label="bathymetry")
        plot([0], [-1], "kx", label="Computed")  ## For legend
        legend()

    plotaxes.afteraxes = plot_topo_xsec

    plotitem = plotaxes.new_plotitem(plot_type="1d_from_2d_data")

    def xsec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import ravel, where

        x = current_data.x
        y = ravel(current_data.y)
        dy = current_data.dy
        q = current_data.q

        ij = where((y <= dy / 2.0) & (y > -dy / 2.0))
        x_slice = ravel(x)[ij]
        eta_slice = ravel(q[-1, :, :])[ij]
        return x_slice, eta_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = "kx"  ## need to be able to set amr_plotstyle
    plotitem.kwargs = {"markersize": 3}
    plotitem.amr_show = [1]  # plot on all levels

    # -----------------------------------------
    # Figure for grids alone
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name="grids", figno=2)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-2, 2]
    plotaxes.ylimits = [-2, 2]
    plotaxes.title = "grids"
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type="2d_patch")
    plotitem.amr_patch_bgcolor = ["#ffeeee", "#eeeeff", "#eeffee"]
    plotitem.amr_celledges_show = [1, 1, 0]
    plotitem.amr_patchedges_show = [1]

    # -----------------------------------------

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True  # print figures
    plotdata.print_format = "png"  # file format
    plotdata.print_framenos = "all"  # list of frames to print
    plotdata.print_gaugenos = []  # list of gauges to print
    plotdata.print_fignos = "all"  # list of figures to print
    plotdata.html = True  # create html files of plots?
    plotdata.html_homelink = "../README.html"  # pointer for top of index
    plotdata.latex = True  # create latex file of plots?
    plotdata.latex_figsperline = 2  # layout of plots
    plotdata.latex_framesperline = 1  # layout of plots
    plotdata.latex_makepdf = False  # also run pdflatex?
    plotdata.parallel = True  # make multiple frame png's at once

    return plotdata

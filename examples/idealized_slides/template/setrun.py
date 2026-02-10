"""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.

"""

import os
import sys

import numpy as np

try:
    CLAW = os.environ["CLAW"]
except:
    raise Exception("*** Must first set CLAW environment variable")

import yaml
from clawpack.amrclaw.data import FlagRegion
from clawpack.geoclaw import fgmax_tools, fgout_tools

with open("params.yml", "r") as file:
    params = yaml.safe_load(file)

# problem geometry from yaml
m0 = params["m0"]
phi = params["phi"]
kr = float(params["kr"])
mc = params["mc"]
alphamethod = params["alphamethod"]
src2method = params["src2method"]
order = params["order"]

amr = params["amr"]

xlower = params["xlower"]
ylower = params["ylower"]
xupper = params["xupper"]
yupper = params["yupper"]
num_cells_x = params["num_cells_x"]
num_cells_y = params["num_cells_y"]


if order == 2:
    transverse = 2
    cfl_desired = 0.85
    cfl_max = 1.0
else:
    transverse = 0
    cfl_desired = 0.45
    cfl_max = 0.5


if amr:
    mxnest = 3
    ratios = [5, 2]
    factor = 1
else:
    mxnest = 1
    ratios = [1]
    factor = 10

bouss = False
if bouss:
    i_eta = 10
    i_hm = 6
    num_eqn = 9
else:
    i_eta = 8
    i_hm = 4
    num_eqn = 7


# ------------------------------
def setrun(claw_pkg="dclaw"):
    # ------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == "dclaw", "Expected claw_pkg = 'dclaw'"

    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # ------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    # ------------------------------------------------------------------

    # probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    # probdata.add_param('variable_eta_init', True)  # now in qinit info

    # ------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    # ------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # Set single grid parameters first.
    # See below for AMR parameters.

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    # Lower and upper edge of computational domain:
    clawdata.lower[0] = xlower
    clawdata.upper[0] = xupper

    clawdata.lower[1] = ylower
    clawdata.upper[1] = yupper

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = num_cells_x * factor
    clawdata.num_cells[1] = num_cells_y * factor

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = num_eqn

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 10

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 0

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0

    # Restart from checkpoint file of a previous run?
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False  # True to restart from prior results
    clawdata.restart_file = ""

    # -------------
    # Output times:
    # --------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = 1

    if clawdata.output_style == 1:
        # Output nout frames at equally spaced times up to tfinal:
        dt = 10
        clawdata.tfinal = 60.0  # 240.
        clawdata.num_output_times = int(clawdata.tfinal / dt)

        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = 100
        clawdata.total_steps = 2000
        clawdata.output_t0 = True

    clawdata.output_format = "ascii"

    clawdata.output_q_components = "all"  # need all
    clawdata.output_aux_components = "none"  # eta=h+B is in q
    clawdata.output_aux_onlyonce = True  # output aux arrays each frame

    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = True

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = 0.0001

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = cfl_desired
    clawdata.cfl_max = cfl_max

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = order

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = "unsplit"

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = transverse

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 5

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = [4, 4, 4, 4, 4]  # TODO VERIFY THAT 4 in old and new are the same

    clawdata.use_fwaves = True  # True ==> use f-wave version of algorithms
    # TODO This is not in old setrun.py

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = "godunov"

    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = "extrap"
    clawdata.bc_upper[0] = "extrap"

    clawdata.bc_lower[1] = "extrap"
    clawdata.bc_upper[1] = "extrap"

    # --------------
    # Checkpointing:
    # --------------

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    # negative checkpoint_style means alternate between aaaaa and bbbbb files
    # so that at most 2 checkpoint files exist at any time, useful when
    # doing frequent checkpoints of large problems.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif clawdata.checkpt_style == 1:
        # Checkpoint only at tfinal.
        pass

    elif abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.
        clawdata.checkpt_times = 3600.0 * np.arange(1, 16, 1)

    elif abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # ---------------
    # AMR parameters:
    # ---------------
    amrdata = rundata.amrdata

    # max number of refinement levels:
    amrdata.amr_levels_max = mxnest

    # List of refinement ratios at each level (length at least mxnest-1)
    # dx = dy = 2', 10", 2", 1/3":
    amrdata.refinement_ratios_x = ratios
    amrdata.refinement_ratios_y = ratios
    amrdata.refinement_ratios_t = ratios

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = [
        "center",
        "center",
        "yleft",
        "center",
        "center",
        "center",
        "center",
        "center",
        "center",
        "center",
    ]

    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False  # use Richardson?
    amrdata.flag2refine = True

    # steps to take on each level L between regriddings of level L+1:
    amrdata.regrid_interval = 3

    # width of buffer zone around flagged points:
    # (typically the same as regrid_interval so waves don't escape):
    amrdata.regrid_buffer_width = 2

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.700000

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = 1

    # ---------------
    # Regions:
    # ---------------
    # rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]
    # NO OLD STYLE REGIONS USED HERE

    # ---------------
    # NEW flagregions
    # ---------------

    flagregions = rundata.flagregiondata.flagregions  # initialized to []

    # now append as many flagregions as desired to this list:

    # ---------------
    # Gauges:
    # ---------------
    # for gauges append lines of the form  [gaugeno, x, y, t1, t2]
    rundata.gaugedata.gauges = []

    # Set GeoClaw specific runtime parameters.

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 1
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = False

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = -9999
    geo_data.dry_tolerance = 1.0e-3
    geo_data.friction_forcing = True  # TODO change?
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1e6

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = True
    refinement_data.wave_tolerance = 0.01

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]
    topofiles.append([3, "basal_topo.tt3"])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data

    # == setqinit.data values ==
    qinitdclaw_data = rundata.qinitdclaw_data  # initialized when rundata instantiated

    etafile = "surface_topo.tt3"
    qinitdclaw_data.qinitfiles.append([3, i_eta, etafile])

    # == setauxinit.data values ==
    # auxinitdclaw_data = rundata.auxinitdclaw_data  # initialized when rundata instantiated

    # == fgmax_grids.data values ==
    # NEW STYLE STARTING IN v5.7.0

    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 5  # Save all quantities
    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # fg = fgmax_tools.FGmaxGrid()

    # fg.point_style = 2  # uniform rectangular x-y grid

    # dx_fg = 20

    # fg.x1 = clawdata.lower[0] + dx_fg / 2
    # fg.x2 = clawdata.upper[0] - dx_fg / 2
    # fg.y1 = clawdata.lower[1] + dx_fg / 2
    # fg.y2 = clawdata.upper[1] - dx_fg / 2
    # fg.dx = dx_fg  # desired resolution of fgmax grid
    # fg.tstart_max = 0.0  # when to start monitoring max values
    # fg.tend_max = 1.0e10  # when to stop monitoring max values
    # fg.dt_check = 1.0  # target time (sec) increment between updating
    # # max values
    # fg.min_level_check = 2  # which levels to monitor max on
    # fg.arrival_tol = 1.0e-2  # tolerance for flagging arrival
    # fg.interp_method = 0  # 0 ==> pw const in cells, recommended
    # fgmax_grids.append(fg)  # written to fgmax_grids.data

    # == setdclaw.data values ==
    dclaw_data = rundata.dclaw_data  # initialized when rundata instantiated

    dclaw_data.rho_f = 1000.0
    dclaw_data.rho_s = 2700.0
    dclaw_data.m_crit = mc
    dclaw_data.m0 = m0
    dclaw_data.mref = 0.6
    dclaw_data.kref = kr
    dclaw_data.phi = phi
    dclaw_data.delta = 0.001
    dclaw_data.mu = 0.005
    dclaw_data.alpha_c = 0.01
    dclaw_data.c1 = 1
    dclaw_data.sigma_0 = 1.0e3

    dclaw_data.src2method = src2method
    dclaw_data.alphamethod = alphamethod

    dclaw_data.segregation = 0
    dclaw_data.beta_seg = 0.0
    dclaw_data.chi0 = 0.5
    dclaw_data.chie = 0.5

    dclaw_data.bed_normal = 0
    dclaw_data.theta_input = 0.0

    dclaw_data.entrainment = 0
    dclaw_data.entrainment_rate = 0.0
    dclaw_data.entrainment_method = 1
    dclaw_data.me = 0.6

    dclaw_data.momlevel = amrdata.amr_levels_max
    dclaw_data.mom_autostop = False

    # == pinitdclaw.data values ==
    pinitdclaw_data = rundata.pinitdclaw_data  # initialized when rundata instantiated

    pinitdclaw_data.init_ptype = 0  # hydrostatic (-1 ==> zero everywhere)

    # == flowgrades.data values ==
    flowgrades_data = rundata.flowgrades_data  # initialized when rundata instantiated

    flowgrades_data.flowgrades = []
    # for using flowgrades for refinement append lines of the form
    # [flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel]
    # where:
    # flowgradevalue: floating point relevant flowgrade value for following measure:
    # flowgradevariable: 1=depth, 2= momentum, 3 = sign(depth)*(depth+topo) (0 at sealevel or dry land).
    # flowgradetype: 1 = norm(flowgradevariable), 2 = norm(grad(flowgradevariable))
    # flowgrademinlevel: refine to at least this level if flowgradevalue is exceeded.

    # flowgrades_data.keep_fine = True
    flowgrades_data.flowgrades.append([1.0e-3, 1, 1, amrdata.amr_levels_max])
    # flowgrades_data.flowgrades.append([1.0e-6, 1, 1, 1])

    # Grid info
    dx_coarse = (clawdata.upper[0] - clawdata.lower[0]) / clawdata.num_cells[0]
    dy_coarse = (clawdata.upper[1] - clawdata.lower[1]) / clawdata.num_cells[1]

    dx_fine = dx_coarse
    dy_fine = dy_coarse

    total_refinement_x = 1
    total_refinement_y = 1

    for i in range(amrdata.amr_levels_max - 1):
        dx_fine /= amrdata.refinement_ratios_x[i]
        dy_fine /= amrdata.refinement_ratios_y[i]
        total_refinement_x *= amrdata.refinement_ratios_x[i]
        total_refinement_y *= amrdata.refinement_ratios_y[i]

    # == fgout_grids.data values ==
    # NEW IN v5.9.0
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2  # will specify a 2d grid of points
    # fgout.output_format = 'binary32'  # 4-byte, float32

    fgout.output_format = "binary32"
    fgout.nx = clawdata.num_cells[0] * total_refinement_x
    fgout.ny = clawdata.num_cells[1] * total_refinement_y
    fgout.x1 = clawdata.lower[0]
    fgout.x2 = clawdata.upper[0]
    fgout.y1 = clawdata.lower[1]
    fgout.y2 = clawdata.upper[1]

    fgout.tstart = 0.0
    fgout.tend = clawdata.tfinal
    fgout.nout = int(clawdata.tfinal) + 1
    fgout.q_out_vars = [1, 2, 3, 4, 5, 8]
    fgout_grids.append(fgout)  # written to fgout_grids.data

    #  ----- For developers -----
    # Toggle debugging print statements:
    amrdata.dprint = False  # print domain flags
    amrdata.eprint = False  # print err est flags
    amrdata.edebug = False  # even more err est flags
    amrdata.gprint = False  # grid bisection/clustering
    amrdata.nprint = False  # proper nesting output
    amrdata.pprint = False  # proj. of tagged points
    amrdata.rprint = False  # print regridding summary
    amrdata.sprint = False  # space/memory output
    amrdata.tprint = False  # time step reporting each level
    amrdata.uprint = False  # update/upbnd reporting

    amrdata.max1d = 300
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata

    # end of function setrun
    # ----------------------


if __name__ == "__main__":
    # Set up run-time parameters and write all data files.
    import sys

    rundata = setrun(*sys.argv[1:])
    rundata.write()

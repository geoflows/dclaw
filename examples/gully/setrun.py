"""
Module to set up run time parameters for D-Claw.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
"""

import os
import sys

import numpy as np
from clawpack.geoclaw import fgmax_tools, fgout_tools

# get information about the extent of the domain and landslide from setinput.py
from setinput import x0, x2, xl1, xl2, y0, y2, yl1, yl2, xr1, xr2, yr1, yr2

try:
    CLAW = os.environ["CLAW"]
except:
    raise Exception("*** Must first set CLAW environment variable")


# This setrun has been set up to make it easy to compare results under a
# variety of numerical options:
# - with/without AMR (amr = True/False)
# - with order = 1 or order = 2
# - with transverse = 1, 2, 3
# - with each type of limiter, limiter = 0, 1, 2, 3, 4
# - multiple amr test options, which control the gridding and timestepping (see
#   below).
# One might create two applications with the same setup (setinput.py,setrun.py,
# setplot.py, and Makefile) that differ only in the numerical options.

amr = True
order = 1
transverse = 0
limiter = 4
amr_test_no = 0


# Use either one level with 25 m grid cells (AMR = False)
# or three levels with 100 -> 50 -> 25 m grid cells.
# keep track of the entire refinement factor (needed for timestep control below)
if amr:
    factor = 1
    mxnest = 3
    refinement_ratios = [2, 2]
else:
    factor = 4
    mxnest = 1
    refinement_ratios = [1]


# amr_test_no = 0
# ---------------
# use standard output and nonuniform time.
# this is the value to use for applications not testing behavior with/without
# amr.

if amr_test_no == 0:
    output_style = 1
    uniform_time = False
    dt_initial = 0.1
    kratio = refinement_ratios
    force_full_region = False
    max1d = 150

# amr_test_no = 1
# ---------------
# this is the simplest test for equivalence between results with and without
# amr.
# - the entire simulation domain is forced to have grids at all levels
# - max1d is large so that each level only has one grid
# - all timesteps on all levels are the same (no refinement in time)

if amr_test_no == 1:
    output_style = 3
    uniform_time = True
    aligned_steps = True
    force_full_region = True
    max1d = 1000


# amr_test_no = 2
# ---------------
# this differs from amr_test_no = 1 in that it relaxes the alignment of all
# timesteps and allows different levels to take different sized timesteps. the
# finest level (level 1 without amr and level 3 with amr) will take equivalent
# timesteps.
# - the entire simulation domain is forced to have grids at all levels
# - max1d is large so that each level only has one grid
# - timesteps on each level differ (refinement in time)

if amr_test_no == 2:
    output_style = 3
    uniform_time = True
    aligned_steps = False
    force_full_region = True
    max1d = 1000


# amr_test_no = 3
# ---------------
# this differs from amr_test_no = 2 in that it does not require that the AMR
# case has a fine grid over the entire domain. It thus adds regridding but
# keeps only one grid on each level.
# - Only flagged regions have refinement.
# - max1d is large so that each level only has one grid
# - timesteps on each level differ (refinement in time)

if amr_test_no == 3:
    output_style = 3
    uniform_time = True
    aligned_steps = False
    force_full_region = False
    max1d = 1000


# amr_test_no = 4
# ---------------
# this differs from amr_test_no = 3 in that max1d is reduced such that there
# are multiple grids on each level.
# - Only flagged regions have refinement.
# - max1d is small so that each level has multiple grids
# - timesteps on each level differ (refinement in time)

if amr_test_no == 4:
    output_style = 3
    uniform_time = True
    aligned_steps = False
    force_full_region = False
    max1d = 50


# set dt variable based on whether using uniform time.
if uniform_time:
    dt_variable = False
else:
    dt_variable = True
variable_dt_refinement_ratios = dt_variable

# Timestep control.
# depending on whether timesteps are aligned, set values for dt_initial, total_steps, output_step_interval, and kratio.
total_time = 100
if amr_test_no > 0:
    if aligned_steps:
        dt_initial = 0.5
        total_steps = int(total_time / dt_initial)
        output_step_interval = 10
        kratio = [1, 1, 1]
    else:
        dt_initial = 1.0 / factor
        total_steps = int(total_time / dt_initial)
        output_step_interval = 5 * factor
        kratio = refinement_ratios

# set target and max CFL conditions based on Order.
if order == 1:
    cfl_desired = 0.45
    cfl_max = 0.50
elif order == 2:
    cfl_desired = 0.75
    cfl_max = 0.85


# ------------------------------
def setrun(claw_pkg="dclaw"):
    # ------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "dclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData
    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == "dclaw", "Expected claw_pkg = 'dclaw'"

    # Number of dimensions
    num_dim = 2
    rundata = data.ClawRunData(claw_pkg, num_dim)

    # ------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    # ------------------------------------------------------------------

    # probdata = rundata.new_UserData(name='probdata',fname='setprob.data')

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
    # x-dimension
    clawdata.lower[0] = x0
    clawdata.upper[0] = x2
    # y-dimension
    clawdata.lower[1] = y0
    clawdata.upper[1] = y2

    # Number of grid cells: Coarsest grid
    clawdata.num_cells[0] = int((x2 - x0) / 100) * factor  # x
    clawdata.num_cells[1] = int((y2 - y0) / 100) * factor  # y

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 9

    # Number of auxiliary variables in the aux array (initialized in setaux)
    # num_aux must be at least the number needed
    clawdata.num_aux = 8

    # Index of aux array corresponding to capacity function, if there is one:
    # typically only used if coordinate_system=2
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

    clawdata.output_style = output_style

    if clawdata.output_style == 1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.num_output_times = 20
        clawdata.tfinal = 100.0
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = output_step_interval
        clawdata.total_steps = total_steps
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
    clawdata.dt_variable = dt_variable

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = dt_initial

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    # D-Claw requires CFL<0.5
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
    clawdata.limiter = [
        limiter,
        limiter,
        limiter,
        limiter,
        limiter,
    ]  # TODO VERIFY THAT 4 in old and new are the same

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
    amrdata.refinement_ratios_x = refinement_ratios
    amrdata.refinement_ratios_y = refinement_ratios
    amrdata.refinement_ratios_t = kratio

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
    amrdata.regrid_buffer_width = 3

    # clustering alg. cutoff for (# flagged pts) / (total # of cells refined)
    # (closer to 1.0 => more small grids may be needed to cover flagged cells)
    amrdata.clustering_cutoff = 0.70

    # print info about each regridding up to this level:
    amrdata.verbosity_regrid = mxnest

    # ---------------
    # Regions:
    # ---------------
    rundata.regiondata.regions = []
    # to specify regions of refinement append lines of the form
    #  [minlevel,maxlevel,t1,t2,x1,x2,y1,y2]

    rundata.regiondata.regions.append(
        [3, 3, 0, 1, xl1 - 200, xl2 + 200, yl1 - 200, yl2 + 200]
    )
    rundata.regiondata.regions.append(
        [3, 3, 0, 1, xr1 - 200, xr2 + 200, yr1 - 200, yr2 + 200]
    )

    # for testing and comparison with/without amr, set mxnest region to entire domain.
    if force_full_region:
        rundata.regiondata.regions.append(
            [
                mxnest,
                mxnest,
                0,
                1e9,
                clawdata.lower[0],
                clawdata.upper[0],
                clawdata.lower[1],
                clawdata.upper[1],
            ]
        )

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
    geo_data.sea_level = -9999.0
    geo_data.dry_tolerance = 1.0e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = 0.025
    geo_data.friction_depth = 1e6
    geo_data.speed_limit = 30

    # Refinement settings
    refinement_data = rundata.refinement_data
    refinement_data.variable_dt_refinement_ratios = variable_dt_refinement_ratios

    # == settopo.data values ==
    topofiles = rundata.topo_data.topofiles
    # for topography, append lines of the form
    #    [topotype, fname]
    topofiles.append([3, "basal_topo.tt3"])

    # == setdtopo.data values ==
    dtopo_data = rundata.dtopo_data

    # == setqinit.data values ==
    qinitdclaw_data = rundata.qinitdclaw_data  # initialized when rundata instantiated

    hfile = "thickness.tt3"
    qinitdclaw_data.qinitfiles.append([3, 1, hfile])

    # == fgmax_grids.data values ==
    # set num_fgmax_val = 1 to save only max depth,
    #                     2 to also save max speed,
    #                     5 to also save max hs,hss,hmin
    rundata.fgmax_data.num_fgmax_val = 5  # Save all quantities
    fgmax_grids = rundata.fgmax_data.fgmax_grids  # empty list to start

    # Now append to this list objects of class fgmax_tools.FGmaxGrid
    # specifying any fgmax grids.
    fg = fgmax_tools.FGmaxGrid()
    fg.point_style = 2  # uniform rectangular x-y grid

    # determine total number of fine cells
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

    fg.x1 = (
        clawdata.lower[0] + dx_fine / 2
    )  # specify pts to align with finite volume cell centers
    fg.x2 = clawdata.upper[0] - dx_fine / 2
    fg.y1 = clawdata.lower[1] + dx_fine / 2
    fg.y2 = clawdata.upper[1] - dx_fine / 2

    fg.dx = dx_fine  # desired resolution of fgmax grid
    fg.tstart_max = 0.0  # when to start monitoring max values
    fg.tend_max = 1.0e10  # when to stop monitoring max values
    fg.dt_check = 1.0  # target time (sec) increment between updating
    # max values
    fg.min_level_check = 1  # which levels to monitor max on
    fg.arrival_tol = 1.0e-2  # tolerance for flagging arrival
    fg.interp_method = 0  # 0 ==> pw const in cells, recommended
    # fgmax_grids.append(fg)  # written to fgmax_grids.data

    # == fgout_grids.data values ==
    # Set rundata.fgout_data.fgout_grids to be a list of
    # objects of class clawpack.geoclaw.fgout_tools.FGoutGrid:
    fgout_grids = rundata.fgout_data.fgout_grids  # empty list initially

    fgout = fgout_tools.FGoutGrid()
    fgout.fgno = 1
    fgout.point_style = 2  # will specify a 2d grid of points
    # fgout.output_format = 'binary32'  # 4-byte, float32
    fgout.output_format = "ascii"  # 4-byte, float32
    fgout.nx = int(clawdata.num_cells[0] * total_refinement_x)
    fgout.ny = int(clawdata.num_cells[1] * total_refinement_y)
    fgout.x1 = clawdata.lower[0]  # specify edges (fgout pts will be cell centers)
    fgout.x2 = clawdata.upper[0]
    fgout.y1 = clawdata.lower[1]
    fgout.y2 = clawdata.upper[1]
    fgout.tstart = 0.0
    fgout.tend = 5
    fgout.nout = 6
    fgout.q_out_vars = [1, 4, 8]
    # fgout_grids.append(fgout)    # written to fgout_grids.data

    # == setauxinit.data values ==
    auxinitdclaw_data = rundata.auxinitdclaw_data  # initialized when rundata instantiated
    dhdtfile = "dhdt.tt3"
    auxinitdclaw_data.auxinitfiles.append([3, 8, dhdtfile]) # with coordinate system 1, dhdt is aux8

    # == fgmax.data values ==
    # fgmax_files = rundata.fgmax_data.fgmax_files
    # for fixed grids append to this list names of any fgmax input files

    # == setdclaw.data values ==
    dclaw_data = rundata.dclaw_data  # initialized when rundata instantiated

    dclaw_data.rho_f = 1000.0
    dclaw_data.rho_s = 2700.0
    dclaw_data.m_crit = 0.64
    dclaw_data.m0 = 0.63
    dclaw_data.mref = 0.63
    dclaw_data.kref = 1.0e-10
    dclaw_data.phi = 32.0
    dclaw_data.delta = 0.001
    dclaw_data.mu = 0.005
    dclaw_data.alpha_c = 0.05
    dclaw_data.c1 = 1.0
    dclaw_data.sigma_0 = 1.0e3
    dclaw_data.dd_manning = True

    dclaw_data.src2method = 2
    dclaw_data.alphamethod = 1

    dclaw_data.segregation = 0
    dclaw_data.beta_seg = 0.0
    dclaw_data.chi0 = 0.0
    dclaw_data.chie = 0.0

    dclaw_data.bed_normal = 0
    dclaw_data.theta_input = 0.0

    dclaw_data.entrainment = 0
    dclaw_data.entrainment_rate = 0.0
    dclaw_data.entrainment_method = 1
    dclaw_data.me = 0.63

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

    flowgrades_data.keep_fine = True
    flowgrades_data.flowgrades.append([1.0e-6, 1, 1, 3])

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

    amrdata.max1d = max1d
    # More AMR parameters can be set -- see the defaults in pyclaw/data.py

    return rundata

    # end of function setrun
    # ----------------------


if __name__ == "__main__":
    # Set up run-time parameters and write all data files.
    import sys

    rundata = setrun(*sys.argv[1:])
    rundata.write()

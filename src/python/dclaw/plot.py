"""
Functions used for plotting D-Claw results using `visclaw <https://www.clawpack.org/plotting.html>`__

When using `visclaw <https://www.clawpack.org/plotting.html>`__
it is common to provide a function to the attribute ``plot_var`` of a
`ClawPlotItem <https://www.clawpack.org/ClawPlotItem.html#attributes>`__
to specify what variable is plotted.

These functions provide the most commonly used functions for D-Claw.
"""

import numpy as np
from numpy import ma as ma

# Indices to elements of q and aux

# q indices
_i_h = 0
_i_hu = 1
_i_hv = 2
_i_hm = 3
_i_pb = 4
_i_hchi = 5
_i_beroded = 6
_i_eta = -1

# aux indices
_i_topo = 0
_i_capax = 1
_i_capay = 2

# if coordinate_system == 1 _i_dig = 1
# else _i_dig = 3

_i_dig = 1
_i_phi = _i_dig
_i_theta = _i_dig + 1
_i_taudir_x = _i_dig + 4
_i_taudir_y = _i_dig + 5

_grav = 9.81
_rho_f = 1000
_rho_s = 2700
_kappita = 1e-8
_mu = 0.001
_m_crit = 0.52
_delta = 0.01
_m0 = 0.6
_c1 = 1.0
_bed_normal = 0
_drytol = 0.001
_phi = 40

# TODO update all elements here to reflect new changes to the data model
# (e.g., mref, kref, etc)


# Gravity, adjusted for bed normal.
def gmod(current_data):
    """Gravity, adjusted for bed normal"""
    if hasattr(current_data.plotdata, "geoclaw_data"):
        grav = current_data.plotdata.geoclaw_data.gravity
    else:
        grav = _grav

    if hasattr(current_data.plotdata, "dclaw_data"):
        bed_normal = current_data.plotdata.dclaw_data.bed_normal
    else:
        bed_normal = _bed_normal

    if bed_normal == 1:
        aux = current_data.aux
        theta = aux[_i_theta, :, :]
        gmod = grav * np.cos(theta)
    else:
        gmod = grav

    return gmod


# Upper surface and topo
def eta(current_data):
    """
    Return eta = h + b
    """
    q = current_data.q
    eta = q[_i_eta, :, :]
    return eta


def topo(current_data):
    """
    Return topography = eta - h.
    """
    q = current_data.q
    h = q[_i_h, :, :]
    eta = q[_i_eta, :, :]
    topo = eta - h
    return topo


def land(current_data):
    """
    Return a masked array containing the surface elevation only in dry cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    eta = q[_i_eta, :, :]
    land = ma.masked_where(h > drytol, eta)
    return land


def surface(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    eta = q[_i_eta, :, :]
    water = ma.masked_where(h <= drytol, eta)
    return water


def surface_solid_frac_lt03(current_data):
    """
    Return a masked array containing the surface elevation only in wet cells.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    eta = q[_i_eta, :, :]
    hm = q[_i_hm, :, :]

    with np.errstate(divide="ignore", invalid="ignore"):
        m = hm / h

    water = ma.masked_where(h <= drytol, eta)
    water = ma.masked_where(m > 0.3, water)

    return water


# Depth
def depth(current_data):
    """
    Return a masked array containing the depth of fluid only in wet cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    depth = ma.masked_where(h <= drytol, h)
    return depth


def surface_or_depth(current_data):
    """
    Return a masked array containing the surface elevation where the topo is
    below sea level or the water depth where the topo is above sea level.
    Mask out dry cells.  Assumes sea level is at topo=0.
    Surface is eta = h+topo, assumed to be output as 4th column of fort.q
    files.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    eta = q[_i_eta, :, :]
    topo = eta - h
    surface = ma.masked_where(h <= drytol, eta)
    depth = ma.masked_where(h <= drytol, h)
    surface_or_depth = np.where(topo < 0, surface, depth)
    return surface_or_depth


# Velocity
def velocity_u(current_data):
    """
    Return a masked array containing velocity u in wet cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hu = q[_i_hu, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = ma.masked_where(h <= drytol, hu / h)
    return u


def velocity_v(current_data):
    """
    Return a masked array containing velocity v in wet cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hv = q[_i_hv, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        v = ma.masked_where(h <= drytol, hv / h)
    return v


def velocity_unorm(current_data):
    """
    Return a masked array containing velocity u in wet cells.
    Normalized by the velocity magnitude.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hu = q[_i_hu, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = ma.masked_where(h <= drytol, hu / h)
    return u / velocity_magnitude(current_data)


def velocity_vnorm(current_data):
    """
    Return a masked array containing velocity v in wet cells.
    Normalized by the velocity magnitude.

    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hv = q[_i_hv, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        v = ma.masked_where(h <= drytol, hv / h)
    return v / velocity_magnitude(current_data)


def velocity(current_data):
    """
    Return a masked array containing a tuple of x and y directed velocity (u,v)
    in wet cells.

    velocity defined as sqrt(u**2 + v**2)
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hu = q[_i_hu, :, :]
    hv = q[_i_hv, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = ma.masked_where(h <= drytol, hu / h)
        v = ma.masked_where(h <= drytol, hv / h)
    return (u, v)


def velocity_magnitude(current_data):
    """
    Return a masked array of the magnitude of velocity at wet cells.

    velocity defined as sqrt(u**2 + v**2)
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hu = q[_i_hu, :, :]
    hv = q[_i_hv, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        u = hu / h
        v = hv / h
    vel = np.sqrt(u**2 + v**2)

    vel = ma.masked_where(h <= drytol, vel)

    return vel


def solid_frac(current_data):
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hm = q[_i_hm, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        m = ma.masked_where(h < drytol, hm / h)
    return m


def solid_frac_gt03(current_data):
    m = solid_frac(current_data)
    return ma.masked_where(m < 0.3, m)


def basalP(current_data):
    # basal pressure.
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    basalP = ma.masked_where(q[_i_h, :, :] < drytol, q[_i_pb, :, :])
    return basalP


def species1_fraction(current_data):
    """
    Return a masked array containing the fraction of species A in wet cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hchi = q[_i_hchi, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi1 = ma.masked_where(h <= drytol, hchi / h)
    return chi1


def species2_fraction(current_data):
    """
    Return a masked array containing the fraction of species B in wet cells.
    """
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    h = q[_i_h, :, :]
    hchi = q[_i_hchi, :, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = ma.masked_where(h <= drytol, 1.0 - (hchi / h))
    return chi2


def b_eroded(current_data):
    # eroded depth
    q = current_data.q
    b_eroded = q[_i_beroded, :, :]
    # if np.any(b_eroded>0):
    #        print(f'b_eroded gt zero, {np.sum(b_eroded>0)}')
    return b_eroded


# Ancillary values


def density(current_data):
    # new segregation might modify.
    m = solid_frac(current_data)

    if hasattr(current_data.plotdata, "dclaw_data"):
        rho_f = current_data.plotdata.dclaw_data.rho_f
        rho_s = current_data.plotdata.dclaw_data.rho_s
    else:
        rho_f = _rho_f
        rho_s = _rho_s

    rho = (1.0 - m) * rho_f + m * rho_s
    return rho


def kperm(current_data):
    # permeability (m^2)

    if hasattr(current_data.plotdata, "dclaw_data"):
        kappita = current_data.plotdata.dclaw_data.kappita
        m0 = current_data.plotdata.dclaw_data.m0
    else:
        kappita = _kappita
        m0 = _m0

    m = solid_frac(current_data)
    return kappita * np.exp(-(m - m0) / (0.04))


def shear(current_data):
    # units of 1/second
    h = depth(current_data)
    vnorm = velocity_magnitude(current_data)
    return 2.0 * vnorm / h  # in code refered to as hbounded (defined as h)


def dilatency(current_data):
    # depth averaged dilatency rate (m/s)

    if hasattr(current_data.plotdata, "dclaw_data"):
        mu = current_data.plotdata.dclaw_data.mu
    else:
        mu = _mu

    h = depth(current_data)
    # Royal Society, Part 2, Eq 2.6
    D = (
        2.0
        * (kperm(current_data) / (mu * h))
        * hydrostatic_minus_basal_pressure(current_data)
    )
    vnorm = velocity_magnitude(current_data)
    D[vnorm <= 0] = 0
    # depth averaged dilatency has units of L/T (this is consistent with Part 1 Eq 4.6)
    # k [=] L**2
    # mu [=] Pa-s = M / (L * T)
    # rho [=] M/L**3
    # g [=] L/T**2
    # h [=] L
    # L**2 / ((ML)/(LT)) * (M/L**3)*L*(L/T**2)
    # (L**2 T / M) * (M/ (LT**2))
    # L/T
    return D


def tanpsi(current_data):
    # tangent of dilation angle (#)

    if hasattr(current_data.plotdata, "dclaw_data"):
        c1 = current_data.plotdata.dclaw_data.c1
    else:
        c1 = _c1

    gamma = shear(current_data)
    # in code, m-meqn is regularized based on shear
    vnorm = velocity_magnitude(current_data)
    tpsi = c1 * m_minus_meqn(current_data) * np.tanh(gamma / 0.1)
    tpsi[vnorm <= 0] = 0
    return tpsi


def psi(current_data):
    return np.arctan(
        m_minus_meqn(current_data)
    )  # maybe this should be arctan of tanpsi (with the regularization as is discussed for tanpsi


# Related to m_eqn
def m_crit(current_data):
    # eventually this may need modification based on segregation (just like kperm)
    if hasattr(current_data.plotdata, "dclaw_data"):
        return current_data.plotdata.dclaw_data.m_crit

    else:
        return _m_crit


def m_eqn(current_data):
    # equilibrium value for m (not currently correct if segregation is used (but segregation may change))
    m_c = m_crit(current_data)
    m_eqn = m_c / (1.0 + np.sqrt(N(current_data)))
    return m_eqn


def meqn_over_mcrit(current_data):
    return m_eqn(current_data) / m_crit(current_data)


def m_minus_meqn(current_data):
    return solid_frac(current_data) - m_eqn(current_data)


def m_minus_mcrit(current_data):
    return solid_frac(current_data) - m_crit(current_data)


# Calculated based on pressure
def basal_pressure_over_hydrostatic(current_data):
    return basalP(current_data) / hydrostaticP(current_data)


def lithostaticP(current_data):
    # lithostatic pressure
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    h = depth(current_data)
    rho = density(current_data)
    var = ma.masked_where(h < drytol, gmod(current_data) * rho * h)
    return var


def hydrostaticP(current_data):
    # hydrostatic pressure
    if hasattr(current_data.plotdata, "dclaw_data"):
        rho_f = current_data.plotdata.dclaw_data.rho_f
    else:
        rho_f = _rho_f

    h = depth(current_data)
    return gmod(current_data) * rho_f * h


def sigma_e(current_data):
    # effective basal pressure (lithostatic less basal pressure)
    se = lithostaticP(current_data) - basalP(current_data)
    se[se < 0.0] = 0.0  # cannot be negative.
    return se


def hydrostatic_minus_basal_pressure(current_data):
    # effective basal pressure (lithostatic less basal pressure)
    se = hydrostaticP(current_data) - basalP(current_data)
    # can be negative.
    return se


def sigma_e_over_hydrostatic(current_data):
    # effective basal pressure over rho_f * g * h
    return sigma_e(current_data) / hydrostaticP(current_data)


def sigma_e_over_lithostatic(current_data):
    # this is the same as the liquifaction ratio.
    return sigma_e(current_data) / lithostaticP(current_data)


def liquefaction_ratio(current_data):
    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    q = current_data.q
    p = ma.masked_where(q[_i_h, :, :] < drytol, q[_i_pb, :, :])
    litho = lithostaticP(current_data)
    ratio = ma.masked_where(q[_i_h, :, :] < drytol, p / litho)
    return ratio


# Nondimentional numbers


def nondimentional_c(current_data):

    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    if hasattr(current_data.plotdata, "dclaw_data"):
        rho_f = current_data.plotdata.dclaw_data.rho_f
        kappita = current_data.plotdata.dclaw_data.kappita
        mu = current_data.plotdata.dclaw_data.mu
    else:
        rho_f = _rho_f
        kappita = _kappita
        mu = _mu

    U = velocity_magnitude(current_data)
    g = gmod(current_data)

    c = np.array((rho_f * g * kappita) / (mu * U), dtype=float)

    q = current_data.q
    h = q[_i_h, :, :]

    # set to zero if nan or inf
    c[np.isnan(c)] = 0
    c[np.isinf(c)] = 0
    c_out = ma.masked_where(h <= drytol, c)
    return c_out


def Iv(current_data):
    # inertial number
    if hasattr(current_data.plotdata, "dclaw_data"):
        mu = current_data.plotdata.dclaw_data.mu
    else:
        mu = _mu

    gamma = shear(current_data)
    return (mu * gamma) / sigma_e(current_data)



def Froude(current_data):
    # Froude number
    if hasattr(current_data.plotdata, "geoclaw_data"):
        grav = current_data.plotdata.geoclaw_data.gravity
    else:
        grav = _grav

    h = depth(current_data)
    u = velocity_magnitude(current_data)
    return u/np.sqrt(grav*h)


def Stokes(current_data):
    # stokes number

    if hasattr(current_data.plotdata, "dclaw_data"):
        mu = current_data.plotdata.dclaw_data.mu
        rho_s = current_data.plotdata.dclaw_data.rho_s
        delta = current_data.plotdata.dclaw_data.delta
    else:
        mu = _mu
        rho_s = _rho_s
        delta = _delta

    gamma = shear(current_data)
    return rho_s * gamma * delta**2 / mu


def N(current_data):  # dimensionless state parameter N

    if hasattr(current_data.plotdata, "dclaw_data"):
        mu = current_data.plotdata.dclaw_data.mu
        rho_s = current_data.plotdata.dclaw_data.rho_s
        delta = current_data.plotdata.dclaw_data.delta
    else:
        mu = _mu
        rho_s = _rho_s
        delta = _delta

    gamma = shear(current_data)
    sigbedc = (rho_s * ((gamma * delta) ** 2.0)) + sigma_e(current_data)
    N = (mu * gamma) / (sigbedc)
    N[sigbedc < 0.0] = 0.0
    # print(N.max())
    return N


# Force balance calculations


def local_slope(current_data):
    h = depth(current_data)

    # # gravitational driving force per unit area.
    # if hasattr(current_data.plotdata, "dclaw_data"):
    #     bed_normal = current_data.plotdata.dclaw_data.bed_normal
    # else:
    #     bed_normal = _bed_normal

    eta = surface(current_data)

    if hasattr(current_data.plotdata, "geoclaw_data"):
        drytol = current_data.plotdata.geoclaw_data.dry_tolerance
    else:
        drytol = _drytol

    # if bed_normal == 1:
    #     q = current_data.q
    #     theta = q[_i_theta, :, :]
    #     sintheta = np.sin(theta)
    # else:
    #     sintheta = 0.0

    dx = current_data.dx
    dy = current_data.dy

    etaL = np.roll(eta.copy(), 1, axis=1)
    etaL[:, 0] = np.nan
    etaR = np.roll(eta.copy(), -1, axis=1)
    etaR[:, -1] = np.nan
    etaB = np.roll(eta.copy(), 1, axis=0)
    etaB[0, :] = np.nan
    etaT = np.roll(eta.copy(), -1, axis=0)
    etaT[-1, :] = np.nan

    row, col = eta.shape

    detadx = np.zeros((row, col, 4))
    detadx[:, :, 0] = np.abs((eta - etaL) / (dx))
    detadx[:, :, 1] = np.abs((eta - etaB) / (dy))
    detadx[:, :, 2] = np.abs((etaR - eta) / (dx))
    detadx[:, :, 3] = np.abs((etaT - eta) / (dy))

    maxdetadx = np.max(detadx, axis=-1)
    slope = np.rad2deg(np.arctan(maxdetadx))

    return ma.masked_where(h <= drytol, slope)


def Fgravitational(current_data):
    # gravitational driving force per unit area.
    if hasattr(current_data.plotdata, "dclaw_data"):
        bed_normal = current_data.plotdata.dclaw_data.bed_normal
    else:
        bed_normal = _bed_normal

    g = gmod(current_data)
    h = depth(current_data)
    eta = surface(current_data)
    rho = density(current_data)

    if bed_normal == 1:
        q = current_data.q
        theta = q[_i_theta, :, :]
        sintheta = np.sin(theta)
    else:
        sintheta = 0.0

    dx = current_data.dx
    dy = current_data.dy

    hL = np.roll(
        h.copy(), 1, axis=1
    )  # roll right on columns so that value at (i, j-1) is at (i,j)
    hL[:, 0] = np.nan  # first column has undefined values
    hR = np.roll(h.copy(), -1, axis=1)
    hR[:, -1] = np.nan
    hB = np.roll(h.copy(), 1, axis=0)
    hB[0, :] = np.nan
    hT = np.roll(h.copy(), -1, axis=0)
    hT[-1, :] = np.nan

    etaL = np.roll(eta.copy(), 1, axis=1)
    etaL[:, 0] = np.nan
    etaR = np.roll(eta.copy(), -1, axis=1)
    etaR[:, -1] = np.nan
    etaB = np.roll(eta.copy(), 1, axis=0)
    etaB[0, :] = np.nan
    etaT = np.roll(eta.copy(), -1, axis=0)
    etaT[-1, :] = np.nan

    FxL = rho * (
        np.abs(
            -g * 0.5 * (h + hL) * (eta - etaL) / (dx) + g * 0.5 * (h + hL) * sintheta
        )
    )
    FyB = rho * (np.abs(-g * 0.5 * (h + hB) * (eta - etaB) / (dy)))

    FxR = rho * (
        -g * 0.5 * (h + hR) * (etaR - eta) / (dx) + g * 0.5 * (h + hR) * sintheta
    )
    FyT = rho * (-g * 0.5 * (h + hT) * (etaT - eta) / (dy))

    # units are M/L**3 * L/T**2 * L
    # = M / (L * T**2) = Pressure  =  Force per unit area (OK)
    different_sign_x = (FxL * FxR) < 0
    different_sign_y = (FyT * FyB) < 0

    Fx = np.abs(FxL)
    FxR_mag_smaller = np.abs(FxR) < np.abs(FxL)
    Fx[FxR_mag_smaller] = np.abs(FxR)[FxR_mag_smaller]
    Fx[different_sign_x] = 0

    Fy = np.abs(FyT)
    FyB_mag_smaller = np.abs(FyB) < np.abs(FyT)
    Fy[FyB_mag_smaller] = np.abs(FyB)[FyB_mag_smaller]
    Fy[different_sign_y] = 0
    Fg = np.sqrt(Fx**2, Fy**2)
    # Royal Proceedings, Part 2, equation 2.4b,c (momentum source terms) first term on RHS
    # deta/dx portion taken from calc_taudir
    return Fg


def Fdrag(current_data):  # units of force per unit area
    # Royal Proceedings, Part 2, equation 2.4b,c (momentum source terms) second term on RHS
    # what is this term?
    # katy asks: driving force due to longitudinal stress gradients?
    # (based on text right before part 2 eq 2.15
    # dave says: I don't know if there's a simple interpretation...it sort of
    # drops out from the derivation and then rearrangement of the equations into
    # conservative form (ie. derivatives on hu not u). I think it might be
    # similar to a drag term that appears on fully two-phase equations for solid
    # and fluid velocity fields.

    if hasattr(current_data.plotdata, "dclaw_data"):
        rho_f = current_data.plotdata.dclaw_data.rho_f
    else:
        rho_f = _rho_f

    vnorm = velocity_magnitude(current_data)
    D = dilatency(current_data)
    rho = density(current_data)
    # units:
    # dilatency has units L/T
    # L/T * L/T * M/L*3
    # L**2 M/(L**3 T**2)
    # M/(L T**2) -> Pressure, Force per unit area, OK.
    return vnorm * D * (rho - rho_f)


def Ffluid(current_data):  # units of force per unit area
    # Resisting force due to fluid.
    if hasattr(current_data.plotdata, "dclaw_data"):
        mu = current_data.plotdata.dclaw_data.mu
    else:
        mu = _mu

    h = depth(current_data)
    m = solid_frac(current_data)
    vnorm = velocity_magnitude(current_data)
    tauf = 2.0 * mu * (1.0 - m) * vnorm / h
    # units
    # mu [=] Pa-s = M / (L * T)
    # M/(LT) * L/T &* 1/L
    # M/(LT**2) --> pressure, force per unit area, OK
    return tauf


def phi(current_data):
    # angle of internal friction (radians)
    if hasattr(current_data.plotdata, "dclaw_data"):
        _p = current_data.plotdata.dclaw_data.phi
    else:
        _p = _phi  # otherwise use default _phi above.

    # results are returned in radians.
    return np.deg2rad(_p)


def Fsolid(current_data):  # units of force per unit area
    # resisting force due to solid.

    tau_s = sigma_e(current_data) * np.tan(phi(current_data) + psi(current_data))
    # tau = dmax1(0.d0,mreg*sigbed*tan(atan(mu_bf)+atan(tanpsi)))
    tau_s[tau_s < 0] = 0
    # units
    # sigma_e * [-]
    # sigma_e is pressure, so OK.
    # where material is static, tau_s must be no greater than Fdriving.
    vnorm = velocity_magnitude(current_data)
    Fd = Fdriving(current_data)
    static = vnorm == 0
    reduce_tau_s = static * (tau_s > Fd)
    tau_s[reduce_tau_s] = Fd[reduce_tau_s]
    return tau_s


def Fdriving(current_data):  # units of force per unit area
    # driving force.
    return Fgravitational(current_data) + Fdrag(current_data)


def Fresisting(current_data):  # units of force per unit area
    # resisting force
    return Ffluid(current_data) + Fsolid(current_data)


def Fnet(current_data):  # units of force per unit area
    # at lowest Fnet is zero because friction will just balance driving force
    return Fdriving(current_data) - Fresisting(current_data)


# Stability related values


def static_angle(current_data):
    se_over_sl = sigma_e_over_lithostatic(current_data)
    # phi is in radians
    # static limit, so choose psi = 0

    deta_dx = se_over_sl * np.tan(0.0 + phi(current_data))
    theta = np.arctan(deta_dx)
    theta_deg = np.rad2deg(theta)
    return theta_deg

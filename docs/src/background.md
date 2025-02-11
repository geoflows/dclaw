# Background

D-Claw is a single-layer, depth-averaged, two-phase model that simulates the evolution of granular-fluid mixtures over arbitrary topography based on conservation of mass and momentum (George and Iverson, 2014; Iverson and George, 2014). The D-Claw software was originally developed as an extension of the GeoClaw software (Berger et al., 2011; George, 2008; LeVeque, 2002; LeVeque et al., 2011; Mandli et al., 2016) for clear-water flows, including coseismic tsunami propagation and inundation, overland flow, and dam-breach outburst flood inundation (George, 2011). In the absence of solid granular material (solid volume fraction equal to zero), the D-Claw equations reduce to the nonlinear shallow water equations—the governing equations for GeoClaw. 

D-Claw is written in fortran and has a python wrapper. It installs and runs as a submodule of the Clawpack software, which must be installed in additon to D-Claw. An interested user of D-Claw should expect to read the Clawpack documentation. 

## Model, Numerics, Application

The D-Claw model is the combination of
- Mathematical governing equations (see Iverson & George, 2014 and [the documentation page on theory](#theory))
- Numerical algortithms (see George and Iverson, 2014 and Clawpack references listed below)
- Software implementation

Any individual simulation run with the D-Claw model requires the specification of parameter values (physical and numerical), and spatially variable initial and boundary conditions. Each simulation is typically called an "application" and is run in an application directory

## Filestructure

Within each application directory, running D-Claw will generate multiple output files. See [the Clawpack documentation](https://www.clawpack.org/contents.html) for information about these output files.

## Visualization

Because D-Claw uses adaptive mesh plotting of results can sometimes be a challenge. The [Visclaw](https://www.clawpack.org/contents.html#visclaw-plotting-and-visualization-tools) plotting and visualization tools are a useful resource for plotting results.

## References 

Berger, M. J., George, D. L., LeVeque, R. J., & Mandli, K. T. (2011). The GeoClaw software for depth-averaged flows with adaptive refinement. Advances in Water Resources, 34(9), 1195–1206, [https://doi.org/10.1016/j.advwatres.2011.02.016](https://doi.org/10.1016/j.advwatres.2011.02.016)

George, D. L. (2008). Augmented Riemann solvers for the shallow water equations over variable topography with steady states and inundation. Journal of Computational Physics, 227(6), 3089–3113. [https://doi.org/10.1016/j.jcp.2007.10.027](https://doi.org/10.1016/j.jcp.2007.10.027)

George, D.L., and Iverson, R.M., 2014, A depth-averaged debris-flow model that includes the effects of evolving dilatancy—II. Numerical predictions and experimental tests: Proceedings of the Royal Society of London. Series A, v. 470, no. 2170, p. 20130820, [https://doi.org/10.1098/rspa.2013.0820](https://doi.org/10.1098/rspa.2013.0820).

George, D.L., 2011, Adaptive finite volume methods with well-balanced Riemann solvers for modeling floods in rugged terrain: Application to the Malpasset dam-break flood (France, 1959). Int. J. Numer. Meth. Fluids, 66: 1000-1018, [https://doi.org/10.1002/fld.2298](https://doi.org/10.1002/fld.2298).

Iverson, R.M., and George, D.L., 2014, A depth-averaged debris-flow model that includes the effects of evolving dilatancy—I. Physical basis: Proceedings of the Royal Society of London. Series A, v. 470, no. 2170, p. 20130819, [https://doi.org/10.1098/rspa.2013.0819](https://doi.org/10.1098/rspa.2013.0819).

LeVeque, R. J., 2002, Finite Volume Methods for Hyperbolic Problems. Cambridge: Cambridge University Press, [https://doi.org/10.1017/CBO9780511791253](https://doi.org/10.1017/CBO9780511791253).

LeVeque, R. J., George, D. L., & Berger, M. J., 2011, Tsunami modelling with adaptively refined finite volume methods. Acta Numerica, 20, 211–289, [https://doi.org/10.1017/S0962492911000043](https://doi.org/10.1017/S0962492911000043).

Mandli, K. T., Ahmadia, A. J., Berger, M., Calhoun, D., George, D. L., Hadjimichael, Y., et al., 2016, Clawpack: building an open source ecosystem for solving hyperbolic PDEs. PeerJ Computer Science, 2, e68 [https://doi.org/10.7717/peerj-cs.68](https://doi.org/10.7717/peerj-cs.68).


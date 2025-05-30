Tests
=====

There is one test located in ``dclaw/tests/bowl-slosh`` and described below. Additional not-yet-ready tests are located in ``dclaw/tests/dev``. The applications within this test development directory are not expected to work.


Sloshing water in a parabolic bowl
----------------------------------

This is a modified version of the `bowl-slosh test in geoclaw <https://github.com/clawpack/geoclaw/tree/master/tests/bowl_slosh>`_ .


Background
^^^^^^^^^^

Waves in a parabolic bowl with a flat surface sloshing around.
An exact analytic solution is known in which the surface stays flat.

In this code, :math:`x` and :math:`y` are in meters (coordinate_system=1
in `setrun.py`).

Topography: :math:`B(x,y) = h_0((x^2 + y^2)/a^2 -1)`,

Depth: :math:`h(x,y,t) = \max\left(0,~~ (\sigma h_0/a^2)(2x\cos(\omega t) + 2y\sin(\omega t) -
\sigma) - B(x,y)\right)`

Velocities:  :math:`u(x,y,t) = -\sigma \omega \sin(\omega t),\qquad
v(x,y,t) = \sigma \omega \cos(\omega t).`

where :math:`\omega = \sqrt{2gh_0} / a`.

The period of oscillation is  :math:`T = 2\pi / \omega`.

The following parameters are currently hardwired several places:

:math:`a = 1, ~~\sigma = 0.5, ~~h = 0.1,~~g = 9.81`

Instructions
^^^^^^^^^^^^

This example uses a custom ``qinit.f90``. First recompile the code from within the ``bowl_slosh`` directory::

	make clean
	make clobber
	make .exe

Create the topo file before running the code::

    make topo
    make .data
    make .output
    make .plots

Alternatively, you may run:

	make all

to run all steps. Examine the plots by opening ``bowl-slosh/_plots/_PlotIndex.html`` in a browser. We suggest the cross-section movies as the best means for simluation-analytical solution comparison. Finally, we suggest that you also run the equivalent test case in geoclaw (``geoclaw/tests/bowl-slosh``) and compare the output.


.. note::

	This should be cleaned up: better to put them in a setprob.data file that is read in where needed.


References
----------

* W. C. Thacker, Some exact solutions to the nonlinear shallow water wave equations,
  J. Fluid Mech. 107 (1981), 499-508.

* J.M. Gallardo, C. Pares, and M. Castro, On a well-balanced high-order
  finite volume scheme for shallow water equations with topography and dry
  areas, J. Comput. Phys. 227(2007) 574-601.

* Y. Xing, X. Zhang and C.-W. Shu, Positivity preserving high order well
  balanced discontinuous Galerkin methods for the shallow water equations ,
  Advances in Water Resources  33 (2010), pp. 1476-1493.

This test problem has been used in several other papers too.

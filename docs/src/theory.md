# Theory

This page describes the equations solved by D-Claw. A user interested in the details of the equation derivation is referred to Iverson and George (2014). Those interested in an explanation of the numerical implementation are referred to George and Iverson (2014).

## State variables and auxiliary variables

D-Claw considers the flow of material under gravity ({math}`\vec{g} = (g_x,g_y,g_z)^\mathrm{T}`) with thickness {math}`h`, x- and y- directed velocities {math}`u` and {math}`v`, solid-volume fraction {math}`m`, and basal pore-fluid pressure {math}`p_b` over an arbitrary surface {math}`b`. Optionally, entrainment of erodible sediment and particle size segregation are implemented. The depth of eroded and entrained material is given by {math}`\Delta b` (positive indicating entrainment has occurred). Segregation considers two species fractions {math}`A` and {math}`B` with {math}`\chi` representing the fraction of species {math}`A`.

The variables which vary in space and time are:

:::{list-table}
:widths: 5 20 10
:header-rows: 1

*   - Name
    - Description
    - Units
*   - {math}`h`
    - Flow depth
    - meters
*   - {math}`hu`
    - Flow x-momentum: depth times x-directed velocity
    - meters squared per second
*   - {math}`hv`
    - Flow y-momentum: depth times y-directed velocity
    - meters squared per second
*   - {math}`hm`
    - Solid volume: depth times solid-volume fraction
    - meters
*   - {math}`p_b`
    - Basal pore-fluid pressure
    - kilograms per meter per time squared
*   - {math}`h\chi`
    - Species A volume: depth times species A fraction
    - meters
*   - {math}`\Delta b`
    - Depth of erosion (material entrained)
    - meters
:::

The variables which vary in space only are

:::{list-table}
:widths: 5 20 10
:header-rows: 1

*   - Name
    - Description
    - Units
*   - {math}`b`
    - Topobathymetric surface elevation
    - meters
*   - {math}`\Theta`
    - Slope angle in x-direction (used only if bed-normal coordinates are enabled)
    - degrees
*   - {math}`h_e`
    - Initial thickness of material that can be entrained.
    - meters
:::


## Core equations

D-Claw solves the following set of partial differential equations. See Iverson and George (2014) and George and Iverson (2014) for derivation, explanation, and references.


```{math}
\frac{\partial h}{\partial t} +
\frac{\partial (hu)}{\partial x} +
\frac{\partial (hv)}{\partial y}=
\varphi_1
```

```{math}
\frac{\partial (hu)}{\partial t} +
\frac{\partial}{\partial x} \left ( hu^2\right ) +
\kappa \frac{\partial}{\partial x}\left (\frac{1}{2} g_z h^2\right ) +
\frac{\partial (huv)}{\partial y} +  
\frac{h(1-\kappa)}{\rho}\frac{\partial p_b}{\partial x} =
\varphi_2
```

```{math}
\frac{\partial (hv)}{\partial t} +
\frac{\partial (huv)}{\partial x} +
\frac{\partial}{\partial y} \left ( hv^2 \right )  +
\kappa\frac{\partial}{\partial y}\left (\frac{1}{2}  g_z h^2\right ) +  
\frac{h(1-\kappa)}{\rho}\frac{\partial p_b}{\partial y} =
\varphi_3
```

```{math}
\frac{\partial (hm)}{\partial t} +
\frac{\partial (hum)}{\partial x} +
\frac{\partial (hvm)}{\partial y} =
\varphi_4
```

```{math}
\frac{\partial p_b}{\partial t} -
\rho g_z\varrho u \frac{\partial h}{\partial x}  +
\rho g_z\varrho \frac{\partial (hu)}{\partial x} +
u\frac{\partial p_b}{\partial x} -
\rho g_z\varrho  v \frac{\partial h}{\partial y}  +
\rho g_z\varrho \frac{\partial (hv)}{\partial y} +
v\frac{\partial p_b}{\partial y}= \varphi_5
```

where {math}`\varphi_1`,  {math}`\varphi_2`,  {math}`\varphi_3`,  {math}`\varphi_4`, and {math}`\varphi_5`  represent the source terms.

The bulk density of the flow, {math}`\rho`, is calculated based on the fluid and solid densities, {math}`\rho_f` and {math}`\rho_s`, respectively.

```{math}
\rho = \rho_s m + (1-m)\rho_f
```

The lateral pressure coefficient {math}`\kappa` is typically set to {math}`\kappa=1`.

The source terms {math}`\varphi_1--\varphi_5` are defined as

```{math}
\varphi_1 = \frac{(\rho-\rho_f)}{\rho} D
```

```{math}
\varphi_2 = hg_x + uD\frac{(\rho-\rho_f)}{\rho}-\frac{(\tau_{s,x} + \tau_{f,x})}{\rho}
```

```{math}
\varphi_3 = hg_y + vD\frac{(\rho-\rho_f)}{\rho}-\frac{(\tau_{s,y} + \tau_{f,y})}{\rho}
```

```{math}
\varphi_4 = -Dm\frac{\rho_f}{\rho}
```

```{math}
\varphi_5 = \zeta D - \frac{3}{\alpha h}||(u,v)^\mathrm{T}||\tan\psi
```

where {math}`\vec{\tau}_s = (\tau_{s,x},\tau_{s,y})^\mathrm{T}` and {math}`\vec{\tau}_f = (\tau_{f,x},\tau_{f,y})^\mathrm{T}` are the shear tractions exerted by the solid phase and fluid phase, respectively, at the base of the flow, {math}`\alpha` is the debris' elastic compressibility, {math}`\psi` is the granular dilatancy angle, and {math}`D` is the granular dilation rate.

For brevity, the variables {math}`\varrho` in and {math}`\varrho` are defined.

```{math}
	\varrho = \frac{(\rho_f + 3\rho)}{4\rho }
```

```{math}
	\zeta = \frac{3}{2\alpha h} + \frac{g_z \rho_f(\rho -\rho_f)}{4\rho}
```

{math}`\alpha` is calculated based on {math}`m`, {math}`h`, and {math}`p_b`

```{math}
\alpha = \frac{a}{m(\rho g_z h - p_b + \sigma_0)}
```

where {math}`a` and {math}`\sigma_0` are constants (typically set to {math}`a=0.01-0.3`, and {math}`\sigma_0=1000` Pa).

An optional definition for {math}`\sigma_0` (see setrun documentation) utilizes, instead of a fixed constant, a diagnostic function:
```{math}
\sigma_0 = \frac{a}{2}\rho_f g_z h \frac{(\rho_s-\rho_f)}{\rho}.
```
Selecting this option for the definition of {math}`\sigma_0` may provide enhanced stability of numerical simulations but is still undergoing beta testing.

The dilation rate,{math}`D`, is postulated to be proportional to the excess fluid pressure (see Iverson and George, 2014)

```{math}
D = -\frac{2k}{h\mu}\left ( p_b - \rho_f g_z h \right )
```
where {math}`k` is the hydraulic permeability (a function of {math}`m`), and {math}`\mu`  is the effective shear viscosity of the pore-fluid.

A form for {math}`k` is
```{math}
k = k_r\exp{\frac{m-m_r}{0.04}}
```
where {math}`k_r` and {math}`m_r` are the reference permeability and solid volume fraction.

The fluid and solid components of basal resistance {math}`\vec{\tau}_f` and {math}`\vec{\tau}_f` are calculated as

```{math}
\vec{\tau}_f = \frac{2\mu\left ( 1-m \right )}{h}(u,v)^\mathrm{T},
```

```{math}
\vec{\tau}_s = \sigma_e \tan(\phi + \psi)\frac{(u,v)^\mathrm{T}}{||(u,v)^\mathrm{T}||}.
```

with {math}`\sigma_e = \rho g_z h - p_b`, the basal effective normal stress.

{math}`\psi` is the granular dilatancy angle:

```{math}
c_1 \tan\psi = m - \frac{m_\mathrm{crit}}{1+\sqrt{N}}
```

:math:`\c_1` is a dimensionless constant (typically 1), {math}`m_\mathrm{crit}` is the critical-state solid fraction, and {math}`N` is a dimensionless variable defined as

```{math}
N = \frac{\mu \dot{\gamma}}{\rho_s \dot{\gamma}^2 \delta^2 + \sigma_e}.
```

where {math}`\dot{\gamma}` is the shear rate and {math}`\delta` is a characteristic length scale associated with grain collisions.


## Additional equations

The elements of D-Claw described in this section are experimental and may change. They are optional and are not enabled by default.

### Segregation

Fully representing the influence of segregation on flow dynamics requires (1) an expression of how flow behavior results in the segregation of different particle species and (2) an expression of how flow behavior is affected by resulting particle species ratios.

To treat the first element, the influence of flow behavior on the segregation of particle species, D-Claw implements a model developed by Gray and Kokelaar (2010). This model considers the segregation of two particle species (species {math}`A` and {math}`B`) in a depth-averaged flow with a linear velocity profile defined as

```{math}
u(z) = (1-\beta)\bar{u} + 2\beta\bar{u}\frac{z}{h}
```
where {math}`\beta` is a constant between 0 (plug flow) and 1 (simple shear) and {math}`\bar{u}` is the depth-averaged velocity.

As the mixture shears, species {math}`A` moves to the surface of the flow, and is preferentially advected by the flow. The lateral transport of {math}`\chi`, the fraction of species {math}`A` is described by the following equation, which, for brevity, shows only x-directed transport.

```{math}
\frac{\partial}{\partial t}(h\chi) +
\frac{\partial}{\partial x}(h \chi u) -
\frac{\partial}{\partial x} \left (\beta h \chi u \left( 1-\chi\right) \right)=0
```

:::{admonition} todo
Reconcile this with
G&K, solver, and plot.py.
is chi the fraction of the faster or slower species?
:::

The representation of the feedback between the value of {math}`\chi` and flow behavior is highly experimental (c.f., Jones et al., 2023).

:::{admonition} todo
Add what is in segeval
:::

### Entrainment

Multiple options for entrainment are present in D-Claw. A user specifies a spatially variable value of the thickness of erodible material that is initially available for entrainment, {math}`h_e`. An entrainment rate, {math}`\frac{\partial h_e}{\partial t}` is calculated if {math}`h_e` is more than the thickness of material that has already been entrained, {math}`\Delta b` (i.e., erodible material remains at a given location). The available erodible material is assumed to have a constant solid volume fraction, {math}`m_e`.

The options for {math}`\frac{\partial h_e}{\partial t}` are:

:::{admonition} todo
Add what is in entrainment
:::

## References

George, D.L., and Iverson, R.M., 2014, A depth-averaged debris-flow model that includes the effects of evolving dilatancy—II. Numerical predictions and experimental tests: Proceedings of the Royal Society of London. Series A, v. 470, no. 2170, p. 20130820,  <https://doi.org/10.1098/rspa.2013.0820>.

Gray, J.M.N.T., and Kokelaar, B.P., 2010, Large particle segregation, transport and accumulation in granular free-surface flows: Journal of Fluid Mechanics, v. 652, p. 105–137 <https://doi.org/10.1017/S002211201000011X>.

Iverson, R.M., and George, D.L., 2014, A depth-averaged debris-flow model that includes the effects of evolving dilatancy—I. Physical basis: Proceedings of the Royal Society of London. Series A, v. 470, no. 2170, p. 20130819, <https://doi.org/10.1098/rspa.2013.0819>.

Jones, R.P., Rengers, F.K., Barnhart, K.R., George, D.L., Staley, D.M., and Kean, J.W., 2023, Simulating Debris Flow and Levee Formation in the 2D Shallow Flow Model D‐Claw: Channelized and Unconfined Flow: Earth and Space Science, v. 10, p. e2022EA002590, <https://doi.org/10.1029/2022EA002590>.

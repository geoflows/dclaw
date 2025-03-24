# setrun.py

## Introduction
The `setrun.py` file is the parameter file used to specify a D-Claw simulation. Many elements of the `setrun.py` are not unique to D-Claw but are inherited from Clawpack. See [this page](https://www.clawpack.org/setrun_geoclaw.html) for a detailed description of these elements.

The data model for the D-Claw specific parts of the setrun are defined by assigning values to attributes of the following five classes:

- [DClawInputData](#dclaw.data.DClawInputData)
- [QinitDClawData](#dclaw.data.QinitDClawData)
- [AuxInitDClawData](#dclaw.data.AuxInitDClawData)
- [PInitDClawInputData](#dclaw.data.PInitDClawInputData)
- [FlowGradesData](#dclaw.data.FlowGradesData)

See the [API documentation of the D-Claw data model](#dclaw.data) for additional details.

## State and auxiliary variable index

D-Claw stores the state variables in an array called **q** and the auxiliary variables in an array called **aux**. There may be a different number of elements of **aux** depending on how D-Claw is specified in the `setrun.py`.

The elements of **q** are defined as in the following table. See the [theory](theory.md) page for further explanation.

:::{list-table}
:widths: 5 20 10 40
:header-rows: 1

*   - Name
    - Description
    - Units
    - Element of **q** (base-1 index)
*   - {math}`h`
    - Flow depth
    - meters
    - 1
*   - {math}`hu`
    - Flow x-momentum: depth times x-directed velocity
    - meters squared per second
    - 2
*   - {math}`hv`
    - Flow y-momentum: depth times y-directed velocity
    - meters squared per second
    - 3
*   - {math}`hm`
    - Flow solid volume: depth times solid-volume fraction
    - meters
    - 4
*   - {math}`p_b`
    - Basal pore-fluid pressure
    - kilograms per meter per time squared
    - 5
*   - {math}`h\chi`
    - Species volume: depth times species A fraction
    - meters
    - 6
*   - {math}`\Delta b`
    - Depth of erosion (material entrained)
    - meters
    - 7
*   - {math}`\eta`
    - Surface elevation, {math}`\eta = b - \Delta b + h`
    - meters
    - 8
:::

The user-facing elements of **aux** are listed in the following table. The numbering depends on the value for ``coordinate_system``, which is a [geoclaw parameter](https://www.clawpack.org/setrun_geoclaw.html#general-geo-parameters) that toggles between horizontal units of meters (`geo_data.coordinate_system = 1`) and latitude-longitude (`geo_data.coordinate_system = 2`).

In addition to the user-facing elements of **aux**, other elements are internally used for the calculation of spatially-variable quantities.

:::{list-table}
:widths: 5 20 10 40
:header-rows: 1

*   - Name
    - Description
    - Units
    - Element of **aux** (base-1 index)
*   - {math}`b`
    - Topobathymetric surface
    - meters
    - 1
*   - {math}`\Theta`
    - Slope angle in x-direction (used only if `bed_normal=1`)
    - degrees
    - 3 (5 for `geo_data.coordinate_system = 2`)
*   - {math}`h_e`
    - Initial thickness of material that can be entrained (used only if `entrainment=1`)
    - meters
    - 7 (9 for `geo_data.coordinate_system = 2`)
:::

## Example `setrun.py`
The example `setrun.py` includes extensive comments intended to support the user in understanding what each part accomplishes.

```{literalinclude} ../../examples/radial_slide/setrun.py
```

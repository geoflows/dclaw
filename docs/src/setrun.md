# setrun.py

The `setrun.py` file is the parameter file used to specify a D-Claw simulation. Many elements of the `setrun.py` are not unique to D-Claw but are inherited from Clawpack. See [this page](https://www.clawpack.org/setrun_geoclaw.html) for a detailed description of these elements.

The data model for the D-Claw specific parts of the setrun are defines by assigning values to attributes of the following five classes:

- [DClawInputData](#dclaw.data.DClawInputData)
- [QinitDClawData](#dclaw.data.QinitDClawData)
- [AuxInitDClawData](#dclaw.data.AuxInitDClawData)
- [PInitDClawInputData](#dclaw.data.PInitDClawInputData)
- [FlowGradesData](#dclaw.data.FlowGradesData)

See the [API documentation of the D-Claw data model](#dclaw.data) for additional details.

The example `setrun.py` includes extensive comments intended to support the user in understanding what each part accomplishes.

```{literalinclude} ../../examples/radial_slide/setrun.py
```

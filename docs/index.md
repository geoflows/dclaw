# D-Claw

D-Claw is a software library for the simulation of dense granular flows over spatially variable topography. D-Claw builds on [clawpack](https://www.clawpack.org/), is written in Fortran, and has a python wrapper. The code is hosted on [code.usgs.gov](https://code.usgs.gov/claw/dclaw) and mirrored on [github](https://github.com/geoflows/dclaw).

D-Claw is supported on Linux/Unix only.

## Using the documentation

The documentation contains multiple elements. Most users should start with the [background](src/background.md) and [theory](src/theory.md) pages, which provide an introduction to what D-Claw is and what equations it solves. Should a user find that they want to install and use D-Claw themselves, they can refer to the [installation instructions](src/installation.md). A few [example applications](src/examples.md) are provided.

Each application of D-Claw is specified through a parameter input file called `setrun.py`. This file includes information about the simulation domain, the initial conditions, and scalar values. Many parameters in setrun are inherited from [clawpack](https://www.clawpack.org/) and a user should expect to read the clawpack documentation to fully understand all options. [This page](src/setrun.md) provides a template `setrun.py`. This template includes extensive comments that describe the nature and format of D-Claw specific parameters.

D-Claw has a data model written in python and a set of utility functions that assist with use of clawpack's [visclaw](https://www.clawpack.org/plotting.html) tools. Refer to the [API](src/autodoc2/index.rst) for the documentation of these python tools. Most users will find the [setrun.py](src/setrun.md) page a more usable description of the D-Claw data model than the API documentation.

:::{warning}
D-Claw is scientific software. A user should expect to gain a high-level of understanding by reading the documentation. For detailed understanding, a user should expect to read the source code. Additionally, some elements of D-Claw are experimental and may change.
:::

D-Claw has been used in [many scientific publications](src/used_by.md). If you use D-Claw in a publication, cite as below and send us a message so we can update the bibliography.

## Citation

If you use D-Claw in a publication, cite it as follows:


``` none
George, D.L., Barnhart, K.R., and Iverson, R.M., 2025,
D-Claw - A library for simulation of granular-fluid flows,
version 1.0.0: U.S. Geological Survey software release,
https://doi.org/10.5066/P13GUXUT.
```

Bibtex:
``` none
@misc{george_2025,
        title = {D-Claw - A library for simulation of
                granular-fluid flows, version 1.0.0},
        author = {George, D.L. and
                  Barnhart, K.R. and
                  Iverson, R.M.},
        url = {https://code.usgs.gov/claw/dclaw},
        year = {2025},
        doi = {10.5066/P13GUXUT}
    }
```


```{toctree}
---
caption: Documentation
hidden: True
---
src/background.md
src/theory.md
src/installation.md
src/examples.md
src/tests.rst
src/setrun.md
src/used_by.md
src/version_control.md
src/build-docs.md
src/autodoc2/index.rst
```

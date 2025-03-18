(install-dclaw)=
# Install D-Claw

Installing D-Claw requires getting the source code for Clawpack and D-Claw, having a fortran compiler, installing a few python packages, and setting environment variables.

## Dependencies

### Platform

D-Claw is supported on Linux/Unix (including Mac OSX). Windows users should expect to use Windows Subsystem for Linux.

### Python libraries

The python wrapper used by D-Claw relies the following python libraries. Install the following in the environment used to run D-Claw.

- [numpy](https://numpy.org/doc/stable/index.html)
- [matplotlib](https://matplotlib.org/)

As a user, you may find other python libraries are helpful in preparing input files or analyzing results.

### Fortran compiler

You need a fortran compiler that works on your machine. On Mac OSX, we recommend using [Homebrew](https://brew.sh/).

```bash
brew install gcc
```

On Linux or WSL this may be accomplished using [apt-get](https://help.ubuntu.com/community/AptGet/Howto) or similar package management tools.

```bash
apt-get gcc
```

## Clone repositories

First, clone the source code. You likely want to do this in a place you store source code.

``` bash
git clone https://github.com/clawpack/clawpack.git
cd clawpack
git submodule init
git submodule update
source pull_all.sh
git clone https://code.usgs.gov/claw/dclaw.git
```
This will place the D-Claw source code in a subfolder (`dclaw/`) of the parent Clawpack directory (`clawpack/`) along with other clawpack submodules (for example, GeoClaw, AMRClaw etc).
## Set environment variables

Running D-Claw, or any Clawpack application requires that the following two environment variables be set.

```
CLAW=/path/to/top/level/clawpack/directory
PYTHONPATH=/path/to/top/level/clawpack/directory
```

If you typically have other elements of your python path, the `CLAW` variable may be appended.

:::{tip}
These environment variables **must** be set for D-Claw to run. Consider adding them to your shell profile (`~/.bashrc` or similar).
:::

## Compile D-Claw

Navigate to one of the example directories (`dclaw/examples/<example-name>`) and compile D-Claw by executing the following:

```bash
make new
```

Execute `make help` to identify all possible commands within the common Clawpack makefile.

```bash
make help
   "make .objs"    to compile object files
   "make .exe"     to create executable
   "make .data"    to create data files using setrun.py
   "make .output"  to run code
   "make output"   to run code with no dependency checking
   "make .plots"   to produce plots
   "make plots"    to produce plots with no dependency checking
   "make .htmls"   to produce html versions of files
   "make .program" to produce single program file
   "make new"      to remove all objs and then make .exe
   "make clean"    to clean up compilation and html files
   "make clobber"  to also clean up output and plot files
   "make help"     to print this message
```

In addition, the standard Clawpack make commands listed above, the standard D-Claw Makefile defines two additional commands, `make input` and `make postprocess`.

## Run an example application

Assuming that compilation has executed successfully, you should have a file called `xdclaw` in the example directory. To run a typical application, you may run three to five commands.

**make input**

The first command executes a file `setinput.py`. This file may not exist in all applications. It is commonly used when preprocessing steps are used to set up initial conditions.

```bash
make input
```

**make .data**

The second command uses the `setrun.py` file to generate required D-Claw input files that end in `.data`. This file is a standard Clawpack input file. Detailed information about its contents may be found [here](https://www.clawpack.org/setrun_geoclaw.html) and in the [D-Claw specific](setrun.md) setrun page.

This step includes checking that any topography or material thickness files needed to run the simulation exist.  

```bash
make .data
```

**make .output**

The third command uses the `.data` files and conducts the D-Claw simulation, generating output (located in the `_output` directory).

```bash
make .output
```

**make .plots**

The fourth command uses a `setplot.py` file to generate plots based on D-Claw output. Similar to the `setrun.py` file, the `setplot.py` file is a standard Clawpack file. It specifies how to plot using the [visclaw tools](https://www.clawpack.org/plotting.html).

```bash
make .plots
```


**make postprocess**

The final command executes the file `setpostprocess.py` and conducts postprocessing analysis with specified within this file. Not all applications will use this command.

```bash
make postprocess
```

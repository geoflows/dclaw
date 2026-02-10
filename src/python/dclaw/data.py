from __future__ import absolute_import, print_function

import os

import clawpack.clawutil.data


class DClawInputData(clawpack.clawutil.data.ClawData):
    r"""Data object describing D-Claw parameters

    Within the ``setrun.py`` a ``DClawInputData`` class is initialized and
    values are assigned to attributes. If a value is not assigned, the default
    value is used.

    .. code-block::

        # First, the rundata object is initialized.
        from clawpack.clawutil import data
        assert claw_pkg.lower() == 'dclaw',  "Expected claw_pkg = 'dclaw'"
        num_dim = 2
        rundata = data.ClawRunData(claw_pkg, num_dim)

        # For a D-Claw model, rundata will have an attribute
        # dclaw_data. This is an instance of the ``DClawInputData``
        # class.
        dclaw_data = rundata.dclaw_data

        # Assign values to dclaw_data attributes to select
        # not-default values.
        dclaw_data.rho_f = 1000.0 # for example, for fluid
        dclaw_data.rho_s = 2700.0 # and solid densities.


    The attributes fall into a few categories. First are the attributes
    that pertain to parameter values in the core set of equations. If
    the parameter is defined in the theory section, a symbol is listed.

    .. list-table::
       :widths: 10 5 10 10 10 10 10
       :header-rows: 1

       * - Attribute Name
         - Symbol
         - Description
         - Default Value
         - Typical Range
         - Type
         - Units
       * - ``rho_s``
         - :math:`\rho_s`
         - Solid density
         - 2700.0
         -
         - float
         - kilograms per cubic meter
       * - ``rho_f``
         - :math:`\rho_s`
         - Fluid density
         - 1000.0
         -
         - float
         - kilograms per cubic meter
       * - ``m_crit``
         - :math:`m_\mathrm{crit}`
         - Critical solid volume fraction
         - 0.62
         -
         - float
         - unitless
       * - ``m0``
         -
         - Initial solid volume fraction. This may be used to set a
           spatially uniform value. Alternatively, use
           :py:class:`dclaw.data.QinitDClawData` to set a spatially
           variable value.
         - 0.52
         -
         - float
         - unitless
       * - ``mref``
         - :math:`m_r`
         - Reference solid volume fraction
         - 0.60
         -
         - float
         - unitless
       * - ``kref``
         - :math:`k_r`
         - Reference permeability
         - 0.0001
         -
         - float
         - square meters
       * - ``phi``
         - :math:`\phi`
         - Friction angle
         - 40.0
         -
         - float
         - degrees
       * - ``delta``
         - :math:`\delta`
         - Characteristic length scale associated with grain collisions
         - 0.01
         -
         - float
         - meters
       * - ``mu``
         - :math:`\mu`
         - Effective shear viscosity of the pore-fluid
         - 0.001
         -
         - float
         - Pa-s
       * - ``alpha_c``
         - :math:`a`
         - Debris compressibility coefficient
         - 0.01
         - 0.01-0.3
         - float
         - unitless
       * - ``c1``
         - :math:`c_1`
         - Granular dilation coefficient
         - 1
         -
         - float
         - unitless
       * - ``sigma_0``
         - :math:`\sigma_0`
         - Reference stress
         - 1000
         -
         - float
         - Pa


    Optionally, D-Claw may use a depth-dependent Manning friction value. The value provided to
    ``geo_data.manning_coefficient`` is used for flows greater than 6 cm. For flows less
    than 4 cm the Manning coefficient will be the value specified by ``manning_max`` if
    ``dd_manning == True`` (a smooth transition is used between the two values for flow depth
    between 4 and 6 cm.

    .. list-table::
       :widths: 10 30 10 10
       :header-rows: 1

       * - Attribute Name
         - Description
         - Default Value
         - Type
       * - ``dd_manning``
         - Flag to indicate whether depth dependent Manning coefficient is used (if
           ``True`` then it is used).
         - ``False``
         - bool
       * - ``manning_max``
         - Maximum Manning coeficient used when :math:`h` approaches the dry tolerance.
         - 0.06
         - float


    Two parameters control the numerical method used for integrating the source
    term (``src2method``) and the equation used to calculate :math:`\alpha`,
    the debris compressibility (``alphamethod``). All combinations except for
    ``src2method=alphamethod=0`` should be considered experimental.


    .. list-table::
       :widths: 10 30 10 10
       :header-rows: 1

       * - Attribute Name
         - Description
         - Default Value
         - Type
       * - ``src2method``
         - Method used to integrate the source term.
         - 0
         - int
       * - ``alphamethod``
         - Method used to calculate :math:`\alpha`, the debris
           compressibility.
         - 0
         - int


    .. list-table::
       :widths: 10 30
       :header-rows: 1

       * - ``src2method``
         - Description
       * - -1
         - Ignore :math:`m` and :math:`p` coevolution. This option is experimental and may
           used for be for shallow water with friction and advection of :math:`m`.
       * - 0
         - Traditional method for source term integration.
       * - 1
         - Experimental method intermediate between methods 0 and 2.
       * - 2
         - New, experimental, method for source term integration. Requires
           ``alphamethod==1``.



    .. list-table::
       :widths: 10 30
       :header-rows: 1

       * - ``alphamethod``
         - Description
       * - 0
         - Traditional method, in which :math:`\sigma_0` is a constant as
           specified in the ``setrun.py``.
       * - 1
         - New method, in which :math:`\sigma_0` is dynamically determined
           as :math:`\sigma_0 = 0.5 \alpha (\rho_s-\rho_f)g_z h/\rho`.


    The following parameters are used to control whether and how bed
    normal coordinates, segregation, and entrainment are implemented.

    .. list-table::
       :widths: 10 30 10 10
       :header-rows: 1

       * - Attribute Name
         - Description
         - Default Value
         - Type
       * - ``bed_normal``
         - Whether to use bed normal coordinates in the x-direction
           (1) or not (0). Should ``bed_normal=1``, a value for the
           x-directed slope, :math:`\theta` should be set.
         - 0
         - int
       * - ``theta_input``
         - Constant value to use for :math:`\theta`. Alternatively, a
           spatially variable :math:`\theta` may be set with
           :py:class:`dclaw.data.QinitDClawData`.
         - 0.0
         - float
       * - ``entrainment``
         - Whether to use entrainment (1) or not (0).
         - 0
         - int
       * - ``entrainment_method``
         - Entrainment method used. See the :ref:`entrainment section <entrainment-theory>` on the theory page for explaination.
         - 1
         - int
       * - ``entrainment_rate``
         - A coefficient for use of entrainment. See documentation for entraiment method used.
         - 0.2
         - float
       * - ``me``
         - Solid volume fraction of entrainable material.
         - 0.65
         - float
       * - ``segregation``
         - Whether to use segregation (1) or not (0).
         - 0
         - int
       * - ``beta_seg``
         - The value of :math:`\beta`. When :math:`\beta>0`, segregation is
           active.
         - 0.0
         - float
       * - ``chi_init_val``
         - The initial value of :math:`\chi`. Only used if :math:`\beta>0`.
         - 0.5
         - float

    Finally, the following parameters control whether the simulation will halt
    when the momentum on a specific adaptive mesh refinement level reaches zero
    and whether to use curvature terms.

    .. list-table::
       :widths: 10 30 10 10
       :header-rows: 1

       * - Attribute Name
         - Description
         - Default Value
         - Type
       * - ``mom_autostop``
         - Whether to halt the simulation when the momentum on ``momlevel``
           is equal to zero.
         - False
         - bool
       * - ``momlevel``
         - The level to consider for determining if flow has no momentum.
         - False
         - bool
       * - ``curvature``
         - Whether use curvature terms (0=No, 1=Yes).
         - 0
         - int

    """

    def __init__(self):
        """Initialize a `DClawInputData`"""
        super(DClawInputData, self).__init__()

        # Set default values:

        self.add_attribute("rho_s", 2700.0)
        self.add_attribute("rho_f", 1000.0)
        self.add_attribute("m_crit", 0.62)
        self.add_attribute("m0", 0.52)
        self.add_attribute("mref", 0.60)
        self.add_attribute("kref", 0.0001)
        self.add_attribute("phi", 40.0)
        self.add_attribute("delta", 0.01)
        self.add_attribute("mu", 0.001)
        self.add_attribute("alpha_c", 0.01)
        self.add_attribute("c1", 1.0)
        self.add_attribute("sigma_0", 1.0e3)

        self.add_attribute("dd_manning", False)
        self.add_attribute("manning_max", 0.06)

        self.add_attribute("src2method", 0)
        self.add_attribute("alphamethod", 0)

        self.add_attribute("bed_normal", 0)
        self.add_attribute("theta_input", 0.0)

        self.add_attribute("entrainment", 0)
        self.add_attribute("entrainment_method", 1)
        self.add_attribute("entrainment_rate", 0.2)
        self.add_attribute("me", 0.65)

        self.add_attribute("segregation", 0)
        self.add_attribute("beta_seg", 0.0)
        self.add_attribute("chi0", 0.5)
        self.add_attribute("chie", 0.5)

        self.add_attribute("mom_autostop", False)
        self.add_attribute("momlevel", 1)
        self.add_attribute("curvature", 0)

    def write(self, out_file="dclaw.data", data_source="setrun.py"):
        """Write the contents of ``DClawInputData`` to a file."""

        # check for any parameter conflicts.
        if self.bed_normal == 1:
            raise ValueError(
                "bed_normal=1 not currently supported because of dx dy not accessible in riemann solver"
            )

        if self.src2method == 2 and self.alphamethod < 1:
            raise ValueError(
                "D-Claw parameter conflict: if src2method == 2 then alphamethod must = 1"
            )

        self.open_data_file(out_file, data_source)

        self.data_write("rho_s", description="solid grain density (kg/m^3)")
        self.data_write("rho_f", description="pore-fluid density  (kg/m^3)")
        self.data_write("m_crit", description="critical state value of m (#)")
        self.data_write("m0", description="initial solid volume fraction (#)")
        self.data_write("mref", description="reference solid volume fraction (#)")
        self.data_write(
            "kref",
            description=" reference permeability",
        )
        self.data_write("phi", description="basal friction angle (degrees)")
        self.data_write("delta", description="characteristic grain diameter (m)")
        self.data_write("mu", description="viscosity of pore-fluid (Pa-s)")
        self.data_write("alpha_c", description="debris compressibility constant (#)")
        self.data_write("c1", description="granular dilatency constant (#)")
        self.data_write(
            "sigma_0", description="baseline stress for definition of compressibility"
        )
        self.data_write(
            "dd_manning", description="Depth dependent Manning flag (False = Not used)"
        )
        self.data_write("manning_max", description="Maximum manning coefficient")

        self.data_write(
            "src2method", description="-1=swe, 0=orig,1=intermediate,2=new"
        )  # DIG: update text
        self.data_write("alphamethod", description="0=,1=,2=")  # DIG: update text

        self.data_write(
            "bed_normal",
            description="use of bed normal coordinates (0=false, 1=true). bed_normal = 1 requires theta in aux for slope in one direction",
        )
        self.data_write("theta_input", description="slope angle (degrees)")

        self.data_write(
            "entrainment", description="flag for entrainment, 0 = no entrainment"
        )
        self.data_write(
            "entrainment_method", description="flag for entrainment method"
        )  # DIG: update text
        self.data_write(
            "entrainment_rate", description="rate of entrainment parameter 0-1"
        )
        self.data_write(
            "me", description="Solid volume fraction of entrainable material"
        )

        self.data_write(
            "segregation",
            description="flag for segregation, 0 = no segregation",
        )
        self.data_write(
            "beta_seg",
            description="coefficient of segregation velocity profile",
        )
        self.data_write(
            "chi0",
            description="initial fraction of species A (#). Between 0-1.",
        )
        self.data_write(
            "chie",
            description="fraction of species A of entrainable material (#). Between 0-1.",
        )

        self.data_write(
            "mom_autostop",
            description="flag for momentum autostop False = no autostop, True = autostop",
        )  # currently only works with ascii output
        self.data_write(
            "momlevel",
            description="level to do momentum calculation IF mom_autostop==True",
        )
        self.data_write(
            "curvature",
            description="flag for curvature correction 0 = not used, 1 = used",
        )

        self.close_data_file()


class QinitDClawData(clawpack.clawutil.data.ClawData):
    r"""
    Data object describing initialization of D-Claw state variables stored in **q**.

    Within the ``setrun.py`` the ``QinitDClawData`` class is initialized and
    lists are appended to the attribute ``qinitfiles``.

    To initialize **q**, the user provides a file indicating the spatially
    variable values of that element of **q**. Multiple file types are permitted
    and are documented in
    `the clawpack documentation <https://www.clawpack.org/topo.html>`__ .

    Most users will use *topotype* = 3, which is equivalent to an
    `ESRI ascii raster <https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm>`__
    file type.

    One list is appended for each file used to initialize **q**. Multiple files
    may be used to initialize each element of **q**. If files for the same
    element of **q** are provided then values from the higher resolution
    file is used. If multiple files for the same element of **q** are provided
    that overlap, the value from the last entry into ``qinitfiles`` for that
    element of **q** is used.

    For each provided file, the user will append a list containing the following
    three elements (in this order) to the attribute ``qinitfiles``.

    * ``qinitftype``: integer indicating the topotype of the provided file (ie: 1, 2 or 3)
    * ``iqinit``: integer indicating the element of **q** associated with this file.
    * ``fname``: string indicating the path to the file

    All elements of **q** except for :math:`\Delta b` are accessible through
    ``setqinit``. Some may not typically be set and some are set using the
    value divided by :math:`h`.

    .. list-table::
       :widths: 10 30 10 10 30
       :header-rows: 1

       *   - Name
           - Description
           - Units
           - Element of **q**
           - Notes
       *   - :math:`h`
           - Flow depth
           - meters
           - 1
           -
       *   - :math:`hu`
           - Flow depth times x-directed velocity
           - meters squared per second
           - 2
           - Specified as :math:`u` not :math:`hu`
       *   - :math:`hv`
           - Flow depth times y-directed velocity
           - meters squared per second
           - 3
           - Specified as :math:`v` not :math:`hv`
       *   - :math:`hm`
           - Flow depth times solid volume fraction
           - meters
           - 4
           - Specified as :math:`m` not :math:`hm`
       *   - :math:`p_b`
           - Basal pore pressure
           - kilograms per meter per time squared
           - 5
           - Specified as :math:`p_b/h` not :math:`p_b`
       *   - :math:`h\chi`
           - Depth times species A fraction
           - meters
           - 6
           - Specified as :math:`\chi` not :math:`h\chi`
       *   - :math:`\Delta b`
           - Depth of material entrained
           - meters
           - 7
           - Not set because this value must start at zero.
       *   - :math:`\eta`
           - Surface elevation :math:`\eta=b+h`
           - meters
           - 8
           - If files are provided for both :math:`h` and :math:`\eta`, the file
             provided last is used. If :math:`\eta` is provided and :math:`eta<b`,
             the value for :math:`h` is set to 0.

    The following is an example of providing files to initialize **q** using
    the ``QinitDClawData`` class.

    .. code-block::

        # First, the rundata object is initialized.
        from clawpack.clawutil import data
        assert claw_pkg.lower() == 'dclaw',  "Expected claw_pkg = 'dclaw'"
        num_dim = 2
        rundata = data.ClawRunData(claw_pkg, num_dim)

        # For a D-Claw model, rundata will have an attribute
        # qinitdclaw_data. This is an instance of the
        # QinitDClawData class.
        qinitdclaw_data = rundata.qinitdclaw_data

        # To set input data files for the q array, similar to setting
        # topo files append lists with the following elements:

        #    [qinitftype,iqinit, fname]
        #
        # where
        #   qinitftype = the file type
        #   iqinit     = the element of q
        #   fname      = the path to the file

        qinitdclaw_data.qinitfiles.append([3, 1, "h.tt3"])

    """

    def __init__(self):
        """Initialize a ``QinitDClawData`` instance."""
        super(QinitDClawData, self).__init__()
        self.add_attribute("qinitfiles", [])
        self.add_attribute("nqinits", None)

    def write(self, data_source="setrun.py", out_file="qinit_dclaw.data"):
        """Write the content of a ``QinitDClawData`` to a file."""

        self.open_data_file(out_file, data_source)
        self.nqinits = len(self.qinitfiles)

        self.data_write("nqinits", description="nqinits")
        self._out_file.write("\n")

        for tfile in self.qinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                raise ValueError(f"*** Error: file not valid string {tfile[-1]}")

            if len(fname) > 150:
                raise ValueError(
                    f"*** Error: file name too long (must be <150)  {tfile[-1]}"
                )

            if not os.path.exists(tfile[-1]):
                raise ValueError(f"*** Error: file not found: {tfile[-1]}")

            self._out_file.write("\n%s  \n" % fname)
            self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))

        self.close_data_file()


class AuxInitDClawData(clawpack.clawutil.data.ClawData):
    r"""
    Data object describing initialization of D-Claw auxiliary variables stored in **aux**.

    Within the ``setrun.py`` the ``AuxInitDClawData`` class is initialized and
    lists are appended to the attribute ``auxinitfiles``.

    To initialize **aux**, the user provides a file indicating the spatially
    variable values of that element of **aux**. Multiple file types are permitted
    and are documented in
    `the clawpack documentation <https://www.clawpack.org/topo.html>`__ .

    Most users will use *topotype* = 3, which is equivalent to an
    `ESRI ascii raster <https://desktop.arcgis.com/en/arcmap/latest/manage-data/raster-and-images/esri-ascii-raster-format.htm>`__
    file type.

    One list is appended for each file used to initialize **aux**. Multiple files
    may be used to initialize each element of **aux**. If files for the same
    element of **aux** are provided then values from the higher resolution
    file is used. If multiple files for the same element of **aux** are provided
    that overlap, the value from the last entry into ``auxinitfiles`` for that
    element of **aux** is used.

    For each provided file, the user will append a list containing the following
    three elements (in this order) to the attribute ``auxinitfiles``.

    * ``auxinitftype``: integer indicating the topotype of the provided file (ie: 1, 2 or 3)
    * ``iauxinit``: integer indicating the element of **aux** associated with this file.
    * ``fname``: string indicating the path to the file

    The elements of aux that are possible to set are:

    .. list-table::
       :widths: 10 30 10 10 30
       :header-rows: 1

       *   - Name
           - Description
           - Units
           - Element of **aux** if ``coordinate_system=1``
           - Element of **aux** if ``coordinate_system=2``
       *   - :math:`\theta`
           - Slope in the x-direction
           - degrees
           - 3
           - 5
       *   - :math:`h_e`
           - Thickness of entrainable material
           - meters
           - 7
           - 9

    The topobathymetric surface :math:`b` is specified using the attribute
    ``rundata.topo_data.topofiles`` (see
    `the clawpack documentation <https://www.clawpack.org/setrun_geoclaw.html#topography-data-file-parameters>`__
    for details)

    The following is an example of providing files to initialize **aux** using
    the ``AuxInitDClawData`` class.

    .. code-block::

        # First, the rundata object is initialized.
        from clawpack.clawutil import data
        assert claw_pkg.lower() == 'dclaw',  "Expected claw_pkg = 'dclaw'"
        num_dim = 2
        rundata = data.ClawRunData(claw_pkg, num_dim)

        # For a D-Claw model, rundata will have an attribute
        # auxinitdclaw_data. This is an instance of the
        # AuxInitDClawData class.
        auxinitdclaw_data = rundata.auxinitdclaw_data

        # To set input data files for the aux array, similar to setting
        # topo files append lists with the following elements:

        #    [auxinitftype, iauxinit, fname]
        #
        # where
        #   auxinitftype = the file type
        #   iauxinit     = the element of aux
        #   fname        = the path to the file

        auxinitdclaw_data.auxinitfiles.append([3, 3, "theta.tt3"])

    """

    def __init__(self):
        """Initialize an ``AuxInitDClawData`` object."""

        super(AuxInitDClawData, self).__init__()
        self.add_attribute("auxinitfiles", [])
        self.add_attribute("nauxinits", None)

    def write(self, data_source="setrun.py", out_file="auxinit_dclaw.data"):
        """Write the contents of an ``AuxInitDClawData`` object to a file."""

        self.open_data_file(out_file, data_source)
        self.nauxinits = len(self.auxinitfiles)

        self.data_write("nauxinits", description="nauxinits")
        self._out_file.write("\n")

        for tfile in self.auxinitfiles:
            try:
                fname = "'%s'" % os.path.abspath(tfile[-1])
            except:
                raise ValueError(f"*** Error: file not valid string {tfile[-1]}")

            if len(fname) > 150:
                raise ValueError(
                    f"*** Error: file name too long (must be <150)  {tfile[-1]}"
                )

            if not os.path.exists(tfile[-1]):
                raise ValueError(f"*** Error: file not found: {tfile[-1]}")

            self._out_file.write("\n%s  \n" % fname)
            self._out_file.write("%3i %3i \n" % tuple(tfile[:-1]))

        self.close_data_file()


class PInitDClawInputData(clawpack.clawutil.data.ClawData):
    r"""
    Data object describing the initialization of the D-Claw pressure field.

    The attributes of ``PInitDClawInputData`` control how the initial values
    of the basal pressure field, :math:`p_b` are initialized.

    Within the setrun, the ``rundata`` object will have an attribute
    ``rundata.pinitdclaw_data``. This object has one attribute, ``init_ptype``
    that controls the pressure initialization.

    .. list-table::
       :widths: 10 30
       :header-rows: 1

       *   - ``init_ptype``
           - Description
       *   - -1
           - :math:`p_b=0`. If a user-defined spatially variable file is provided
             through :py:class:`dclaw.data.QinitDClawData`, use this option.
       *   - 0
           - Hydrostatic pressure, :math:`p_b=\rho_f g_z h`
       *   - 1
           - Based on the values of :math:`h` and :math:`b`, a static force
             balance is conducted at each mesh cell to determine the failure
             pressure, :math:`p_f`, at which the granular material is not stable.
             At each cell, the failure pressure ratio,
             :math:`R_f=p_f/(\rho_f g_z h)` is calculated.
             The minimum value of :math:`R_f` is used to calculate the initial
             basal pressure, :math:`p_b = \mathrm{min} (R_f) \rho_f g_z h`.
       *   - 2
           - Same as case 1, except that the average failure pressure is used.
             :math:`p_b = \mathrm{mean} (R_f) \rho_f g_z h`.
       *   - 3
           - Hydrostatic pressure times a specified coefficient, :math:`C_p`
             given by ``init_pratio``, :math:`p_b=C_p \rho_f g_z h`.
             Cannot exceed lithostatic pressure, :math:`p_b=((m\rho_s + (1-m)\rho_f) g_z h`
             such that :math:`C_p\leq((m\rho_s + (1-m)\rho_f)/\rho_f``.

    """

    def __init__(self):
        """Initialize a ``PInitDClawInputData`` object"""
        super(PInitDClawInputData, self).__init__()

        # Set default values:
        self.add_attribute("init_ptype", 0)
        self.add_attribute("init_pratio", 1)


    def write(self, out_file="pinit_dclaw.data", data_source="setrun.py"):
        """Write the contents of a ``PInitDClawInputData`` to a file."""
        self.open_data_file(out_file, data_source)

        # open file and write a warning header:
        self.data_write(
            "init_ptype",
            description="-1 = zero pressure or user defined files in qinit, 0 = hydrostatic, 1,2 = failure pressure (1=min, 2=avg), 3=hydrostatic*init_pratio",
        )
        self.data_write(
            "init_pratio",
            description="Initialization pressure ratio, used if init_ptype==3.",
        )
        self.close_data_file()


class FlowGradesData(clawpack.clawutil.data.ClawData):
    r"""
    Data object describing the flowgrades used to control refinement.

    In D-Claw refinement will be flagged within a mesh cell on an adaptive
    mesh refinement level if the characteristics of **q** meet user-defined
    criteria.

    A flagged cell will refine to the minimum level specified only if refinement
    to this level is permitted based on refinement regions. See the
    `clawpack documentation <https://www.clawpack.org/setrun_geoclaw.html#amr-refinement-region-parameters>`__
    for more details.

    Note that the
    `geoclaw-specific refinement criteria <https://www.clawpack.org/setrun_geoclaw.html#additional-amr-parameters>`__
    ``wave_tolerance`` and ``speed_tolerance`` are not supported in D-Claw.

    .. code-block::

        # First, the rundata object is initialized.
        from clawpack.clawutil import data
        assert claw_pkg.lower() == 'dclaw',  "Expected claw_pkg = 'dclaw'"
        num_dim = 2
        rundata = data.ClawRunData(claw_pkg, num_dim)

        # For a D-Claw model, rundata will have an attribute
        # flowgrades_data. This is an instance of the
        # FlowGradesData class.
        flowgrades_data = rundata.flowgrades_data

        flowgrades_data.flowgrades = []

        # for using flowgrades for refinement append lines of the form
        # [flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel]
        # where:
        #  flowgradevalue:
        #      floating point relevant flowgrade value for following measure:
        #  flowgradevariable:
        #      1 = depth
        #      2 = velocity, sqrt(u**2+v**2)
        #      3 = sign(depth)*(depth+topo) (0 at sealevel or dry land)
        #  flowgradetype:
        #      1 = norm(flowgradevariable)
        #      2 = norm(grad(flowgradevariable))
        #  flowgrademinlevel:
        #      refine to at least this level if flowgradevalue is exceeded

        flowgrades_data.flowgrades.append([1.0e-6, 2, 1, 3])

        # multiple flowgrades may be used by appending multiple lists.
        # refinement will occur if *any* of the criteria are met.

    Traditionally, refinement in clawpack adaptive mesh codes only evaluates
    whether refinement should be flagged up to one minus the maximum level.
    In confined topography, it is possible to have a fine grid that meets the
    refinement criteria but the coarse grid does not. This may result in
    unexpected derefinement.

    To address this issue, use

    .. code-block::

        flowgradesdata.keep_fine = True
        flowgrades.fine_level = 3

    If the finest grid level meets the flowgrades refinement criteria,
    a coarser grid containing it that does not itself meet the refinement
    criteria will remain refined.

    """

    def __init__(self):
        """Initialize a ``FlowGradesData`` object"""

        super(FlowGradesData, self).__init__()
        self.add_attribute("flowgrades", [])
        self.add_attribute("keep_fine", False)
        self.add_attribute("fine_level", -1)
        self.add_attribute("nflowgrades", None)

    def write(self, out_file="flowgrades.data", data_source="setrun.py"):
        """Write the contents of a ``FlowGradesData`` object to a file"""
        self.nflowgrades = len(self.flowgrades)
        self.open_data_file(out_file, data_source)
        self.data_write("nflowgrades", description="nflowgrades")
        self._out_file.write("\n")
        for flowgrade in self.flowgrades:
            self._out_file.write(4 * "%g  " % tuple(flowgrade) + "\n")
        self.data_write("keep_fine", description="keep_fine")
        self.data_write("fine_level", description="fine_level")

        self.close_data_file()

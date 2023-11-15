#!/usr/bin/env python

"""

Classes representing parameters for D-Claw runs

:Classes:

 - DClawInputData
 - QinitDClawData
 - AuxInitDClawData
 - PInitDClawInputData
 - FlowGradesData


:Constants:

 - Rearth - Radius of earth in meters
 - DEG2RAD factor to convert degrees to radians
 - RAD2DEG factor to convert radians to degrees
 - LAT2METER factor to convert degrees in latitude to meters
"""

from __future__ import absolute_import
from __future__ import print_function
import os
import numpy
import clawpack.clawutil.data
import warnings


# Radius of earth in meters.
# For consistency, should always use this value when needed, e.g.
# in setrun.py or topotools:
Rearth = 6367.5e3  # average of polar and equatorial radii

DEG2RAD = numpy.pi / 180.0
RAD2DEG = 180.0 / numpy.pi
LAT2METER = Rearth * DEG2RAD


class DClawInputData(clawpack.clawutil.data.ClawData):
    r"""
    D-Claw data object

    See description text for explaination.
    """

    def __init__(self):
        super(DClawInputData, self).__init__()

        # Set default values:
        self.add_attribute("rho_s", 2700.0)
        self.add_attribute("rho_f", 1000.0)
        self.add_attribute("phi_bed", 40.0)
        self.add_attribute("theta_input", 0.0)
        self.add_attribute("delta", 0.01)
        self.add_attribute("kappita", 0.0001)
        self.add_attribute("mu", 0.001)
        self.add_attribute("alpha_c", 1.0) # DIG: is this the right default value?
        self.add_attribute("m_crit", 0.62)
        self.add_attribute("c1", 1.0)
        self.add_attribute("m0", 0.52)
        self.add_attribute("sigma_0", 1.0e3)
        self.add_attribute("alpha_seg", 0.0)
        self.add_attribute("bed_normal", 0)
        self.add_attribute("phi_seg_coeff",0.0) # Not currently used
        self.add_attribute("entrainment", 0)
        self.add_attribute("entrainment_rate", 0.2)

        # DIG: currently unused self.add_attribute("curvature", 1)
        # DIG: currently unused self.add_attribute("momlevel", 1)
        # DIG: currently unused self.add_attribute("mom_autostop", False)

        # DIG: currently unused self.add_attribute("mom_perc", 0.05)
        # DIG: currently unused self.add_attribute("src_ftn_num", 0)
        # DIG: currently unused self.add_attribute("fric_offset_val", 0.0)
        # DIG: currently unused self.add_attribute("fric_star_val", 0.0)
        # DIG: currently unused self.add_attribute("chi_init_val", 0.0)
        # DIG: currently unused self.add_attribute("kappita_diff", 1.0)
        # DIG: currently unused self.add_attribute("outaux", False)


        # DIG: KRB would prefer to add back momlevel & mom autostop. It could go 
        # elsewhere. 


    def write(self,out_file='setdclaw.data',data_source='setrun.py'):
        self.open_data_file(out_file,data_source)

        self.data_write("rho_s", description="solid grain density (kg/m^3)")
        self.data_write("rho_f", description="pore-fluid density  (kg/m^3)")
        self.data_write("phi_bed", description="basal friction angle (degrees)")
        self.data_write("theta_input", description="slope angle (degrees)")
        self.data_write("delta", description= "characteristic grain diameter (m)")
        self.data_write("kappita", description="permeability at m=setdig.m0 (m^2), k0 in G&I eq 2.7 if m0 is 0.6")
        self.data_write("mu", description="viscosity of pore-fluid (Pa-s)")
        self.data_write("alpha_c", description="debris compressibility constant (#)")
        self.data_write("m_crit", description="critical state value of m (#)")
        self.data_write("c1", description="dilation regularization coefficient 1 (#)")
        self.data_write("m0", description="initial solid volume fraction (#)")
        self.data_write("sigma_0", description="baseline stress for definition of compressibility")
        self.data_write("alpha_seg", description="coefficient of segregation velocity profile. When alpha_seg = 0, no segregation occurs")
        self.data_write("bed_normal", description="use of bed normal coordinates (0=false, 1=true). bed_normal = 1 requires theta in aux for slope in one direction")
        self.data_write("phi_seg_coeff", description="adjustment to friction coefficient based on segregation") # not currently used.)
        self.data_write("entrainment", description="flag for entrainment, 0 = no entrainment")
        self.data_write("entrainment_rate", description="rate of entrainment parameter 0-1")
        # DIG: These currently unused: self.data_write("mom_autostop", description= "flag for momentum autostop False = no autostop, True = autostop")
        # DIG: These currently unused: self.data_write("mom_perc", description="percentage of max momentum for autostop, default is 0.05 (5%)")
        # DIG: These currently unused: self.data_write("src_ftn_num", description="number of in-domain sources, if used the file 'sethydrographs.data' is required")
        # DIG: These currently unused: self.data_write("fric_offset_val", description="start/stop friction offset (degrees). if this value is >0, then hysteretic friction is used (Rocha, Johnson, Gray, 2019)")
        # DIG: These currently unused: self.data_write("fric_star_val", description="deep friction offset (degrees). only used when fric_offset_val > 0 (Rocha, Johnson, Gray, 2019)")
        # DIG: These currently unused: self.data_write("chi_init_val", description="initial fraction of species 1, (#). Between 0-1.")
        # DIG: These currently unused: self.data_write("kappita_diff", description="permeability multiplier for different size species. Only used when alpha_seg>0. kappita is used for species1, kappita*kappita_diff used for species2")
        # DIG: These currently unused: self.data_write("outaux", description="flag for writing aux to output F = not written, T = written")
        # DIG: These currently unused: self.data_write("curvature", description="flag for curvature correction 0 = not used, 1 = used")
        # DIG: These currently unused: self.data_write("momlevel", description="level to do momentum calculation IF mom_autostop==True")

        self.close_data_file()

class QinitDClawData(clawpack.clawutil.data.ClawData):
    r"""
    Qinit data class for D-Claw

    To set input data files for the q array, similar to setting 
    topo files append lists with the following elements:
    
        [qinitftype,iqinit, minlev, maxlev, fname]

    where:
    
    - qinitftype: 
        file-type, same as topo files, ie: 1, 2 or 3
    - iqinit: 
        The following values are allowed for iqinit:
            n=1,meqn perturbation of q(i,j,n)
            except for n=7 (b_eroded)
            n=meqn+1: surface elevation eta is defined by 
            the file and results in h=max(eta-b,0)

            TODO, if h and eta are both defined, which wins
    
    - minlev:
    - maxlev: 
    - fname: 

    TODO: will this way of setting flagregions work win 5.x?

    The elements of q accessible through setqinit are:
    - q1, h: 
        depth (most common value to be set this way)
    - q2, hu: 
        depth * x-directed velocity (unlikley to be set 
        in this way)
        If specified in this way, provide u rather than h*u
        If not provided, hu is initialized to 0
    - q3, hv: 
        depth * y-directed velocity (unlikley to be set 
        in this way)
        If specified in this way, provide v rather than h*v
        If not provided, hu is initialized to 0
    - q4, hm: 
        depth * solid volume fraction (provide m to 
        setqinit, not h*m)
        If not specified, hm is initialized as h*m0
    - q5, pb: 
        basal pressure
        Provided as pb/h rather than the absolute value
        of pb
    - q6, hchi: 
        depth * chi, the fraction of species 1
        If specified this way, provide chi rather than h*chi
        TODO: should be chi_init_val rather than 0.5 by default

    - q7, b_eroded: depth of material removed by
      entrainment (cannot be set in this way as it
      must start as zero).
    - q8, eta: h+b 


    """
    def __init__(self):
        super(QinitDClawData,self).__init__()
        self.add_attribute('qinitfiles',[])
        self.add_attribute('nqinits',None)

    def write(self,data_source='setrun.py', out_file='setqinit_dclaw.data'):

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
                raise ValueError(f"*** Error: file name too long (must be <150)  {tfile[-1]}")

            if not os.path.exists(tfile[-1]):
                raise ValueError(f"*** Error: file not found: {tfile[-1]}")
            
            self._out_file.write("\n%s  \n" % fname)
            self._out_file.write("%3i %3i %3i %3i \n" % tuple(tfile[:-1]))

        self.close_data_file()

class AuxInitDClawData(clawpack.clawutil.data.ClawData):
    r"""
    AuxInit data class for D-Claw

    To set input data files for the aux array, similar to setting 
    topo files append lists with the following elements:
    
        [auxinitftype,iauxinit, minlev, maxlev, fname]

    where:
    
    - auxinitftype: 
        file-type, same as topo files, ie: 1, 2 or 3
    - iauxinit: 
        The following values are allowed for iauxinit:
            n=1,maux perturbation of aux(i,j,n)
    - minlev:
    - maxlev: 
    - fname: 

    The elements of aux that are accessible in this way are:

    - aux1, topographic elevation, b
    - aux2 and aux3, capacity arrays for handling latitude and longitude.
    - aux4, phi, friction angle (TODO remove)
    - aux5, depth of erodible material OR theta, slope angle in X direction (TODO change so these do not share)
    - aux6, fs, factor of safety (TODO remove?)
    - aux7,fsphi, (TODO what is this?)
    - aux8, taudir_x, friction direction
    - aux9, taudir_y, friction direction

    """
    def __init__(self):
        super(AuxInitDClawData,self).__init__()
        self.add_attribute('auxinitfiles',[])
        self.add_attribute('nauxinits',None)

    def write(self,data_source='setrun.py', out_file='setauxinit_dclaw.data'):

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
                raise ValueError(f"*** Error: file name too long (must be <150)  {tfile[-1]}")

            if not os.path.exists(tfile[-1]):
                raise ValueError(f"*** Error: file not found: {tfile[-1]}")
            
            self._out_file.write("\n%s  \n" % fname)
            self._out_file.write("%3i %3i %3i %3i \n" % tuple(tfile[:-1]))

        self.close_data_file()



class PInitDClawInputData(clawpack.clawutil.data.ClawData):
    r"""
    D-Claw pressure initialization data object
    """

    def __init__(self):
        super(PInitDClawInputData, self).__init__()

        # Set default values:
        self.add_attribute("init_ptype", 0)
        self.add_attribute("init_pmax_ratio", 1.0)
        self.add_attribute("init_ptf", 1.0)
        self.add_attribute("init_ptf2", 0.0)

    def write(self,out_file='setpinit_dclaw.data',data_source='setrun.py'):
        self.open_data_file(out_file,data_source)   

        # open file and write a warning header:
        self.data_write("init_ptype", description="-1 = zero pressure or user defined files in qinit, 0 = hydrostatic, 1,2 = failure pressure (1=min, 2=avg), 3,4 = rising pressure (3=min, 4=avg)")
        self.data_write("init_pmax_ratio", description="p(init_ptf2)= hydro*init_pmax_ratio: pressure will rise to hydrostatic *init_pmax_ratio")
        self.data_write("init_ptf", description="p(init_ptf) = failure, pressure will rise until t = init_ptf without dilatancy")
        self.data_write("init_ptf2", description="p(init_ptf2)= hydro*init_pmax_ratio, pressure will rise until t = init_ptf2")

        self.close_data_file()

class FlowGradesData(clawpack.clawutil.data.ClawData):
    r"""
    Flowgrades data object.

    For using flowgrades for refinement append lines of the form
    
    [flowgradevalue, flowgradevariable, flowgradetype, flowgrademinlevel]
    
    where:
    
    - flowgradevalue: 
        floating point relevant flowgrade value for following measure:
    - flowgradevariable: 
        1=depth, 
        2= momentum, 
        3 = sign(depth)*(depth+topo) (0 at sealevel or dry land).
    - flowgradetype: 
        1 = norm(flowgradevariable), 
        2 = norm(grad(flowgradevariable))
    - flowgrademinlevel: 
        refine to at least this level if flowgradevalue is exceeded.

    flowgradesdata.flowgrades = []
    flowgradesdata.flowgrades.append([1.0e-6, 2, 1, 1])

    Additionally, to turn on "keep fine" refinement, set

    flowgradesdata.keep_fine = True

    If the finest grid level meets the flowgrades refinement criteria, 
    a coarser grid containing it that does not itself meet the refinement
    criteria will stay refined.
    
    """
    def __init__(self):
        super(FlowGradesData,self).__init__()
        self.add_attribute("flowgrades", [])
        self.add_attribute("keep_fine", False)
        self.add_attribute('nflowgrades',None)


    def write(self,out_file='setflowgrades.data',data_source='setrun.py'):
        self.nflowgrades = len(self.flowgrades)
        self.open_data_file(out_file,data_source)
        self.data_write("nflowgrades", description="nflowgrades")        
        self._out_file.write("\n")
        for flowgrade in self.flowgrades:
            self._out_file.write(4 * "%g  " % tuple(flowgrade) + "\n")
        self.data_write("keep_fine", description="keep_fine")        
        self.close_data_file()

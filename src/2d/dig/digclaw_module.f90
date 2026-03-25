! ============================================================================
!  Module for initialization D-Claw parameters and core functions
! ============================================================================


module digclaw_module

   use geoclaw_module, only: dry_tolerance, grav, pi, coordinate_system

   implicit none

    ! ========================================================================
    ! General digclaw parameters
    ! ========================================================================
    double precision :: rho_s,rho_f,m_crit,m0,mref,kref,phi,delta,mu,alpha_c
    double precision :: c1,sigma_0,kappa
    double precision :: theta_input,entrainment_rate,me,beta_seg,chi0,chie
    double precision :: manning_max
    logical :: dd_manning

    integer :: src2method,alphamethod,bed_normal,entrainment,entrainment_method
    integer :: riemann_method
    logical :: rarecorrectortest
    integer :: segregation,curvature,init_ptype
    double precision :: init_pmin_ratio,init_pratio
    double precision :: grad_eta_max,cohesion_max,grad_eta_ave,eta_cell_count

    ! momentum autostop
    logical :: mom_autostop
    integer :: momlevel
    logical :: amidoneyet = .FALSE.

    ! indicies of q
    integer ::  i_h
    integer ::  i_hu
    integer ::  i_hv
    integer ::  i_hm
    integer ::  i_pb
    integer ::  i_hchi
    integer ::  i_bdif

    ! flag for boussinesq
    integer ::  i_bouss

    ! indicies of aux
    integer ::  i_dig
    integer ::  i_phi
    integer ::  i_theta
    integer ::  i_fsphi
    integer ::  i_taudir_x
    integer ::  i_taudir_y
    integer ::  i_ent

    integer, parameter ::  DIG_PARM_UNIT = 78
    integer, parameter ::  PINIT_PARM_UNIT = 79


contains

    ! ========================================================================
    !  set_dig(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_dig(naux, fname)

         implicit none

         ! Input
         integer, intent(in) :: naux
         character*25, intent(in), optional :: fname

         ! Locals
         double precision :: deg2rad
         integer, parameter :: iunit = 127
         character*25 :: file_name
         logical :: found_file


         deg2rad = pi/180.d0

         !Set i_dig based on coordinate system
         if (coordinate_system == 2) then
            i_dig = 4
         else
            i_dig = 2
         endif

         ! flag for whether boussinesq equations are used
         ! (not yet implemented)
         ! if bouss is active (based on reading in a data file)
         ! this will be 2
         i_bouss = 0

         ! hard code q indicies.
         i_h=1
         i_hu=2
         i_hv=3
         i_hm=4 + i_bouss
         i_pb=5 + i_bouss
         i_hchi=6 + i_bouss
         i_bdif=7 + i_bouss

         ! set aux index values based on coordinate system
         i_phi    = i_dig
         i_theta  = i_dig + 1
         i_fsphi  = i_dig + 2
         i_taudir_x = i_dig + 3
         i_taudir_y = i_dig + 4
         i_ent = i_dig + 5

         !kappa: earth pressure coefficient
         kappa = 1.d0

         ! Read user parameters from setgeo.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'dclaw.data'
         endif
         inquire(file=file_name,exist=found_file)
         if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
         endif

         call opendatafile(iunit, file_name)

         read(iunit,*) rho_s
         read(iunit,*) rho_f
         read(iunit,*) m_crit
         read(iunit,*) m0
         read(iunit,*) mref
         read(iunit,*) kref
         read(iunit,*) phi
         phi = deg2rad*phi
         read(iunit,*) delta
         read(iunit,*) mu
         read(iunit,*) alpha_c
         read(iunit,*) c1
         read(iunit,*) sigma_0
         read(iunit,*) dd_manning
         read(iunit,*) manning_max

         read(iunit,*) src2method
         read(iunit,*) alphamethod
         read(iunit,*) rarecorrectortest
         read(iunit,*) riemann_method

         read(iunit,*) bed_normal
         read(iunit,*) theta_input
         theta_input = deg2rad*theta_input

         read(iunit,*) entrainment
         read(iunit,*) entrainment_method
         read(iunit,*) entrainment_rate
         read(iunit,*) me

         read(iunit,*) segregation
         read(iunit,*) beta_seg
         read(iunit,*) chi0
         read(iunit,*) chie

         read(iunit,*) mom_autostop
         read(iunit,*) momlevel
         read(iunit,*) curvature
         close(iunit)

         ! test that naux is large enough
         if (i_dig + 4 + entrainment > naux) then
            write(*,*)
            write(*,*) "*******************************************************"
            write(*,*) "**************** AUX ERROR"
            write(*,*) "**************** Number of aux variables set for D-Claw"
            write(*,*) "**************** inconsistent with coordinate_system"
            write(*,*) "**************** and entrainment."
            write(*,*) "**************** Setrun num_aux    = ", naux
            write(*,*) "**************** Aux required      = ", i_dig + 4 + entrainment
            write(*,*) "**************** coordinate_system = ", coordinate_system
            write(*,*) "**************** entrainment       = ", entrainment
            write(*,*) "*******************************************************"
            write(*,*)
            stop
         endif

         open(unit=DIG_PARM_UNIT,file='fort.dclaw',status="unknown",action="write")

         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETDIG:'
         write(DIG_PARM_UNIT,*) '---------'
         write(DIG_PARM_UNIT,*) '    rho_s:', rho_s
         write(DIG_PARM_UNIT,*) '    rho_f:', rho_f
         write(DIG_PARM_UNIT,*) '    m_crit:', m_crit
         write(DIG_PARM_UNIT,*) '    m0:', m0
         write(DIG_PARM_UNIT,*) '    mref:', mref
         write(DIG_PARM_UNIT,*) '    kref:', kref
         write(DIG_PARM_UNIT,*) '    phi:', phi/deg2rad
         write(DIG_PARM_UNIT,*) '    delta:', delta
         write(DIG_PARM_UNIT,*) '    mu:', mu
         write(DIG_PARM_UNIT,*) '    alpha_c:', alpha_c
         write(DIG_PARM_UNIT,*) '    c1:', c1
         write(DIG_PARM_UNIT,*) '    sigma_0:', sigma_0
         write(DIG_PARM_UNIT,*) '    dd_manning:', dd_manning
         write(DIG_PARM_UNIT,*) '    manning_max:', manning_max
         write(DIG_PARM_UNIT,*) '    src2method:', src2method
         write(DIG_PARM_UNIT,*) '    alphamethod:', alphamethod
         write(DIG_PARM_UNIT,*) '    rarecorrectortest:', rarecorrectortest
         write(DIG_PARM_UNIT,*) '    riemann_method:', riemann_method
         write(DIG_PARM_UNIT,*) '    bed_normal:', bed_normal
         write(DIG_PARM_UNIT,*) '    theta_input:', theta_input/deg2rad
         write(DIG_PARM_UNIT,*) '    entrainment:', entrainment
         write(DIG_PARM_UNIT,*) '    entrainment_method:', entrainment_method
         write(DIG_PARM_UNIT,*) '    entrainment_rate:', entrainment_rate
         write(DIG_PARM_UNIT,*) '    me:', me
         write(DIG_PARM_UNIT,*) '    segregation:', segregation
         write(DIG_PARM_UNIT,*) '    beta_seg:', beta_seg
         write(DIG_PARM_UNIT,*) '    chi0:', chi0
         write(DIG_PARM_UNIT,*) '    chie:', chie
         write(DIG_PARM_UNIT,*) '    mom_autostop:', mom_autostop
         write(DIG_PARM_UNIT,*) '    momlevel:', momlevel
         write(DIG_PARM_UNIT,*) '    curvature:', curvature
         close(DIG_PARM_UNIT)


   end subroutine set_dig

    ! ========================================================================
    !  set_pinit(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
   subroutine set_pinit(fname)

        implicit none

        ! Input
        character*25, intent(in), optional :: fname

        ! Locals
        integer, parameter :: iunit = 127
        character*25 :: file_name
        logical :: found_file

         ! Read user parameters from setgeo.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'pinit_dclaw.data'
         endif
         inquire(file=file_name,exist=found_file)
         if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
         endif

         call opendatafile(iunit, file_name)
         read(iunit,*) init_ptype
         read(iunit,*) init_pratio

         close(unit=iunit)

         init_pmin_ratio = 1.d16
         grad_eta_max = 0.d0
         cohesion_max = 0.d0
         grad_eta_ave = 0.d0
         eta_cell_count = 1.e-6

         open(unit=PINIT_PARM_UNIT,file='fort.pinit',status="unknown",action="write")
         write(PINIT_PARM_UNIT,*) ' '
         write(PINIT_PARM_UNIT,*) '--------------------------------------------'
         write(PINIT_PARM_UNIT,*) 'SETPINIT:'
         write(PINIT_PARM_UNIT,*) '---------'
         write(PINIT_PARM_UNIT,*) '    init_ptype:',init_ptype
         write(PINIT_PARM_UNIT,*) '    init_pratio:',init_pratio
         close(PINIT_PARM_UNIT)

   end subroutine set_pinit

   !====================================================================
   !subroutine qfix
   !accept solution q, return admissible q and primitive vars: u,v,m,rho,chi
   !====================================================================

   subroutine qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

      implicit none

      !i/o
      double precision, intent(in) :: gz
      double precision, intent(inout) :: h,hu,hv,hm,p,hchi
      double precision, intent(inout) :: u,v,m,rho,chi

      !Locals
      double precision :: mmin,mmax,chimin,chimax,pmax,pmin
      logical :: debug

      debug = .false.

      if (h.le.dry_tolerance) then
         if (debug.and.h<0.d0) write(*,*) '****WARNING******** QFIX: dry state encountered, h:', h, hm
         h =  0.d0
         hu = 0.d0
         hv = 0.d0
         hm = 0.d0
         p  = 0.d0
         hchi = 0.d0
         u = 0.d0
         v = 0.d0
         m = 0.d0 ! m and chi of dry state should not matter because in
         chi = 0.d0 ! riemann solver (and elsewhere) the neighbor of a dry
         ! cell should be used instead. (nowhere should a physically undefined
         ! variable be used).
         rho = 0.d0
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h
      chi = hchi/h

      !mlo = 1.d-3
      if (debug.and.hm<-dry_tolerance) write(*,*) '****WARNING******** QFIX: negative solid h,hm,m:', h, hm, m
      mmin = 0.0d0
      mmax = 1.d0
      m = max(m,mmin)
      m = min(m,mmax)
      hm = h*m

      chimin = 0.0d0
      chimax = 1.d0
      chi = max(chi,chimin)
      chi = min(chi,chimax)
      hchi = h*chi

      rho = rho_s*m + (1.d0-m)*rho_f
      pmax = rho*gz*h
      pmin = 0.d0
      p = max(0.d0,p)
      p = min(pmax,p)

      return

   end subroutine qfix

   !====================================================================
   !subroutine qfix_cmass
   !find physically admissible values of h,m,p,rho while holding rhoh,u,v const.
   ! 'cmass' = constant mass (constant rho*h)
   !====================================================================

   subroutine qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)

      implicit none

      !i/o
      double precision, intent(in) :: gz
      double precision, intent(in) :: u,v,rhoh
      double precision, intent(inout) :: h,m,rho,p_exc,p
      double precision, intent(out):: hu,hv,hm
      !Locals
      double precision :: mmin,mmax,p_exc_max,p_exc_min

      mmin = 0.0d0
      mmax = 1.d0
      m = max(m,mmin)
      m = min(m,mmax)

      rho = m*(rho_s - rho_f) + rho_f
      h = rhoh/rho

      !p_exc = p - rho_f*gz*h
      p_exc_max = rhoh*gz-rho_f*h*gz
      p_exc_min = 0.d0-rho_f*h*gz
      p_exc = max(p_exc_min,p_exc)
      p_exc = min(p_exc_max,p_exc)

      p = rho_f*gz*h + p_exc

      hu = h*u
      hv = h*v
      hm = h*m

      return

   end subroutine qfix_cmass

   !====================================================================
   ! subroutine setvars: evaluates needed variable parameters that are
   !                     functions of the solution vector q
   !====================================================================

   subroutine setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)

      implicit none

      !i/o
      real(kind=8), intent(in)  :: h,u,v,m,p,chi,gz
      real(kind=8), intent(out) :: tau,m_eq,tanpsi,kperm,alphainv
      real(kind=8), intent(out) :: rho

      !local
      real(kind=8) :: vnorm,shear,Nden,Nnum,psi,delta_kr_order,kr_chi
      real(kind=8) :: sig_0,sig_eff
      real(kind=8) :: alpha_meq ! tanh function for setting meqn to m when m is less than ~0.3


      !calculate rho ! DIG: Cleaner approach would be for rho to be in, and
      ! always have been calculated beforehand.
      rho = m*(rho_s-rho_f) + rho_f

      !determine kperm
      !segregation effect on kperm
      ! DIG: Eventually extent and create a 'segregation.f90' file with
      ! segeval(...) subroutine that handles all aspects of segregation-
      ! induced changes (similar to how entrainment options are handled).

      kr_chi = 1.0d0
      if (segregation==1) then
         !how many orders of magnitude should kperm change with chi=[0,1]
         delta_kr_order = 2.d0
         !DIG: check whether chi = 1 is coarse or fine
         kr_chi = 10.d0**(delta_kr_order*2.d0*(chi-0.5d0))
      endif
      kperm = kr_chi*kref*exp(-(m-mref)/(0.04d0))

      !determine vars related to m_eq and compressibility
      sig_eff = max(0.d0,rho*gz*h - p)

      !  determine 1/alpha
      select case (alphamethod)
      case(0)
         ! case 0: old d-claw, constant sig_0, set as user defined sigma_0
         sig_0 = sigma_0
      case(1)
         ! case 1: new method that ensures pressure always decays to hydrostatic
         ! sig_0 has a theoretical lower bound such that alphainv.ge.0.
         ! define sigma_0 so that it is just above the bound (use 0.51 instead
         ! of 0.50).
         sig_0 = 0.50d0*alpha_c*rho_f*(rho_s-rho_f)*gz*h/rho
      end select
      alphainv = m*(sig_eff + sig_0)/alpha_c

      ! Calculate velocity magnitude
      vnorm = sqrt(u**2 + v**2)

      ! Calculate the shear rate (gamma dot in the equations)
      if (h.gt.dry_tolerance) then
        shear = 2.d0*vnorm/h
      else
        shear = 0.d0
      endif

      !determine m_eq
      !m_eq = m_crit* 1/(1 + sqrt(Nnum/Nden))
      !m_eq = m_crit* sqrt(Nden)/(sqrt(Nden)+sqrt(Nnum))
      Nden = rho_s*(shear*delta)**2 + sig_eff
      Nnum = mu*shear

      alpha_meq = min(max(0.5d0*(1-tanh(60.d0*(m-0.35d0))), 0.d0), 1.d0)
      if (m.gt.0.4) alpha_meq = 0.d0
      if (m.lt.0.3) alpha_meq = 1.d0

      if (Nnum<=0.d0) then
         m_eq = alpha_meq*m + (1-alpha_meq)*m_crit
      else
         m_eq = alpha_meq*m + (1-alpha_meq)*m_crit*(sqrt(Nden)/(sqrt(Nden)+ sqrt(Nnum)))
      endif

      !Note: c1 is an adjustable parameter that we have traditionally set to 1.0.
      !For problems where one wants to avoid the complications of dilatancy and use
      !     a simpler rheological model, can be set to 0 or a small value.
      tanpsi = c1*(m-m_eq)
      psi = atan(c1*(m-m_eq))

      !calculate coulomb friction tau
      ! Note: for v=0, bounds on tau for static friction are determined in Riemann solver
      !        because the bounds are due to gradients in q. This routine determines vars
      !        from pointwise cell-centered value q(i).
      tau = sig_eff*tan(phi + psi)

      ! Note:  for water (m=0) sig_eff=0.0 and so tau=0.
      !        However, for dilute suspensions,
      !        we make o(m) not O(m))
      tau = tau*0.5d0*(1.d0 + tanh(100.d0*(m-0.1d0)))

      return

end subroutine setvars

   ! ========================================================================
   !  calc_taudir_riemann0
   ! ========================================================================
   !  Determines the resistive force vector for static cells
   !  outputs direction cosines at each interface
   !
   !  This method is used when riemann_method = 0
   ! ========================================================================

subroutine calc_taudir_riemann0(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc,meqn,maux

      !Locals
      double precision :: gz
      double precision :: h,hu,hv,hm,p,hchi,b,eta
      double precision :: hL,huL,hvL,hmL,pL,hchiL,bL,etaL
      double precision :: hR,huR,hvR,hmR,pR,hchiR,bR,etaR
      double precision :: hB,huB,hvB,hmB,pB,hchiB,bB,etaB
      double precision :: hT,huT,hvT,hmT,pT,hchiT,bT,etaT

      double precision :: u,v,m,chi
      double precision :: uL,vL,mL,chiL,rhoL
      double precision :: uR,vR,mR,chiR,rhoR
      double precision :: uB,vB,mB,chiB,rhoB
      double precision :: uT,vT,mT,chiT,rhoT
      double precision :: theta
      double precision :: tau,rho,alphainv
      double precision :: tanpsi,kperm,m_eq
      double precision :: Fx,Fy,FxL,FxR,FyL,FyR,FyC,FxC,dot,hvnorm,Fproj

      integer :: i,j

      gz = grav

      do j=2-mbc,my+mbc-1
         do i=2-mbc,mx+mbc-1

            if (bed_normal.eq.1) then
              theta = aux(i_theta,i,j)
              gz = grav*dcos(theta)
            else
              theta = 0.d0
            endif

            ! Get h, hu, hv, hm, p, hchi at cell center (no suffix)
            ! (L)eft (R)ight (T)op and (B)ottom cells of the 3x3
            ! cell stencil surrounding cell i,j

            h = q(i_h,i,j)
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p  = q(i_pb,i,j)
            hchi = q(i_hchi,i,j)

            call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

            b = aux(1,i,j)-q(i_bdif,i,j)
            eta = h+b
            phi = aux(i_phi,i,j)

            hL = q(i_h,i-1,j)
            huL= q(i_hu,i-1,j)
            hvL= q(i_hv,i-1,j)
            hmL = q(i_hm,i-1,j)
            pL  = q(i_pb,i-1,j)
            hchiL = q(i_hchi,i-1,j)

            call qfix(hL,huL,hvL,hmL,pL,hchiL,uL,vL,mL,chiL,rhoL,gz)

            bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
            etaL= hL+bL
            if (hL<dry_tolerance) then
               etaL = min(etaL,eta)
            endif

            hR = q(i_h,i+1,j)
            huR= q(i_hu,i+1,j)
            hvR= q(i_hv,i+1,j)
            hmR = q(i_hm,i+1,j)
            pR  = q(i_pb,i+1,j)
            hchiR = q(i_hchi,i+1,j)

            call qfix(hR,huR,hvR,hmR,pR,hchiR,uR,vR,mR,chiR,rhoR,gz)

            bR = aux(1,i+1,j)-q(i_bdif,i+1,j)
            etaR= hR+bR
            if (hR<dry_tolerance) then
               etaR = min(etaR,eta)
            endif

            hB = q(i_h,i,j-1)
            huB= q(i_hu,i,j-1)
            hvB= q(i_hv,i,j-1)
            hmB = q(i_hm,i,j-1)
            pB  = q(i_pb,i,j-1)
            hchiB = q(i_hchi,i,j-1)

            call qfix(hB,huB,hvB,hmB,pB,hchiB,uB,vB,mB,chiB,rhoB,gz)

            bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
            etaB= hB+bB
            if (hB<dry_tolerance) then
               etaB = min(etaB,eta)
            endif

            hT = q(i_h,i,j+1)
            huT= q(i_hu,i,j+1)
            hvT= q(i_hv,i,j+1)
            hmT = q(i_hm,i,j+1)
            pT  = q(i_pb,i,j+1)
            hchiT = q(i_hchi,i,j+1)

            call qfix(hT,huT,hvT,hmT,pT,hchiT,uT,vT,mT,chiT,rhoT,gz)

            bT = aux(1,i,j+1)-q(i_bdif,i,j+1)
            etaT= hT+bT
            if (hT<dry_tolerance) then
               etaT = min(etaT,eta)
            endif

            ! If cell center is dry, set value used for
            ! cell center here based on left/right cell
            ! values.

            if (h<dry_tolerance) then
               eta = min(etaL,eta)
               eta = min(etaB,eta)
            endif

            ! If all cells in the center are dry, return
            ! zeros for taudir_x, _y, and fsphi

            if ((h+hL+hB+hR+hT)<dry_tolerance) then
               aux(i_taudir_x,i,j) = 0.d0
               aux(i_taudir_y,i,j) = 0.d0
               aux(i_fsphi,i,j) = 0.d0
               cycle
            endif

            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)

            !minmod gradients
            FxC = -gz*h*(EtaR-EtaL)/(2.d0*dx) + gz*h*dsin(theta)
            FyC = -gz*h*(EtaT-EtaB)/(2.d0*dy)

            FxL = -gz*0.5d0*(h+hL)*(Eta-EtaL)/(dx) + gz*0.5d0*(h+hL)*dsin(theta)
            FyL = -gz*0.5d0*(h+hB)*(Eta-EtaB)/(dy)

            FxR = -gz*0.5d0*(h+hR)*(EtaR-Eta)/(dx) + gz*0.5d0*(h+hR)*dsin(theta)
            FyR = -gz*0.5d0*(h+hT)*(EtaT-Eta)/(dy)

            if (FxL*FxR>0.d0) then
               Fx = sign(min(abs(FxL),abs(FxR)),FxL)
            else
               Fx = 0.d0
            endif

            if (FyL*FyR>0.d0) then
               Fy = sign(min(abs(FyL),abs(FyR)),FyL)
            else
               Fy = 0.d0
            endif

            ! Calculate the magnitude of momentum.
            hvnorm = sqrt(hu**2 + hv**2)
            if (hvnorm>0.d0) then ! If moving.

               ! In D-Claw 4 dx and dy were accessible in the riemann solver.
               ! To fix in the transition, dx and dy are multiplied by taudir_x and y
               ! here. This works for everything except for bed normal and produces
               ! equivalent results (1/13/24)
               aux(i_taudir_x,i,j) = -dx*hu/hvnorm
               aux(i_taudir_y,i,j) = -dy*hv/hvnorm
               dot = min(max(0.d0,Fx*hu) , max(0.d0,Fy*hv))
               if (dot>0.d0) then
                  !friction should oppose direction of velocity
                  !if net force is in same direction, split friction source term
                  !splitting is useful for small velocities and nearly balanced forces
                  !only split amount up to maximum net force for large velocities
                  !aux has cell centered interpretation in Riemann solver
                  Fproj = dot/hvnorm
                  aux(i_fsphi,i,j) = min(1.d0,Fproj*rho/max(tau,1.d-16))

                  ! DIG: ensure that very low m does not contribute to static friction.
                  ! KRB thinks that this might be the right place to handle.

               else
                  !net force is in same direction as friction
                  !if nearly balanced steady state not due to friction
                  !no splitting, integrate friction in src
                  aux(i_fsphi,i,j) = 0.d0
               endif

            else
               !aux now have cell edge interpretation in Riemann solver
               !friction should oppose net force. resolve in Riemann solver
               if ((FxL**2+Fy**2)>0.d0) then
                  aux(i_taudir_x,i,j) = -dx*FxL/sqrt(FxL**2+Fy**2)
               else
                  aux(i_taudir_x,i,j) = dx*1.d0
               endif

               if ((Fx**2+FyL**2)>0.d0) then
                  aux(i_taudir_y,i,j) = -dy*FyL/sqrt(Fx**2+FyL**2)
               else
                  !There is no motion or net force. Resolve in src after Riemann
                  aux(i_taudir_y,i,j) = dy*1.d0
               endif

               ! Calculate factor of safety
               if ((aux(i_taudir_y,i,j)**2 + aux(i_taudir_x,i,j)**2)>0.d0) then
                  aux(i_fsphi,i,j) = 1.d0
               else
                  aux(i_fsphi,i,j) = 0.d0
               endif
            endif

         enddo
      enddo

end subroutine calc_taudir_riemann0

   ! ========================================================================
   !  calc_taudir_riemann1
   ! ========================================================================
   !  Determines the force balance for static cell interface
   !  using estimates of the total force direction
   !  from the 4 transverse forces connected to interface
   !  if any of the 4 configurations fail, the interface is allowed to
   ! fail and dynamic friction will be integrated in the RS
   !  aux(i_taudir_(x,y),i,j) = dx,dy if failure allowed, otherwise 0. to be
   ! interpretted correctly in RS. Note that if failure is allowed, the normal
   ! RS may not necessarily lead to waves (e.g., steep slope in transverse direction
   ! but uniform in normal direction, like a flume).
   ! for interfaces with motion or connected to an interface with motion, failure allowed
   !
   ! This method is used for riemann_method = 1
   ! ========================================================================

subroutine calc_taudir_riemann1(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc,meqn,maux

      !Locals
      double precision :: FdyV(4),FdxV(4),wdy(4),wdx(4)
      double precision :: gz
      double precision :: h,hu,hv,hm,p,hchi,b,eta
      double precision :: hL,huL,hvL,hmL,pL,hchiL,bL,etaL
      double precision :: hR,huR,hvR,hmR,pR,hchiR,bR,etaR
      double precision :: hB,huB,hvB,hmB,pB,hchiB,bB,etaB
      double precision :: hT,huT,hvT,hmT,pT,hchiT,bT,etaT
      double precision :: hTL,huTL,hvTL,hmTL,pTL,hchiTL,bTL,etaTL
      double precision :: hBL,huBL,hvBL,hmBL,pBL,hchiBL,bBL,etaBL
      double precision :: hBR,huBR,hvBR,hmBR,pBR,hchiBR,bBR,etaBR

      double precision :: u,v,m,chi,veltol
      double precision :: uL,vL,mL,chiL,rhoL
      double precision :: uR,vR,mR,chiR,rhoR
      double precision :: uB,vB,mB,chiB,rhoB
      double precision :: uT,vT,mT,chiT,rhoT
      double precision :: uTL,vTL,mTL,chiTL,rhoTL
      double precision :: uBL,vBL,mBL,chiBL,rhoBL
      double precision :: uBR,vBR,mBR,chiBR,rhoBR
      double precision :: theta
      double precision :: tau,rho,alphainv,tauB,tauL
      double precision :: tanpsi,kperm,m_eq
      double precision :: Fdxmax,Fdymax,Fmag,Fresist
      double precision :: FdyBL,FdyB,FdyTL,FdyT,FdxL,FdxBL,FdxR,FdxBR
      double precision :: Fdx,Fdy,nx,ny,dirnorm

      integer :: i,j

      logical reproject_cosines


      gz = grav
      veltol = 1.d-10

      !reproject to correct if force direction is inconsistent
      !between x and y stencils for 2 edges (left and lower)
      !one should expect ||taudir_x,taudir_y||=||(dx,dy)||
      reproject_cosines = .true.

      !first for cell edge between i,j and i-1,j for x direction
      do j=2-mbc,my+mbc-1
         do i=2-mbc,mx+mbc-1

            if (bed_normal.eq.1) then
              theta = aux(i_theta,i,j)
              gz = grav*dcos(theta)
            else
              theta = 0.d0
            endif
            aux(i_fsphi,i,j) = 10.d0
            ! Get h, hu, hv, hm, p, hchi at cell center (no suffix)
            ! (L)eft (R)ight (T)op and (B)ottom cells of the 3x3
            ! cell stencil surrounding cell i,j

            h = q(i_h,i,j)
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p  = q(i_pb,i,j)
            hchi = q(i_hchi,i,j)
            b = aux(1,i,j)-q(i_bdif,i,j)
            eta = h+b
            phi = aux(i_phi,i,j)
            call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

            hL = q(i_h,i-1,j)
            huL= q(i_hu,i-1,j)
            hvL= q(i_hv,i-1,j)
            hmL = q(i_hm,i-1,j)
            pL  = q(i_pb,i-1,j)
            hchiL = q(i_hchi,i-1,j)
            bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
            etaL= hL+bL
            call qfix(hL,huL,hvL,hmL,pL,hchiL,uL,vL,mL,chiL,rhoL,gz)

            ! If both cells are dry, no RS anyway, move on
            if (max(h,hL)<dry_tolerance) then
               aux(i_taudir_x,i,j) = dx
               cycle
            endif

            hB = q(i_h,i,j-1)
            huB= q(i_hu,i,j-1)
            hvB= q(i_hv,i,j-1)
            hmB = q(i_hm,i,j-1)
            pB  = q(i_pb,i,j-1)
            hchiB = q(i_hchi,i,j-1)
            bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
            etaB= hB+bB
            call qfix(hB,huB,hvB,hmB,pB,hchiB,uB,vB,mB,chiB,rhoB,gz)

            hT = q(i_h,i,j+1)
            huT= q(i_hu,i,j+1)
            hvT= q(i_hv,i,j+1)
            hmT = q(i_hm,i,j+1)
            pT  = q(i_pb,i,j+1)
            hchiT = q(i_hchi,i,j+1)
            bT = aux(1,i,j+1)-q(i_bdif,i,j+1)
            etaT= hT+bT
            call qfix(hT,huT,hvT,hmT,pT,hchiT,uT,vT,mT,chiT,rhoT,gz)

            hTL = q(i_h,i-1,j+1)
            huTL= q(i_hu,i-1,j+1)
            hvTL= q(i_hv,i-1,j+1)
            hmTL = q(i_hm,i-1,j+1)
            pTL  = q(i_pb,i-1,j+1)
            hchiTL = q(i_hchi,i-1,j+1)
            bTL = aux(1,i-1,j+1)-q(i_bdif,i-1,j+1)
            etaTL= hTL+bTL
            call qfix(hTL,huTL,hvTL,hmTL,pTL,hchiTL,uTL,vTL,mTL,chiTL,rhoTL,gz)

            hBL = q(i_h,i-1,j-1)
            huBL= q(i_hu,i-1,j-1)
            hvBL= q(i_hv,i-1,j-1)
            hmBL = q(i_hm,i-1,j-1)
            pBL  = q(i_pb,i-1,j-1)
            hchiBL = q(i_hchi,i-1,j-1)
            bBL = aux(1,i-1,j-1)-q(i_bdif,i-1,j-1)
            etaBL= hBL+bBL
            call qfix(hBL,huBL,hvBL,hmBL,pBL,hchiBL,uBL,vBL,mBL,chiBL,rhoBL,gz)

            ! if any cells in stencil have motion allow failure (determined in RS)
            if ((hu**2+huL**2+huB**2+huT**2+huBL**2+huTL**2+hv**2+hvL**2+hvB**2+hvT**2+hvBL**2+hvTL**2)>veltol) then
               aux(i_taudir_x,i,j) = dx
               cycle
            endif

            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)
            call setvars(hL,uL,vL,mL,pL,chiL,gz,rhoL,kperm,alphainv,m_eq,tanpsi,tauL)
            ! if one cell is dry, use resistance of wet cell
            if (h<dry_tolerance) then
                tau = tauL
                rho = rhoL
            elseif (hL<dry_tolerance) then
                tauL = tau
                rhoL = rho
            endif

            !Find if failure for i'th cell edge (x-direction, j row)
            !Force in normal direction (x)
            call calc_interface_force(Fdx,gz,h,hL,b,bL,dx,theta)

            !Force in y-direction at 4 edges connected to i'th cell edge
            call calc_interface_force(FdyBL,gz,hL,hBL,bL,bBL,dy,0.d0)
            call calc_interface_force(FdyB,gz,h,hB,b,bB,dy,0.d0)
            call calc_interface_force(FdyTL,gz,hTL,hT,bTL,bT,dy,0.d0)
            call calc_interface_force(FdyT,gz,hT,h,bT,b,dy,0.d0)
            FdyV = [FdyT,FdyTL,FdyBL,FdyB]

            !find largest force in y skip below if ~0
            Fdymax = maxval(abs(FdyV))
            if (Fdymax<1.d-6) then
                aux(i_taudir_x,i,j) = dx
               cycle
            endif
            !weight the forces in y to determine Fdy
            !average forces gives 2nd order for smooth gradients
            !weight strategy in case there is a sharp discontinuity at one interface that dominates local slope (reduce Fdx)
            wdy = abs(FdyV)/sum(abs(FdyV))
            Fdy = dot_product(FdyV,wdy)
            Fmag = sqrt(Fdy**2 + Fdx**2)
            !find resistance at cell edge center
            Fresist = 0.5d0*(tauL/rhoL + tau/rho)
            if (Fmag.gt.Fresist.and.Fmag.gt.1.d-16) then
                aux(i_taudir_x,i,j) = dx*abs(Fdx)/Fmag
                !aux(i_taudir_x,i,j) = dx
                aux(i_fsphi,i,j) = min(aux(i_fsphi,i,j), Fresist/Fmag)
            elseif (Fmag.gt.1.d-16) then
                aux(i_taudir_x,i,j) = dx*abs(Fdx)/Fmag
                aux(i_fsphi,i,j) = min(aux(i_fsphi,i,j), Fresist/Fmag)
            else
                aux(i_taudir_x,i,j) = dx
                aux(i_fsphi,i,j) = 0.d0
            endif

         enddo
      enddo

      !repeat for cell edge between i,j and i,j-1 for y direction
      ! note: stencil to check for transverse forces (x-direction) is different
      ! (left 2/3 of stencil above and bottom 2/3 of stencil below)
      ! so new loop is performed because of the cycle instructions in first loop

      do j=2-mbc,my+mbc-1
         do i=2-mbc,mx+mbc-1

            if (bed_normal.eq.1) then
              theta = aux(i_theta,i,j)
              gz = grav*dcos(theta)
            else
              theta = 0.d0
            endif

            ! Get h, hu, hv, hm, p, hchi at cell center (no suffix)
            ! (L)eft (R)ight (T)op and (B)ottom cells of the 3x3
            ! cell stencil surrounding cell i,j

            h = q(i_h,i,j)
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p  = q(i_pb,i,j)
            hchi = q(i_hchi,i,j)
            b = aux(1,i,j)-q(i_bdif,i,j)
            eta = h+b
            phi = aux(i_phi,i,j)
            call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

            hB = q(i_h,i,j-1)
            huB= q(i_hu,i,j-1)
            hvB= q(i_hv,i,j-1)
            hmB = q(i_hm,i,j-1)
            pB  = q(i_pb,i,j-1)
            hchiB = q(i_hchi,i,j-1)
            bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
            etaB= hB+bB
            call qfix(hB,huB,hvB,hmB,pB,hchiB,uB,vB,mB,chiB,rhoB,gz)

            ! If both cells are dry, no RS anyway, move on
            if (max(h,hB)<dry_tolerance) then
               aux(i_taudir_y,i,j) = dy
               cycle
            endif

            hL = q(i_h,i-1,j)
            huL= q(i_hu,i-1,j)
            hvL= q(i_hv,i-1,j)
            hmL = q(i_hm,i-1,j)
            pL  = q(i_pb,i-1,j)
            hchiL = q(i_hchi,i-1,j)
            bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
            etaL= hL+bL
            call qfix(hL,huL,hvL,hmL,pL,hchiL,uL,vL,mL,chiL,rhoL,gz)

            hBL = q(i_h,i-1,j-1)
            huBL= q(i_hu,i-1,j-1)
            hvBL= q(i_hv,i-1,j-1)
            hmBL = q(i_hm,i-1,j-1)
            pBL  = q(i_pb,i-1,j-1)
            hchiBL = q(i_hchi,i-1,j-1)
            bBL = aux(1,i-1,j-1)-q(i_bdif,i-1,j-1)
            etaBL= hBL+bBL
            call qfix(hBL,huBL,hvBL,hmBL,pBL,hchiBL,uBL,vBL,mBL,chiBL,rhoBL,gz)

            hBR = q(i_h,i+1,j-1)
            huBR= q(i_hu,i+1,j-1)
            hvBR= q(i_hv,i+1,j-1)
            hmBR = q(i_hm,i+1,j-1)
            pBR  = q(i_pb,i+1,j-1)
            hchiBR = q(i_hchi,i+1,j-1)
            bBR = aux(1,i+1,j-1)-q(i_bdif,i+1,j-1)
            etaBR= hBR+bBR
            call qfix(hBR,huBR,hvBR,hmBR,pBR,hchiBR,uBR,vBR,mBR,chiBR,rhoBR,gz)

            hR = q(i_h,i+1,j)
            huR= q(i_hu,i+1,j)
            hvR= q(i_hv,i+1,j)
            hmR = q(i_hm,i+1,j)
            pR  = q(i_pb,i+1,j)
            hchiR = q(i_hchi,i+1,j)
            bR = aux(1,i+1,j)-q(i_bdif,i+1,j)
            etaR= hR+bR
            call qfix(hR,huR,hvR,hmR,pR,hchiR,uR,vR,mR,chiR,rhoR,gz)

            ! if any cells in stencil have motion allow failure
            if ((hu**2+huL**2+huR**2+huBL**2+huB**2+huBR**2+hv**2+hvL**2+hvR**2+hvBL**2+hvB**2+hvBR**2)>veltol) then
               aux(i_taudir_y,i,j) = dy
               cycle
            endif

            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)
            call setvars(hB,uB,vB,mB,pB,chiB,gz,rhoB,kperm,alphainv,m_eq,tanpsi,tauB)
            ! if one cell is dry, use resistance of wet cell
            if (h<dry_tolerance) then
                tau = tauB
                rho = rhoB
            elseif (hB<dry_tolerance) then
                tauB = tau
                rhoB = rho
            endif


            !Find if failure for j'th cell edge (y-direction, i column)
            !Force in normal direction (y)
            call calc_interface_force(Fdy,gz,h,hB,b,bB,dx,0.d0)

            !Force in y-direction at 4 edges connected to i'th cell edge
            call calc_interface_force(FdxL,gz,h,hL,b,bL,dx,theta)
            call calc_interface_force(FdxBL,gz,hBL,hB,bBL,bB,dx,theta)
            call calc_interface_force(FdxR,gz,hR,h,bR,b,dx,theta)
            call calc_interface_force(FdxBR,gz,hBR,hB,bBR,bB,dx,theta)
            FdxV = [FdxL,FdxBL,FdxR,FdxBR]

            !find largest force in y skip below if ~0
            Fdxmax = maxval(abs(FdxV))
            if (Fdxmax<1.d-6) then
                aux(i_taudir_y,i,j) = dy
               cycle
            endif
            !weight the forces in y to determine Fdy
            !average forces gives 2nd order for smooth gradients
            !weight strategy in case there is a sharp discontinuity at one interface that dominates local slope (reduce Fdx)
            wdx = abs(FdxV)/sum(abs(FdxV))
            Fdx = dot_product(FdxV,wdx)
            Fmag = sqrt(Fdy**2 + Fdx**2)
            !find resistance at cell edge center
            Fresist = 0.5d0*(tauB/rhoB + tau/rho)
            if (Fmag.gt.Fresist.and.Fmag.gt.1.d-16) then
                aux(i_taudir_y,i,j) = dy*abs(Fdy)/Fmag
                !aux(i_taudir_y,i,j) = dy
                aux(i_fsphi,i,j) = min(aux(i_fsphi,i,j), Fresist/Fmag)
            elseif (Fmag.gt.1.d-16) then
                aux(i_taudir_y,i,j) = dy*abs(Fdy)/Fmag
                aux(i_fsphi,i,j) = min(aux(i_fsphi,i,j), Fresist/(Fmag+1.d-16))
            else
                aux(i_taudir_y,i,j) = dy
                aux(i_fsphi,i,j) = 0.d0
            endif

         enddo
      enddo

      if (reproject_cosines) then
      !fix if new dx,dy ~= nx*dx + ny*dy not ~ |(dx,dy)|
      !rescale both so direction cosines at cell corner
        do j=2-mbc,my+mbc-1
            do i=2-mbc,mx+mbc-1
                h = q(i_h,i,j)
                nx = aux(i_taudir_x,i,j)/dx
                ny = aux(i_taudir_y,i,j)/dy
                dirnorm = max(sqrt(nx**2 + ny**2),1.d-12)

                aux(i_taudir_x,i,j) = aux(i_taudir_x,i,j)/(dirnorm)
                aux(i_taudir_y,i,j) = aux(i_taudir_y,i,j)/(dirnorm)
            enddo
        enddo
      endif

end subroutine calc_taudir_riemann1


   ! ========================================================================
   !  calc_pmin
   ! ========================================================================
   !  Determines minimum pore pressure for mobilization
   !  Determines factor of safety and cohesion for static states
   ! ========================================================================

subroutine calc_pmin(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc,meqn,maux

      !Locals
      double precision :: h,hL,hR,hu,hv,b,bL,bR,bT,bB,hT,hB
      double precision :: phi,theta,rho
      double precision :: gmod,dry_tol
      double precision :: EtaL,EtaR,EtaT,EtaB,Eta
      double precision :: detadx,detadxL,detadxR,detady,detadyT,detadyB
      double precision :: grad_eta

      integer :: i,j

      dry_tol = dry_tolerance
      gmod = grav
      rho = m0*rho_s + (1.d0-m0)*rho_f

      do j=1,my
         do i=1,mx


            ! Get values of q from the 3x3 cell stencil

            h = q(i_h,i,j)
            hL = q(i_h,i-1,j)
            hR = q(i_h,i+1,j)
            if (h<dry_tol) then
               cycle
            endif

            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)

            if ((hu**2+hv**2)>0.d0) then
               cycle
            endif

            b = aux(1,i,j)
            bR = aux(1,i+1,j)
            bL = aux(1,i-1,j)
            phi = aux(i_phi,i,j)

            if ((phi)==0.d0) then
               init_pmin_ratio = 0.d0
               cycle
            endif

            hT = q(i_h,i,j+1)
            bT = aux(1,i,j+1)
            hB = q(i_h,i,j-1)
            bB = aux(1,i,j-1)

            if (bed_normal.eq.1) then
               theta = aux(i_theta,i,j)
               gmod = grav*cos(theta)
            else
               theta = 0.d0
            endif

            Eta  = h+b

            ! calculate the max deta/dx
            EtaR = hR+bR
            EtaL = hL+bL
            if (hR<=dry_tol) then
               EtaR = min(Eta,bR)
            endif
            if (hL<=dry_tol) then
               EtaL = min(Eta,bL)
            endif
            detadxR = (EtaR-Eta)/dx -tan(theta)
            detadxL = (Eta-EtaL)/dx -tan(theta)
            if (detadxR*detadxL<=0.d0) then
               detadx = 0.d0
            elseif (abs(detadxR)>abs(detadxL)) then
               detadx = detadxL
            else
               detadx = detadxR
            endif

            ! calculate the max deta/dy
            EtaT = hT+bT
            EtaB = hB+bB
            if (hT<=dry_tol) then
               EtaT = min(Eta,bT)
            endif
            if (hB<=dry_tol) then
               EtaB = min(Eta,bB)
            endif
            detadyT = (EtaT-Eta)/dy
            detadyB = (Eta-EtaB)/dy
            if (detadyT*detadyB<=0.d0) then
               detady = 0.d0
            elseif (abs(detadyT)>abs(detadyB)) then
               detady = detadyB
            else
               detady = detadyT
            endif

            grad_eta = sqrt(detadx**2 + detady**2)
            grad_eta_ave = grad_eta_ave + grad_eta/tan(phi)
            eta_cell_count = eta_cell_count + 1.d0

            grad_eta_max = max(grad_eta_max,grad_eta/tan(phi))

            init_pmin_ratio = min(init_pmin_ratio, 1.d0-grad_eta/tan(phi))

         enddo
      enddo

      if (init_ptype==2) then
         init_pmin_ratio = 1.d0-(grad_eta_ave/eta_cell_count)
      endif
      if (init_ptype>0) then
         write(*,*) '--------------------------------------------'
         write(*,*) 'hydrostatic liquefaction ratio:', rho_f/rho
         write(*,*) 'initiation liquefaction  ratio:',init_pmin_ratio, grad_eta_ave
         write(*,*) 'maximum surface slope angle:',180.d0*atan(tan(phi)*grad_eta_max)/3.14d0, grad_eta_max
         write(*,*) 'average failure liquefaction ratio:', 1.d0-(grad_eta_ave/eta_cell_count) , eta_cell_count
         write(*,*) '--------------------------------------------'
      endif
   end subroutine calc_pmin

    !========================================
    !calc_interface_force
    !========================================
    !calculate the driving force at a cell edge between arbitrary
    !(static) left and right states
    !=======================================

    subroutine calc_interface_force(Fx,gz,hR,hL,bR,bL,dx,theta)

        implicit none

        !Input/Output
        double precision :: gz,hR,hL,bR,bL,dx,theta,Fx

        !Locals
        double precision :: bLdx,bRdx,etaL,etaR

        etaL = hL+bL
        etaR = hR+bR
        !adjust b for db/dx term if wet/dry data
        if (hL<dry_tolerance) then
            bLdx = min(etaL,etaR)
        else
            bLdx = bL
        endif
        if (hR<dry_tolerance) then
            bRdx = min(etaL,etaR)
        else
            bRdx = bR
        endif
        !Fx = (-.5d0*gz*(hR**2-hL**2)-gz*0.5*(hL+hR)*(bRdx-bLdx))/dx-(gz*0.5*(hL+hR)*tan(theta))*dx
        Fx = (-.5d0*gz*(hR**2-hL**2)-gz*0.5*(hL+hR)*(bRdx-bLdx))/dx-0.d0*(gz*0.5*(hL+hR)*tan(theta))*dx


    end subroutine calc_interface_force


end module digclaw_module

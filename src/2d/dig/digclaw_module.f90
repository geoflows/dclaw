
module digclaw_module

   use geoclaw_module, only: dry_tolerance, grav, pi, coordinate_system

   implicit none

    ! ========================================================================
    ! General digclaw parameters
    ! ========================================================================
    double precision :: rho_s,rho_f,phi_bed,theta_input,delta,kappita
    double precision :: mu,alpha,m_crit,c1,m0,beta_seg,sigma_0,entrainment_rate

    integer :: init_ptype,bed_normal,entrainment,curvature,alphamethod,src2method 
    integer :: segregation
    double precision :: init_pmax_ratio,init_ptf2,init_ptf,init_pmin_ratio
    double precision :: grad_eta_max,cohesion_max,grad_eta_ave,eta_cell_count
    double precision :: chi_init_val

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

    ! indicies of aux
    integer ::  i_dig
    integer ::  i_phi
    integer ::  i_theta
    integer ::  i_fsphi
    integer ::  i_taudir_x
    integer ::  i_taudir_y
    integer ::  i_ent

    integer, parameter ::  DIG_PARM_UNIT = 78


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

         ! hard code q indicies.
         i_h=1
         i_hu=2
         i_hv=3
         i_hm=4
         i_pb=5
         i_hchi=6
         i_bdif=7

         ! set aux index values based on coordinate system
         i_phi    = i_dig
         i_theta  = i_dig + 1
         i_fsphi  = i_dig + 2
         i_taudir_x = i_dig + 3
         i_taudir_y = i_dig + 4
         i_ent = i_dig + 5

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
         read(iunit,*) phi_bed
         phi_bed = deg2rad*phi_bed
         read(iunit,*) theta_input
         theta_input = deg2rad*theta_input
         read(iunit,*) delta
         read(iunit,*) kappita
         read(iunit,*) mu
         read(iunit,*) alpha
         read(iunit,*) m_crit
         read(iunit,*) c1
         read(iunit,*) m0
         read(iunit,*) sigma_0
         read(iunit,*) beta_seg
         read(iunit,*) bed_normal
         read(iunit,*) entrainment
         read(iunit,*) entrainment_rate
         read(iunit,*) chi_init_val
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
         write(DIG_PARM_UNIT,*) '    rho_s:',rho_s
         write(DIG_PARM_UNIT,*) '    rho_f:',rho_f
         write(DIG_PARM_UNIT,*) '    phi_bed:', phi_bed/deg2rad
         write(DIG_PARM_UNIT,*) '    theta_input:', theta_input/deg2rad
         write(DIG_PARM_UNIT,*) '    delta:', delta
         write(DIG_PARM_UNIT,*) '    kappita:', kappita
         write(DIG_PARM_UNIT,*) '    mu:', mu
         write(DIG_PARM_UNIT,*) '    alpha:', alpha
         write(DIG_PARM_UNIT,*) '    m_crit:', m_crit
         write(DIG_PARM_UNIT,*) '    c1:', c1
         write(DIG_PARM_UNIT,*) '    m0:', m0
         write(DIG_PARM_UNIT,*) '    sigma_0:', sigma_0
         write(DIG_PARM_UNIT,*) '    beta_seg:', beta_seg
         write(DIG_PARM_UNIT,*) '    bed_normal:', bed_normal
         write(DIG_PARM_UNIT,*) '    entrainment:', entrainment
         write(DIG_PARM_UNIT,*) '    entrainment_rate:', entrainment_rate
         write(DIG_PARM_UNIT,*) '    chi_init_val:', chi_init_val
         write(DIG_PARM_UNIT,*) '    mom_autostop:', mom_autostop
         write(DIG_PARM_UNIT,*) '    momlevel:', momlevel
         write(DIG_PARM_UNIT,*) '    curvature:', curvature


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
         read(iunit,*) init_pmax_ratio
         read(iunit,*) init_ptf ! DIG - some of these parameters are no longer needed
         read(iunit,*) init_ptf2
         close(unit=iunit)

         init_pmin_ratio = 1.d16
         grad_eta_max = 0.d0
         cohesion_max = 0.d0
         grad_eta_ave = 0.d0
         eta_cell_count = 1.e-6


         write(DIG_PARM_UNIT,*) ' '
         write(DIG_PARM_UNIT,*) '--------------------------------------------'
         write(DIG_PARM_UNIT,*) 'SETPINIT:'
         write(DIG_PARM_UNIT,*) '---------'
         write(DIG_PARM_UNIT,*) '    init_ptype:',init_ptype
         write(DIG_PARM_UNIT,*) '    init_pmax_ratio:',init_pmax_ratio
         write(DIG_PARM_UNIT,*) '    init_ptf:',init_ptf
         close(DIG_PARM_UNIT)

   end subroutine set_pinit

   !====================================================================
   !subroutine qfix
   !accept solution q, return admissible q and primitive vars: u,v,m,rho,chi
   !====================================================================

   subroutine qfix(h,hu,hv,hm,hchi,p,u,v,m,chi,rho,gz)

      implicit none

      !i/o
      double precision, intent(in) :: gz
      double precision, intent(inout) :: h,hu,hv,hm,p,hchi
      double precision, intent(inout) :: u,v,m,rho,chi

      !Locals
      double precision :: mmin,mmax,chimin,chimax,pmax,pmin

      if (h.le.drytolerance) then
         h =  0.d0
         hu = 0.d0
         hv = 0.d0
         hm = 0.d0
         p  = 0.d0 
         u = 0.d0
         v = 0.d0
         m = 0.d0
         chi = 0.d0
         rho = 0.d0
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h
      chi = hchi/h

      !mlo = 1.d-3
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
   !subroutine admissibleq
   !accept solution q, return q in admissible space
   !====================================================================

   subroutine admissibleq(h,hu,hv,hm,p,u,v,m,theta)

      implicit none

      !i/o
      double precision, intent(in) :: theta
      double precision, intent(inout) :: h,hu,hv,hm,p
      double precision, intent(out) :: u,v,m

      !Locals
      double precision :: mlo,mhi,hlo,pmax,phi,plo,rho,dry_tol,m_min,gmod

      gmod = grav
      dry_tol = dry_tolerance
      if (bed_normal.eq.1) gmod = grav*dcos(theta)

      if (h.le.dry_tol) then
         h =  0.d0
         hu = 0.d0
         hv = 0.d0
         hm = 0.d0
         p  = 0.d0 !h*gmod*rho_f
         u = 0.d0
         v = 0.d0
         m = 0.d0
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h

      mlo = 0.0d0
      mhi = 1.d0 - mlo

      if (m.lt.mlo) then
         m = dmax1(m,mlo)
         hm = h*m
      elseif (m.gt.mhi) then
         m = dmin1(m,1.d0)
         hm = h*m
      endif

      rho = rho_s*m + (1.d0-m)*rho_f
      pmax = rho*gmod*h
      p = dmin1(pmax,p)
      p = dmax1(0.d0,p)

      if (m < 1d-5) then
         ! reset p to hydrostatic pressure in pure water
         ! (need to figure out proper tolerance in test)
         ! (better way? E.g. relaxation toward hydrostatic in src2?)
         p = h*gmod*rho_f
      endif

      return

   end subroutine admissibleq

   !====================================================================
   ! subroutine setvars: evaluates needed variable parameters that are
   !                     functions of the solution vector q
   !====================================================================

   subroutine setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)

      implicit none

      !i/o
      real(kind=8), intent(inout) :: rho
      real(kind=8), intent(in)  :: h,u,v,m,p,gz
      real(kind=8), intent(out) :: alphainv,sig_0,sig_eff
      real(kind=8), intent(out) :: tau,m_eq,tanpsi,kperm

      !local
      real(kind=8) :: vnorm,shear,Nden,Nnum,psi


      !check rho
      rho = m*(rho_s-rho_f) + rho_f

      !determine kperm
      kperm = kappita*exp(-(m-m0)/(0.04d0))

      !determine vars related to m_eq and compressibility
      sig_eff = max(0.d0,rho*gz*h - p)

      alphamethod = 2
      select case (alphamethod) 
      case(0:1)
         sig_0 = sigma_0
      case(2)
         !sig_0 = alpha*(rho_s-rho_f)*gz*h
         sig_0 = 0.5d0*alpha*(rho_s-rho_f)*gz*h/rho
         alphainv = m*(sig_eff + sig_0)/alpha
      end select

      vnorm = sqrt(u**2 + v**2)
      shear = 2.d0*vnorm/h
      

      !determine m_eq
      !m_eq = m_crit* 1/(1 + sqrt(Nnum/Nden))
      !m_eq = m_crit* sqrt(Nden)/(sqrt(Nden)+sqrt(Nnum))
      Nden = rho_s*(shear*delta)**2 + sig_eff
      Nnum = mu*shear
      if (Nnum<=0.d0) then
         m_eq = m_crit
      else
         m_eq = m_crit*(sqrt(Nden)/(sqrt(Nden)+ sqrt(Nnum)))
      endif
      !Note: c1 is an adjustable parameter that we have traditionally set to 1.0.
      !For problems where one wants to avoid the complications of dilatancy and use
      !     a simpler rheological model, can be set to 0 or a small value.
      tanpsi = c1*(m-m_eq)
      psi = atan(c1*(m-m_eq)) !does atan return the correct angle always?

      !calculate coulomb friction tau
      ! Note: for v=0, bounds on tau for static friction are determined in Riemann solver
      !        because the bounds are due to gradients in q. This routine determines vars
      !        from pointwise cell-centered value q(i).
      tau = sig_eff*tan(phi_bed + 0.d0*psi)
      ! Note:  for water (m=0) sig_eff=0.0 and so tau=0.         
      !        However, for dilute suspensions, 
      !        we make o(m) not O(m))
      !if (m<0.55d0) then
      tau = tau*0.5d0*(1.d0 + tanh(50.d0*(m-0.40d0))) !+ 100.d0*shear))
      !endif

      return

      end subroutine setvars
   !====================================================================
   ! subroutine auxeval: evaluates the auxiliary variables as functions
   !                     of the solution vector q
   !====================================================================

   subroutine auxeval(h,u,v,m,p,phi_bed,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

      implicit none

      !i/o
      double precision, intent(inout) :: chi
      double precision, intent(in)  :: h,u,v,m,p,phi_bed,theta
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: sigbed,kperm,compress

      !local
      double precision :: m_eqn,vnorm,gmod,sigbedc,shear,tanphi,rho_fp
      double precision :: seg,pmtanh01,m_crit_m,m_crit_pm

      if (h.lt.dry_tolerance) return

      gmod=grav
      chi = max(0.0d0,chi)
      chi = min(1.0d0,chi)

      if (bed_normal.eq.1) gmod=grav*dcos(theta)
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_fp*(1.d0-m)
      shear = 2.0d0*vnorm/h
      sigbed = dmax1(0.d0,rho*gmod*h - p)
      sigbedc = rho_s*(shear*delta)**2 + sigbed
      if (sigbedc.gt.0.0d0) then
         S = (mu*shear/(sigbedc))
      else
         S = 0.d0
      endif
      !Note: m_eqn = m_crit/(1+sqrt(S))
      !From Boyer et. al

      kperm = kappita*exp(-(m-m0)/(0.04d0))!*(10**(pmtanh01))
      m_crit_pm =  0.d0

      m_crit_pm = pmtanh01*0.09d0
      m_crit_m = m_crit - m_crit_pm
      m_eqn = m_crit_m/(1.d0 + sqrt(S))
      tanpsi = c1*(m-m_eqn)*tanh(shear/0.1d0)

      if (m.le.1.d-16) then
         compress = 1.d16
         kperm = 0.0d0
         tanpsi = 0.0d0
         sigbed=0.0d0
      else
         compress = alpha/(m*(sigbed +  sigma_0))
      endif

      if (vnorm.le.0.d0) then
         tanpsi = 0.d0
         D = 0.d0
      elseif (h*mu.gt.0.d0) then
         D = 2.0d0*(kperm/(mu*h))*(rho_fp*gmod*h - p)
      else
         D = 0.d0
      endif

      tanphi = dtan(phi_bed + datan(tanpsi))

      tau = dmax1(0.d0,sigbed*tanphi)*((tanh(100.0*(m-0.05))+1.0)*0.5)

      !kappa: earth pressure coefficient
      kappa = 1.d0

   end subroutine auxeval


   ! ========================================================================
   !  calc_taudir
   ! ========================================================================
   !  Determines the resistive force vector for static cells
   !  outputs direction cosines at each interface
   ! ========================================================================

subroutine calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc,meqn,maux

      !Locals
      double precision :: gmod,dry_tol
      double precision :: h,hu,hv,hm,p,b,eta
      double precision :: hL,huL,hvL,hmL,pL,bL,etaL
      double precision :: hR,huR,hvR,hmR,pR,bR,etaR
      double precision :: hB,huB,hvB,hmB,pB,bB,etaB
      double precision :: hT,huT,hvT,hmT,pT,bT,etaT

      double precision :: u,v,m
      double precision :: uL,vL,mL
      double precision :: uR,vR,mR
      double precision :: uB,vB,mB
      double precision :: uT,vT,mT
      double precision :: thetaL,thetaB,theta
      double precision :: tau,tauL,tauR,tauB,tauT,rho,rhoL,rhoR,rhoT,rhoB
      double precision :: phi,kappa,S,tanpsi,D,sigbed,kperm,compress,pm
      double precision :: Fx,Fy,FxL,FxR,FyL,FyR,FyC,FxC,dot,vnorm,Fproj

      integer :: i,j

      dry_tol = dry_tolerance
      gmod = grav

      do j=2-mbc,my+mbc-1
         do i=2-mbc,mx+mbc-1

            if (bed_normal.eq.1) then
              theta = aux(i_theta,i,j)
              thetaL = aux(i_theta,i-1,j)
              thetaB = aux(i_theta,i,j-1)
              gmod = grav*cos(theta)
            else
              theta = 0.d0
              thetaL = 0.d0
              thetaB = 0.d0
            endif

            h = q(i_h,i,j)
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p  = q(i_pb,i,j)
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            b = aux(1,i,j)-q(i_bdif,i,j)
            eta = h+b
            phi = aux(i_phi,i,j)

            hL = q(i_h,i-1,j)
            huL= q(i_hu,i-1,j)
            hvL= q(i_hv,i-1,j)
            hmL = q(i_hm,i-1,j)
            pL  = q(i_pb,i-1,j)
            call admissibleq(hL,huL,hvL,hmL,pL,uL,vL,mL,theta)
            bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
            etaL= hL+bL
            if (hL<dry_tol) then
               etaL = min(etaL,eta)
            endif

            hR = q(i_h,i+1,j)
            huR= q(i_hu,i+1,j)
            hvR= q(i_hv,i+1,j)
            hmR = q(i_hm,i+1,j)
            pR  = q(i_pb,i+1,j)
            call admissibleq(hR,huR,hvR,hmR,pR,uR,vR,mR,theta)
            bR = aux(1,i+1,j)-q(i_bdif,i+1,j)
            etaR= hR+bR
            if (hR<dry_tol) then
               etaR = min(etaR,eta)
            endif

            hB = q(i_h,i,j-1)
            huB= q(i_hu,i,j-1)
            hvB= q(i_hv,i,j-1)
            hmB = q(i_hm,i,j-1)
            pB  = q(i_pb,i,j-1)
            call admissibleq(hB,huB,hvB,hmL,pB,uB,vB,mB,theta)
            bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
            etaB= hB+bB
            if (hB<dry_tol) then
               etaB = min(etaB,eta)
            endif

            hT = q(i_h,i,j+1)
            huT= q(i_hu,i,j+1)
            hvT= q(i_hv,i,j+1)
            hmT = q(i_hm,i,j+1)
            pT  = q(i_pb,i,j+1)
            call admissibleq(hT,huT,hvT,hmT,pT,uT,vT,mT,theta)
            bT = aux(1,i,j+1)-q(i_bdif,i,j+1)
            etaT= hT+bT
            if (hT<dry_tol) then
               etaT = min(etaT,eta)
            endif

            if (h<dry_tol) then
               eta = min(etaL,eta)
               eta = min(etaB,eta)
            endif

            if ((h+hL+hB+hR+hT)<dry_tol) then
               aux(i_taudir_x,i,j) = 0.d0
               aux(i_taudir_y,i,j) = 0.d0
               aux(i_fsphi,i,j) = 0.d0
               cycle
            endif


            pm = 0.5d0 !does not effect tau. only need tau in different cells
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            call auxeval(hL,uL,vL,mL,pL,phi,theta,kappa,S,rhoL,tanpsi,D,tauL,sigbed,kperm,compress,pm)
            call auxeval(hR,uR,vR,mR,pR,phi,theta,kappa,S,rhoR,tanpsi,D,tauR,sigbed,kperm,compress,pm)
            call auxeval(hB,uB,vB,mB,pB,phi,theta,kappa,S,rhoB,tanpsi,D,tauB,sigbed,kperm,compress,pm)
            call auxeval(hT,uT,vT,mT,pT,phi,theta,kappa,S,rhoT,tanpsi,D,tauT,sigbed,kperm,compress,pm)

            !minmod gradients
            FxC = -gmod*h*(EtaR-EtaL)/(2.d0*dx) + gmod*h*sin(theta)
            FyC = -gmod*h*(EtaT-EtaB)/(2.d0*dy)

            FxL = -gmod*0.5d0*(h+hL)*(Eta-EtaL)/(dx) + gmod*0.5d0*(h+hL)*sin(theta)
            FyL = -gmod*0.5d0*(h+hB)*(Eta-EtaB)/(dy)

            FxR = -gmod*0.5d0*(h+hR)*(EtaR-Eta)/(dx) + gmod*0.5d0*(h+hR)*sin(theta)
            FyR = -gmod*0.5d0*(h+hT)*(EtaT-Eta)/(dy)

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

            vnorm = sqrt(hu**2 + hv**2)
            if (vnorm>0.d0) then ! If moving.

               ! In D-Claw 4 dx and dy were accessible in the riemann solver.
               ! To fix in the transition, dx and dy are multiplied by taudir_x and y
               ! here. This works for everything except for bed normal and produces
               ! equivalent results (1/13/24)
               aux(i_taudir_x,i,j) = -dx*hu/vnorm
               aux(i_taudir_y,i,j) = -dy*hv/vnorm
               dot = min(max(0.d0,Fx*hu) , max(0.d0,Fy*hv))
               if (dot>0.d0) then
                  !friction should oppose direction of velocity
                  !if net force is in same direction, split friction source term
                  !splitting is useful for small velocities and nearly balanced forces
                  !only split amount up to maximum net force for large velocities
                  !aux has cell centered interpretation in Riemann solver
                  Fproj = dot/vnorm
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

end subroutine calc_taudir


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
            !---------max deta/dx-------------------
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


            !---------max deta/dy-------------------
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

      if (init_ptype==2.or.init_ptype==4) then
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

   !====================================================================
   !subroutine admissibleq
   !accept solution q, return q in admissible space
   !====================================================================

subroutine calc_pmtanh(pm,seg,pmtanh)

      implicit none

      !i/o
      double precision, intent(in) :: pm,seg
      double precision, intent(out) :: pmtanh

      !Locals


      pmtanh = seg*(0.5d0*(tanh(40.d0*(pm-0.9d0))+1.d0))
      !pmtanh = 0.8d0*seg*(0.5d0*(tanh(40.d0*(pm-0.98d0))+1.d0))

      return

   end subroutine calc_pmtanh


end module digclaw_module

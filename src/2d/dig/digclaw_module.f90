
module digclaw_module

   use geoclaw_module, only: dry_tolerance, grav, pi

   implicit none

    ! ========================================================================
    ! General digclaw parameters
    ! ========================================================================
    double precision :: rho_s,rho_f,phi_bed,theta_input,delta,kappita
    double precision :: mu,alpha,m_crit,c1,m0,alpha_seg,sigma_0,phi_seg_coeff,entrainment_rate

    integer :: init_ptype,p_initialized,bed_normal,entrainment
    double precision :: init_pmax_ratio,init_ptf2,init_ptf,init_pmin_ratio
    double precision :: grad_eta_max,cohesion_max,grad_eta_ave,eta_cell_count

    integer, parameter ::  i_dig    = 4 !Start of digclaw aux variables
    integer, parameter ::  i_phi    = i_dig
    integer, parameter ::  i_theta  = i_dig + 1
    integer, parameter ::  i_fs     = i_dig + 2
    integer, parameter ::  i_fsphi  = i_dig + 3
    integer, parameter ::  i_taudir_x = i_dig + 4
    integer, parameter ::  i_taudir_y = i_dig + 5
    integer, parameter ::  DIG_PARM_UNIT = 78


contains

    ! ========================================================================
    !  set_dig(fname)
    ! ========================================================================
    !  Reads in user parameters from the given file name if provided
    ! ========================================================================
    subroutine set_dig(fname)

         implicit none

         ! Input
         character*25, intent(in), optional :: fname

         ! Locals
         double precision :: deg2rad
         integer, parameter :: iunit = 127
         character*25 :: file_name
         logical :: found_file


         deg2rad = pi/180.d0

         ! Read user parameters from setgeo.data
         if (present(fname)) then
            file_name = fname
         else
            file_name = 'setdclaw.data'
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
         read(iunit,*) alpha_seg
         read(iunit,*) bed_normal
         read(iunit,*) phi_seg_coeff
         read(iunit,*) entrainment
         read(iunit,*) entrainment_rate

         close(iunit)
         alpha_seg = 1.0 - alpha_seg

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
         write(DIG_PARM_UNIT,*) '    alpha_seg:', alpha_seg
         write(DIG_PARM_UNIT,*) '    bed_normal:', bed_normal
         write(DIG_PARM_UNIT,*) '    phi_seg_coeff:', phi_seg_coeff
         write(DIG_PARM_UNIT,*) '    entrainment:', entrainment
         write(DIG_PARM_UNIT,*) '    entrainment_rate:', entrainment_rate


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
            file_name = 'setpinit_dclaw.data'
         endif
         inquire(file=file_name,exist=found_file)
         if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            stop
         endif

         call opendatafile(iunit, file_name)
         read(iunit,*) init_ptype
         read(iunit,*) init_pmax_ratio
         read(iunit,*) init_ptf
         read(iunit,*) init_ptf2
         close(unit=iunit)

         p_initialized = 0
         init_pmin_ratio = 1.d16
         grad_eta_max = 0.0
         cohesion_max = 0.0
         grad_eta_ave = 0.0
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
         hm = h*m0
         p  = h*gmod*rho_f
         u = 0.d0
         v = 0.d0
         m = 0.d0
         return
      endif

      u = hu/h
      v = hv/h
      m = hm/h

      !mlo = 1.d-3
      mlo = 0.0d0
      mhi = 1.d0 - mlo

      if (m.lt.mlo) then
         m = dmax1(m,mlo)
         !m = (m**2 + mlo**2)/(2.d0*mlo)
         hm = h*m
      elseif (m.gt.mhi) then
         m = dmin1(m,1.d0)
         !m = 1.d0 - ((1.d0-mhi)**2 + (1.d0-m)**2)/(2.d0*(1.d0-mhi))
         hm = h*m
      endif

      rho = rho_s*m + (1.d0-m)*rho_f
      pmax = rho*gmod*h
      plo = rho_f*dry_tol*gmod*dry_tol
      phi = pmax - plo
      if (p.lt.plo) then
         if ((u**2+v**2)>0.d0) then
            p = dmax1(0.d0,p)
            ! WHY not when static? 
         endif
         !p = dmax1(-5.0*pmax,p)
         !p = (p**2 + plo**2)/(2.d0*plo)
      elseif (p.gt.phi) then
         p = dmin1(pmax,p)
         !p = pmax - ((pmax-p)**2+ (pmax-phi)**2)/(2.d0*(pmax-phi))
      endif

      return

   end subroutine admissibleq

   !====================================================================
   ! subroutine auxeval: evaluates the auxiliary variables as functions
   !                     of the solution vector q
   !====================================================================

   subroutine auxeval(h,u,v,m,p,phi_bed,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

      implicit none

      !i/o
      double precision, intent(inout) :: pm
      double precision, intent(in)  :: h,u,v,m,p,phi_bed,theta
      double precision, intent(out) :: S,rho,tanpsi,D,tau,kappa
      double precision, intent(out) :: sigbed,kperm,compress

      !local
      double precision :: m_eqn,vnorm,gmod,sigbedc,hbounded,shear,tanphi,rho_fp
      double precision :: seg,pmtanh01,m_crit_m,m_crit_pm

      if (h.lt.dry_tolerance) return

      hbounded = h!max(h,0.1)
      gmod=grav
      pm = max(0.0d0,pm)
      pm = min(1.0d0,pm)
      if (dabs(alpha_seg-1.0d0)<1.d-6) then
         seg = 0.0d0
         rho_fp = rho_f
         pmtanh01=0.0d0
      else
         seg = 1.0d0
         call calc_pmtanh(pm,seg,pmtanh01)
         rho_fp = (1.0d0-pmtanh01)*rho_f
      endif
      !pmtanh01 = seg*(0.5*(tanh(20.0*(pm-0.80))+1.0))
      !pmtanh01 = seg*(0.5*(tanh(40.0*(pm-0.90))+1.0))

      if (bed_normal.eq.1) gmod=grav*dcos(theta)
      vnorm = dsqrt(u**2 + v**2)
      rho = rho_s*m + rho_fp*(1.d0-m)
      shear = 2.0d0*vnorm/hbounded
      sigbed = dmax1(0.d0,rho*gmod*h - p)
      sigbedc = rho_s*(shear*delta)**2 + sigbed
      if (sigbedc.gt.0.0d0) then
         S = (mu*shear/(sigbedc))
      else
         S = 0.d0
      endif
      !Note: m_eqn = m_crit/(1+sqrt(S))
      !From Boyer et. al
      !S = 0.0
      !m_eqn = m_crit/(1.d0 + sqrt(S))
      !if (m.gt.m_eqn) write(*,*) 'm,m_eqn,S:',m,m_eqn,S,sigbed,shear
      !tanpsi = c1*(m-m_eqn)*tanh(shear/0.1)
      !pmlin = seg*2.0*(pm-0.5)
      !pmtan = seg*0.06*(tan(3.*(pm-0.5)))
      !pmtanh = seg*tanh(3.*pmlin)
      !pmtanh01 = seg*0.5*(tanh(8.0*(pm-0.75))+1.0)
      !pmtanh01 = seg*0.5*(tanh(20.0*(pm-0.80))+1.0)
      !pmtanh01s = seg*4.0*(tanh(8.0*(pm-0.95))+1.0)

   
      kperm = kappita*exp(-(m-m0)/(0.04d0))!*(10**(pmtanh01))
      !m_crit_pm - max(pm-0.5,0.0)*(0.15/0.5) - max(0.5-pm,0.0)*(0.15/0.5)
      !m_crit_pm =  max(pm-0.7,0.0)*((m_crit- 0.55)/0.5) + max(0.3-pm,0.0)*((m_crit-0.55)/0.5)
      m_crit_pm =  0.d0! max(pm-0.6,0.0)*((m_crit- 0.55)/0.4) + max(0.3-pm,0.0)*((m_crit-0.55)/0.3)
      !m_crit_pm = max(pm-0.9,0.0)*((m_crit- 0.55)/0.1) + max(0.1-pm,0.0)*((m_crit-0.55)/0.1);

      m_crit_pm = pmtanh01*0.09d0
      m_crit_m = m_crit - m_crit_pm
      m_eqn = m_crit_m/(1.d0 + sqrt(S))
      tanpsi = c1*(m-m_eqn)*tanh(shear/0.1)

      !kperm = kperm + 1.0*pm*kappita
      !compress = alpha/(sigbed + 1.d5)

      !if (m.le.0.04) then
          !eliminate coulomb friction for hyperconcentrated/fluid problems
          !klugey, but won't effect debris-flow problems
         !sigbed = sigbed*0.5d0*(tanh(400.d0*(m-0.02d0)) + 1.0d0)
      !endif

      if (m.le.1.d-16) then
         compress = 1.d16
         kperm = 0.0d0
         tanpsi = 0.0d0
         sigbed=0.0d0
      else
         compress = alpha/(m*(sigbed +  sigma_0))
      endif

      !if (m.le.0.4) then
      !   kperm = tanh(10.d0*m)*kperm
      !   tanpsi = tanh(10.d0*m)*tanpsi
      !   sigbed= tanh(10.d0*m)*sigbed
      !endif

      if (p_initialized.eq.0.and.vnorm.le.0.d0) then
      !if (vnorm.le.0.d0) then
         tanpsi = 0.d0
         D = 0.d0
      elseif (h*mu.gt.0.d0) then
         D = 2.0d0*(kperm/(mu*h))*(rho_fp*gmod*h - p)
      else
         D = 0.d0
      endif

      tanphi = dtan(phi_bed + datan(tanpsi))! + phi_seg_coeff*pmtanh01*dtan(phi_bed)
      !if (S.gt.0.0) then
      !   tanphi = tanphi + 0.38*mu*shear/(shear + 0.005*sigbedc)
      !endif

      tau = dmax1(0.d0,sigbed*tanphi)

      !tau = (grav/gmod)*dmax1(0.d0,sigbed*tanphi)
      !kappa: earth pressure coefficient
      !if (phi_int.eq.phi_bed) then
      !   sqrtarg = 0.d0
      !else
      !   sqrtarg = 1.d0-(dcos(phi_int)**2)*(1.d0 + dtan(phi_bed)**2)
      !endif

      !kappa = (2.d0 - pm*2.d0*dsqrt(sqrtarg))/(dcos(phi_int)**2)
      !kappa = kappa - 1.d0
      kappa = 1.d0
      !kappa = 0.4d0

   end subroutine auxeval


   !====================================================================
   !subroutine psieval: evaluate the source term
   !====================================================================

   subroutine psieval(tau,rho,D,tanpsi,kperm,compress,h,u,m,psi)

      implicit none

      !i/o
      double precision, intent(out) :: psi(4)
      double precision, intent(in)  :: tau,rho,D,tanpsi,kperm,compress
      double precision, intent(in)  :: h,u,m

      !local
      double precision :: taushear,drytol,vnorm

      drytol = dry_tolerance

      taushear = (tau/rho)*dsign(1.d0,u)
      vnorm = dabs(u)
      if (h.lt.drytol.or..true.) then
         psi(1) = 0.d0
         psi(2) = 0.d0
         psi(3) = 0.d0
         psi(4) = 0.d0
      else
         psi(1) =  D*(rho-rho_f)/rho
         psi(2) =  u*D*(rho-rho_f)/rho
         psi(3) = -D*m*(rho_f/rho)
         psi(4) = 0.d0
      endif

   end subroutine psieval

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


      do i=2-mbc,mx+mbc-1
         do j=2-mbc,my+mbc-1


            h = q(1,i,j)
            hu = q(2,i,j)
            hv = q(3,i,j)
            if (h<dry_tol) then
               hu=0.d0
               hv=0.d0
            endif
            hm = q(4,i,j)
            p  = q(5,i,j)
            b = aux(1,i,j)
            eta = h+b
            phi = aux(i_phi,i,j)

            hL = q(1,i-1,j)
            huL= q(2,i-1,j)
            hvL= q(3,i-1,j)
            hmL = q(4,i-1,j)
            pL  = q(5,i-1,j)
            bL = aux(1,i-1,j)
            etaL= hL+bL
            if (hL<dry_tol) then
               etaL = min(etaL,eta)
            endif

            hR = q(1,i+1,j)
            huR= q(2,i+1,j)
            hvR= q(3,i+1,j)
            hmR = q(4,i+1,j)
            pR  = q(5,i+1,j)
            bR = aux(1,i+1,j)
            etaR= hR+bR
            if (hR<dry_tol) then
               etaR = min(etaR,eta)
            endif

            hB = q(1,i,j-1)
            huB= q(2,i,j-1)
            hvB= q(3,i,j-1)
            hmB = q(4,i,j-1)
            pB  = q(5,i,j-1)
            bB = aux(1,i,j-1)
            etaB= hB+bB
            if (hB<dry_tol) then
               etaB = min(etaB,eta)
            endif

            hT = q(1,i,j+1)
            huT= q(2,i,j+1)
            hvT= q(3,i,j+1)
            hmT = q(4,i,j+1)
            pT  = q(5,i,j+1)
            bT = aux(1,i,j+1)
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

            pm = 0.50d0 !does not effect tau. only need tau in different cells
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call admissibleq(hL,huL,hvL,hmL,pL,uL,vL,mL,theta)
            call admissibleq(hB,huB,hvB,hmL,pB,uB,vB,mB,theta)
            call admissibleq(hR,huR,hvR,hmR,pR,uR,vR,mR,theta)
            call admissibleq(hT,huT,hvT,hmT,pT,uT,vT,mT,theta)


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
            if (vnorm>0.d0) then
               !DIG: check to see if adding dx and dy works everywhere. added in 4.x-->5.0
               aux(i_taudir_x,i,j) = -dx*hu/sqrt(hv**2+hu**2)
               aux(i_taudir_y,i,j) = -dy*hv/sqrt(hv**2+hu**2)

               dot = min(max(0.d0,Fx*hu) , max(0.d0,Fy*hv))
               if (dot>0.0) then
                  !friction should oppose direction of velocity
                  !if net force is in same direction, split friction source term
                  !splitting is useful for small velocities and nearly balanced forces
                  !only split amount up to maximum net force for large velocities
                  !aux has cell centered interpretation in Riemann solver
                  Fproj = dot/vnorm
                  aux(i_fsphi,i,j) = min(1.d0,Fproj*rho/max(tau,1.d-16))
               else
                  !net force is in same direction as friction
                  !if nearly balanced steady state not due to friction
                  !no splitting, integrate friction in src
                  aux(i_fsphi,i,j) = 0.0
               endif


            else
               !aux now have cell edge interpretation in Riemann solver
               !friction should oppose net force. resolve in Riemann solver
               !DIG: check to see if adding dx and dy works everywhere. added in 4.x-->5.0
               if ((FxL**2+Fy**2)>0.d0) then
                  aux(i_taudir_x,i,j) = -dx*FxL/sqrt(FxL**2+Fy**2)
               else
                  aux(i_taudir_x,i,j) = dx*1.d0
               endif

               if ((Fx**2+FyL**2)>0.0) then
                  aux(i_taudir_y,i,j) = -dy*FyL/sqrt(Fx**2+FyL**2)
               else
                  !there is no motion or net force. resolve in src after Riemann
                  aux(i_taudir_y,i,j) = dy*1.d0
               endif

               if ((aux(i_taudir_y,i,j)**2 + aux(i_taudir_x,i,j)**2)>0.0) then
                  aux(i_fsphi,i,j) = 1.d0
               else
                  aux(i_fsphi,i,j) = 0.d0
               endif
            endif

         enddo
      enddo

end subroutine calc_taudir

   ! ========================================================================
   !  calc_tausplit
   ! ========================================================================
   !  Determines splitting of tau for rp vs. src.
   ! ========================================================================

subroutine calc_tausplit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)


      implicit none

      !Input
      double precision :: dx,dy,xlower,ylower
      double precision :: q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      integer :: mx,my,mbc,meqn,maux

      !Locals
      double precision :: h,hL,hR,hu,hv,hm,p,b,bL,bR,bT,bB,hT,hB,u,v,m
      double precision :: phi,theta,rho,kappa,S,tanpsi,D,tau,sigbed,kperm,compress,pm
      double precision :: gmod,dry_tol
      double precision :: EtaL,EtaR,EtaT,EtaB,Eta
      double precision :: detadx,detadxL,detadxR,detady,detadyT,detadyB
      double precision :: grad_eta


      integer :: i,j

      dry_tol = dry_tolerance
      gmod = grav
      rho = m0*rho_s + (1.d0-m0)*rho_f

      do i=2-mbc,mx+mbc-1
         do j=2-mbc,my+mbc-1

            h = q(1,i,j)
            hL = q(1,i-1,j)
            hR = q(1,i+1,j)
            if (h<dry_tol) then
               aux(i_fsphi,i,j) = 1.d0
               cycle
            endif

            hu = q(2,i,j)
            hv = q(3,i,j)
            hm = q(4,i,j)
            p  = q(5,i,j)

            if ((hu**2 + hv**2)==0.d0) then
               aux(i_fsphi,i,j) = 1.0
               !cycle
            endif

            b = aux(1,i,j)
            bR = aux(1,i+1,j)
            bL = aux(1,i-1,j)
            phi = aux(i_phi,i,j)

            hT = q(1,i,j+1)
            bT = aux(1,i,j+1)
            hB = q(1,i,j-1)
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

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            pm = q(6,i,j)/q(1,i,j)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            if (tau>0.0) then
               aux(i_fsphi,i,j) = min(1.d0,grad_eta*rho*gmod*h/tau)
            else
               aux(i_fsphi,i,j) = 1.d0
            endif
         enddo
      enddo

   end subroutine calc_tausplit

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

      do i=1,mx
         do j=1,my

            h = q(1,i,j)
            hL = q(1,i-1,j)
            hR = q(1,i+1,j)
            if (h<dry_tol) then
               aux(i_fs,i,j) = 10.d0
               cycle
            endif

            hu = q(2,i,j)
            hv = q(3,i,j)

            if ((hu**2+hv**2)>0.d0) then
               aux(i_fs,i,j) = 0.d0
               cycle
            endif

            b = aux(1,i,j)
            bR = aux(1,i+1,j)
            bL = aux(1,i-1,j)
            phi = aux(i_phi,i,j)

            if ((phi)==0.0) then
               aux(i_fs,i,j) = 0.d0
               init_pmin_ratio = 0.d0
               cycle
            endif

            hT = q(1,i,j+1)
            bT = aux(1,i,j+1)
            hB = q(1,i,j-1)
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

            if (grad_eta>0.d0) then
               aux(i_fs,i,j) = tan(phi)/grad_eta
            else
               aux(i_fs,i,j) = 10.d0
            endif
         enddo
      enddo

      if (init_ptype==2.or.init_ptype==4) then
         init_pmin_ratio = 1.d0-(grad_eta_ave/eta_cell_count)
      endif
      if (init_ptype>0) then
         write(*,*) '--------------------------------------------'
         write(*,*) 'hydrostatic liquefaction ratio:', rho_f/rho
         write(*,*) 'initiation liquefaction  ratio:',init_pmin_ratio, grad_eta_ave
         write(*,*) 'maximum surface slope angle:',180.*atan(tan(phi)*grad_eta_max)/3.14, grad_eta_max
         write(*,*) 'average failure liquefaction ratio:', 1.0-(grad_eta_ave/eta_cell_count) , eta_cell_count
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

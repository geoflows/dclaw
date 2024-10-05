

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module, only: grav, dry_tolerance,deg2rad,friction_depth
      use geoclaw_module, only: manning_coefficient,friction_forcing

      use digclaw_module, only: alpha_c,beta_seg,bed_normal,curvature
      use digclaw_module, only: entrainment,entrainment_rate,entrainment_method
      use digclaw_module, only: src2method
      use digclaw_module, only: phi,segregation
      use digclaw_module, only: i_ent,i_fsphi,i_phi,i_theta
      use digclaw_module, only: mu,rho_f,rho_s, sigma_0
      use digclaw_module, only: admissibleq,auxeval
      use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_bdif
      use digclaw_module, only: qfix,setvars
      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

      ! Output
      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gz,gx,h,hu,hv,hm,u,v,m,p
      real(kind=8) :: rhoh
      real(kind=8) :: b,bR,bL,bT,bB,bTR,bTL,bBR,bBL
      real(kind=8) :: kappa,S,rho,tanpsi
      real(kind=8) :: D,tau,sigbed,kperm,compress,chi,coeffmanning
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,taucf,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg,m_eq
      real(kind=8) :: b_xx,b_yy,b_xy,hchi,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining
      real(kind=8) :: gamma,zeta,krate,p_eq,dgamma,sig_0,sig_eff
      real(kind=8) :: dtk,dtremaining,alphainv,gacc


      ! TO REMOVE
      double precision :: alpha,kappita,phi_bed,alpha_seg
      ! END TO REMOVE

      integer :: i,j,ii,jj,icount,itercount,itercountmax

      logical :: debug
      debug = .false.


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

      if (friction_forcing) then
         coeffmanning = manning_coefficient(1)
      else
         coeffmanning = 0.d0
      endif

      ! Current implementation of friction has manning as an array
      ! take the first element for now. If only one value is
      ! provided to geo_data.manning_coefficient
      ! it will be manning_coefficient(1)
      ! DIG: Decide if this should be handled in some other way.


      gz = grav  !needed later for bed-normal direction gravity
      gx = 0.d0
      theta=0.d0

      do j=1,my
         do i=1,mx
         ! DIG: 1/12/24: KRB and MJB notice that here we are looping over ghost cells.
         ! These ghost cells have not been updated by the riemann solver? Are they used
         ! meaningfully (e.g., theta is diffed below.)
         ! 1/30/2024 - Leaving this as is for the moment, this is something to evaluate later.
         ! DIG: 10/3/24: DLG is eliminating ghost cells from loop.

            ! Get state variable
            h = q(i_h,i,j)
            if (h<=dry_tolerance) cycle
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p =  q(i_pb,i,j)
            hchi = q(i_hchi,i,j)
            rhoh = hm*rho_s + (h-hm)*rho_f
            call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gz = grav*cos(theta)
               gx = grav*sin(theta)
            endif
            if (curvature==1.or.segregation==1) then
               b = aux(1,i,j)-q(i_bdif,i,j)
               bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
               bR = aux(1,i+1,j)-q(i_bdif,i+1,j)
               bT = aux(1,i,j+1)-q(i_bdif,i,j+1)
               bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
               bTR = aux(1,i+1,j+1)-q(i_bdif,i+1,j+1)
               bTL = aux(1,i-1,j+1)-q(i_bdif,i-1,j+1)
               bBR = aux(1,i+1,j-1)-q(i_bdif,i+1,j-1)
               bBL = aux(1,i-1,j-1)-q(i_bdif,i-1,j-1)
               b_x = (bR-bL)/2.d0*dx
               b_y = (bT-bB)/2.d0*dy
               b_xx=(bR - 2.d0*b + bL)/(dx**2)
               b_yy=(bT - 2.d0*b + bB)/(dy**2)
               b_xy=((bTR-bTL) - (bBR-bBL))/(4.0*dx*dy)
            endif
            if (curvature==1) then
               dtheta = -(aux(i_theta,i+1,j) - aux(i_theta,i-1,j))/(2.d0*dx)
               gacc = max(u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy + u**2*dtheta,0.d0)!max:currently only consider enhancement not reduction of gz (ie. basin not a hump)
               gz = gz + gacc
            endif

            !Manning friction
            if ((friction_forcing).and.(coeffmanning>0.d0)) then
               if (h<=friction_depth) then
                  beta = 1.0d0-m
                  gamma= beta*sqrt(hu**2 + hv**2)*(gz*coeffmanning**2)/(h**(7.d0/3.d0))
                  dgamma=1.d0 + dt*gamma
                  hu= hu/dgamma
                  hv= hv/dgamma
                  !new u,v below
                  call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)
                  !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
               endif
            endif

!-----------!integrate momentum source term------------------------
            ! need tau:
            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)

!-----------! changes only hu,hv,u,v ------------------------------
            !integrate momentum source term
            hvnorm0 = sqrt(hu**2 + hv**2)
            vnorm = hvnorm0/h

            if (hvnorm0>0.d0) then
               !integrate dynamic friction !DIG: TO DO - move dynamic friction to Riemann solver
               vnorm = max(0.d0,vnorm - dt*tau/rhoh) !exact solution for Coulomb friction
               vnorm = vnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(h*rhoh)) !exact solution (prior to h change) for effective viscous friction
               ! velocity determined, calculate directions etc. from vnorm
               hvnorm = h*vnorm
               hu = hvnorm*hu/hvnorm0 + gx*h*dt !gx=0 unless bed-normal !DIG: last term should ultimately be in Riemann solver
               hv = hvnorm*hv/hvnorm0
               u = hu/h
               v = hv/h
               vnorm = sqrt(u**2 + v**2)
               ! velocity now constant for remainder of src2. hu,hv adjusted due to change in h
            endif
!-------------------------------------------------------------------------
            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
!----------- ! integrate p  & m-------------------------------------------

            select case (src2method)

            case(0:1)
               call mp_update_relax_Dclaw4(dt,h,u,v,m,p,chi,rhoh,gz)
               hu = h*u
               hv = h*v
               hm = h*m
               hchi = h*chi
            case(2)
            ! changes only p,m,hm,& h. hrho constant
            ! takes in general multiple interior timesteps, dtk
            ! sum(dtk) = dt
            ! explicit integration for each dt
               dtremaining = dt
               itercountmax=10
               itercount=0

               do while (dtremaining>1.d-99)
                  call mp_update_FE_4quad(dtremaining,h,u,v,m,p,chi,rhoh,gz,dtk)
                  dtremaining = dtremaining-dtk
                  itercount = itercount + 1
                  if (dtk==0.d0.and.(.not.debug)) exit
                  if (h<dry_tolerance) exit
                  if (itercount>=itercountmax) then
                     exit
                  endif
                  if (debug.and.(itercount>50)) then
                     write(*,*) 'WARNING SRC2: update iteration'
                     write(*,*) 'itercount,dt,dtremaining:',itercount,dt,dtremaining
                  endif
               enddo

               hu = h*u
               hv = h*v
               hm = h*m
               hchi = h*chi
               call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)
               if (h<=dry_tolerance) then
                  cycle
               endif

            end select
            !========================== end src integration ======================

            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)

            !======================mass entrainment===========================
            if (entrainment==1) then
               b_eroded = q(i_bdif,i,j)
               b_remaining = aux(i_ent,i,j)-b_eroded
               select case(entrainment_method)
               case(0)
                  call ent_dclaw4(dt,h,u,v,m,p,rho,hchi,gz,b_x,b_y,b_eroded,b_remaining)
                  q(i_bdif,i,j) = b_eroded
               case(1)
                  !do nothing yet
               end select
            endif

            !===================================================================
            ! end of entrainment, put state variables back in q.

            q(i_h,i,j) = h
            q(i_hu,i,j) = hu
            q(i_hv,i,j) = hv
            q(i_hm,i,j) = hm
            q(i_pb,i,j) = p
            q(i_hchi,i,j) = hchi

         enddo
      enddo

      return
      end

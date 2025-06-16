!================================================================
! Integrate right-hand side (source term) of the D-Claw equations
! on one grid at one AMR level.
!
! This subroutine is called by stepgrid to integrate the right-
! hand side of the equations. It is called after the Riemann solver
! integrates the left-hand side.
!
! Variables
! ---------
! meqn, number of equations
! mbc, number of ghost cells
! mx, number of cells in x-direction
! my, number of cells in y-direction
! xlower, x-coordinate of lower left corner
! ylower, y-coordinate of lower left corner
! dx, grid cell size in x-direction
! dy, grid cell size in y-direction
! q, solution state array
! maux, number of auxilary variables
! aux, auxilary array
! t, current time
! dt, time step
!=================================================================

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================

      use geoclaw_module, only: grav, dry_tolerance,friction_depth
      use geoclaw_module, only: manning_coefficient,friction_forcing

      use digclaw_module, only: bed_normal,curvature
      use digclaw_module, only: entrainment,entrainment_method
      use digclaw_module, only: src2method
      use digclaw_module, only: i_ent,i_theta,i_dhdt
      use digclaw_module, only: mu,rho_f,rho_s
      use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_hs,i_hf
      use digclaw_module, only: qfix,setvars
      use digclaw_module, only: dd_manning,manning_max
      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

      ! Output
      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gz,gx,h,hu,hv,hm,u,v,m,p,chi,hs,hf
      real(kind=8) :: b,bR,bL,bT,bB,bTR,bTL,bBR,bBL

      real(kind=8) :: rhoh,hchi
      real(kind=8) :: rho,tanpsi
      real(kind=8) :: tau,kperm
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,hvnorm0
      real(kind=8) :: m_eq
      real(kind=8) :: b_x,b_y,b_xx,b_yy,b_xy
      real(kind=8) :: beta,coeffmanning
      real(kind=8) :: gamma,dgamma
      real(kind=8) :: dtk,dtremaining,alphainv,gacc

      integer :: i,j,itercount,itercountmax

      logical :: debug
      debug = .false.


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

      ! Current implementation of friction has manning as an array
      ! take the first element for now. If only one value is
      ! provided to geo_data.manning_coefficient
      ! it will be manning_coefficient(1)
      ! DIG: Decide if this should be handled in some other way.

      if (friction_forcing) then
         coeffmanning = manning_coefficient(1)
         ! if depth dependent manning, will reset coefficient later based on depth.
      else
         coeffmanning = 0.d0
      endif

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

            ! Get state variables and evaluate rain.
            h = q(i_h,i,j)

            if (aux(i_dhdt,i,j).gt.0.d0) then
               ! add rain to hf
               write(*,*) 'src2: raining', dt*aux(i_dhdt,i,j)
               q(i_hf,i,j) = q(i_hf,i,j) + (dt*aux(i_dhdt,i,j))
            endif

            hf = q(i_hf,i,j)
            hs = q(i_hs,i,j)

            ! if h + (hf-hs) > 2*drytol, move hf-hs to h
            ! hm does not need adjusting as rain has m=0
            if ((h+(hf-hs)>2.0d0*dry_tolerance).and.(hf-hs>0.d0)) then

               write(*,*) 'src2: moving hf to h: hf=', hf,'h=', h

               h = h + (hf-hs)
               hf = hs
               q(i_hf,i,j) = hf
            endif

            if (h<=dry_tolerance) cycle
            hu = q(i_hu,i,j)
            hv = q(i_hv,i,j)
            hm = q(i_hm,i,j)
            p =  q(i_pb,i,j)
            hchi = q(i_hchi,i,j)
            rhoh = hm*rho_s + (h-hm)*rho_f

            if (hm.lt.0.d0) then
              !write(*,*) 'SRC2: warning: hm negative, h, hm, rhoh: ',h, hm, rhoh
              ! This most likely occurs using second order corrections.
              !stop
              cycle
            endif

            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gz = grav*dcos(theta)
               gx = grav*dsin(theta)
            endif
            if (curvature==1.or.entrainment==1) then
               b = aux(1,i,j)       -aux(i_ent,i,j)     +q(i_hs,i,j)
               bL = aux(1,i-1,j)    -aux(i_ent,i-1,j)   +q(i_hs,i-1,j)
               bR = aux(1,i+1,j)    -aux(i_ent,i+1,j)   +q(i_hs,i+1,j)
               bT = aux(1,i,j+1)    -aux(i_ent,i,j+1)   +q(i_hs,i,j+1)
               bB = aux(1,i,j-1)    -aux(i_ent,i,j-1)   +q(i_hs,i,j-1)
               bTR = aux(1,i+1,j+1) -aux(i_ent,i+1,j+1) +q(i_hs,i+1,j+1)
               bTL = aux(1,i-1,j+1) -aux(i_ent,i-1,j+1) +q(i_hs,i-1,j+1)
               bBR = aux(1,i+1,j-1) -aux(i_ent,i+1,j-1) +q(i_hs,i+1,j-1)
               bBL = aux(1,i-1,j-1) -aux(i_ent,i-1,j-1) +q(i_hs,i-1,j-1)
               b_x = (bR-bL)/2.d0*dx
               b_y = (bT-bB)/2.d0*dy
               b_xx=(bR - 2.d0*b + bL)/(dx**2)
               b_yy=(bT - 2.d0*b + bB)/(dy**2)
               b_xy=((bTR-bTL) - (bBR-bBL))/(4.0*dx*dy)
               ! entrainment only needs b_x and b_y.
               ! curvature also needs b_xx,b_yy,and b_xy
            endif

            if (curvature==1) then
               u = hu/h
               v = hv/h
               dtheta = -(aux(i_theta,i+1,j) - aux(i_theta,i-1,j))/(2.d0*dx)
               gacc = max(u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy + u**2*dtheta,0.d0)!max:currently only consider enhancement not reduction of gz (ie. basin not a hump)
               gz = gz + gacc
            endif

            call qfix(h,hu,hv,hm,p,hchi,hf,u,v,m,chi,rho,gz)

!-----------!integrate momentum source term------------------------
            ! need tau:
            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)

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
            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)
!----------- ! integrate p  & m-------------------------------------------

            select case (src2method)

            case(-1)
            !ignore source integration of m and p
            !solve shallow water equations with possible friction and passive advection of m
            !could be used for sediment transport/entrainment with shallow water
               hu = h*u
               hv = h*v
               hm = h*m
               hchi = h*chi
            case(0:1)
               call mp_update_relax_Dclaw4(dt,h,u,v,m,p,chi,hf,rhoh,gz)
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
               itercountmax=100
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
               call qfix(h,hu,hv,hm,p,hchi,hf,u,v,m,chi,rho,gz)

            end select

            if (h<=dry_tolerance) then
            ! put state variables back in q.
               q(i_h,i,j) = h
               q(i_hu,i,j) = hu
               q(i_hv,i,j) = hv
               q(i_hm,i,j) = hm
               q(i_pb,i,j) = p
               q(i_hchi,i,j) = hchi
               q(i_hs,i,j) = hs
               q(i_hf,i,j) = hf
               cycle
            endif

            !========================== end src integration ======================

            call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)


            !======================mass entrainment===========================
            if (entrainment==1) then
               hs = q(i_hs,i,j)
               select case(entrainment_method)
               case(0)
                  call ent_dclaw4(dt,h,u,v,m,p,rho,hchi,gz,tau,b_x,b_y,hs)
                  q(i_hs,i,j) = hs
               case(1)
                  !do nothing yet
               end select
            endif

            !===================================================================
            ! end of entrainment


            !===========Manning friction========================================
            if ((friction_forcing).and.(coeffmanning>0.d0)) then

               if (h<=friction_depth) then

                  if (dd_manning) then
                     coeffmanning = manning_max + (manning_coefficient(1)-manning_max)*0.5D0*(1.d0+tanh(300.d0*(h-0.05d0)))
                  endif 


                  !beta = 1.d0-m  ! reduced friction led to high velocities
                  beta = 1.d0     ! use full Manning friction
                  gamma= beta*sqrt(hu**2 + hv**2)*(gz*coeffmanning**2)/(h**(7.d0/3.d0))
                  dgamma=1.d0 + dt*gamma
                  hu= hu/dgamma
                  hv= hv/dgamma
                  !new u,v below
                  call qfix(h,hu,hv,hm,p,hchi,hf,u,v,m,chi,rho,gz)
               endif
            endif

            !===================================================================
            ! end of Manning friction

            ! put state variables back in q.
            q(i_h,i,j) = h
            q(i_hu,i,j) = hu
            q(i_hv,i,j) = hv
            q(i_hm,i,j) = hm
            q(i_pb,i,j) = p
            q(i_hchi,i,j) = hchi
            q(i_hs,i,j) = hs
            q(i_hf,i,j) = hf

         enddo
      enddo

      return
      end

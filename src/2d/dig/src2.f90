

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module, only: grav, dry_tolerance,deg2rad,friction_depth
      use geoclaw_module, only: manning_coefficient,friction_forcing

      use digclaw_module, only: alpha,alpha_seg,bed_normal,curvature
      use digclaw_module, only: entrainment,entrainment_rate,phi_bed
      use digclaw_module, only: i_ent,i_fsphi,i_phi,i_theta
      use digclaw_module, only: mu,rho_f,rho_s, sigma_0
      use digclaw_module, only: admissibleq,auxeval
      use digclaw_module, only: calc_pmtanh
      use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_bdif

      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt

      ! Output
      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gz,gx,h,hu,hv,hm,u,v,m,p
      real(kind=8) :: b,bR,bL,bT,bB,bTR,bTL,bBR,bBL
      real(kind=8) :: phi,kappa,S,rho,tanpsi
      real(kind=8) :: D,tau,sigbed,kperm,compress,chi,coeffmanning
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,taucf,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,hchi,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining
      real(kind=8) :: gamma,zeta,krate,p_eq,dgamma

      integer :: i,j,ii,jj,icount
      

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
      phi = phi_bed

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
            call qfix(h,hu,hv,hm,p,u,v,m,rho,gz)
            ! DIG: 10/3/24: DLG: not sure need chi check below...should be checked somewhere
            ! else where hchi or chi is computed. Or in the physical check of q (qfix or admissible q)
            !chi = max(0.0d0,chi)
            !chi = min(1.0d0,chi)
           
            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gz = grav*cos(theta)
               gx = grav*sin(theta)
            endif
            if (curvature==1) then
               b = aux(1,i,j)-q(i_bdif,i,j)
               bL = aux(1,i-1,j)-q(i_bdif,i-1,j)
               bR = aux(1,i+1,j)-q(i_bdif,i+1,j)
               bT = aux(1,i,j+1)-q(i_bdif,i,j+1)
               bB = aux(1,i,j-1)-q(i_bdif,i,j-1)
               bTR = aux(1,i+1,j+1)-q(i_bdif,i+1,j+1)
               bTL = aux(1,i-1,j+1)-q(i_bdif,i-1,j+1)
               bBR = aux(1,i+1,j-1)-q(i_bdif,i+1,j-1)
               bBL = aux(1,i-1,j-1)-q(i_bdif,i-1,j-1)
               b_xx=(bR - 2.d0*b + bL)/(dx**2)
               b_yy=(bT - 2.d0*b + bB)/(dy**2)
               b_xy=((bTR-bTL) - (bBR-bBL))/(4.0*dx*dy)
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
                  call qfix(h,hu,hv,hm,p,u,v,m,rho,gz)
                  !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
               endif
            endif

c-----------!integrate momentum source term------------------------
c-----------! changes only hu,hv,u,v ------------------------------
            ! need tau:
            call setvars(h,u,v,m,p,gz,rho,kperm,alphainv,sig_0,sig_eff,m_eq,tanpsi,tau)
           
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
c-------------------------------------------------------------------------

c----------- ! integrate p  & m-------------------------------------------
            select case (src2method)

            case(0)
               call mp_update_Dclaw4()

            case(1)
               call mp_update_relax()

            case(2)
            ! changes only p,m,hm,& h. hrho constant
            ! takes in general multiple interior timesteps, dtk
            ! sum(dtk) = dt
            ! explicit integration for each dt
               dtremaining = dt
               itercountmax=10
               itercount=0

               do while (dtremaining>1.d-99)
                  call mp_update_FE_4quad(dtremaining,h,u,v,m,p,rhoh,gz,dtk)
                  dtremaining = dtremaining-dtk
                  itercount = itercount + 1
                  if (dtk==0.d0.and.(.not.debug)) exit
                  if (h<drytolerance) exit
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
               call qfix(h,hu,hv,hm,p,u,v,m,rho,gz)
               if (h<=drytolerance) then
                  cycle
               endif

            case(2)

            end select
            
c------------------------------------------------------------------------

c
            ! call admissible q and auxeval before moving on to shear induced
            ! dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

            ! calculate velocity
            vnorm = sqrt(u**2 + v**2)

            !integrate shear-induced dilatancy
            shear = 2.d0*vnorm/h
            krate = 1.5d0*shear*m*tanpsi/alpha
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dt*3.d0*vnorm*tanpsi/(h*compress)
            endif

            ! update pmtanh01 and rho_fp for segregation
            ! DIG: if segregation is compartmentalized
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               pmtanh01=0.d0
      		else
         		seg = 1.d0
               call calc_pmtanh(chi,seg,pmtanh01)
               rho_fp = max(0.d0,(1.d0-pmtanh01))*rho_f
      		endif

            ! integrate pressure relaxation
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.d0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
            p_eq = h*rho_fp*gmod
            p = p_eq + (p-p_eq)*exp(krate*dt)

            ! call admissible q and auxeval before moving on to dilatancy.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

            ! calculate rate of change
            krate = D*(rho-rho_fp)/rho

            ! integrate hu, hv, hm, and h.
            hu = hu*exp(dt*krate/h)
            hv = hv*exp(dt*krate/h)
            hm = hm*exp(-dt*D*rho_fp/(h*rho))
            h = h + krate*dt


            !======================mass entrainment===========================

            ! call admissible q and auxeval before moving on to mass entrainment.
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

            vnorm = sqrt(u**2 + v**2)
            vlow = 0.1d0 ! minimum velocity for entrainment to occur. ! DIG: should this be a user
            ! specified variable.

            if (ent.and.vnorm.gt.vlow) then
               if (aux(i_ent,i,j)>0.d0) then

                  ! calculate the basal surface (aux(1)-q(i_bdif)) gradients in x and y
                  ! to determine dbdv, the slope in the direction of flow by taking the dot
                  ! product of the slope vector and the unit vector in the direction of flow.

                  b_x = (aux(1,i+1,j)-q(i_bdif,i+1,j)-aux(1,i-1,j)+q(i_bdif,i-1,j))/(2.d0*dx)
                  b_y = (aux(1,i,j+1)-q(i_bdif,i,j+1)-aux(1,i,j-1)+q(i_bdif,i,j-1))/(2.d0*dy)
                  dbdv = (u*b_x+v*b_y)/vnorm
                  slopebound = 1.d10
                  b_eroded = q(i_bdif,i,j)

                  ! erode if material to erode is still available and the slope is
                  ! less than a critical slope value. !DIG: Improve critical slope value.

                  if (dbdv<slopebound.and.b_eroded<aux(i_ent,i,j)) then

                     ! calculate remaining material that may be eroded.
                     b_remaining = aux(i_ent,i,j)-b_eroded

                     ! value for m for entrained material
                     m2 = 0.6d0
                     ! DIG: eventually make this a user defined variable or an aux value.
                     ! DIG: should there also be a substrate value for chi?

                     ! calculate entrained material density, using hard coded values
                     ! for rho_f and rho_s
                     rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0

                     beta2 = 0.66d0

                     ! calculate top and bottom shear stress.
                     t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                     beta = 1.d0-m
                     gamma= rho*beta2*(vnorm**2)*(beta*gmod*coeff**2)/(tanh(h+1.d-2)**(1.0d0/3.0d0))

                     t1bot = t1bot + gamma
                     t1bot = t1bot + tau

                     t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))

                     ! calculate pressure ratio
                     prat = p/(rho*h)

                     ! calculate dh
                     dh = entrainment_rate*dt*(t1bot-t2top)/(rho2*beta2*vnorm)
                     dh = min(dh,b_remaining)

                     ! increment h based on dh
                     h = h + dh
                     hm = hm + dh*m2

                     ! store amount eroded in q7
                     q(i_bdif,i,j) = q(i_bdif,i,j) + dh

                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                     call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

                     ! update pressure based on prior pressure ratio.
                     p = prat*rho*h

                     ! DIG: should hchi (pm) be updated here, just as there is a m2, should there be a pm2?

                     call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  endif
               endif
            endif

            !===================================================================
            ! end of entrainment, put state variables back in q.

            q(i_h,i,j) = h
            q(i_hu,i,j) = hu
            q(i_hv,i,j) = hv
            q(i_hm,i,j) = hm
            q(i_pb,i,j) = p
            q(i_hchi,i,j) = chi*h

         enddo
      enddo



      ! Manning friction------------------------------------------------
      if (friction_forcing) then
      if (coeff>0.d0.and.friction_depth>0.d0) then

         do j=1,my
            do i=1,mx

               if (bed_normal==1) gmod = grav*cos(aux(i_theta,i,j))
                  h=q(i_h, i,j)
               if (h<=friction_depth) then
                 !# apply friction source term only in shallower water
                  hu=q(i_hu,i,j)
                  hv=q(i_hv,i,j)
                  hm = q(i_hm,i,j)
                  p =  q(i_pb,i,j)
                  phi = aux(i_phi,i,j)
                  theta = aux(i_theta,i,j)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<dry_tolerance) cycle
                  pm = q(i_hchi,i,j)/h
                  pm = max(0.0d0,pm)
                  pm = min(1.0d0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                  if (h.lt.dry_tolerance) then
                     q(i_h,i,j)=0.d0
                     q(i_hu,i,j)=0.d0
                     q(i_hv,i,j)=0.d0
                  else
                     !beta = 1.d0-m  ! reduced friction led to high velocities
                     beta = 1.d0     ! use full Manning friction
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*coeff**2)/(h**(7.0d0/3.0d0))
                     dgamma=1.d0 + dt*gamma
                     q(i_hu,i,j)= q(i_hu,i,j)/dgamma
                     q(i_hv,i,j)= q(i_hv,i,j)/dgamma
                  endif
               endif
            enddo
         enddo
         endif
      endif
     ! ----------------------------------------------------------------

      return
      end

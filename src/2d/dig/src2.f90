

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module, only: grav, dry_tolerance,deg2rad,friction_depth,manning_coefficient,friction_forcing
      use digclaw_module

      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
      double precision, intent(in) :: xlower,ylower,dx,dy,t,dt
      
      ! Output
      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gmod,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi,dti
      real(kind=8) :: D,tau,sigbed,kperm,compress,pm,coeff,tol
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,taucf,fsphi,hvnorm0
      real(kind=8) :: shear,sigebar,pmtanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,chi,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining,dtcoeff
      real(kind=8) :: gamma,zeta,krate,p_eq,p_litho,p_hydro,dry_tol,dgamma

      integer :: i,j,ii,jj,jjend,icount,curvature
      logical :: ent


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

      gmod=grav
      coeff = manning_coefficient(1) ! Current implementation of friction has manning as an array 
      ! take the first element for now. If only one value is provided to geo_data.manning_coefficient 
      ! it will be manning_coefficient(1)
      ! DIG: FIX.

      tol = dry_tolerance !# to prevent divide by zero in gamma
      curvature = 0 !add friction due to curvature acceleration KRB: why is this hardcoded to 0?
      !write(*,*) 'src:init,value',p_initialized,init_pmin_ratio
      if (entrainment>0) then
         ent = .true.
      else
         ent = .false.
      endif

      do i=1-mbc+1,mx+mbc-1
         do j=1-mbc+1,my+mbc-1
            theta = 0.d0
            dtheta = 0.d0
            if (bed_normal==1) then
               theta = aux(i_theta,i,j)
               gmod = grav*cos(theta)
               dtheta = -(aux(i_theta,i+1,j) - theta)/dx
            endif

            !call admissibleq(q(1,i,j),q(2,i,j),q(3,i,j),q(4,i,j),q(5,i,j),u,v,m,theta)
            h = q(1,i,j)
            if (h<=dry_tolerance) cycle
            hu = q(2,i,j)
            hv = q(3,i,j)
            hm = q(4,i,j)
            p =  q(5,i,j)
            phi = aux(i_phi,i,j)
            pm = q(6,i,j)/h
            pm = max(0.0,pm)
            pm = min(1.0,pm)
            fsphi = aux(i_fsphi,i,j)

            jjend = 1
            dti = dt!/real(jjend,kind=8)
            !do jj=1,jjend

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !integrate momentum source term
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            !tau = max(tau*(1.0-fsphi),0.0)

            vnorm = sqrt(u**2 + v**2)
            hvnorm = sqrt(hu**2 + hv**2)
            hvnorm0 = hvnorm

            !integrate friction
            hvnorm = dmax1(0.d0,hvnorm - dti*tau/rho)
            hvnorm = hvnorm*exp(-(1.d0-m)*2.0d0*mu*dti/(rho*h**2))
            !hvnorm = hvnorm*exp(-(1.d0-m)*2.0d0*0.1*dti/(rho*h**2.0))
            if (hvnorm<1.d-16) hvnorm = 0.d0


            if (hvnorm>0.d0.and.curvature==1) then
               b_xx=(aux(1,i+1,j)-2.d0*aux(1,i,j)+aux(1,i-1,j))/(dx**2)
               b_yy=(aux(1,i,j+1)-2.d0*aux(1,i,j)+aux(1,i,j-1))/(dy**2)
               b_xy=(aux(1,i+1,j+1)-aux(1,i-1,j+1) -aux(1,i+1,j-1)+aux(1,i-1,j-1))/(4.0*dx*dy)
               chi = (u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy)/gmod
               chi = max(chi,-1.d0)
               taucf = chi*tau
               hvnorm = dmax1(0.d0,hvnorm - dti*taucf/rho)
               taucf = u**2*dtheta*tau/gmod
               hvnorm = dmax1(0.d0,hvnorm - dti*taucf/rho)
            endif

            if (hvnorm0>0.d0) then
               hu = hvnorm*hu/hvnorm0
               hv = hvnorm*hv/hvnorm0
            endif

            if (p_initialized==0) cycle

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            
            vnorm = sqrt(u**2 + v**2)

            !integrate shear-induced dilatancy
            sigebar = rho*gmod*h - p + sigma_0
            shear = 2.d0*vnorm/h
            krate = 1.5d0*shear*m*tanpsi/alpha
            sigebar = sigebar*exp(krate*dti)
            !p = rho*gmod*h + sigma_0 - sigebar
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dti*3.d0*vnorm*tanpsi/(h*compress)
            endif

            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               pmtanh01=0.d0

      		else
         		seg = 1.d0
               call calc_pmtanh(pm,seg,pmtanh01)
               rho_fp = max(0.d0,(1.d0-pmtanh01))*rho_f
      		endif
            !pmtanh01 = seg*(0.5*(tanh(20.0*(pm-0.80))+1.0))
            !pmtanh01 = seg*(0.5*(tanh(40.0*(pm-0.90))+1.0))
            
            !integrate pressure relaxation
            !if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
            !   zeta = 3.d0/(compress*h*2.0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !else
            !   zeta = (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            !endif
            zeta = ((m*(sigbed +  sigma_0))/alpha)*3.d0/(h*2.d0)  + (rho-rho_fp)*rho_fp*gmod/(4.d0*rho)
            krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
            p_hydro = h*rho_fp*gmod
            p_litho = (rho_s*m + (1.d0-m)*rho_fp)*gmod*h

            !if (abs(compress*krate)>0.0) then
            !   p_eq = p_hydro + 3.0*vnorm*tanpsi/(compress*h*krate)
            !else
            !   p_eq = p_hydro
            !endif
            !if (abs(pm-.5)>.49) then
            !pmtanh01 = 0.5*(tanh(20.0*(pm-0.80))+1.0)
            p_eq = p_hydro !*(1.0-pmtanh01)
            !p_eq = max(p_eq,0.0)
            !p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dti)


            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
            

            krate = D*(rho-rho_fp)/rho
            hu = hu*exp(dti*krate/h)
            hv = hv*exp(dti*krate/h)
            hm = hm*exp(-dti*D*rho_fp/(h*rho))
            h = h + krate*dti

            !enddo

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

            !======================mass entrainment===========================
            
            vnorm = sqrt(u**2 + v**2)
            vlow = 0.1d0

            if (ent.and.vnorm.gt.vlow.and.(aux(i_theta,i,j)>0.d0)) then
               b_x = (aux(1,i+1,j)+q(7,i+1,j)-aux(1,i-1,j)-q(7,i-1,j))/(2.d0*dx)
               b_y = (aux(1,i,j+1)+q(7,i,j+1)-aux(1,i,j-1)-q(7,i,j-1))/(2.d0*dy)
               dbdv = (u*b_x+v*b_y)/vnorm
               slopebound = 1.d10
               b_eroded = q(7,i,j)
               if (dbdv<slopebound.and.b_eroded<aux(i_theta,i,j)) then
                  b_remaining = aux(i_theta,i,j)-b_eroded
                  m2 = 0.6d0
                  rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0
                  beta2 = 0.66d0
                  t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                  !write(*,*) '------------'
                  !write(*,*) 'vu',t1bot
                  beta = 1.d0-m!tanh(10.d0*m) !tan(1.5*p/(rho*gmod*h))/14.0
                  gamma= rho*beta2*(vnorm**2)*(beta*gmod*coeff**2)/(tanh(h+1.d-2)**(1.0/3.0))
                  !write(*,*) 'gamma', gamma
                  t1bot = t1bot + gamma
                  t1bot = t1bot + tau!+p*tan(phi)
                  !write(*,*) 'tau',tau
                  t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))
                  !write(*,*) 't2top',t2top
                  prat = p/(rho*h)
                  !dh = dti*(t1bot-t2top)/(beta2*tanh(vnorm+1.d-2)*rho2)
                  vreg = ((vnorm-vlow)**2/((vnorm-vlow)**2+1.d0))
                  dtcoeff = entrainment_rate*dti*vreg/(beta2*(vnorm+vlow)*rho2)
                  !dh = dtcoeff*t1bot/(1.d0 + dtcoeff*tan(phi))
                  dh = dtcoeff*(t1bot-t2top)
                  dh = entrainment_rate*dti*(t1bot-t2top)/(rho2*beta2*vnorm)
                  !write(*,*) 'dh',dh
                  !write(*,*) 'dh/dti', dh/dti
                  dh = min(dh,b_remaining)
                  h = h + dh
                  hm = hm + dh*m2
                  q(7,i,j) = q(7,i,j) + dh

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)
                  p = prat*rho*h
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
               endif
            endif
            !===================================================================

            q(1,i,j) = h
            q(2,i,j) = hu
            q(3,i,j) = hv
            q(4,i,j) = hm
            q(5,i,j) = p
            q(6,i,j) = pm*h

         enddo
      enddo



      ! Manning friction------------------------------------------------
      if (friction_forcing) return
      if (coeff>0.d0.and.friction_depth>0.d0) then
         do i=1,mx
            do j=1,my

               if (bed_normal==1) gmod = grav*cos(aux(i_theta,i,j))
               h=q(i,j,1)
               if (h<=friction_depth) then
                 !# apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)
                  hm = q(4,i,j)
                  p =  q(5,i,j)
                  phi = aux(i_phi,i,j)
                  theta = aux(i_theta,i,j)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<dry_tol) cycle
                  pm = q(6,i,j)/h
                  pm = max(0.0,pm)
                  pm = min(1.0,pm)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,pm)

                  if (h.lt.tol) then
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     beta = 1.d0-m !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*coeff**2)/(h**(7.0/3.0))
                     dgamma=1.d0 + dt*gamma
                     q(2,i,j)= q(2,i,j)/dgamma
                     q(3,i,j)= q(3,i,j)/dgamma
                  endif
               endif
            enddo
         enddo
      endif
     ! ----------------------------------------------------------------

      return
      end

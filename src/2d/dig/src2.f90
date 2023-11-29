

   !=========================================================
      subroutine src2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt)
   !=========================================================
      use geoclaw_module, only: grav, dry_tolerance,deg2rad,friction_depth,manning_coefficient,friction_forcing
      use digclaw_module

      implicit none

      ! Input parameters
      integer, intent(in) :: meqn,mbc,mx,my,maux
      real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
      
      ! Output
      real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      !local
      real(kind=8) :: gz,gx,h,hu,hv,hm,u,v,m,p,phi,kappa,S,rho,tanpsi
      real(kind=8) :: rhoh,manning_n
      real(kind=8) :: D,tau,sigbed,kperm,compress,chi,tol
      real(kind=8) :: vnorm,hvnorm,theta,dtheta,w,hvnorm0
      real(kind=8) :: shear,sigebar,chitanh01,rho_fp,seg
      real(kind=8) :: b_xx,b_yy,b_xy,gacc,beta
      real(kind=8) :: t1bot,t2top,beta2,dh,rho2,prat,b_x,b_y,dbdv
      real(kind=8) :: vlow,m2,vreg,slopebound
      real(kind=8) :: b_eroded,b_remaining,dtcoeff
      real(kind=8) :: gamma,zeta,krate,p_eq,p_litho,p_hydro,dgamma

      integer :: i,j,ii,jj,jjend,icount


      ! check for NANs in solution:
      call check4nans(meqn,mbc,mx,my,q,t,2)

      manning_n = manning_coefficient(1) ! Current implementation of friction has manning as an array 
      ! take the first element for now. If only one value is provided to geo_data.manning_coefficient 
      ! it will be manning_coefficient(1)
      ! DIG: FIX.

      tol = dry_tolerance !# to prevent divide by zero in gamma
      
      gz = grav  !needed later for bed-normal direction gravity
      gx = 0.d0
      theta=0.d0 

      do i=1-mbc+1,mx+mbc-1
         do j=1-mbc+1,my+mbc-1
            
            if (q(1,i,j)<=dry_tolerance) cycle
            
            h = q(1,i,j)
            hu = q(2,i,j)
            hv = q(3,i,j)
            hm = q(4,i,j)
            p =  q(5,i,j)
            phi = aux(ia_phi,i,j)
            chi = q(iq_seg,i,j)/h
            chi = max(0.0,chi)
            chi = min(1.0,chi)
            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)

            !modified gravity: bed-normal weight and acceleration
            if (bed_normal==1) then
               theta = aux(ia_theta,i,j)
               gz = grav*cos(theta)
               gx = grav*sin(theta)
            endif
            if (curvature==1) then
               b_xx=(aux(1,i+1,j)-2.d0*aux(1,i,j)+aux(1,i-1,j))/(dx**2)
               b_yy=(aux(1,i,j+1)-2.d0*aux(1,i,j)+aux(1,i,j-1))/(dy**2)
               b_xy=(aux(1,i+1,j+1)-aux(1,i-1,j+1) -aux(1,i+1,j-1)+aux(1,i-1,j-1))/(4.0*dx*dy)
               dtheta = -(aux(ia_theta,i+1,j) - theta)/dx
               gacc = max((u**2*b_xx + v**2*b_yy + 2.0*u*v*b_xy + u**2*dtheta,0.d0)!max:currently only consider enhancement not reduction of gz (ie. basin not a hump)
               gz = gz + gacc
            endif

            call vareval(h,u,v,m,p,rho,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation)

            rhoh = h*rho !this is invariant in src and always >0 below
            hvnorm0 = sqrt(hu**2 + hv**2)
            vnorm = hvnorm0/h

            if (hvnorm0>0.d0) then
               !integrate dynamic friction !DIG: TO DO - move dynamic friction to Riemann solver
               vnorm = dmax1(0.d0,vnorm - dt*tau/rhoh) !exact solution for Coulomb friction
               vnorm = vnorm*exp(-(1.d0-m)*2.0d0*mu*dt/(h*rhoh)) !exact solution (prior to h change) for effective viscous friction
               ! velocity determined, calculate directions etc. from vnorm
               hvnorm = h*vnorm
               hu = hvnorm*hu/hvnorm0 + gx*h*dt !gx=0 unless bed-normal !DIG: last term should ultimately be in Riemann solver
               hv = hvnorm*hv/hvnorm0
               u = hu/h
               v = hv/h
               ! velocity now constant for remainder of src2
            endif

            if (p_initialized==0) cycle !DIG: deprecate?

            !determine m and p evolution from multiple integration steps maintaining physically admissable m,p. h = rhoh/rho
            ! Note: if m = 0, it cannot increase as dm/dt - m D. In that case, rho=rho_f ==> D = 0 for all t, and dp/dt =0.
            if (m>0.d0) then
               dt_remaining = dt
               if (D==0.d0 & (vnorm==0.d0)) then !at a critical point already, rhs = 0
                     dt_remaining = 0.d0
               endif
               do while (dt_remaining>0.d0)
                  call integrate_mp(h,u,v,m,p,rho,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation,dt_remaining,dt_taken)
                  dt_remaining = max(0.d0,dt_remaining-dt_taken)
               enddo
            endif





            !integrate shear-induced dilatancy
            sigebar = rho*gmod*h - p + sigma_0
            shear = 2.d0*vnorm/h
            krate = 1.5d0*shear*m*tanpsi/alpha
            !sigebar = sigebar*exp(krate*dt)
            !p = rho*gmod*h + sigma_0 - sigebar
            if (compress<1.d15) then !elasticity is = 0.0 but compress is given 1d16 in auxeval
               p = p - dt*3.d0*vnorm*tanpsi/(h*compress)
            endif

            !call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            !call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
            if (dabs(alpha_seg-1.d0)<1.d-6) then
         		seg = 0.d0
               rho_fp = rho_f
               chitanh01=0.d0

      		else
         		seg = 1.d0
               call calc_chitanh(chi,seg,chitanh01)
               rho_fp = max(0.d0,(1.d0-chitanh01))*rho_f
      		endif
            !chitanh01 = seg*(0.5*(tanh(20.0*(chi-0.80))+1.0))
            !chitanh01 = seg*(0.5*(tanh(40.0*(chi-0.90))+1.0))
            
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
            !if (abs(chi-.5)>.49) then
            !chitanh01 = 0.5*(tanh(20.0*(chi-0.80))+1.0)
            p_eq = p_hydro !*(1.0-chitanh01)
            !p_eq = max(p_eq,0.0)
            !p_eq = min(p_eq,p_litho)

            p = p_eq + (p-p_eq)*exp(krate*dt)


            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
            

            krate = D*(rho-rho_fp)/rho
            hu = hu*exp(dt*krate/h)
            hv = hv*exp(dt*krate/h)
            hm = hm*exp(-dt*D*rho_fp/(h*rho))
            h = h + krate*dt

            !enddo

            call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
            call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

            !======================mass entrainment===========================
            
            vnorm = sqrt(u**2 + v**2)
            vlow = 0.1d0

            if (ent.and.vnorm.gt.vlow.and.(aux(ia_theta,i,j)>0.d0)) then
               b_x = (aux(1,i+1,j)+q(7,i+1,j)-aux(1,i-1,j)-q(7,i-1,j))/(2.d0*dx)
               b_y = (aux(1,i,j+1)+q(7,i,j+1)-aux(1,i,j-1)-q(7,i,j-1))/(2.d0*dy)
               dbdv = (u*b_x+v*b_y)/vnorm
               slopebound = 1.d10
               b_eroded = q(7,i,j)
               if (dbdv<slopebound.and.b_eroded<aux(ia_theta,i,j)) then
                  b_remaining = aux(ia_theta,i,j)-b_eroded
                  m2 = 0.6d0
                  rho2 = m2*2700.d0 + (1.d0-m2)*1000.d0
                  beta2 = 0.66d0
                  t1bot = beta2*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                  !write(*,*) '------------'
                  !write(*,*) 'vu',t1bot
                  beta = 1.d0-m!tanh(10.d0*m) !tan(1.5*p/(rho*gmod*h))/14.0
                  gamma= rho*beta2*(vnorm**2)*(beta*gmod*manning_n**2)/(tanh(h+1.d-2)**(1.0/3.0))
                  !write(*,*) 'gamma', gamma
                  t1bot = t1bot + gamma
                  t1bot = t1bot + tau!+p*tan(phi)
                  !write(*,*) 'tau',tau
                  t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))
                  !write(*,*) 't2top',t2top
                  prat = p/(rho*h)
                  !dh = dt*(t1bot-t2top)/(beta2*tanh(vnorm+1.d-2)*rho2)
                  vreg = ((vnorm-vlow)**2/((vnorm-vlow)**2+1.d0))
                  dtcoeff = entrainment_rate*dt*vreg/(beta2*(vnorm+vlow)*rho2)
                  !dh = dtcoeff*t1bot/(1.d0 + dtcoeff*tan(phi))
                  dh = dtcoeff*(t1bot-t2top)
                  dh = entrainment_rate*dt*(t1bot-t2top)/(rho2*beta2*vnorm)
                  !write(*,*) 'dh',dh
                  !write(*,*) 'dh/dt', dh/dt
                  dh = min(dh,b_remaining)
                  h = h + dh
                  hm = hm + dh*m2
                  q(7,i,j) = q(7,i,j) + dh

                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)
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
            q(6,i,j) = chi*h

         enddo
      enddo



      ! Manning friction------------------------------------------------
      if (friction_forcing) then
      if (manning_n>0.d0.and.friction_depth>0.d0) then
         do i=1,mx
            do j=1,my
                  h=q(1, i,j)
               if (h<=friction_depth) then
                 !# apply friction source term only in shallower water
                  hu=q(2,i,j)
                  hv=q(3,i,j)
                  hm = q(4,i,j)
                  p =  q(5,i,j)
                  phi = aux(ia_phi,i,j)
                  theta = aux(ia_theta,i,j)
                  call admissibleq(h,hu,hv,hm,p,u,v,m,theta)
                  if (h<dry_tolerance) cycle
                  chi = q(6,i,j)/h
                  chi = max(0.0,chi)
                  chi = min(1.0,chi)
                  call auxeval(h,u,v,m,p,phi,theta,kappa,S,rho,tanpsi,D,tau,sigbed,kperm,compress,chi)

                  if (h.lt.dry_tolerance) then
                     q(2,i,j)=0.d0
                     q(3,i,j)=0.d0
                  else
                     beta = 1.d0-m !tan(1.5*p/(rho*gmod*h))/14.0
                     gamma= beta*dsqrt(hu**2 + hv**2)*(gmod*manning_n**2)/(h**(7.0/3.0))
                     dgamma=1.d0 + dt*gamma
                     q(2,i,j)= q(2,i,j)/dgamma
                     q(3,i,j)= q(3,i,j)/dgamma
                  endif
               endif
            enddo
         enddo
         endif
      endif
     ! ----------------------------------------------------------------

      return
      end subroutine src2

   !====================================================================
   ! subroutine integrate_mp: integrate portion of rhs for m and p
   !     dm/dt ~ -Dm
   !     dp/dt ~ D + (m-m_eqn)
   !
   !     integrated for p_tilde = p - p_eq
   !     dm/dt ~ -p_tilde m
   !     dp_tilde/dt ~ p_tilde + (m-m_eqn)
   !     note: p_eq is h(t) dependent. dp_tilde/dt = dp/dt - rho_f gz *dh/dt ~ dp/dt + p_tilde
   !====================================================================

   subroutine integrate_mp(h,u,v,m,p,rhoh,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation,dt_remaining,dt)

      use geoclaw_module, only: grav, dry_tolerance
      use digclaw_module, only: rho_f, rho_s, m_crit
      implicit none

      use geoclaw_module, only: grav, dry_tolerance

      !i/o
      real(kind=8), intent(in) :: vnorm,rhoh,phi_bed,gz
      real(kind=8), intent(inout)  :: h,m,p,kperm,dt_remaining,dt_taken
      real(kind=8), intent(out) ::

      !local
      real(kind=8) :: p_eq, p_tilde
      real(kind=8) :: 

      call vareval(h,u,v,m,p,rhoh,sigma,sigma_e,tanphi,tanpsi,D,tau,kperm,compress,chi,gz,M_saturation)
      shear = vnorm/h !

      p_eq = rho_f*gz*h*M_saturation !Note: M_saturation =1, unless experimenting with segregation models
      p_tilde = p - p_eq
      m_eq = !DIG: WIP - maybe this and above calculated in vareval

      if (D==0.d0 & (sigma_e==0.d0|shear==0.d0|tanpsi==0.d0)) then
         dt = dt_remaining
         return
      endif

      if (D==0.d0) then !integrate only shear induced dilatancy. no change in h or m.
         rhs_p = -(3.d0/alpha_c)*m*sigma_e*shear*tanpsi
         if (rhs_p>0.d0) then !m<m_eq => increasing p
            dt = min(abs((rhoh*gz-p)/rhs_p),dt_remaining) !enforce dt st p1<= rho g h
            !forward euler. rhs_p is nonlinear function of p through tanpsi(m_eq)
            !DIG: - update to RK. 
            p = p + rhs_p*dt 
            return
         elseif (rhs_p<0.d0) then
            dt = min(abs(p/rhs_p),dt_remaining) !enforce dt st p1>=0
            !forward euler. rhs_p is nonlinear function of p through tanpsi(m_eq)
            !DIG: - update to RK. 
            p = p + rhs_p*dt 
            return
         else
            dt = dt_remaining
            return
         endif 
      elseif (D>0.d0) then
         !note: by integrating m (=> rho and hence h), dh/dt terms on rhs of p are avoided by using new p_eq(h).

         rhs_ptilde =  -(3.d0*alphainv/h)*(vnorm*tanpsi +0.5d0*

         
      else

      endif

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

   end subroutine integrate_mp
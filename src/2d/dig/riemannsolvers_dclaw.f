
c-----------------------------------------------------------------------
      subroutine riemann_dig2_aug_sswave_ez(ixy,meqn,mwaves,hL,hR,
     &         huL,huR,hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         thetaL,thetaR,phiL,phiR,sw,fw,wallprob,taudirL,
     &         taudirR,chiL,chiR,fsL,fsR,ilook)

      !-----------------------------------------------------------------
      ! solve the D-Claw Riemann problem at each interface for debris flow eqn
      ! this is for 2d version
      !
      ! This solver is an extension of GeoClaw's solver described in (George, 2008).
      ! See (George & Iverson, 2014) for an overview of this D-Claw extension.
      ! 
      ! The Riemann solver solves the left-hand side of the PDEs, neglecting source terms
      ! except for the topographic source terms, -ghb_x and -ghb_y.
      ! In the case of static material that does not fail (deform)
      ! topographic source terms and static friction balance, accounted for here.
      ! For moving material, frictional resistance is handled with the other 
      ! source terms (the right hand side) accounted for in a fractional step
      ! implemented in src2 subroutine
      !-----------------------------------------------------------------

      use geoclaw_module, only: grav, dry_tolerance
      use digclaw_module, only: beta_seg, rho_f, kappa

      implicit none

*     !i/o
      integer ixy,meqn,mwaves,ilook

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uL,uR,vL,vR,mL,mR,chiL,chiR,seg_L,seg_R
      double precision thetaL,thetaR,phiL,phiR
      double precision taudirL,taudirR,fsL,fsR
      logical wallprob


      double precision fw(meqn,mwaves)
      double precision sw(mwaves)

*     !local
      integer m,mw,k,cwavetype
      !double precision h,u,v
      double precision det1,det2,det3,determinant
      double precision R(0:2,1:3),del(0:4) 
      double precision beta(3)
      double precision pratL,pratR,tan_phi_max,source2dxf
      double precision phiL_effective, phiR_effective,phi_eff
      double precision rhoL,rhoR,tauL,tauR,rho_bar,rhoedge
      double precision tanpsi, delbf, deldelhf,huedge
      double precision kperm,m_eq,alphainv,s1s2_denom,tauratL,tauratR
      double precision theta,gamma,eps,taudirUfrac,hustarHLLn
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision delb,s1m,s2m,hm,criticaltol,criticaltol_2
      double precision s1s2bar,s1s2tilde,hbar,source2dx,veltol1,veltol2
      double precision s1s2_ratio,ss_delta,dels,hustarHLL,hmin
      double precision hstarHLL,deldelh,drytol,gz,geps,tausource
      double precision raremin,raremax,rare1st,rare2st,sdelta
      double precision gammaL,gammaR,theta1,theta2,theta3,vnorm
      double precision alpha_seg,a,b,c,ctn,hsmallest,uedge,vedge
      logical sonic,rare1,rare2
      logical rarecorrectortest,rarecorrector

      veltol1=1.d-6
      veltol2=0.d0
      !criticaltol=1.d-6
      drytol = dry_tolerance
      criticaltol = 1.d-3! max(drytol*grav, 1d-6)
      criticaltol_2 = sqrt(criticaltol)

      theta = 0.5d0*(thetaL + thetaR)
      gz = grav*dcos(theta)


      ! alpha seg here reflects alpha in Gray and Kokelaar (2010)
      ! alpha_seg 0 = simple shear and 1 = plug flow
      ! we have defined beta_seg as 1-alpha_seg
      alpha_seg=1.0d0-beta_seg

         call setvars(hL,uL,vL,mL,pL,chiL,gz,rhoL,kperm,alphainv,m_eq,
     &        tanpsi,tauL)

         call setvars(hR,uR,vR,mR,pR,chiR,gz,rhoR,kperm,alphainv,m_eq,
     &        tanpsi,tauR)


      hbar = 0.5d0*(hL + hR)
      rho_bar = 0.5d0*(rhoL + rhoR)
      
      !tau = 0.5d0*(tauL + tauR)
      select case (src2method)
      case(-1)
         ! gamma is
         gamma = rho_f/rho_bar
         gammaL = rho_f/rhoL
         gammaR = rho_f/rhoR

      case(0:1)
         gamma = 0.25d0*(rho_f + 3.0d0*rho)/rho_bar
         gammaL = 0.25d0*(rho_f + 3.0d0*rhoL)/rhoL
         gammaR = 0.25d0*(rho_f + 3.0d0*rhoR)/rhoR
      end select

      eps = kappa + (1.d0-kappa)*gamma
      geps = gz*eps

      !determine wave speeds
      sL=uL-dsqrt(geps*hL) ! 1 wave speed of left state
      sR=uR+dsqrt(geps*hR) ! 2 wave speed of right state
      uhat=(dsqrt(hL)*uL + dsqrt(hR)*uR)/(dsqrt(hR)+dsqrt(hL)) ! Roe average
      chat=dsqrt(geps*0.5d0*(hR+hL)) ! Roe average
      sRoe1=uhat-chat ! Roe wave speed 1 wave
      sRoe2=uhat+chat ! Roe wave speed 2 wave

      sE1 = min(sL,sRoe1) ! Eindfeldt speed 1 wave
      sE2 = max(sR,sRoe2) ! Eindfeldt speed 2 wave
      sw(1) = sE1
      sw(3) = sE2
      !u = uhat

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &                                          1,drytol,geps)
      sw(1)= min(sw(1),s2m) !Modified Einfeldt speed
      sw(3)= max(sw(3),s1m) !Modified Einfeldt speed
      sw(2) = 0.5d0*(sw(3)+sw(1))
      dels = sw(3)-sw(1)
      !hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! single middle state in an HLL solve
      hstarHLL = max((huL-huR+sw(3)*hR-sw(1)*hL)/(sw(3)-sw(1)),drytol) ! single middle state in an HLL solve
      hustarHLL = (dels*huL + sw(1)*sw(3)*(hR-hL)+sw(1)*(huL-huR))/dels
c     !determine the middle entropy corrector wave------------------------
      rarecorrectortest = .false.
      rarecorrector=.false.
      if (rarecorrectortest) then
         sdelta=sw(3)-sw(1)
         raremin = 0.5d0
         raremax = 0.9d0
         if (rare1.and.sE1*s1m.lt.0.d0) raremin=0.2d0
         if (rare2.and.sE2*s2m.lt.0.d0) raremin=0.2d0
         if (rare1.or.rare2) then
            !see which rarefaction is larger
            rare1st=3.d0*(dsqrt(geps*hL)-dsqrt(geps*hm))
            rare2st=3.d0*(dsqrt(geps*hR)-dsqrt(geps*hm))
            if (max(rare1st,rare2st).gt.raremin*sdelta.and.
     &         max(rare1st,rare2st).lt.raremax*sdelta) then
                  rarecorrector=.true.
               if (rare1st.gt.rare2st) then
                  sw(2)=s1m
               elseif (rare2st.gt.rare1st) then
                  sw(2)=s2m
               else
                  sw(2)=0.5d0*(s1m+s2m)
               endif
            endif
         endif
         if (hstarHLL.lt.min(hL,hR)/5.d0) rarecorrector=.false.
      endif

      !determine cell edge conditions for friction 

      delb=(bR-bL)!#kappa

      !determine ss-wave
      s1s2bar = 0.25d0*(uL+uR)**2- gz*hbar
      s1s2tilde= max(0.d0,uL*uR) - gz*hbar

c     !find if sonic problem (or very far from steady state)
      !sonic=.false.
      
      !if (s1s2bar*s1s2tilde.le.criticaltol**2) then
      !   sonic=.true.
      !elseif (dabs(s1s2bar).le.criticaltol) then
      !   sonic=.true.
      !elseif (uL*uR.lt.0.d0) then 
      !   sonic =.true.
      !elseif (s1s2bar*sw(1)*sw(3).le.criticaltol**2) then
      !   sonic = .true.
      !elseif (min(dabs(sE1),dabs(sE2)).lt.criticaltol_2) then
      !   sonic=.true.
      !elseif (sE1.lt.criticaltol_2.and.s1m.gt.-criticaltol_2) then
      !   sonic = .true.
      !elseif (sE2.gt.-criticaltol_2.and.s2m.lt.criticaltol_2) then
      !   sonic = .true.
      !elseif ((uL+dsqrt(geps*hL))*(uR+dsqrt(geps*hR)).lt.0.d0) then
      !   sonic=.true.
      !elseif ((uL-dsqrt(geps*hL))*(uR-dsqrt(geps*hR)).lt.0.d0) then
      !   sonic=.true.
      !endif

      ! find effective friction at interface
      ! note: this is a workaround for pressure/effective stress
      ! a more rigorous approach is under development but would involve
      ! finding exact steady state solutions when d(p/rho*g*h)/dx neq 0.d0
      if (hR>drytol) then
        pratR = max(min(pR/(rhoR*gz*hR),1.d0),0.d0)
        tauratR = tauR/(rhoR*gz*hR)
      else
        pratR = max(min(pL/(rhoL*gz*hL),1.d0),0.d0)
        tauratR = tauL/(rhoL*gz*hL)
      endif

      if (hL>drytol) then
        pratL = max(min(pL/(rhoL*gz*hL),1.d0),0.d0)
        tauratL = tauL/(rhoL*gz*hL)
      else
        pratL = max(min(pR/(rhoR*gz*hR),1.d0),0.d0)
        tauratL = tauR/(rhoR*gz*hR)
      endif
      phiL_effective = atan(tauratR)
      phiR_effective = atan(tauratL)
        !write(*,*) 'pratL,pratR',pratL,pratR
       !phiL_effective = atan(max(0.d0,(1.d0-pratL))*tan(phiL))
       !phiR_effective = atan(max(0.d0,(1.d0 - pratR))*tan(phiR))
       !phiL_effective = atan(max(0.d0,(tauratL+pratL))*tan(phiL))
       !phiR_effective = atan(max(0.d0,(tauratR+pratR))*tan(phiR))
       phi_eff = max(0.5d0*(phiL_effective + phiR_effective),0.d0)
       if (phi_eff.lt.1.d-6) then
          phi_eff=0.d0
        endif
        
      ! determine if steady state Riemann invariants are close or far
      ! if far, critical excess ratio, s1s2_ratio, evaluated using
      ! far left and right states (qR,qL) are not a good approx at interface
      ! because cell edge values are not necessarily near qL and qR 
      ! ss_delta = 0.0 => steady state data. ss_delta = 1.0 => far from steady state
      ! replaces approach where near sonic states are tested. better generally I think
      ! invariants for steady states are hu and 0.5u**2 + g*eta + x*sgn(u)tanphi
      if ((dabs(huL)+dabs(huR)).gt.0.d0) then
        ss_delta = dabs(huR - huL)/(dabs(huL)+dabs(huR)) !1 if opposite sign (non-steady)
      else
        ss_delta = 0.d0
      endif
      ! 2nd Riemann invariant
      ss_delta = max(ss_delta,
     &  dabs(0.5*uR**2 + gz*hR +gz*bR - 0.5*uL**2 - gz*hL -gz*bL 
     &   + gz*taudirR*tan(phi_eff)*0.5d0*
     &     (dsign(1.d0,uR)+dsign(1.d0,uL)))/
     &   (dabs(0.5*uR**2 +gz*hR)+dabs(0.5*uL**2 +gz*hL)+dabs(bR-bL)
     &    + dabs(gz*taudirR*tan(phi_eff))))
      ! transcritical (metric either 1 or zero)
       ss_delta = max(ss_delta,
     &  dabs(dsign(1.d0,(uR**2-gz*hR))-dsign(1.d0,(uL**2-gz*hL))))/2.d0
        !fix rounding error if any ss_delta in [0,1]
       ss_delta = max(0.d0, min(ss_delta,1.d0))
       !ss_delta = 1.d0
      ! bound jump in h at interface, positivity constraint, also constrains source term
      sonic = .false.
      ctn = criticaltol/dsqrt(gz*hbar) !normalize the tolerance to depth (very small depths not necessarily near critical)
      if (sw(1).gt.ctn) then 
        if (hL.gt.drytol) then
            s1s2bar = max(delb*gz*hbar*sw(1)/(hstarHLL*dels),
     &         max(s1s2bar,gz*hbar*delb/-hL))
        else
            s1s2bar = max(s1s2bar,0.d0)
        endif
      elseif (sw(3).lt.-ctn) then
        if (hR.gt.drytol) then
            s1s2bar = max(gz*hbar*delb*sw(3)/(hstarHLL*dels),
     &       max(s1s2bar,gz*hbar*delb/hR))    
        else      
            s1s2bar = max(s1s2bar,0.d0)
        endif 
      elseif (sw(1).lt.-ctn*0.d0.and.sw(3).gt.
     &                          0.d0*ctn) then
         if (hstarHLL.gt.drytol) then
            s1s2bar=min(sw(1)*gz*hbar*delb/(hstarHLL*dels),
     &           min(s1s2bar,sw(3)*gz*hbar*delb/(hstarHLL*dels)))
         else 
            !s1s2bar = min(s1s2bar,0.d0)
            s1s2bar = min(s1s2bar,min(uL**2-gz*hL,uR**2-gz*hR)) !for either shock or rare we want subcrit. estimate
            ss_delta = 1.d0
         endif
      else
        s1s2bar = min(uL**2-gz*hL,uR**2-gz*hR) !for either shock or rare we want subcrit. estimate
        sonic = .true.
        ss_delta = 1.d0
      endif

      !if (s1s2bar*sw(1)*sw(3).le.criticaltol**4) then
      !   s1s2bar = dsign(max(abs(s1s2bar),abs(sw(1)*sw(3))),sw(1)*sw(3))
      !   sonic=.true.
      !   ss_delta = 1.d0
      !endif

      !if (s1s2bar*s1s2tilde.le.criticaltol**4) then
       !  s1s2bar = dsign(max(abs(s1s2bar),abs(sw(1)*sw(3))),sw(1)*sw(3))
       !  sonic=.true.
       !  ss_delta = 1.d0
      !endif

      ! if actual ubar**2-ghbar at interface = 0
      ! then cant be a true solution if delb !=0
      ! true solution in this case could be stationary shock
      ! at interface, but only if delb!=0.
      ! but s1s2bar is only approximation using qR,qL
      ! if close to zero (near critical flow) use standard averages
      ! following set of s1s2bar just prevents divide by 0
      ! s1s2_ratio won't use s1s2bar because ss_delta=1=> s1s2_ratio = 1
      if (dabs(s1s2bar).le.ctn**2) then
        s1s2bar = dsign(ctn**2,s1s2bar)
        !s1s2bar = min(uL**2-gz*hL,uR**2-gz*hR)
        !sonic=.true.
        ss_delta = 1.d0
      endif
      !if (sonic) then
      !   s1s2_ratio = 1.d0
      !   deldelh =  -delb
         !source2dx = -gz*hbar*delb
      !else
      s1s2_ratio = max(0.d0,s1s2tilde/s1s2bar)
      s1s2_ratio = min(max(hL/hbar,hR/hbar),
     &     max(s1s2_ratio,min(hL/hbar,hR/hbar))) !bounds |source term| in (g*min(h)delb,g*max(h)delb)
         !deldelh = delb*((1.d0-ss_delta)*(gz*hbar/s1s2bar)-ss_delta)

      s1s2_ratio = ((1.d0-ss_delta)*s1s2_ratio + ss_delta)
      s1s2_denom = ((1.d0-ss_delta)*(1.d0/s1s2bar) 
     &        + ss_delta*(1.d0/(-gz*hbar)))
      deldelh = delb*gz*hbar*s1s2_denom !jump in h at interface from bathy not friction
      source2dx = -gz*hbar*delb*s1s2_ratio !source from bathy no friction yet
      !source2dx=min(source2dx,gz*max(-hL*delb,-hR*delb)) 
      !source2dx=max(source2dx,gz*min(-hL*delb,-hR*delb))
      
      vnorm = sqrt(uR**2 + uL**2 + vR**2 + vL**2)
      !hu at interface not considering friction. friction should oppose this velocity
      !if this is a static problem (vnorm=0) then friction opposes net force and determined further below
      hustarHLLn = (dels*huL + sw(1)*sw(3)*(hR-hL-deldelh)
     &    +sw(1)*(huL-huR+source2dx))/dels
      !if (vnorm.gt.veltol1) then
        if (sw(1).gt.0.d0) then
          uedge = uL
          huedge = huL
        elseif (sw(3).lt.0.d0) then
          uedge = uR
          huedge = huR
        else
          uedge = hustarHLLn/hstarHLL
          huedge = hustarHLLn
        endif
        if (uedge >0.d0) then
          vedge = vL
          rhoedge = rhoL
        elseif (uedge<0.d0) then
          vedge = vR
          rhoedge = rhoR
        else
          vedge = 0.5d0*(vL+vR)
          rhoedge = 0.5d0*(rhoL + rhoR)
        endif
        if (abs(uedge).gt.1.d-6) then
          taudirUfrac = taudirR*uedge/(sqrt(uedge**2 + vedge**2))
        else
          taudirUfrac = 0.d0
        endif

        tan_phi_max = dabs(huedge)/
     &       max(abs(gz*hbar*s1s2_ratio*taudirUfrac),1.d-12)
        tan_phi_max = 0.99999d0*tan_phi_max
        !write(*,*) 'tan_phi_max,tan(phi_eff)',tan_phi_max,tan(phi_eff)
        

      if (vnorm.gt.veltol1) then
        delbf = taudirUfrac*min(tan(phi_eff),tan_phi_max)
        deldelhf = delbf*gz*hbar*s1s2_denom
        if (abs(deldelhf).eq.0.d0) then
            delbf = 0.d0
        endif
      else
        delbf = 0.d0
        deldelhf = 0.d0
      endif
        !write(*,*) 'delbf,deldelhf', delbf,deldelhf
        !write(*,*) 'ss_delta,taudir', ss_delta,taudirUfrac,taudirR
        deldelh = deldelh + deldelhf
        source2dxf = - gz*hbar*delbf*s1s2_ratio
        if (abs(source2dxf).gt.abs(huedge)) then
            write(*,*) 'source2dxf,huedge',source2dxf,huedge
            write(*,*) 's1s2_ratio',s1s2_ratio
            write(*,*) 'ratio:', abs(source2dxf)-abs(huedge)
        endif
        !source2dx = source2dx + source2dxf

      ! if (dabs(u).le.veltol2) then
      !   source2dx=-hbar*gz*delb
      !endif

c     !find bounds in case of critical state resonance, or negative states
c     !find jump in h, deldelh

      !deldelh = delb*gz*hbar/s1s2bar
      !if (sonic) then
      !   deldelh =  -delb
      !else
      !   deldelh = delb*((1.d0-ss_delta)*(gz*hbar/s1s2bar)-ss_delta)
      !endif
c     !find bounds on deltah at interface based on depth positivity constraint and energy
      !constraint ensures energy does not increase accross interface
      
      if (sE1.lt.-ctn.and.sE2.gt.ctn) then
        hsmallest = 0.d0*((hustarHLL)**2/(2.d0*gz))**(1.d0/3.d0) 
        hmin = max(hstarHLL -hsmallest,0.d0)!small depth for energy constraint
        deldelh = min(deldelh,hmin*(sE2-sE1)/sE2)
        deldelh = max(deldelh,hmin*(sE2-sE1)/sE1)
      elseif (sE1.ge.ctn) then
        deldelh = min(deldelh,hmin*(sE2-sE1)/sE1)
        deldelh = max(deldelh,-hL)
            
      elseif (sE2.le.-ctn) then
        deldelh = min(deldelh,hR)
        deldelh = max(deldelh,hmin*(sE2-sE1)/sE2)
      endif


*     !determine R
      R(0,2) = 0.d0
      R(1,2) = 0.d0
      R(2,2) = 1.d0
      cwavetype = 1

      R(0,1) = 1.d0
      R(1,1) = sw(1)
      R(2,1) = sw(1)**2

      R(0,3) = 1.d0
      R(1,3) = sw(3)
      R(2,3) = sw(3)**2

      if (rarecorrector) then
         R(0,2) = 1.d0
         R(1,2) = sw(2)
         R(2,2) = sw(2)**2
        cwavetype = 2
      endif

      !determine del
      del(0) = hR- hL 
      del(1) = huR - huL
      del(2) = hR*uR**2 + 0.5d0*kappa*gz*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*gz*hL**2)
      del(2) = del(2) + (1.d0-kappa)*hbar*(pR-pL)/rho_bar
      del(3) = pR - pL - gamma*rho_bar*gz*deldelh
      del(4) = -gamma*rho_bar*gz*uhat*(hR-hL) 
     &    + gamma*rho_bar*gz*del(1) + uhat*(pR-pL)

*     !determine the source term

      !if (ixy.eq.1) then
         ! DIG: note that theta = 0.0 unless bed_normal is true. For now, assume bed_normal is false. Resolve if dx is needed later.
         !source2dx = source2dx !+ dx*hbar*grav*dsin(theta)
         ! DIG: this is the only place dx is needed
         ! until fixed, bed_normal = 1 yields error in make .data (1/30/24)
      !endif

      if (vnorm.le.veltol1) then 
        ! static problem
        ! 2 possibilities: friction should not be larger than net force in x-direction. 
        ! taudirR has already been reduced to dx Fx/|F| in module
        ! friction can be less than Fx but it cannot be larger than Fx
        ! otherwise Riemann problem would fail in the wrong direction
        !if (abs(hbar*taudirR*gz*tan(phi_eff))
        if (abs(hbar*gz*taudirR*min(tan(phi_eff),tan_phi_max))
     &            .lt.abs(del(2)-source2dx)) then
            !failure and friction opposes failure
            delbf = dsign(1.d0,del(2)-source2dx)*
     &       taudirR*min(tan(phi_eff),tan_phi_max)
            source2dxf = - gz*hbar*delbf
            deldelhf = -delbf
            deldelh = deldelh + deldelhf
        else 
            !friction opposes net force, no waves
            del(0) = deldelh
            del(1) = 0.0d0
            source2dx = del(2)
            source2dxf = 0.d0
            del(4) = 0.0d0
        endif
      endif

      !if (wallprob) then
      !   tausource = 0.0d0
      !endif
      del(0) = del(0) - deldelh 
      del(2) = del(2) - source2dx -source2dxf

      !--------theta--------------------
      if (sw(1).ge.0.d0) then
         theta1 = thetaR
         theta2 = thetaR
         theta3 = thetaR
      elseif (sw(3).le.0.d0) then
         theta1 = thetaL
         theta2 = thetaL
         theta3 = thetaL
      elseif (sw(2).ge.0.d0) then
         theta1 = thetaL
         theta2 = thetaR
         theta3 = thetaR
      else
         theta1 = thetaL
         theta2 = thetaL
         theta3 = thetaR
      endif


       !R beta = del
        a = sw(1)
        b = sw(2)
        c = sw(3)
  
        !solve for beta = Rinv*delta
        if (cwavetype==1) then
          !r2 is (0,0,1)
          beta(1) = (c*del(0) - del(1))/(c-a)
          beta(2) = a*c*del(0) - (a+c)*del(1) + del(2)
          beta(3) = (del(1)-a*del(0))/(c-a)
        elseif (cwavetype==2) then
          !r2 is (1, s2, s2**2)
          beta(1) = (b*c*del(0) - (b+c)*del(1) +del(2))/
     &          (a**2- a*b - a*c + b*c)
          beta(2) = (-a*c*del(0) + (a+c)*del(1) -del(2))/
     &      (a*b - a*c- b**2 + b*c)
          beta(3) = (a*b*del(0) -(a+b)*del(1) + del(2))/
     &      (a*b - a*c - b*c + c**2)
        endif

      do mw=1,3
         do m=1,2
            fw(m,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

      !split the jump in phi (mom. flux) in middle wave into outer waves
      fw(2,1) = fw(2,1) + 0.5d0*fw(2,2)
      fw(2,3) = fw(2,3) + 0.5d0*fw(2,2)
      fw(2,2) = 0.d0
      !waves and fwaves for delta hum
      fw(4,1) = fw(1,1)*mL
      fw(4,3) = fw(1,3)*mR
      fw(4,2) = hmR*uR-hmL*uL - fw(4,1)- fw(4,3)

      !waves and fwaves for delta huv
      fw(3,1) = fw(1,1)*vL
      fw(3,3) = fw(1,3)*vR
      fw(3,2) = hvR*uR-hvL*uL -fw(3,1) -fw(3,3)

      !fwaves for delta p
      fw(5,1) = fw(1,1)*gammaL*rhoL*grav*dcos(theta1)
      fw(5,3) = fw(1,3)*gammaR*rhoR*grav*dcos(theta3)
      fw(5,2) = del(4) - fw(5,3) - fw(5,1)


      !fwaves for segregation
      seg_L = chiL*hL*uL*(1.0d0+(1.0d0-alpha_seg)*(1.0d0-chiL))
      seg_R = chiR*hR*uR*(1.0d0+(1.0d0-alpha_seg)*(1.0d0-chiR))
      fw(6,1) = fw(1,1)*chiL*(1.0+(1.0d0-alpha_seg)*(1.0d0-chiL))
      fw(6,3) = fw(1,3)*chiR*(1.0+(1.0d0-alpha_seg)*(1.0d0-chiR))
      fw(6,2) = seg_R - seg_L - fw(6,1) - fw(6,3)

        ! NaN / Inf diagnostic block
        if (.not.(hstarHLL == hstarHLL)) then
            write(*,*) 'NaN detected: hstarHLL'
        endif
        if (.not.(hustarHLL == hustarHLL)) then
            write(*,*) 'NaN detected: hustarHLL'
        endif
        if (.not.(hustarHLLn == hustarHLLn)) then
            write(*,*) 'NaN detected: hustarHLLn'
        endif
        if (.not.(deldelh == deldelh)) then
            write(*,*) 'NaN detected: deldelh'
        endif
        if (.not.(deldelhf == deldelhf)) then
            write(*,*) 'NaN detected: deldelhf'
        endif
        if (.not.(source2dx == source2dx)) then
            write(*,*) 'NaN detected: source2dx'
        endif
        if (.not.(s1s2bar == s1s2bar)) then
            write(*,*) 'NaN detected: s1s2bar'
        endif
        if (.not.(s1s2_denom == s1s2_denom)) then
            write(*,*) 'NaN detected: s1s2_denom'
        endif
        if (.not.(phiL_effective == phiL_effective)) then
            write(*,*) 'NaN detected: phiL_effective'
        endif
        if (.not.(phiR_effective == phiR_effective)) then
            write(*,*) 'NaN detected: phiR_effective'
        endif
        if (.not.(phi_eff == phi_eff)) then
            write(*,*) 'NaN detected: phi_eff'
        endif
        if (.not.(taudirUfrac == taudirUfrac)) then
            write(*,*) 'NaN detected: taudirUfrac'
        endif
        if (.not.(delbf == delbf)) then
            write(*,*) 'NaN detected: delbf'
        endif
        if (.not.(beta(1) == beta(1))) then
            write(*,*) 'NaN detected: beta(1)'
        endif
        if (.not.(beta(2) == beta(2))) then
            write(*,*) 'NaN detected: beta(2)'
        endif
        if (.not.(beta(3) == beta(3))) then
            write(*,*) 'NaN detected: beta(3)'
            write(*,*) 'del0,del1', del(0),del(1)
            write(*,*) 'a,b,c',a,b,c
        endif
        if (.not.(taudirR == taudirR)) then
            write(*,*) 'NaN detected: taudir'
            stop
        endif
        if (.not.(source2dx == source2dx)) then
            write(*,*) 'NaN detected: src'
        endif
        if (abs(source2dx) > 1d30) then
            write(*,*) 'Large value (Inf?) detected: src =', source2dx
            write(*,*) 'source2dxf=',source2dxf
            write(*,*) 's1s2ratio,bar', s1s2_ratio,s1s2bar
            write(*,*) 'a,b,c',a,b,c
            write(*,*) 'ss_delta', ss_delta
            write(*,*) 'hustarHLLn',hustarHLLn
            write(*,*) 'delb,delbf,deldelh,deldelhf',delb,delbf,
     &          deldelh,deldelhf
            write(*,*) 'taudirUfrac,phi_eff,tan_phi_max',
     &            taudirUfrac,phi_eff,tan_phi_max
            write(*,*) 'uedge,vedge', uedge,vedge,hL,hR
            write(*,*) 'taudirR',taudirR
            stop
        endif
        if (abs(source2dxf) > 1d30) then
            write(*,*) 'Large value (Inf?) detected: srcf =', source2dxf
            write(*,*) 'source2dx=',source2dx
            write(*,*) 's1s2ratio,bar', s1s2_ratio,s1s2bar
            write(*,*) 'a,b,c',a,b,c
            write(*,*) 'ss_delta', ss_delta
            write(*,*) 'hstarHLLn',hustarHLLn
            stop
        endif
        if (abs(dels) < 1d-30) then
            write(*,*) 'Large value (Inf?) detected: dels =', dels
        endif
        if (.not.(del(1) == del(1))) then
            write(*,*) 'NaN detected: del(1)'
        endif
        if (.not.(del(2) == del(2))) then
            write(*,*) 'NaN detected: del(2)'
        endif
        if (.not.(del(3) == del(3))) then
            write(*,*) 'NaN detected: del(3)'
        endif


      return
      end !subroutine riemann_dig2_aug_sswave_ez

c-----------------------------------------------------------------------

c=============================================================================
      subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &             maxiter,drytol,g)

      !determine the Riemann structure (wave-type in each family)


      implicit none

      !input
      double precision hL,hR,uL,uR,drytol,g
      integer maxiter

      !output
      double precision s1m,s2m
      logical rare1,rare2

      !local
      double precision hm,u1m,u2m,um,delu
      double precision h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
      integer iter



c     !Test for Riemann structure

      h_min=min(hR,hL)
      h_max=max(hR,hL)
      delu=uR-uL

      if (h_min.le.drytol) then
         hm=0.d0
         um=0.d0
         s1m=uR+uL-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hL)
         s2m=uR+uL-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hL)
         if (hL.le.0.d0) then
            rare2=.true.
            rare1=.false.
         else
            rare1=.true.
            rare2=.false.
         endif

      else
         F_min= delu+2.d0*(dsqrt(g*h_min)-dsqrt(g*h_max))
         F_max= delu +
     &         (h_max-h_min)*(dsqrt(.5d0*g*(h_max+h_min)/(h_max*h_min)))

         if (F_min.gt.0.d0) then !2-rarefactions

            hm=(1.d0/(16.d0*g))*
     &               max(0.d0,-delu+2.d0*(dsqrt(g*hL)+dsqrt(g*hR)))**2
            um=dsign(1.d0,hm)*(uL+2.d0*(dsqrt(g*hL)-dsqrt(g*hm)))

            s1m=uL+2.d0*dsqrt(g*hL)-3.d0*dsqrt(g*hm)
            s2m=uR-2.d0*dsqrt(g*hR)+3.d0*dsqrt(g*hm)

            rare1=.true.
            rare2=.true.

         elseif (F_max.le.0.d0) then !2 shocks

c           !root finding using a Newton iteration on dsqrt(h)===
            h0=h_max
            do iter=1,maxiter
               gL=dsqrt(.5d0*g*(1/h0 + 1/hL))
               gR=dsqrt(.5d0*g*(1/h0 + 1/hR))
               F0=delu+(h0-hL)*gL + (h0-hR)*gR
               dfdh=gL-g*(h0-hL)/(4.d0*(h0**2)*gL)+
     &                   gR-g*(h0-hR)/(4.d0*(h0**2)*gR)
               slope=2.d0*dsqrt(h0)*dfdh
               h0=(dsqrt(h0)-F0/slope)**2
            enddo
               hm=h0
               u1m=uL-(hm-hL)*dsqrt((.5d0*g)*(1/hm + 1/hL))
               u2m=uR+(hm-hR)*dsqrt((.5d0*g)*(1/hm + 1/hR))
               um=.5d0*(u1m+u2m)

               s1m=u1m-dsqrt(g*hm)
               s2m=u2m+dsqrt(g*hm)
               rare1=.false.
               rare2=.false.

         else !one shock one rarefaction
            h0=h_min

            do iter=1,maxiter
               F0=delu + 2.d0*(dsqrt(g*h0)-dsqrt(g*h_max))
     &                  + (h0-h_min)*dsqrt(.5d0*g*(1.d0/h0+1.d0/h_min))
               slope=(F_max-F0)/(h_max-h_min)
               h0=h0-F0/slope
            enddo

            hm=h0
            if (hL.gt.hR) then
               um=uL+2.d0*dsqrt(g*hL)-2.d0*dsqrt(g*hm)
               s1m=uL+2.d0*dsqrt(g*hL)-3.d0*dsqrt(g*hm)
               s2m=uL+2.d0*dsqrt(g*hL)-dsqrt(g*hm)
               rare1=.true.
               rare2=.false.
            else
               s2m=uR-2.d0*dsqrt(g*hR)+3.d0*dsqrt(g*hm)
               s1m=uR-2.d0*dsqrt(g*hR)+dsqrt(g*hm)
               um=uR-2.d0*dsqrt(g*hR)+2.d0*dsqrt(g*hm)
               rare2=.true.
               rare1=.false.
            endif
         endif
      endif

      return

      end ! subroutine riemanntype----------------------------------------------------------------


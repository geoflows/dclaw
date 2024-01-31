
c-----------------------------------------------------------------------
      subroutine riemann_dig2_aug_sswave_ez(ixy,meqn,mwaves,hL,hR,
     &         huL,huR,hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         thetaL,thetaR,phiL,phiR,sw,fw,w,wallprob,taudirL,
     &         taudirR,chiL,chiR,fsL,fsR,ilook)

      !-----------------------------------------------------------------
      ! solve the dig Riemann problem for debris flow eqn
      ! this is for 2d version
      !
      ! This solver is an extension of that described in:
      ! J. Comput. Phys. (6): 3089-3113, March 2008
      ! Augmented Riemann Solvers for the Shallow Equations,
      !                   with Steady States and Inundation
      !-----------------------------------------------------------------

      use geoclaw_module, only: grav, dry_tolerance
      use digclaw_module ! DIG: specify which variables.

      implicit none

*     !i/o
      integer ixy,meqn,mwaves,ilook

      double precision hL,hR,huL,huR,hvL,hvR,hmL,hmR,pL,pR
      double precision bL,bR,uL,uR,vL,vR,mL,mR,chiL,chiR,seg_L,seg_R
      double precision thetaL,thetaR,phiL,phiR
      double precision taudirL,taudirR,fsL,fsR
      logical wallprob


      double precision fw(meqn,mwaves),w(meqn,mwaves)
      double precision sw(mwaves)
      double precision psi(4)

*     !local
      integer m,mw,ks,mp
      double precision h,u,v,mbar
      double precision det1,det2,det3,detR
      double precision R(0:2,1:3),A(3,3),del(0:4)
      double precision beta(3)
      double precision rho,rhoL,rhoR,tauL,tauR,tau
      double precision kappa,kperm,compress,tanpsi,D,SN,sigbed,pmL,pmR
      double precision theta,gamma,eps,pm
      double precision sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
      double precision delb,s1m,s2m,hm,heL,heR,criticaltol
      double precision s1s2bar,s1s2tilde,hbar,source2dx,veltol1,veltol2
      double precision hstarHLL,deldelh,drytol,gmod,geps,tausource
      double precision raremin,raremax,rare1st,rare2st,sdelta,rho_fp,seg
      double precision gammaL,gammaR,theta1,theta2,theta3,vnorm,pmtanh01
      logical sonic,rare1,rare2
      logical rarecorrectortest,rarecorrector

      veltol1=1.d-6
      veltol2=0.d0
      criticaltol=1.d-6
      drytol = dry_tolerance

      do m=1,4
         psi(m) = 0.d0
      enddo


      ! DIG: move segregation into segeval() that
      ! adjusts rho_fp, etc.

      pmL = chiL
      pmR = chiR
      pm = 0.5*(chiL + chiR)
      pm = min(1.0d0,pm)
      pm = max(0.0d0,pm)
      if (alpha_seg==1.0d0) then
         seg = 0.0d0
      else
         seg = 1.0d0
      endif
      call calc_pmtanh(pm,seg,pmtanh01)
      rho_fp = (1.0d0-pmtanh01)*rho_f


      if (hL.ge.drytol.and.hR.ge.drytol) then
         call auxeval(hL,uL,vL,mL,pL,phiL,thetaL,
     &        kappa,SN,rhoL,tanpsi,D,tauL,sigbed,kperm,compress,pmL)
         call auxeval(hR,uR,vR,mR,pR,phiR,thetaR,
     &        kappa,SN,rhoR,tanpsi,D,tauR,sigbed,kperm,compress,pmR)
         h = 0.5d0*(hL + hR)
         v = 0.5d0*(vL + vR)
         mbar = 0.5d0*(mL + mR)
      elseif (hL.ge.drytol) then
         call auxeval(hL,uL,vL,mL,pL,phiL,thetaL,
     &        kappa,SN,rhoL,tanpsi,D,tauL,sigbed,kperm,compress,pmL)
         call auxeval(hL,uL,vL,mL,pL,phiL,thetaL,
     &        kappa,SN,rhoR,tanpsi,D,tauR,sigbed,kperm,compress,pmL)
         tauR=0.5d0*tauL
         h = 0.5d0*hL
         v = vL
         mbar = mL
      else
         call auxeval(hR,uR,vR,mR,pR,phiR,thetaR,
     &        kappa,SN,rhoL,tanpsi,D,tauL,sigbed,kperm,compress,pmR)
         call auxeval(hR,uR,vR,mR,pR,phiR,thetaR,
     &        kappa,SN,rhoR,tanpsi,D,tauR,sigbed,kperm,compress,pmR)
         tauL=0.5d0*tauR
         h = 0.5d0*hR
         v = vR
         mbar = mR
      endif

      tauL = fsL*tauL
      tauR = fsR*tauR
      rho = 0.5d0*(rhoL + rhoR)
      tau = 0.5d0*(tauL + tauR)
      theta = 0.5d0*(thetaL + thetaR)
      gamma = 0.25d0*(rho_fp + 3.0d0*rho)/rho
      gammaL = 0.25d0*(rho_fp + 3.0d0*rhoL)/rhoL
      gammaR = 0.25d0*(rho_fp + 3.0d0*rhoR)/rhoR

      eps = kappa + (1.d0-kappa)*gamma

      gmod = grav*dcos(theta)
      geps = gmod*eps

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
      u = uhat

      call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,
     &                                          1,drytol,geps)
      sw(1)= min(sw(1),s2m) !Modified Einfeldt speed
      sw(3)= max(sw(3),s1m) !Modified Einfeldt speed
      sw(2) = 0.5d0*(sw(3)+sw(1))

      hstarHLL = max((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.d0) ! middle state in an HLL solve
c     !determine the middle entropy corrector wave------------------------
      rarecorrectortest = .true.
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

      delb=(bR-bL)!#kappa

      !determine ss-wave
      hbar =  0.5d0*(hL+hR)
      s1s2bar = 0.25d0*(uL+uR)**2- gmod*hbar
      s1s2tilde= max(0.d0,uL*uR) - gmod*hbar

c     !find if sonic problem
      sonic=.false.
      if (dabs(s1s2bar).le.criticaltol) sonic=.true.
      if (s1s2bar*s1s2tilde.le.criticaltol) sonic=.true.
      if (s1s2bar*sE1*sE2.le.criticaltol) sonic = .true.
      if (min(dabs(sE1),dabs(sE2)).lt.criticaltol) sonic=.true.
      if (sE1.lt.0.d0.and.s1m.gt.0.d0) sonic = .true.
      if (sE2.gt.0.d0.and.s2m.lt.0.d0) sonic = .true.
      if ((uL+dsqrt(geps*hL))*(uR+dsqrt(geps*hR)).lt.0.d0) sonic=.true.
      if ((uL-dsqrt(geps*hL))*(uR-dsqrt(geps*hR)).lt.0.d0) sonic=.true.

      if (sonic) then
         source2dx = -gmod*hbar*delb
      else
         source2dx = -delb*gmod*hbar*s1s2tilde/s1s2bar
      endif

      source2dx=min(source2dx,gmod*max(-hL*delb,-hR*delb))
      source2dx=max(source2dx,gmod*min(-hL*delb,-hR*delb))

      if (dabs(u).le.veltol2) then
         source2dx=-hbar*gmod*delb
      endif

c     !find bounds in case of critical state resonance, or negative states
c     !find jump in h, deldelh

      if (sonic) then
         deldelh =  -delb
      else
         deldelh = delb*gmod*hbar/s1s2bar
      endif
c     !find bounds in case of critical state resonance, or negative states
      if (sE1.lt.-criticaltol.and.sE2.gt.criticaltol) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE2)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE1)
      elseif (sE1.ge.criticaltol) then
         deldelh = min(deldelh,hstarHLL*(sE2-sE1)/sE1)
         deldelh = max(deldelh,-hL)
      elseif (sE2.le.-criticaltol) then
         deldelh = min(deldelh,hR)
         deldelh = max(deldelh,hstarHLL*(sE2-sE1)/sE2)
      endif


*     !determine R
      R(0,2) = 0.d0
      R(1,2) = 0.d0
      R(2,2) = 1.d0

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
      endif

      !determine del
      del(0) = hR- hL - deldelh
      del(1) = huR - huL
      del(2) = hR*uR**2 + 0.5d0*kappa*gmod*hR**2 -
     &      (hL*uL**2 + 0.5d0*kappa*gmod*hL**2)
      del(2) = del(2) + (1.d0-kappa)*h*(pR-pL)/rho
      del(3) = pR - pL - gamma*rho*gmod*deldelh
      del(4) = -gamma*rho*gmod*u*(hR-hL) + gamma*rho*gmod*del(1)
     &         + u*(pR-pL)

*     !determine the source term

      if (ixy.eq.1) then
         ! DIG: note that theta = 0.0 unless bed_normal is true. For now, assume bed_normal is false. Resolve if dx is needed later.
         source2dx = source2dx !+ dx*hbar*grav*dsin(theta)
         ! DIG: this is the only place dx is needed
         ! until fixed, bed_normal = 1 yields error in make .data (1/30/24)
      endif

      vnorm = sqrt(uR**2 + uL**2 + vR**2 + vL**2)
      if (vnorm>0.0d0) then

         tausource =  0.0d0 ! if vnorm>0 then src2 handles friction.

      elseif (0.5d0*abs(taudirR*tauR/rhoR + tauL*taudirR/rhoL)
     &      .gt.abs(del(2) - source2dx)) then


!       DIG Symmetry: should this be RRR, LLL? It is symmetric in the line above
!       that is commented out. Leaving as is b/c it is the same in dclaw4 and dclaw5
!       KRB&MJB - 1/12/24

         ! no failure of static material
         tausource = del(2) - source2dx
         del(1) = 0.0d0
         del(0) = 0.0d0
         del(4) = 0.0d0
      else
         ! failure of static material
         tausource = 0.5d0*((taudirR*tauR/rhoR)+(tauL*taudirR/rhoL))!*dx
         tausource = dsign(tausource,del(2)-source2dx)
      endif

      if (wallprob) then
         tausource = 0.0d0
      endif

      del(2) = del(2) - source2dx  - tausource

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


*     !R beta = del
*     !gauss routine replaces del with beta and R with it's inverse
      !want to keep R, so replacing with A
      do mw=0,2
         beta(mw+1) = del(mw)
         do m=0,2
            A(m+1,mw+1)=R(m,mw+1)
         enddo
      enddo

      call gaussj(A,3,3,beta,1,1)

      do mw=1,3
         do m=1,2
            fw(m,mw) = beta(mw)*R(m,mw)
         enddo
      enddo

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

c=======================================================================
c  file: gaussj
c  routine: gaussj
c  solves linear system ax=b directly using Guassian Elimination
c  a(1:n,1:n) is the input matrix of size np by np
c  b(1:n,1:m) is the right hand side (can solve for m vectors)
c  On output a is replaced by it's inverse, and b replaced by x.

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(np,np),b(np,mp)
      PARAMETER (NMAX=50)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      DOUBLE PRECISION big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.d0
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (dabs(a(j,k)).ge.big)then
                  big=dabs(a(j,k))
                  irow=j
                  icol=k
                endif

              else if (ipiv(k).gt.1) then
                write(*,*) 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.d0) write(*,*) 'singular matrix in gaussj'

        pivinv=1./a(icol,icol)
        a(icol,icol)=1.d0
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.d0
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then

          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

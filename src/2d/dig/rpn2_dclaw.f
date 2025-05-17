
c D-Claw specific core file
c This file is a modified version of
c clawpack/geoclaw/src/2d/shallow/rpn2.f 
c
c
c rpn2: solve normal Riemann problems at cell interfaces (x- or y- oriented)
c       this subroutine is given an entire row or column of cells 
c       from a grid. 
c
c       The standard clawpack version is primarily empty and intended
c       to be user-written depending on model equations.
c       See (Leveque, 2002) (D-Claw bibliography) for details.
c
c       This version is similar to the geoclaw version in that it
c       deals with a depth-averaged flow over topography by testing 
c       wall-boundary ghost Riemann problems at wet-dry interface.
c       For related details see (George, 2008;2011) in D-Claw reference list  
c
c       this routine calls the subroutine riemann_dig2_aug_sswave_ez
c       found in the file riemannsolvers_dclaw.f to solve the actual
c       Riemann problem at each interface, after the wall-boundary tests.
c======================================================================
       subroutine rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,
     &                 ql,qr,auxl,auxr,fwave,s,amdq,apdq)
c======================================================================
c
c GeoClaw version solves normal Riemann problems for the 2D SHALLOW WATER equations
c     with topography:
c     #        h_t + (hu)_x + (hv)_y = 0                           #
c     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
c     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #

c DClaw versions solves related depth-averaged system, see (Iverson & George, 2014)

c On input, ql contains the state vector at the left edge of each cell
c     qr contains the state vector at the right edge of each cell
c
c This data is along a slice in the x-direction if ixy=1
c     or the y-direction if ixy=2.

c  Note that the i'th Riemann problem has left state qr(i-1,:)
c     and right state ql(i,:)
c  From the basic clawpack routines, this routine is called with
c     ql = qr
c
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      # This Riemann solver is for debris flow equations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use geoclaw_module, only: grav, dry_tolerance
      use geoclaw_module, only: earth_radius, deg2rad
      use amr_module, only: mcapa

      use digclaw_module, only: bed_normal,i_theta,qfix
      use digclaw_module, only: i_fsphi,i_phi,i_taudir_x,i_taudir_y
      use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_hf
      use digclaw_module, only: i_ent,i_hs

      implicit none

      !input
      integer maxm,meqn,maux,mwaves,mbc,mx,ixy

      double precision  fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      double precision  s(mwaves, 1-mbc:maxm+mbc)
      double precision  ql(meqn, 1-mbc:maxm+mbc)
      double precision  qr(meqn, 1-mbc:maxm+mbc)
      double precision  apdq(meqn,1-mbc:maxm+mbc)
      double precision  amdq(meqn,1-mbc:maxm+mbc)
      double precision  auxl(maux,1-mbc:maxm+mbc)
      double precision  auxr(maux,1-mbc:maxm+mbc)

      !local
      integer m,i,mw,maxiter,mhu,nhv,waves
      !double precision dtcom,dxcom,dycom,tcom
      double precision wall(3),fw(6,3),sw(3),wave(6,3)
      !logical entropy(5)
      logical rare1,rare2,wallprob,drystate
      !logical entropycorr1,entropycorr2

      double precision drytol,gmod,veltol
      double precision hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL
      double precision pR,pL,hmL,hmR,mL,mR,phi_bedL,phi_bedR
      double precision hstar,hstartest,s1m,s2m,bL,bR
      double precision dxdc,taudirL,taudirR,dx
      double precision theta,thetaL,thetaR
      double precision chiHL,chiHR,chiL,chiR,fsL,fsR
      double precision gz,gzL,gzR,rhoR,rhoL

      gmod=grav
      veltol = 1.d-3
      waves = 3
      drytol = dry_tolerance

      !loop through Riemann problems at each grid cell
      do i=2-mbc,mx+mbc

      !-----------------------Initializing-----------------------------------
         !inform of a bad Riemann problem from the start
         !if((qr(i_h,i-1).lt.0.d0).or.(ql(i_h,i) .lt. 0.d0)) then
         !   write(*,*) 'Negative input: hl,hr,i=',qr(i-1,1),ql(i,1),i
         !endif

         !Initialize Riemann problem for grid interface

         do mw=1,mwaves
            s(mw,i)=0.d0
            !entropy(mw)=.false. ! DIG: this is a difference from old-new rpn2, but it does not seem to be used.
            do m=1,meqn
               fwave(m,mw,i)=0.d0
            enddo
         enddo
         do mw=1,waves
            sw(mw) = 0.d0
            do m=1,6
               wave(m,mw) = 0.d0
               fw(m,mw) = 0.d0
            enddo
         enddo

         !skip problem if in a completely dry area
         if (qr(i_h,i-1).le.drytol.and.ql(i_h,i).le.drytol) then
            go to 30
         endif

c        !set normal direction
         if (ixy.eq.1) then
            mhu=i_hu
            nhv=i_hv
            !dx = dxcom
            taudirR = auxl(i_taudir_x,i)
            taudirL = auxr(i_taudir_x,i-1)
         else
            mhu=i_hv
            nhv=i_hu
            !dx = dycom
            taudirR = auxl(i_taudir_y,i)
            taudirL = auxr(i_taudir_y,i-1)
         endif

         fsL = auxr(i_fsphi,i-1)
         fsR = auxl(i_fsphi,i)

         if (bed_normal.eq.1) then
            thetaL = auxr(i_theta,i-1)
            thetaR = auxl(i_theta,i)
            theta = 0.5d0*(thetaL+thetaR)
            gz=grav*dcos(0.5d0*(thetaL+thetaR))
            gzL=grav*dcos(thetaL)
            gzR=grav*dcos(thetaR)
         else
            thetaL = 0.d0
            thetaR = 0.d0
            theta = 0.d0
            gz = grav
            gzL = grav
            gzR = grav
         endif

         !zero (small) negative values if they exist and set velocities
         call qfix(ql(i_h,i),ql(mhu,i),ql(nhv,i),ql(i_hm,i),
     &             ql(i_pb,i),ql(i_hchi,i),ql(i_hf,i),
     &             uR,vR,mR,chiR,rhoR,gzR)

         call qfix(qr(i_h,i-1),qr(mhu,i-1),qr(nhv,i-1),qr(i_hm,i-1),
     &             qr(i_pb,i-1),qr(i_hchi,i-1),qr(i_hf,i-1),
     &             uL,vL,mL,chiL,rhoL,gzL)

         !Riemann problem variables
         hL = qr(i_h,i-1)
         hR = ql(i_h,i)
         huL = qr(mhu,i-1)
         huR = ql(mhu,i)
         hvL=qr(nhv,i-1)
         hvR=ql(nhv,i)
         hmL = qr(i_hm,i-1)
         hmR = ql(i_hm,i)
         pL = qr(i_pb,i-1)
         pR = ql(i_pb,i)
         bL = auxr(1,i-1) - auxr(i_ent, i-1) + qr(i_hs,i-1)
         bR = auxl(1,i)   - auxl(i_ent, i)   + ql(i_hs,i)
         phi_bedL = auxr(i_phi,i-1)
         phi_bedR = auxl(i_phi,i)
         chiHL = qr(i_hchi,i-1)
         chiHR = ql(i_hchi,i)

         !test for wall problem vs. inundation problem
         do mw=1,waves
            wall(mw) = 1.d0
         enddo
         drystate=.false.
         wallprob = .false.
         if (hR.le.drytol) then
            hR = 0.d0
            pR = 0.d0
            mR=mL ! KRB added this back 1/11/2024
            chiR = chiL
            drystate=.true.
            call riemanntype(hL,hL,uL,-uL,hstar,s1m,s2m,
     &                                 rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hL,hstar)
            if (hstartest+bL.lt.bR) then !right state should become ghost values that mirror left for wall problem
c                bR=hstartest+bL
               wall(2)=0.d0
               wall(3)=0.d0
               wallprob=.true.
               hR=hL
               huR=-huL
               hvR=hvL
               hmR=hmL
               bR=bL
               uR=-uL
               vR=vL
               mR=mL
               pR=pL
               chiHR=chiHL
               chiR = chiL
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hL+bL.lt.bR) then
               !bR=hL+bL
            endif
         elseif (hL.le.drytol) then ! right surface is lower than left topo
            hL = 0.d0
            pL = 0.d0
            mL = mR ! KRB added this back 1/11/2024
            chiL= chiR
            drystate=.true.
            call riemanntype(hR,hR,-uR,uR,hstar,s1m,s2m,
     &                                  rare1,rare2,1,drytol,gmod)
            hstartest=dmax1(hR,hstar)
            if (hstartest+bR.lt.bL) then  !left state should become ghost values that mirror right
c               bL=hstartest+bR
               wall(1)=0.d0
               wall(2)=0.d0
               wallprob=.true.
               hL=hR
               huL=-huR
               hvL=hvR
               hmL=hmR
               mL = mR
               bL=bR
               uL=-uR
               vL=vR
               pL=pR
               chiHL=chiHR
               chiL = chiR
               !thetaL = 0.d0
               !thetaR = 0.d0
            !elseif (hR+bR.lt.bL) then
               !bL=hR+bR
            endif
         endif
         !--------------------end initializing...finally----------
         !solve Riemann problem.

         maxiter = 1

          ! current dclaw Riemann solver
          call riemann_dig2_aug_sswave_ez(ixy,6,3,hL,hR,huL,huR,
     &         hvL,hvR,hmL,hmR,pL,pR,bL,bR,uL,uR,vL,vR,mL,mR,
     &         thetaL,thetaR,phi_bedL,phi_bedR,sw,fw,wallprob,
     &         taudirL,taudirR,chiL,chiR,fsL,fsR,i)

c        !eliminate ghost fluxes for wall
         do mw=1,3
            sw(mw)=sw(mw)*wall(mw)
            do m=1,6
               fw(m,mw)=fw(m,mw)*wall(mw)
            enddo
         enddo

c============segregation================================================

         s(1,i) = sw(1)
         s(2,i) = sw(2)
         s(3,i) = sw(2)
         s(4,i) = sw(2)
         s(5,i) = sw(3)

         fwave(1,1,i) =   fw(1,1)
         fwave(mhu,1,i) = fw(2,1)
         fwave(nhv,1,i) = fw(3,1)
         fwave(4,1,i)   = fw(4,1)
         fwave(5,1,i) =   fw(5,1)
         fwave(6,1,i) =   fw(6,1)

         fwave(1,5,i) =   fw(1,3)
         fwave(mhu,5,i) = fw(2,3)
         fwave(nhv,5,i) = fw(3,3)
         fwave(4,5,i)   = fw(4,3)
         fwave(5,5,i) =   fw(5,3)
         fwave(6,5,i) =   fw(6,3)

         fwave(1,2,i) =   fw(1,2)
         fwave(mhu,2,i) = fw(2,2)
         fwave(nhv,2,i) = fw(3,2)
         fwave(4,2,i)   = 0.0d0
         fwave(5,2,i) =  0.0d0
         fwave(6,2,i) = fw(6,2)

         fwave(1,3,i) =   0.0d0
         fwave(mhu,3,i) = 0.0d0
         fwave(nhv,3,i) = 0.0d0
         fwave(4,3,i)   = fw(4,2)
         fwave(5,3,i) =  0.0d0
         fwave(6,3,i) =  0.0d0

         fwave(1,4,i) =   0.0d0
         fwave(mhu,4,i) = 0.0d0
         fwave(nhv,4,i) = 0.0d0
         fwave(4,4,i)   = 0.0d0
         fwave(5,4,i) =  fw(5,2)
         fwave(6,4,i) =  0.0d0

 30      continue
      enddo


c==========Capacity for mapping from latitude longitude to physical space====
        if (mcapa.gt.0) then
         do i=2-mbc,mx+mbc
          if (ixy.eq.1) then
             dxdc=(earth_radius*deg2rad)
          else
             dxdc=earth_radius*cos(auxl(3,i))*deg2rad
          endif

          do mw=1,mwaves
c             if (s(mw,i) .gt. 316.d0) then
c               # shouldn't happen unless h > 10 km!
c                write(6,*) 'speed > 316: i,mw,s(mw,i): ',i,mw,s(mw,i)
c                endif
               s(mw,i)=dxdc*s(mw,i)
               do m=1,meqn
                  fwave(m,mw,i)=dxdc*fwave(m,mw,i)
               enddo
          enddo
         enddo
        endif

c===============================================================================


c============= compute fluctuations=============================================
         amdq(1:meqn,:) = 0.d0
         apdq(1:meqn,:) = 0.d0
         do i=2-mbc,mx+mbc
            do  mw=1,mwaves
               if (s(mw,i) < 0.d0) then
                     amdq(1:meqn,i) = amdq(1:meqn,i)
     &                              + fwave(1:meqn,mw,i)
               else if (s(mw,i) > 0.d0) then
                  apdq(1:meqn,i)  = apdq(1:meqn,i)
     &                          + fwave(1:meqn,mw,i)
               else


! DIG: 1/11/2024: KRB & MJB close comparing dclaw4 and dclaw5. These next four
! lines were commented out in dclaw4. We probably want them because in
! geoclaw5 they are used. Keeping commented for now to do a debug ensuring
! identical behavior of dclaw4 and dclaw5.

! if there is a wave speed of zero (s(mw,i)) and there is a jump in the
! flux (fwave(1:meqn,mw,i)>0) then split the fwave value equally between the
! two directions (plus and minus)

! Because of how the source term is handled, DLG does not think we should add
! this back in. (1/30/2024) Also unclear how often s(mw,i)=0 and
! fwave(1:meqn,mw,i)>0 occurs.

!                 amdq(1:meqn,i) = amdq(1:meqn,i)
!     &                              + 0.5d0 * fwave(1:meqn,mw,i)
!                 apdq(1:meqn,i) = apdq(1:meqn,i)
!     &                              + 0.5d0 * fwave(1:meqn,mw,i)
               endif
            enddo
         enddo


!--       do i=2-mbc,mx+mbc
!--            do m=1,meqn
!--                write(51,151) m,i,amdq(m,i),apdq(m,i)
!--                write(51,152) fwave(m,1,i),fwave(m,2,i),fwave(m,3,i)
!--151             format("++3 ampdq ",2i4,2e25.15)
!--152             format("++3 fwave ",8x,3e25.15)
!--            enddo
!--        enddo

      return
      end subroutine

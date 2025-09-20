! ============================================================================
!  Subroutines to integration a portion of the source term
! ============================================================================

    ! This file contains routines for integrating part of source term for m, p
    ! and h. It does not affect u,v, but hu,hv will change due to change in h.
    !
    ! Depending on the value of srcmethod one of the below routines will be
    ! called from src2() with timestep dtk.
    !
    ! Each of the below subroutines integrates:
    !   dm/dt = f_1(p,m)
    !   dp/dt = f_2(p,m)
    ! invariant (rho h)_N+1 = (rho h)_N is exactly maintained
    ! depending on routine.
    ! Note that as m is updated, h = C/rho(m) => h is updated meaning that all variables (h,hu,hv,hm,p) are affected.
    !
    ! three alternatives are provided:
    !   if srcmethod = 0
    !       relax rho h constraint to more easily integrate p (new D)
    !       update h,hm,hu,hv using the same D which only approx maintains rhoh
    !       method that was implemented in 'old dclaw'.
    !       implemented in subroutine mp_update_relax_Dclaw4 below
    !   if srcmethod = 1:
    !       satisfy exactly rho h = constant 
    !       integrate hm and p same as srcmethod = 0.
    !       Then redefine h(m) to give new hu,hv
    !       (this method is experimental)
    !       implemented in subroutine mp_update_relax_Dclaw4 below
    !   if srcmethod = 2:
    !       stay on rho h = constant manifold, integrate stiff
    !       terms for m and p. Then redefine h(m) to give new hu,hv,hm
    !       (this method is experimental)
    !       implemented in subroutine mp_update_FE_4quad below
    !


subroutine mp_update_FE_4quad(dt,h,u,v,m,p,chi,rhoh,gz,dtk)
    ! used if srcmethod = 2
    !====================================================================
    ! subroutine mp_update_FE_4quad: integrate dp_exc/dt,dm/dt by a hybrid
    ! explicit integration that depends on the initial quadrant of phase space.
    ! Basic idea which conforms to phase space vector field: for each quadrant
    ! of phase space (divided by p=p_eq(m_c) and m=m_c, where m_c is m_eqn(p_eq))
    ! there is only one open
    ! interior (physically admissible) boundary that can be crossed. These boundaries
    ! are such that a physically admissible solution must proceed clockwise (m horizontal
    ! axis, p vertical) if it changes quadrants. 1 (UL)-> 2(UR)-> 3 (LR)-> 4(LL).
    ! An open boundary in 1 quadrant is a closed boundary in then next, ie. m=m_c
    ! belongs to quads 2 & 4, p=p_eq belongs to 1 & 3. Explicit solution to
    ! non-homogeneous exponential ode d/dt(m,p_exc) is used with dtk<=dt  such that only
    ! the open/admissible boundary can be crossed in a given substep. For a given
    ! substep dt, either (a) sol remains interior to quadrant (dtk=dt), (b) sol.
    ! reaches physically inadmissible boundary dtk<=dt (c) solution reaches open
    ! quadrant boundary (dtk<=dt) that belongs to next quadrant for next substep.
    ! Quadrants 1 & 3 contain the p and m nullclines. Substeps are taken based on
    ! the initial position wrt these nullclines, because solution is most prone
    ! to oscillations if a substep overshoots a nullcline significantly. Solution
    ! converges to p-nullcline envelope when it should, except when above and below
    ! the p-nullcline on left and right of m_eq respectively,the most problematic regions.
    ! m-nullcline is simply p=p_eq(m). p-nullcline is complicated but = m-nullcline if vnorm = 0.
    ! Otherwise rises toward (above?) lithostatic pressure in left-half plane and 0 (negative?)
    ! in right-half plane as v increases and kperm decreases.


    ! analytical solution of p_exc (t) is
    ! p_exc (t) = [(p_exc0 - (cd0/kp0)) exp (-kp0*t)] + (cd0/kp0)
    ! at long time p_exc relaxes to cd0/kp0


    ! nullclines:
    ! p = p_eq is the m-nullcline
    ! m = m_eq is not the p-nullcline. rather p-nullcline is where tanpsi
    ! pressure change is balanced by relaxation (-kp*p_exc + c_d = 0)
    ! p-nullcline located in UL (q1) and LR (q3) in between p=p_eq and m=m_eq


    ! what is limiting in each region


!     region  |  p limit | m limit
!     ----------------------------
!     1A |  |
!     1B |  |
!     1C |  |
!     2  |  |
!     3A |  |
!     3B |  |
!     3C |  |
!     4  |  |



    !====================================================================

       use digclaw_module, only: rho_f,rho_s,mu,setvars,qfix,qfix_cmass,m_crit,delta

       implicit none

       !i/o
       real(kind=8), intent(inout) :: h,m,p
       real(kind=8), intent(in)  :: u,v,rhoh,dt,chi
       real(kind=8), intent(in)  :: gz
       real(kind=8), intent(out) :: dtk

       !local
       real(kind=8) :: h0,p0,m_0,p_eq0,p_exc0,vnorm,m_eq
       real(kind=8) :: rho,rho0,tanpsi,tau,kperm
       real(kind=8) :: alphainv,p_exc,dtr,dtm,dtp,dts,p_exc_ave
       real(kind=8) :: km0,kp0,c_d0
       real(kind=8) :: hu,hv,hm
       real(kind=8) :: rho_c,h_c,p_eq_c,m_c,m_lower,m_upper
       real(kind=8) :: m_c0,p_eq_c0,sig_c,Nd,Nn,normc,shear,convtol,m_eq1

       integer :: debugloop,iter,itermax,quad0,quad1
       logical :: outquad,debug

       debugloop = 0
       debug = .false.
       outquad = .true.
       vnorm = sqrt(u**2 + v**2)

       !explicit integration (hybrid FE and explicit exponential solution)---------------------------------------------------
       ! q1 = q0 + dtk*f(q0)

       if (debug) then
          write(*,*) '------------SRC MPUPDATE STARTING CONDITIONS ---------->>>>>>>>'
          write(*,*) 'h:', h
          write(*,*) 'u:', u
          write(*,*) 'v:', v
          write(*,*) 'm:', m
          write(*,*) 'p:', p

          rho0 = m*(rho_s-rho_f)+rho_f
          write(*,*) 'rho: ', rho0
          write(*,*) 'p/hydro: ', p/(rho_f*gz*h)
          write(*,*) 'p/litho: ', p/(rho0*gz*h)

          write(*,*) 'dtr: ', dt ! lower this will be set
          write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
       endif


       call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)

       ! Save starting conditions
       h0 = h
       m_0 = m
       rho0 = m_0*(rho_s-rho_f)+rho_f
       p0 = p
       p_eq0 = rho_f*gz*h0
       p_exc0 = p0 - p_eq0

      if (debug) then
          write(*,*) '------------SRC MPUPDATE initial setvars results ---------->>>>>>>>'
          write(*,*) 'kperm:', kperm
          write(*,*) 'alphainv:', alphainv
          write(*,*) 'tanpsi:', tanpsi
          write(*,*) 'm_eq:', m_eq
          write(*,*) 'tau:', tau
          write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
       endif

       ! initial recording of timesteps.
       ! definition of different timesteps used
       ! dtk is : time step used by this call of mpupdate (calculated internally)
       ! dtr is : dt remaining at the beginning of the call
       ! other dts used below
       ! dtp : upper bound for dtk based on pressure evolution equation
       ! dtm : upper bound for dtk based on m evolution equation

       dtk = 0.d0
       dtr = dt

       quad0=0
       quad1=0

       !determine coefficients for update
       ! dm/dt = (2*k*rho^2/mu(rhoh)^2)*p_exc*m = km0 * p_exc * m
       ! dp_p_exc/dt = -kp0*p_exc + c_d see George & Iverson 2014
       ! equations from 2014 recast in terms of p_exc.

       km0 = ((2.d0*kperm*rho0**2)/(mu*rhoh**2))
       c_d0 = -3.d0*vnorm*(alphainv*rho0/(rhoh))*tanpsi
       kp0 = (kperm/(h0*mu))*(3.d0*alphainv*rho0/rhoh - 1.5d0*rho_f*gz*h0*(rho0-rho_f)/rhoh)


      if (debug) then
          write(*,*) '------------SRC MPUPDATE initial coefficients ---------->>>>>>>>'
          write(*,*) 'p_exc0: ', p_exc0
          write(*,*) 'km0:', km0
          write(*,*) 'c_d0:', c_d0
          write(*,*) 'kp0:', kp0
          write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
       endif

       ! kp0 is expected to be greater or equal to zero if alphamethod = 1
       ! but not if alphamethod = 0 (old style).
       ! 9/15/2025 - kp0 observed as -1e-13 with alphamethod 1

       if (debug.and.kp0<0.d0) then
          write(*,*) '------------SRC WARNING: kp0<0 ---------->>>>>>>>'
          write(*,*) 'kp0<0:', kp0
          write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
       endif

       kp0 = max(kp0,1.d-3) ! shouldn't happen with alphamethod=1 but prevent
       ! small rounding error

       ! if m == 0 relax to hydrostatic, call qfix_cmass and return
       if (m==0.d0) then
          p_exc = p_exc0*exp(-kp0*dtr)
          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dtr
          return
       endif

       ! if at critical point dq/dt = 0, return
       if ((abs(p_exc0/(rho_f*gz*h0))<1.d-12).and.(abs(m-m_eq)<1.d-12)) then
          dtk = dtr
          return
       endif

       ! two cases, described within loop.
       if ((vnorm==0.d0).or.(m==m_eq)) then

          ! in all cases, pressure relaxes.
          p_exc = p_exc0*exp(-kp0*dtm)

          ! if vnorm is zero -> pressure relaxes toward equilibrum
          ! (no dilatancy effects on pressure.
          ! m evolves per dmdt. calculate both, then return

          ! if m is equal to m_eq and m is less than 0.3 (value used)
          ! in setvars as the threshold of setting m_eq to m, then
          ! m should not change. Just relax pressure and return

          if ((m==m_eq).and.(m<0.3d0)) then
            ! keep m the same, do nothing
          else
            ! covers both vnorm == 0 AND m==meqn with m>0.3
            m = m_0*exp(km0*p_exc*dtm)
          endif

          ! fix q vars, update time and return
          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dtm

          return
       endif

       !find critical point (mc,peq_c) which defines quadrant origin
       !requires fixed-point iteration, should converge quickly
       m_c = m_crit
       rho_c = m_c*(rho_s-rho_f) + rho_f
       h_c = rhoh/rho_c
       p_eq_c = rho_f*gz*h_c
       shear = 2.d0*vnorm/h0

       ! shear should be guaranteed to be positive b/c of c_d0 check above.

       convtol = 1.d-16
       itermax = 100
       do iter=1,itermax
          m_c0 = m_c
          p_eq_c0 = p_eq_c
          rho_c = m_c*(rho_s-rho_f) + rho_f
          h_c = rhoh/rho_c
          p_eq_c = rho_f*gz*h_c
          sig_c = max(rhoh*gz-p_eq_c,0.d0)
          Nd = rho_s*(shear*delta)**2 + sig_c
          Nn = mu*shear
          m_c = m_crit*(sqrt(Nd)/(sqrt(Nd)+sqrt(Nn)))
          normc = (1.d2*(m_c - m_c0))**2 + (1.d-3*(p_eq_c-p_eq_c0))**2
          if (normc<convtol) then
             exit
          endif
          if (iter>10.and.debug) then
             write(*,*) '--------------------------------------------------'
             write(*,*) 'SRC2 WARNING: fixed point iterations:', iter
             write(*,*) 'norm for convergence:', normc
          endif
       enddo
       if ((iter>=itermax).and.(normc>1.d0)) then
             m_c = m_crit
             rho_c = m_c*(rho_s-rho_f) + rho_f
             h_c = rhoh/rho_c
             p_eq_c = rho_f*gz*h_c
       else !update from m_c from final iteration for consistent quadrant definition
            rho_c = m_c*(rho_s-rho_f) + rho_f
            h_c = rhoh/rho_c
            p_eq_c = rho_f*gz*h_c
       endif
       !

       !determine quadrant of initial solution in state space

       if ((m<m_c).and.(p>=p_eq_c)) then
          quad0=1
       elseif ((m>=m_c).and.(p>p_eq_c)) then 
          quad0=2
       elseif  ((m>m_c).and.(p<=p_eq_c)) then
          quad0=3
       elseif ((m<=m_c).and.(p<p_eq_c)) then
          quad0=4
       endif

       ! dig: div0
       ! at this point it is not guaranteed that p_exc is nonzero.
       ! it may also not be guaranteed that kp0 is nonzero

       select case (quad0)

       case(1) !UL quadrant, variable material and p_exc
          debugloop = 1000
          !loose
          !integration should only cross right boundary (m=m_c)
          !p can increase or decrease if below/above p-nullcline (which is left of m=m_eq, above p=p_eq)
          !m can increase or decrease if above/below m-nullcline (p=p_eq)

          ! first consider region 1A - in UL quadrant, and just below the m-nullcline
          ! in this region, p is increasing and m is decreasing
          if (p_exc0<-1.d-3*rhoh*gz) then
             debugloop=debugloop + 100

             ! set dtp so that p_exc stays below zero.
             if (c_d0>0.d0) then !don't exceed p_exc = 0 until next substep
                dtp = -(1.d0/kp0)*log(c_d0/(c_d0-kp0*p_exc0))
             else
                dtp = dtr
             endif

             ! dtm limit that ensures m does not change more than
             ! 2% or go below 0.0
             ! it uses exponential decay so dtm ensures no change greater than 2%
             dtm = log(0.98d0)/(kp0*p_exc0)

             dts = min(min(dtr,dtp), dtm)

             !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
             dts = max(dts,0.d0)
             p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
             m = m_0*exp(dts*km0*p_exc)

          ! next consider region 1B - in UL quadrant, below the p-nullcline
          ! this includes the portion of UL that is the m-nullcline.
          ! in this region, p increases and m increases
          elseif ((-kp0*p_exc0+c_d0)>0.d0) then
             debugloop=debugloop + 200
             !bound time step
             ! to max allowed (lithostatic,rhoh*g) if nullcline exceeds it
             if ((c_d0/kp0)>(rhoh*gz-p_eq_c)) then
                debugloop=debugloop + 10
                dtp = min(dtr,-(1.d0/kp0)*log((-kp0*rhoh*gz+kp0*p_eq_c+c_d0)/(-kp0*p_exc0+c_d0)))
             !bound time step so not to exceed nullcline
             else
                debugloop=debugloop + 20
                dtp = dtr
             endif

             p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
             p_exc_ave = 0.5d0*(p_exc+p_exc0)
             if (p_exc_ave<=0.d0) then !should only occur if p_exc0<0 and dt small
                debugloop=debugloop + 1
                dtm = dtr
                m = m_0*exp(km0*0.5d0*(p_exc+p_exc0)*dtm)
             else
                debugloop=debugloop + 2
                m_upper = m_eq - max(0.d0,(1.d0/3.d0)*kp0*h0*(p_exc_ave)/(vnorm*alphainv))

                if ((km0*p_exc_ave*m_0).lt.1.d-99) then
                  dtm = dtp
               else
                  dtm = min(dtp,(m_eq-m_0)/(km0*p_exc_ave*m_0))
               endif

             endif
             m = m_0 + km0*m_0*p_exc_ave*dtm
             dts = dtm

          ! finally, consider region 1C - to the right and including the
          ! p-nullcline in UL.
          ! m will increase and p will decrease except at the nullcline
          ! break this region up into two portions, depending on whether to the
          ! right or left of m=m_eq.

          else !-kp0*p_exc0+c_d0)<0
             debugloop=debugloop + 300
             !p_exc is decreasing/above nullcline or static on nullcline
             ! m strictly increasing
             if (c_d0>0.d0) then !p_exc remains in quad1 for any dt
                debugloop=debugloop + 10
                ! only m can take solution beyond nullcline
                dtp = dtr

                if ((km0*p_exc0*m_0)>1.d-99) then
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                else
                  dtm = dtp
                endif

                dts = min(dtm,dtp)
                dts = max(dts,0.d0)
                m = m_0 + km0*m_0*p_exc0*dts
                p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)

             else !p_exc is to the right of m_eq
                 debugloop=debugloop + 20
                ! true solution should remain right of m_eq
                ! FE gives lowest slope of dp/dm
                dts = dtr
                dtp = dtr
                if ((km0*p_exc0*m_0)>1.d-99) then
                  dts = min(dtr,(m_c-m_0)/(km0*p_exc0*m_0))
                endif

                ! dtp based on not crossing p_eq_c

                if ((-kp0*p_eq0 - c_d0)/(kp0*p_exc0 - c_d0) >1.d-99) then
                  dtp = min(dtr,-(1.d0/kp0)*log( (c_d0)/(c_d0 - kp0*p_exc0) ))
                endif

                !if (dtp<dts) write(*,*) 'ERROR: QUAD 1: WRONG DIRECTION ACROSS m=m_eq'
                dts = min(dts,dtp)
                dts = max(dts,0.d0)
                m = m_0 + km0*m_0*p_exc0*dts
                p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
             endif
          endif
          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dts
          dtr = dtr-dts

       case(2) !UR quadrant, dense or equilibrium material, p>p_eq_c>p_eq
          !integration can only cross lower boundary (p=peq_c>p_eq)
          !m is unbounded because at start of integration, m>m_c>m_eq
          !p is strictly decreasing, m strictly increasing

          debugloop = 2000
          dtp = dtr

          ! dtp based on not crossing p_eq_c
          if ( (c_d0 - kp0 * (p_eq_c - p_eq0))/(c_d0-kp0*p_exc0)  >1.d-99) then
            dtp = min(dtr,-(1.d0/kp0)*log( (c_d0 - kp0 * (p_eq_c - p_eq0))/(c_d0-kp0*p_exc0) ))
          endif

          dtm = dtp
          ! ideally, we would bound dtm based on m not going above 1.0
          ! however, based on debugging there are cases with thin (0.001 m)
          ! and fast (>>1m/s) flow that result in dm changing by more than
          ! a few percent (e.g. up to 1.0). enforce that m cannot increase
          ! by more than 0.02
          if (km0*p_exc0*m_0>1.d-99) then
            dtm = min(dtp,(min(1.0d0-m_0,0.02d0)/(km0*p_exc0*m_0))) ! dig div0
          endif

          dts = min(dtp,dtm)
          
          p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)

          m = m_0 + km0*m_0*p_exc*dts

          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dts

          dtr = dtr-dts

       case(3) !LR quadrant, dilative material, p<=p_eq
          !integration can only cross left boundary (m=m_c)
          debugloop = 3000

          ! Region 3A above the m-nullcline and below the horizontal line
          ! through p_eq.
          if (p_exc0>1.d-3*rhoh*gz) then !p decreasing, m increasing
             debugloop=debugloop + 100
             !don't descend below p_exc = 0 until next substep
             !if don't let m increase by more than 0.02 or hit m=1.d0. if m=1.d0,
             ! follow boundary

             dtp = dtr
             if ((c_d0/(-kp0*p_exc0+c_d0))>1.d-99) then

               ! p_exc = (p_exc0 - c_d0/kp0) * exp( - kp0 * dt) + c_d0/kp0 > 0
               ! because here we are solving for the dt at which p_exc = 0,
               ! the dt bound has a different functional form that  in region 2,
               ! where we are solving for the dt at which p_exc = p_excess_critical

               dtp = min(dtr,-(1.d0/kp0)*log(c_d0/(c_d0-kp0*p_exc0)))

             endif
             dtm = dtp

             if ((km0*p_exc0*m_0)>1.d-99) then
               dtm = min(dtp,(min(1.0d0-m_0,0.02d0))/(km0*p_exc0*m_0))
             endif

             !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
             p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
             m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dtm
             dts = dtp


          ! Region 3B below the m-nullcline and above the p-nullcline
          ! p is decreasing and m is decreasing.
          elseif ((-kp0*p_exc0+c_d0)<0.d0) then !p_exc is decreasing above/right of nullcline
             debugloop=debugloop + 200
             if (c_d0/kp0<-p_eq0) then !assymptotic limit of p<0
                debugloop=debugloop + 10
                dtp = dtr
                if (((c_d0+kp0*p_eq0)/(-kp0*p_exc0+c_d0))>1.d-99) then
                  dtp = min(dtr,-(1.d0/kp0)*log((c_d0+kp0*p_eq0)/(-kp0*p_exc0+c_d0)))
                endif
             else
                debugloop=debugloop + 20
                dtp = dtr
             endif
             !p_exc decrease bounded by p=0 or nullcline if dtp>=dtr
             p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dtp) +c_d0)
             p_exc_ave = 0.5d0*(p_exc+p_exc0)
             if (p_exc_ave>0.d0) then
                debugloop=debugloop + 1
                !p_exc still above m nullcline (p=p_eq) even though below p_eq_c, m increasing
                !happens only if still near p_exc=0.
                dtm = min(dtr,(0.02d0)/(km0*p_exc_ave*m_0))
             else
                debugloop=debugloop + 2
                m_lower = m_eq + max(0.d0,(1.d0/3.d0)*kp0*h0*(-p_exc_ave)/(vnorm*alphainv))
                dtm = min(dtr,(m_lower-m_0)/(km0*p_exc_ave*m_0))
             endif
             !dts = min(dts,dtr)
             !m = m_0 + km0*m_0*0.5d0*(p_exc0+p_exc)*dts
             m = m_0 + km0*m_0*p_exc_ave*dtm
             dts = dtm


          ! Region 3C - p is increaseing, m is increasing.
          else !p_exc is increasing below/left of nullcline
             debugloop=debugloop + 300
             if (c_d0<0.d0) then !between nullcline and m_eq
                debugloop=debugloop + 10
                !no limit needed for pressure timestep (p_exc0 < p_exc-->c_d/kp <0 ), limit based on m>m_c
                dtp = dtr
                dtm = dtr
                if (-(km0*p_exc0*m_0)>1.d-99) then
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                endif

                dts = min(dtp,dtm)
                dts = max(dts,0.d0)
                p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
                m = m_0 + km0*m_0*p_exc0*dts
             else !left of m_eq, p_exc0 < p_exc-->c_d/kp > 0, bound by 0 with dtp
                debugloop=debugloop + 20
                dtm = dtr
                dtp = dtr
                if ((km0*p_exc0*m_0)>1.d-99) then
                  dtm = (m_c-m_0)/(km0*p_exc0*m_0)
                endif
                if ((-c_d0/(kp0*p_exc0-c_d0))>1.d-99) then
                  dtp = -(1.d0/kp0)*log(c_d0/(c_d0-kp0*p_exc0))
                endif
                dts = min(dtp,dtm)
                dts = max(dts,0.d0)
                !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
                p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
                m = m_0 + km0*m_0*p_exc0*dts
             endif
          endif
          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dts
          dtr = dtr-dts

       case(4) !LL quadrant
          ! m is strictly decreasing, p strictly increasing
          debugloop = 4000
          !don't exceed p_exc = 0 until next substep
          dtp = dtr



          ! bound dt based on p not crossing p_critical point
          if ( (c_d0 - kp0 * (p_eq_c - p_eq0))/(c_d0-kp0*p_exc0)  >1.d-99) then
            dtp = min(dtr,-(1.d0/kp0)*log( (c_d0 - kp0 * (p_eq_c - p_eq0))/(c_d0-kp0*p_exc0) ))
          endif
          
          dtm = log(0.98d0)/(kp0*p_exc0)

          dts = min(dtr, min(dtp,dtm))

          !integrate dp_exc/dt = -kp0 p_exc + c_d with coefficients at t=0.
          p_exc = (1.d0/kp0)*((kp0*p_exc0 - c_d0)*exp(-kp0*dts) +c_d0)
          m = m_0*exp(dts*km0*0.5d0*(p_exc0+p_exc))
          dts = dtp
          call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
          dtk = dts
          dtr = dtr-dts
       end select

       if (debug) then
          if ((m<m_c).and.(p>=p_eq_c)) then
             quad1=1
          elseif ((m>=m_c).and.(p>p_eq_c)) then
             quad1=2
          elseif  ((m>m_c).and.(p<=p_eq_c)) then
             quad1=3
          elseif ((m<=m_c).and.(p<p_eq_c)) then
             quad1=4
          endif
       endif

       if (debug) then
          call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq1,tanpsi,tau)
          if (dtr>0.d0) then
             if (dtk/dtr<1.d-6.and.quad0==quad1) then
                write(*,*) '---------------SRC WARNING: SMALL TIMESTEP---------->>>>>>'
                write(*,*) 'dtk,dtr:', dtk,dtr
                write(*,*) 'dtm,dtp:', dtm,dtp
                write(*,*) 'debugloop:', debugloop
                write(*,*) ' quad0,quad1', quad0,quad1
                write(*,*) 'm_0, m_c:', m_0, m_c
                write(*,*) 'm_0 - m_c:', m_0 - m_c
                write(*,*) 'm_0 - m_eq0:', m_0 - m_eq
                write(*,*) 'meq0 - m_c:', m_eq - m_c
                write(*,*) 'm1, m_eq1', m, m_eq1
                write(*,*) 'm1 - m_c:', m - m_c
                write(*,*) 'm1 - m_eq1:', m - m_eq1
                write(*,*) 'meq1 - m_c:', m_eq1 - m_c
                write(*,*) 'p0,p1:', p0,p
                write(*,*) 'p_exc0,p_exc:', p_exc0,p_exc
                write(*,*) 'h0,h_c,diff', h0,h_c,h0-h_c
                write(*,*)  'p,p_eq0,diff', p0,p_eq0,p-p_eq0
                write(*,*)  'p,p_eq_c,diff', p0,p_eq_c,p-p_eq_c
                write(*,*)  'p_eq0,p_eq_c,diff', p_eq0,p_eq_c,p_eq0-p_eq_c
                !stop
                write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                !stop
             endif
          endif

          if ((quad1-quad0)==-1) then
             write(*,*) '------------SRC WARNING: COUNTERCLOCKWISE ---------->>>>>>>>'
             write(*,*) ' quad0,quad1',quad0,quad1
             write(*,*) 'debugloop:',debugloop
             write(*,*) 'm_0,m,m_eq,m_c',m_0,m,m_eq,m_c
             write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
             !stop
          endif

          if ((dtm<0.d0).or.(dtp<0.d0)) then
             write(*,*) '------------SRC WARNING: negative dt ---------->>>>>>>>'
             write(*,*) 'dtk,dtr:', dtk,dtr
             write(*,*) 'dtm,dtp:', dtm,dtp
             write(*,*) ' quad0,quad1',quad0,quad1
             write(*,*) 'debugloop:',debugloop
             write(*,*) 'm_0,m,m_eq,m_c',m_0,m,m_eq,m_c
             write(*,*) 'p_exc0,-1.d-3*rhoh*gz', p_exc0,-1.d-3*rhoh*gz
             write(*,*) 'kp0,c_d0,-kp0*p_exc0+c_d0',kp0,c_d0,-kp0*p_exc0+c_d0
             write(*,*) 'h0,h_c,diff', h0,h_c,h0-h_c
             write(*,*)  'p,p_eq0,diff', p0,p_eq0,p-p_eq0
             write(*,*)  'p,p_eq_c,diff', p0,p_eq_c,p-p_eq_c
             write(*,*)  'p_eq0,p_eq_c,diff', p_eq0,p_eq_c,p_eq0-p_eq_c
             write(*,*)  'if small rounding error use max(dt,0)'
             write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
             stop
          endif

         if (abs(m_0-m)>0.05d0) then
                write(*,*) 'm_0-m', m_0-m
   !             write(*,*) '------------SRC WARNING: large change in m ---------->>>>>>>>'
   !             write(*,*) 'dtk,dtr:', dtk,dtr
   !             write(*,*) 'dtm,dtp:', dtm,dtp
   !             write(*,*) ' quad0,quad1',quad0,quad1
   !             write(*,*) 'debugloop:',debugloop
   !             write(*,*) 'm_0,m,m_eq,m_c',m_0,m,m_eq,m_c
   !             write(*,*) "vnorm, tanpsi", vnorm, tanpsi
   !             write(*,*) 'p_exc0,p_exc,-1.d-3*rhoh*gz', p_exc0,p_exc,-1.d-3*rhoh*gz
   !             write(*,*) 'kp0,c_d0,',kp0,c_d0,-kp0*p_exc0+c_d0
   !             write(*,*) '-kp0*p_exc0+c_d0', -kp0*p_exc0+c_d0
   !             write(*,*) 'km0*p_exc0*m_0', km0*p_exc0*m_0
   !             write(*,*) 'h0,h_c,diff', h0,h_c,h0-h_c
   !             write(*,*)  'p,p_eq0,diff', p0,p_eq0,p-p_eq0
   !             write(*,*)  'p,p_eq_c,diff', p0,p_eq_c,p-p_eq_c
   !             write(*,*)  'p_eq0,p_eq_c,diff', p_eq0,p_eq_c,p_eq0-p_eq_c
   !             write(*,*) 'p_exc0/(rho0*gz*h0)', p_exc0/(rho0*gz*h0)
   !             write(*,*) 'p_exc/(rho*gz*h0)', p_exc0/(rho0*gz*h0)
   !             write(*,*)  'if small rounding error use max(dt,0)'
   !             write(*,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
   !             !stop
             endif

       endif



       return
       end subroutine mp_update_FE_4quad

      !========================================================================
subroutine mp_update_relax_Dclaw4(dt,h,u,v,m,p,chi,rhoh,gz)
   !====================================================================
   !used if srcmethod = 0:1
   !subroutine mp_update_relax_Dclaw4: D-Claw 4.x method
   !integrate p, relaxing rho h constant (leave manifold)
   !for easier integration of p.
   !update m,h based on new value of D. rho h not exact/approx.

      use digclaw_module, only: rho_f,rho_s,mu,setvars,qfix,qfix_cmass
      use digclaw_module, only: src2method


      implicit none

      !i/o
      real(kind=8), intent(inout) :: h,u,v,m,p,chi
      real(kind=8), intent(in)  :: rhoh,dt
      real(kind=8), intent(in)  :: gz

      !local
      real(kind=8) :: p_eq,p_exc,vnorm,m_eq
      real(kind=8) :: rho,tanpsi,D,tau,kperm,alphainv
      real(kind=8) :: krate,zeta,hchi,hu,hv,hm

      call setvars(h,u,v,m,p,chi,gz,rho,kperm,alphainv,m_eq,tanpsi,tau)

      hu=h*u
      hv=h*v
      hm=h*m
      hchi=h*chi

      ! integrate pressure response to dilation/contraction
      vnorm = sqrt(hu**2 + hv**2)/h
      p = p - dt*3.d0*vnorm*tanpsi*alphainv/h
      call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)

      ! integrate pressure relaxation to hydrostatic
      zeta = 3.d0*alphainv/(h*2.d0)  + (rho-rho_f)*rho_f*gz/(4.d0*rho)
      krate=-zeta*2.d0*kperm/(h*max(mu,1.d-16))
      p_eq = rho_f*h*gz

      p = p_eq + (p-p_eq)*exp(krate*dt)

      !integrate changes in hm
      p_exc = p - rho_f*gz*h
      !define conservative variabls: (note u,v,m passed in, no hu,hv,hm yet).
      call qfix_cmass(h,m,p,rho,p_exc,hu,hv,hm,u,v,rhoh,gz)
      D = -2.0d0*(kperm/(mu*h))*p_exc
      krate = D*(rho-rho_f)/rhoh
      hm = hm*exp(-dt*D*rho_f/(rhoh))

      select case (src2method)
      !NOTE: at this point hm has changed.
      ! can redefine h by constant rho h and hm, or update h directly
      case(0)
         ! case 0: old method.
         ! integrate p explicitly (gives new D)
         ! integrate hu, hv, hm, and h with same D
         ! rationale is if m is poor, not to base h,hu,hv on it
         ! even if rho h is varied
         h = h + h*krate*dt
         hchi = h*chi
         hu = hu*exp(dt*krate)
         hv = hv*exp(dt*krate)
         call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)
      case(1)
         ! case 1: half-way new method.
         ! integrate hm and p as above.
         ! redefine h by hm and rhoh=constant, set hu,hv by new h
         ! rationale is to maintain rhoh = constant
         h = (rhoh - hm*(rho_s-rho_f))/rho_f
         hu = h*u
         hv = h*v
         hchi = h*chi
         call qfix(h,hu,hv,hm,p,hchi,u,v,m,chi,rho,gz)
      end select

      ! qfix will set prim with conserved

      end subroutine mp_update_relax_Dclaw4

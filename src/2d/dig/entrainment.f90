! ============================================================================
!  Entrain material at a single computational mesh cell
! ============================================================================
!
! Entrainment routines update the solution state based on entrainment and are
! called from within src2().
! If entrainment_method == 0, then ent_dclaw4() is called.
!
! Variables
! ---------
! dt, timestep duration
! h, depth of flow
! u, x-directed velocity
! v, y-directed velocity
! m, solid volume fraction
! p, basal pressure
! rho, density
! hchi, h*species 1 fraction
! gz, z-oriented gravitational acceleration
! tau, shear stress
! b_x, x-directed gradient in basal topography
! b_y, y-directed gradient in basal topography
! b_eroded, total amount of material that has been eroded from the erodible layer
! b_remaining, amount of material remaining in the erodible layer
!--------------------------------------------------------------------------------

        subroutine ent_dclaw4(dt,h,u,v,m,p,rho,hchi,gz,tau,b_x,b_y,b_eroded,b_remaining)

        use digclaw_module, only: rho_f,rho_s,mu,chie,entrainment_rate
        use digclaw_module, only: me,setvars,qfix,qfix_cmass
        use geoclaw_module, only: manning_coefficient,friction_forcing

        implicit none

        !i/o
        real(kind=8), intent(inout) :: h,m,p,u,v,rho,b_eroded,hchi
        real(kind=8), intent(in)  :: b_x,b_y,dt,b_remaining,tau
        real(kind=8), intent(in)  :: gz

        !local
        real(kind=8) :: vlow,vnorm,rhoe,coeff,gamma,visc,beta,t1bot,t2top
        real(kind=8) :: dh,hrho,prat,hm

        if (friction_forcing) then
            coeff = manning_coefficient(1)
        else
            coeff = 0.d0
        endif

        vnorm = sqrt(u**2 + v**2)
        ! minimum velocity for entrainment to occur. ! DIG: should this be a user
        ! specified variable.
        vlow = 0.1d0

        if (vnorm.gt.vlow) then

            if (b_remaining.gt.0.d0) then


                ! calculate entrained material density
                rhoe = me*rho_s + (1.d0-me)*rho_f

                ! calculate top and bottom shear stress. note rho*psi_2

                beta = max(1.d0-m, 0.1d0)

                visc = beta*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                gamma= rho*beta*(vnorm**2)*(beta*gz*coeff**2)/(tanh(h+1.d-2)**(1.0d0/3.0d0))
                t1bot = visc + gamma + tau

                t2top = min(t1bot,(1.d0-beta*entrainment_rate)*tau)

                ! calculate pressure ratio
                prat = p/(rho*h)

                ! calculate dh

                dh = entrainment_rate*dt*(t1bot-t2top)/(rhoe*beta*vnorm)
                dh = min(dh,b_remaining)

                ! increment h based on dh
                h = h + dh
                hm = h*m + dh*me
                m = hm/h

                ! Update values for hchi, hrho, and rho based on entrainment
                hchi = hchi + dh*chie
                hrho = hm*rho_s + (h-hm)*rho_f
                rho = hrho/h

                ! update b_eroded
                b_eroded = b_eroded + dh

                ! maintain pressure as the pressure ratio prior to entrainment
                p = prat*hrho

            endif
        endif

        end subroutine ent_dclaw4

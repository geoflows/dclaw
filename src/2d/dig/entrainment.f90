        
        subroutine ent_dclaw4(dt,h,u,v,m,p,rho,hchi,gz,tau,b_x,b_y,b_eroded,b_remaining)
        
        use digclaw_module, only: rho_f,rho_s,sigma_0,mu,chie
        use digclaw_module, only: me,alpha,setvars,qfix,qfix_cmass,phi,m_crit,delta
        use geoclaw_module, only: grav,dry_tolerance
        use geoclaw_module, only: manning_coefficient,friction_forcing
 
        implicit none
 
        !i/o
        real(kind=8), intent(inout) :: h,m,p,u,v,rho,b_eroded,hchi
        real(kind=8), intent(in)  :: b_x,b_y,dt,b_remaining
        real(kind=8), intent(in)  :: gz

        !local
        real(kind=8) :: vlow,vnorm,dbdv,rhoe,coeff,gamma,visc,beta,t1bot,t2top

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
               
            ! determine dbdv, the slope in the direction of flow by taking the dot
            ! product of the slope vector and the unit vector in the direction of flow.
            dbdv = (u*b_x+v*b_y)/vnorm

            ! calculate entrained material density
            rhoe = me*rho_s + (1.d0-me)*rho_f

            ! calculate top and bottom shear stress. note rho*psi_2
            beta = 1.d0-m 
            visc = beta*vnorm*2.d0*mu*(1.d0-m)/h
            gamma= rho*beta*(vnorm)*(beta*gz*coeff**2)/(h**(4.d0/3.d0))
            t1bot = visc + gamma + tau

            !DIG: this needs to be sorted out...don't remember what's going on here
            !DIG: fix
            t2top = min(t1bot,(1.d0-beta*entrainment_rate)*(tau))

            ! calculate pressure ratio
            prat = p/(rho*h)

            ! calculate dh
            dh = entrainment_rate*dt*(t1bot-t2top)/(rho2*beta*vnorm)
            dh = min(dh,b_remaining)

            ! increment h based on dh
            h = h + dh
            hm = hm + dh*me
            hchi = hchi + dh*chie
            hrho = hm*rho_s + (h-hm)*rho_f
            rho = hrho/h
            b_eroded = b_eroded + dh
            p = prat*hrho


            endif

        end subroutine ent_dclaw4
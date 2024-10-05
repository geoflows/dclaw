        
        subroutine ent_dclaw4(dt,h,u,v,m,p,rho,chi,gz,b_x,b_y,b_eroded,b_remaining)
        
        use digclaw_module, only: rho_f,rho_s,sigma_0,mu,chie
        use digclaw_module, only: me,alpha,setvars,qfix,qfix_cmass,phi,m_crit,delta
        use geoclaw_module, only: grav,drytolerance
 
        implicit none
 
        !i/o
        real(kind=8), intent(inout) :: h,m,p,u,v,rho,b_eroded,b_remaining,chi
        real(kind=8), intent(in)  :: b_x,b_y,dt
        real(kind=8), intent(in)  :: gz

        !local
        real(kind=8) :: vlow,vnorm,dbdv,rhoe

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

            ! calculate top and bottom shear stress.
            t1bot = (1.d0-m)*vnorm*2.d0*mu*(1.d0-m)/(tanh(h+1.d-2))
                     
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

        end subroutine ent_dclaw4
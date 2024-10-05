        
        subroutine ent_dclaw4(dt,h,u,v,m,p,chi,b_x,b_y,b_eroded,b_remaining)
        
        use digclaw_module, only: rho_f,rho_s,sigma_0,mu,alpha,setvars,qfix,qfix_cmass,phi_bed,m_crit,delta
        use geoclaw_module, only: grav,drytolerance
 
        implicit none
 
        !i/o
        real(kind=8), intent(inout) :: dt,h,m,p,u,v,b_eroded,b_remaining
        real(kind=8), intent(in)  :: b_x,b_y
        real(kind=8), intent(in)  :: gz



        end subroutine ent_dclaw4
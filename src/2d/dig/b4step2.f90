! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,actualstep)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET

    use digclaw_module, only: i_theta,bed_normal

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    logical, intent (in) :: actualstep

    ! Local storage
    integer :: index,i,j,k,dummy
    real(kind=8) :: h,u,v,sv,theta

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set other variables appropriately for these cells
    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            theta = 0.d0
            if (bed_normal.eq.1) theta=aux(i_theta,i,j)
            call admissibleq(q(i,j,1),q(i,j,2),q(i,j,3),q(i,j,4),q(i,j,5),u,v,sv,theta)
            q(i,j,6) = min(q(i,j,6),q(i,j,1))
            q(i,j,6) = max(q(i,j,6),0.0)
            q(i,j,7) = max(q(i,j,7),0.0)
        enddo
    enddo


    if (aux_finalized < 2 .and. actualstep) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    ! find factor of safety ratios
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine b4step2

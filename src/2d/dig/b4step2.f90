!! D-Claw specific core file
!! This file is a modified version of
!! clawpack/geoclaw/src/2d/shallow/b4step2.f90 
!!

! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux,actualstep)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(i_h,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(i_h,i,j) is a depth.

!

    use topo_module, only: aux_finalized
    use geoclaw_module, only: grav

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET

    use digclaw_module, only: i_theta,bed_normal,qfix,calc_taudir
    use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_bdif

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
    real(kind=8) :: h,u,v,m,chi,rho,gz

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set other variables appropriately for these cells

    gz=grav
    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (bed_normal.eq.1) gz=grav*dcos(aux(i_theta,i,j))

            call qfix(q(i_h,i,j),q(i_hu,i,j),q(i_hv,i,j),q(i_hm,i,j),q(i_pb,i,j),q(i_hchi,i,j),u,v,m,chi,rho,gz)

            q(i_bdif,i,j) = max(q(i_bdif,i,j),0.0d0) ! DIG: if we store deposition in i_bdif, remove.
        enddo
    enddo

    if (aux_finalized < 2 .and. actualstep) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    ! find factor of safety ratios and friction orientation
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine b4step2

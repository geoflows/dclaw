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
    use auxt_module, only: auxtfill_finalized
    use geoclaw_module, only: grav
    use geoclaw_module, only: speed_limit

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET
    use amr_module, only: outunit

    use digclaw_module, only: i_theta,bed_normal,qfix,calc_taudir
    use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi,i_hs,i_ent,i_hf

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
    real(kind=8) :: sratio, xs,ys, s

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set other variables appropriately for these cells

    gz=grav
    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (bed_normal.eq.1) gz=grav*dcos(aux(i_theta,i,j))

            call qfix(q(i_h,i,j),q(i_hu,i,j),q(i_hv,i,j),q(i_hm,i,j),q(i_pb,i,j),q(i_hchi,i,j),q(i_hf,i,j),u,v,m,chi,rho,gz)

            ! q(i_hs,i,j) = max(q(i_hs,i,j),0.0d0) ! DIG: if we store deposition in i_hs, remove.
        enddo
    enddo


    ! Check for fluid speed sqrt(u**2 + v**2) > speed_limit
    ! and reset by scaling (u,v) down to this value (preserving direction)
    ! Note: similar check is done in getmaxspeed
    ! This helps avoid too many dt reductions when flow off very steep topo
    ! with delta B larger than fluid depth gives big speeds in Riemann solution
    ! (shallow water equations aren't valid for flow off a cliff)

    do j=1-mbc,my+mbc
        do i=1-mbc,mx+mbc
            if (q(1,i,j) > 0.d0) then
                s = sqrt((q(2,i,j)**2 + q(3,i,j)**2)) / q(1,i,j)
                if (s > speed_limit) then
                    sratio = speed_limit / s
                    q(2,i,j) = q(2,i,j) * sratio
                    q(3,i,j) = q(3,i,j) * sratio

                    if (.false.) then
                        ! write out info useful for investigating topo:
                        xs = xlower + (i-0.5d0)*dx
                        ys = ylower + (j-0.5d0)*dy
                        write(6,604) t, i,j, mx,my, dx
                        write(outunit,604) t, i,j, mx,my, dx
                        write(6,603) s,q(1,i,j),aux(1,i,j),xs,ys
                        write(outunit,603) s,q(1,i,j),aux(1,i,j),xs,ys

 604                    format('b4step2 at t =',f10.2, '  i,j,mx,my:',4i4, &
                               '  dx = ',f10.7)
 603                    format('     reset s =',f9.2,'  h=',e11.3, ' B=',f8.2,&
                               '  x,y = ', f11.6,',',f10.6)
                    endif
               endif
            endif
        enddo
    enddo


    if ((aux_finalized < 2 .and. actualstep) .or. (auxtfill_finalized.eqv..false.)) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)

        ! update_auxt() called from within setaux() to ensure it is allways
        ! called when setaux() is called.
    endif

    ! find factor of safety ratios and friction orientation
      call calc_taudir(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

end subroutine b4step2

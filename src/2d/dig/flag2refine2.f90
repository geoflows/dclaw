! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for GeoClaw for tsunami applications and related problems
!
!
!    q   = grid values including ghost cells (bndry vals at specified
!          time have already been set, so can use ghost cell values too)
!
!  aux   = aux array on this grid patch
!
! amrflags  = array to be flagged with either the value
!             DONTFLAG (no refinement needed)  or
!             DOFLAG   (refinement desired)
!
! This routine only checks cells in which amrflags(i,j) = UNSET,
! If the value is already DONTFLAG or DOFLAG then this we determined by
! flagregions or other criteria and is not modified by this routine.
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine flag2refine2(mx,my,mbc,mbuff,meqn,maux,xlower,ylower,dx,dy,t,level, &
                       tolsp,q,aux,amrflags)

    use amr_module, only: mxnest, t0, DOFLAG, UNSET
    use amr_module, only: lfine, lstart, node, rnode  ! DIG
    use amr_module, only: cornxlo, cornylo, levelptr, mxnest, ndihi, ndilo
    use amr_module, only: ndjhi, ndjlo, store1, store2, storeaux
    use amr_module, only: alloc, hxposs, hyposs

    use geoclaw_module, only:dry_tolerance, sea_level
    use geoclaw_module, only: spherical_distance, coordinate_system


    use storm_module, only: storm_specification_type, wind_refine, R_refine
    use storm_module, only: storm_location, wind_forcing, wind_index, wind_refine

    use regions_module, only: num_regions, regions
    use refinement_module, only: wave_tolerance, speed_tolerance
    use refinement_module, only: flowgradevalue, iflowgradevariable, &
        iflowgradetype, iflowgrademinlevel, mflowgrades, keep_fine

    use adjoint_module, only: totnum_adjoints,innerprod_index, &
                              adjoint_flagging,select_snapshots
    use adjointsup_module, only: calculate_innerproduct

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: mx,my,mbc,meqn,maux,level,mbuff
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,tolsp

    real(kind=8), intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Flagging
    real(kind=8), intent(in out) :: amrflags(1-mbuff:mx+mbuff,1-mbuff:my+mbuff)

    ! Generic locals
    integer :: i,j,m,r
    real(kind=8) :: x_c,y_c,x_low,y_low,x_hi,y_hi
    real(kind=8) :: speed, eta, ds

    !!DIG:  Need to declare more flowgrade/keep_fine variables?
    real(kind=8) :: flowgradenorm, flowgradegrad, depth, momentum, surface
    real(kind=8) :: xlow,xhi,ylow,yhi,xxlow,xxhi,yylow,yyhi
    real(kind=8) :: flowgrademeasure, h,hu,hv
    integer :: nx,ny,loc,locaux,mitot,mjtot,mptr,iflow,ii,jj
    

    ! Loop over interior points on this grid
    ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell center at (x_c,y_c)
    y_loop: do j=1,my
        y_low = ylower + (j - 1) * dy
        y_c = ylower + (j - 0.5d0) * dy
        y_hi = ylower + j * dy

        x_loop: do i = 1,mx
            x_low = xlower + (i - 1) * dx
            x_c = xlower + (i - 0.5d0) * dx
            x_hi = xlower + i * dx
            
            if (amrflags(i,j) .ne. UNSET) then
                ! do not consider this cell further
                cycle x_loop
            endif

            ! DIG:  Adapted from D-Claw...

            ! determine if flowgrades are used
            if (mflowgrades > 0) then
                depth= q(1,i,j)
                momentum = sqrt(q(2,i,j)**2 + q(3,i,j)**2)
                surface = q(1,i,j) + aux(1,i,j)
                do iflow=1,mflowgrades
                      if (iflowgradevariable(iflow).eq.1) then
                        flowgradenorm=depth
                        flowgradegrad=depth

                      elseif (iflowgradevariable(iflow).eq.2) then
                        flowgradenorm=momentum
                        flowgradegrad=momentum

                      elseif (iflowgradevariable(iflow).eq.3) then

                        if (depth.gt.dry_tolerance) then
                          flowgradenorm=dabs(surface)
                          flowgradegrad=dabs(surface)
                        else
                          flowgradenorm=0.d0
                          flowgradegrad=0.d0
                        endif
                      endif

                      if (iflowgradetype(iflow).eq.1) then
                        flowgrademeasure=flowgradenorm
                      else
                        write(*,*) 'only flowgradetype = 1 supported'
                        stop
                        flowgrademeasure=flowgradegrad
                      endif

                      if (flowgrademeasure.gt.flowgradevalue(iflow) &
                        .and.level.lt.iflowgrademinlevel(iflow)) then
                        amrflags(i,j)=DOFLAG
                        cycle x_loop
                      endif
                enddo
            endif ! mflowgrades > 0
            
            ! keep fine added by KRB 2022/12/28
            if (keep_fine.and.mflowgrades.gt.0) then
              ! if level is lower than lfine a grid exists here on lfine,
              ! enforce refinement here.
              ! here, ignore ghost cells on the fine grid (e.g., no adjustment
              ! of extent)
              if (level .lt. lfine) then

                 ! loop through lfine grids (code modified by valout)
                 mptr = lstart(lfine)
                 do while (mptr > 0)

                    ! calculate extent of fine grid
                    nx      = node(ndihi,mptr) - node(ndilo,mptr) + 1
                    ny      = node(ndjhi,mptr) - node(ndjlo,mptr) + 1
                    loc     = node(store1, mptr)
                    locaux  = node(storeaux,mptr)
                    mitot   = nx + 2*mbc
                    mjtot   = ny + 2*mbc

                    xlow = rnode(cornxlo,mptr)
                    ylow = rnode(cornylo,mptr)
                    xhi = xlow + nx*hxposs(lfine)
                    yhi = ylow + ny*hyposs(lfine)

                    ! if there is overlap between the fine grid and this
                    ! location on the coarse grid, flag.
                    ! (i,j) grid cell is [x1,x2] x [y1,y2].
                    ! fine grid is [xlow,xhi] x [ylow,yhi]
                    if (x_hi.gt.xlow.and.x_low.lt.xhi.and. &
                                y_hi.gt.ylow.and.y_low.lt.yhi) then

                       ! this loop includes ghosts, update extent and
                       mitot   = nx + 2*mbc
                       mjtot   = ny + 2*mbc
                       xlow = rnode(cornxlo,mptr)-mbc*hxposs(lfine)
                       ylow = rnode(cornylo,mptr)-mbc*hyposs(lfine)
                       xhi = xlow + mitot*hxposs(lfine)
                       yhi = ylow + mjtot*hyposs(lfine)

                       ! get pointers
                       loc     = node(store1, mptr)
                       locaux  = node(storeaux,mptr)

                       ! loop through fine grid cells.
                       do jj = mbc+1, mjtot-mbc
                          do ii = mbc+1, mitot-mbc
                            ! check overlap between fine and coarse cells
                            ! ignore ghost cells
                            xxlow = xlow + hxposs(lfine)*ii
                            xxhi = xxlow + hxposs(lfine)
                            yylow = ylow + hyposs(lfine)*jj
                            yyhi = yylow + hyposs(lfine)

                            ! if fine grid cell is inside of coarse grid cell
                            ! calculate flowgrade values
                            if (x_hi.gt.xxlow.and.x_low.lt.xxhi.and. &
                                   y_hi.gt.yylow.and.y_low.lt.yyhi) then
                              h = alloc(iadd(1,ii,jj))
                              hu = alloc(iadd(2,ii,jj))
                              hv = alloc(iadd(3,ii, jj))
                              momentum = sqrt((hu**2)+(hv**2))
                              surface = h + alloc(iaddaux(1,ii,jj))

                    ! check flowgrade values on fine grid.
                    do iflow=1,mflowgrades
                      if (iflowgradevariable(iflow).eq.1) then
                        flowgradenorm=depth
                        flowgradegrad=depth
                      elseif (iflowgradevariable(iflow).eq.2) then
                        flowgradenorm=momentum
                        flowgradegrad=momentum
                      elseif (iflowgradevariable(iflow).eq.3) then
                        if (depth.gt.dry_tolerance) then
                          flowgradenorm=dabs(surface)
                          flowgradegrad=dabs(surface)
                        else
                          flowgradenorm=0.d0
                          flowgradegrad=0.d0
                        endif
                      endif
                      if (iflowgradetype(iflow).eq.1) then
                        flowgrademeasure=flowgradenorm
                      else
                        write(*,*) 'flowgradetype not supported'
                        stop
                        flowgrademeasure=flowgradegrad
                      endif
                      if (flowgrademeasure.gt.flowgradevalue(iflow) &
                          .and.level.lt.iflowgrademinlevel(iflow)) then
                        amrflags(i,j)=DOFLAG
                        cycle x_loop
                      endif
                    enddo ! end flowgrade loop

                            endif ! end if coarse and fine cell overlap
                          enddo ! endif fine cell loop ii
                       enddo ! endif fine cell loop jj

                    endif ! endif coarse cell overlaps fine grid
                    mptr = node(levelptr, mptr)
                enddo ! end while (mptr > 0)

              endif ! (end level .lt. lfine)
            endif ! (end if keep_fine)
            ! end keep fine


            ! ********* criteria applied only to wet cells:

            if (q(1,i,j) > dry_tolerance) then

! DIG: D-claw not yet compatible with adjoint
!                if(adjoint_flagging) then
!                    if(aux(innerprod_index,i,j) > tolsp) then
!                        amrflags(i,j) = DOFLAG
!                        cycle x_loop
!                    endif
!                else

                ! Check wave criteria
                eta = q(1,i,j) + aux(1,i,j)
                if (abs(eta - sea_level) > wave_tolerance) then
                    amrflags(i,j) = DOFLAG
                    cycle x_loop
                endif

                ! DIG: Speed criteria comparable to flowgrades and keep fine?
                ! Check speed criteria (distinct values for each level), note that it might be useful to
                ! also have a per layer criteria since this is not
                ! gradient based
                speed = sqrt(q(2,i,j)**2 + q(3,i,j)**2) / q(1,i,j)
                do m=1,min(size(speed_tolerance),mxnest)
                    if (speed > speed_tolerance(m) .and. level <= m) then
                        amrflags(i,j) = DOFLAG
                        cycle x_loop
                    endif
                enddo
!                endif
            endif


        enddo x_loop
    enddo y_loop
    
contains

    ! Index into q array
    pure integer function iadd(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iadd = loc + m - 1 + meqn * ((j - 1) * (mitot + 2 * mbc) + i - 1)
    end function iadd

    ! Index into aux array
    pure integer function iaddaux(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iaddaux = locaux + m - 1 + maux * (i - 1) + maux * (mitot + 2 * mbc) * (j - 1)
    end function iaddaux
    
end subroutine flag2refine2

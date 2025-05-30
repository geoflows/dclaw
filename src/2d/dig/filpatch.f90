!! D-Claw specific core file
!! This file is a modified version of
!! clawpack/geoclaw/src/2d/shallow/fipatch.f90 
!!

! :::::::::::::::::::::::::::: FILPATCH :::::::::::::::::::::::::;
!
!  fill the portion of valbig from rows  nrowst
!                             and  cols  ncolst
!  the patch can also be described by the corners (xlp,ybp) by (xrp,ytp).
!  vals are needed at time t, and level level,
!
!  first fill with  values obtainable from the level level
!  grids. if any left unfilled, then enlarge remaining rectangle of
!  unfilled values by 1 (for later linear interp), and recusively
!  obtain the remaining values from  coarser levels.
!
! :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::;
recursive subroutine filrecur(level,nvar,valbig,aux,naux,t,mx,my, &
                              nrowst,ncolst,ilo,ihi,jlo,jhi,patchOnly,msrc,  &
                              do_aux_copy)

    use amr_module, only: nghost, xlower, ylower, xupper, yupper, outunit
    use amr_module, only: xperdom, yperdom, spheredom, hxposs, hyposs
    use amr_module, only: intratx, intraty, iregsz, jregsz
    use amr_module, only: NEEDS_TO_BE_SET
    use digclaw_module, only: m0, rho_f

    use geoclaw_module, only: sea_level, dry_tolerance, grav
    use topo_module, only: topo_finalized
    use qinit_module, only: variable_eta_init
    use qinit_module, only: force_dry,use_force_dry,mx_fdry, my_fdry
    use qinit_module, only: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry
    use qinit_module, only: dx_fdry, dy_fdry
    use qinit_module, only: tend_force_dry

    implicit none

    ! Input
    integer, intent(in) :: level, nvar, naux, mx, my, nrowst, ncolst
    integer, intent(in) :: ilo,ihi,jlo,jhi,msrc
    real(kind=8), intent(in) :: t
    logical  :: patchOnly, do_aux_copy, yes_do_aux_copy

    ! Output
    real(kind=8), intent(in out) :: valbig(nvar,mx,my)
    real(kind=8), intent(in out) :: aux(naux,mx,my)

    ! Local storage
    integer  :: iplo, iphi, jplo, jphi

    ! Flagging of set cells
    logical :: set
    integer :: i, i_coarse, j_coarse, i_fine, j_fine, n, k
    integer :: mx_coarse, my_coarse, mx_patch, my_patch
    integer :: unset_indices(4), coarse_indices(4)
    integer :: refinement_ratio_x, refinement_ratio_y
    real(kind=8) :: dx_fine, dy_fine, dx_coarse, dy_coarse
    real(kind=8) :: xlow_coarse,ylow_coarse, xlow_fine, ylow_fine, xhi_fine,yhi_fine
    real(kind=8) :: h, b, eta_fine, eta1, eta2, up_slope, down_slope
    real(kind=8) :: hv_fine, v_fine, v_new, divide_mass
    real(kind=8) :: h_fine_average, h_fine, h_count, h_coarse
    real(kind=8)::  xcent_fine,xcent_coarse,ycent_fine,ycent_coarse,ratio_x,ratio_y
    integer(kind=1) :: flaguse(ihi-ilo+1, jhi-jlo+1)

    real(kind=8) :: eta1old, eta2old
    real(kind=8) :: veta_init_c
    integer :: ii,jj

    ! Scratch arrays for interpolation
    logical :: fine_flag(nvar, ihi-ilo+2,jhi-jlo + 2)

    logical :: reloop


    ! these are dimensioned at fine size since coarse grid size cant be larger (incl. the +3 that is )
    real(kind=8) ::  fine_mass(ihi-ilo + 3, jhi-jlo + 3)
    real(kind=8) :: eta_coarse(ihi-ilo + 3, jhi-jlo + 3)
    real(kind=8) ::    vel_max(ihi-ilo + 3, jhi-jlo + 3)
    real(kind=8) ::    vel_min(ihi-ilo + 3, jhi-jlo + 3)
    real(kind=8) ::   slope(2, ihi-ilo + 3, jhi-jlo + 3)
    integer ::   fine_cell_count(ihi-ilo+3, jhi-jlo + 3)

    integer :: nghost_patch, lencrse
    logical :: use_force_dry_this_level
    real(kind=8) :: ddxy

    ! Stack storage
    !  use stack-based scratch arrays instead of alloc, since dont really
    !  need to save beyond these routines, and to allow dynami_coarse memory resizing
    !  use 1d scratch arrays that are potentially the same size as
    !  current grid, since may not coarsen.
    !  need to make it 1d instead of 2 and do own indexing, since
    !  when pass it in to subroutines they treat it as having di_fineerent
    !  dimensions than the max size need to allocate here
    ! the +2 is to expand on coarse grid to enclose fine
    real(kind=8) :: valcrse((ihi-ilo+3) * (jhi-jlo+3) * nvar)   ! NB this is a 1D array
    real(kind=8) :: auxcrse((ihi-ilo+3) * (jhi-jlo+3) * naux)
    real(kind=8) :: vetac((ihi-ilo+3) * (jhi-jlo+3))

    yes_do_aux_copy = .true.

    mx_patch = ihi-ilo + 1 ! nrowp
    my_patch = jhi-jlo + 1

    dx_fine     = hxposs(level)
    dy_fine     = hyposs(level)

    use_force_dry_this_level = use_force_dry
    if (use_force_dry) then
        ! check if force_dry resolution the same as this level:
        ddxy = max(abs(dx_fine-dx_fdry), abs(dy_fine-dy_fdry))
        use_force_dry_this_level = (ddxy < 0.01d0*min(dx_fdry,dy_fdry))
        endif

    !if (use_force_dry_this_level .and. (t <= tend_force_dry)) then
    !    write(6,*) '+++ using force_dry in filval, t = ',t
    !    endif

    ! Coordinates of edges of patch (xlp,xrp,ybp,ytp)
    xlow_fine = xlower + ilo * dx_fine
    ylow_fine = ylower + jlo * dy_fine
    xhi_fine = xlower + (ihi + 1) * dx_fine
    yhi_fine = ylower + (jhi + 1) * dy_fine

    ! Fill in the patch as much as possible using values at this level
    ! note that if only a patch, msrc = -1, otherwise a real grid and intfil
    ! uses its boundary list
    ! msrc either -1 (for a patch) or the real grid number
    call intfil(valbig,aux,mx,my,t,flaguse,nrowst,ncolst, ilo,  &
                ihi,jlo,jhi,level,nvar,naux,msrc,do_aux_copy)

    ! Trimbd returns set = true if all of the entries are filled (=1.).
    ! set = false, otherwise. If set = true, then no other levels are
    ! are required to interpolate, and we return.
    !
    ! Note that the used array is filled entirely in intfil, i.e. the
    ! marking done there also takes  into account the points filled by
    ! the boundary conditions. bc2amr will be called later, after all 4
    ! boundary pieces filled.
    call trimbd(flaguse,mx_patch,my_patch,set,unset_indices)
    ! il,ir,jb,jt = unset_indices(4)

    ! If set is .true. then all cells have been set and we can skip to setting
    ! the remaining boundary cells.  If it is .false. we need to interpolate
    ! some values from coarser levels, possibly calling this routine
    ! recursively.
    if (.not.set) then

        ! Error check
        if (level == 1) then
            write(outunit,*)" error in filrecur - level 1 not set"
            write(outunit,'("start at row: ",i4," col ",i4)') nrowst,ncolst
            print *," error in filrecur - level 1 not set"
            print *," should not need more recursion "
            print *," to set patch boundaries"
            print '("start at row: ",i4," col ",i4)', nrowst,ncolst
            stop
        endif

        ! We begin by initializing the level level arrays, so that we can use
        ! purely recursive formulation for interpolating.
        dx_coarse  = hxposs(level - 1)
        dy_coarse  = hyposs(level - 1)

        ! Adjust unset_indices to account for the patch
        ! isl, isr, jsb, jst
        unset_indices(1) = unset_indices(1) + ilo - 1
        unset_indices(2) = unset_indices(2) + ilo - 1
        unset_indices(3) = unset_indices(3) + jlo - 1
        unset_indices(4) = unset_indices(4) + jlo - 1

        ! Coarsened geometry
        refinement_ratio_x = intratx(level - 1)
        refinement_ratio_y = intraty(level - 1)

        ! New patch rectangle (after we have partially filled it in) but in the
        ! coarse patches [iplo,iphi,jplo,jphi]

        iplo = (unset_indices(1) - refinement_ratio_x + nghost * refinement_ratio_x) &
                                                / refinement_ratio_x - nghost
        iphi = (unset_indices(2) + refinement_ratio_x) / refinement_ratio_x
        jplo = (unset_indices(3) - refinement_ratio_y + nghost * refinement_ratio_y) &
                                                / refinement_ratio_y - nghost
        jphi = (unset_indices(4) + refinement_ratio_y) / refinement_ratio_y


        xlow_coarse = xlower + iplo * dx_coarse
        ylow_coarse = ylower + jplo * dy_coarse


        ! Coarse grid number of spatial points (nrowc,ncolc)
        mx_coarse   =  iphi - iplo  + 1
        my_coarse   =  jphi - jplo  + 1

        ! Check to make sure we created big enough scratch arrays
        if (mx_coarse > ihi - ilo + 3 .or. &
            my_coarse > jhi - jlo + 3) then

            print *," did not make big enough work space in filrecur "
            print *," need coarse space with mx_coarse,my_coarse ",mx_coarse,my_coarse
            print *," made space for ilo,ihi,jlo,jhi ",ilo,ihi,jlo,jhi
            stop
        endif

        ! Set the aux array values for the coarse grid, this could be done
        ! instead in intfil using possibly already available bathy data from the
        ! grids

        if (naux > 0) then
            nghost_patch = 0
            lencrse = (ihi-ilo+3)*(jhi-jlo+3)*naux ! set 1 component, not all naux
            ! new system checks initialization before setting aux vals
            ! set here and check later if not all set by copying in intfil
            ! after some more testing shoudl move thiese auxcrse initialization
            ! into the branches that need it (topo_finalized dand periodic
            ! but out of the most used filrecur pathway
            do k = 1, lencrse, naux
              auxcrse(k) = NEEDS_TO_BE_SET
            end do

            ! update topography if needed
            !if ((num_dtopo>0).and.(topo_finalized.eqv..false.)) then
            !   if ((minval(topotime)<maxval(tfdtopo)).and.(t>=minval(t0dtopo))) then
            if (.not. topo_finalized) then
                call topo_update(t)
                call setaux(nghost_patch, mx_coarse, my_coarse,       &
                            xlow_coarse, ylow_coarse,                 &
                            dx_coarse,dy_coarse,naux,auxcrse)
                yes_do_aux_copy = .false. !grids may not have latest topo
            endif
        endif

        ! Fill in the edges of the coarse grid
        if ((xperdom .or. (yperdom .or. spheredom)) .and. sticksout(iplo,iphi,jplo,jphi)) then
            call setaux(nghost_patch, mx_coarse, my_coarse,       &
                        xlow_coarse, ylow_coarse,                 &
                        dx_coarse,dy_coarse,naux,auxcrse)
            call prefilrecur(level-1,nvar,valcrse,auxcrse,naux,t, &
                             mx_coarse,my_coarse,1,1,iplo,iphi,jplo,jphi,iplo,iphi,jplo,jphi,.true.)
        else ! when going to coarser patch, no source grid (for now at least) hence -1
            call filrecur(level-1,nvar,valcrse,auxcrse,naux,t,  &
                          mx_coarse,my_coarse,1,1,iplo,iphi,jplo,jphi,.true.,-1,    &
                          yes_do_aux_copy)
            !call setaux(nghost_patch, mx_coarse, my_coarse,       &
            !            xlow_coarse, ylow_coarse,                 &
            !            dx_coarse,dy_coarse,naux,auxcrse)
        endif

        ! loop through coarse cells determining intepolation slopes
        ! these will be saved for fine grid loop
        ! prevents finding the same slope possibly lratiox*lratioy times
        ! all fine gid depths will be found before any momentum
        reloop = .false.
        fine_cell_count = 0
        fine_flag = .false.
        fine_mass = 0.d0
        slope = 0.d0

        if (variable_eta_init) then
            call set_eta_init(nghost_patch, mx_coarse, my_coarse,  &
                 xlow_coarse, ylow_coarse,dx_coarse,dy_coarse,t,vetac)
          endif

        ! Calculate surface elevation eta using dry limiting
        veta_init_c = sea_level  ! if not variable_eta_init
        do j_coarse = 1, my_coarse
            do i_coarse = 1, mx_coarse
                h = valcrse(ivalc(1,i_coarse,j_coarse))
                b = auxcrse(iauxc(i_coarse,j_coarse))

                if (h < dry_tolerance) then
                    if (variable_eta_init) &
                        veta_init_c = vetac(ivetac(i_coarse,j_coarse))
                    eta_coarse(i_coarse,j_coarse) = veta_init_c
                else
                    eta_coarse(i_coarse,j_coarse) = h + b
                endif
            enddo
        enddo

        !!DIG: Something left out here that sets valcrse for sea level changes
        !! KRB: I think this is handled above by veta_init_c (L272) unless you
        !! think it should be set to sea level in this case if h>drytol.

        ! Calculate limited gradients of coarse grid eta
        do j_coarse = 2, my_coarse - 1
           do i_coarse = 2, mx_coarse - 1

                ! X-Direction
                down_slope = eta_coarse(i_coarse,j_coarse) - eta_coarse(i_coarse-1,j_coarse)
                up_slope = eta_coarse(i_coarse+1,j_coarse) - eta_coarse(i_coarse,j_coarse)
                if (up_slope * down_slope > 0.d0) then
                    slope(1,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope)) &
                        * sign(1.d0,eta_coarse(i_coarse+1,j_coarse) - eta_coarse(i_coarse-1,j_coarse))
                endif

                ! Y-Direction
                down_slope = eta_coarse(i_coarse,j_coarse) - eta_coarse(i_coarse,j_coarse-1)
                up_slope = eta_coarse(i_coarse,j_coarse+1) - eta_coarse(i_coarse,j_coarse)
                if (up_slope * down_slope > 0.d0) then
                    slope(2,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope)) &
                        * sign(1.d0,eta_coarse(i_coarse,j_coarse+1) - eta_coarse(i_coarse,j_coarse-1))
                endif
            enddo
        enddo

        ratio_y = real(refinement_ratio_y,kind=8)  ! needs to be real for "floor" call below
        ratio_x = real(refinement_ratio_x,kind=8)
        ! Loop through patch to be filled, includes multiple coarse cells
        do j_fine = 1, my_patch
            j_coarse     = floor((j_fine + jlo - 1) / ratio_y) - jplo + 1
            ycent_coarse = ylow_coarse + (j_coarse-.5d0)*dy_coarse
            ycent_fine   = ylower + (j_fine-1+jlo + .5d0)*dy_fine
            eta2         = (ycent_fine-ycent_coarse)/dy_coarse
            if (abs(eta2) .gt. .5d0) then
                write(*,*)" filpatch y indexing error: eta2 = ",eta2
            endif
            !eta2old = (-0.5d0 + real(mod(j_fine - 1, refinement_ratio_y),kind=8)) &
            !                        / real(refinement_ratio_y,kind=8)
            do i_fine = 1, mx_patch
                i_coarse     = floor((i_fine+ilo-1) / ratio_x) - iplo + 1
                xcent_coarse = xlow_coarse + (i_coarse-.5d0)*dx_coarse
                xcent_fine   =  xlower + (i_fine-1+ilo + .5d0)*dx_fine
                eta1         = (xcent_fine-xcent_coarse)/dx_coarse
                if (abs(eta1) .gt. .5d0) then
                   write(*,*)" filpatch x indexing error: eta1 = ",eta1
                endif
                !eta1old = (-0.5d0 + real(mod(i_fine - 1, refinement_ratio_x),kind=8)) &
                !                / real(refinement_ratio_x,kind=8)

                if (flaguse(i_fine,j_fine) == 0) then
                    ! Interpolate from coarse cells to fine grid for surface
                    fine_cell_count(i_coarse,j_coarse) = fine_cell_count(i_coarse,j_coarse) + 1
                    eta_fine = eta_coarse(i_coarse,j_coarse) + eta1 * slope(1,i_coarse,j_coarse) &
                                                             + eta2 * slope(2,i_coarse,j_coarse)
                    h_fine = max(eta_fine - aux(1,i_fine + nrowst - 1, j_fine + ncolst - 1), 0.d0)

                    if (variable_eta_init) then
                        veta_init_c = vetac(ivetac(i_coarse,j_coarse))
                        ! else it has been set to the constant sea_level
                    endif

                    if (use_force_dry_this_level &
                            .and. (h_fine > 0) &
                            .and. (t <= tend_force_dry)) then
                        ! check if in force_dry region
                        ii = int((xcent_fine - xlow_fdry + 1d-7) / dx_fdry) + 1
                        jj = int((ycent_fine - ylow_fdry + 1d-7) / dy_fdry) + 1
                        jj = my_fdry - jj  ! since index 1 corresponds to north edge
                        if ((ii>=1) .and. (ii<=mx_fdry) .and. &
                            (jj>=1) .and. (jj<=my_fdry)) then
                            ! grid cell lies in region covered by force_dry,
                            ! check if this cell is forced to be dry
                            ! Otherwise don't change value set above:
                            if (force_dry(ii,jj) == 1) then
                                h_fine = 0.d0
                                endif
                            endif

                        endif


                    valbig(1,i_fine+nrowst-1, j_fine+ncolst-1) = h_fine
                    fine_mass(i_coarse,j_coarse) = fine_mass(i_coarse,j_coarse) + h_fine

                    ! Flag the corresponding coarse cell as needing relimiting
                    ! if one of the fine cells ends up being dry
                    if (h_fine < dry_tolerance) then
                        fine_flag(1,i_coarse,j_coarse) = .true.
                        reloop = .true.
                    endif
                endif
            enddo
        enddo

        ! Momentum Interpolation
        do n = 2, nvar
            slope = 0.d0  ! reinitialize for each var
            do j_coarse = 2, my_coarse - 1
                do i_coarse = 2, mx_coarse - 1

                    ! Determine slopes for interpolation
                    down_slope = (valcrse(ivalc(n,i_coarse,j_coarse)) - valcrse(ivalc(n,i_coarse-1,j_coarse)))
                    up_slope   = (valcrse(ivalc(n,i_coarse+1,j_coarse)) - valcrse(ivalc(n,i_coarse,j_coarse)))
                    if (up_slope * down_slope > 0.d0) then
                        slope(1,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope)) *    &
                                          sign(1.d0, valcrse(ivalc(n,i_coarse+1,j_coarse))    &
                                           - valcrse(ivalc(n,i_coarse-1,j_coarse)))
                    endif

                    down_slope = (valcrse(ivalc(n,i_coarse,j_coarse)) - valcrse(ivalc(n,i_coarse,j_coarse-1)))
                    up_slope   = (valcrse(ivalc(n,i_coarse,j_coarse+1)) - valcrse(ivalc(n,i_coarse,j_coarse)))
                    if (up_slope * down_slope > 0.d0) then
                        slope(2,i_coarse,j_coarse) = min(abs(up_slope), abs(down_slope))       &
                                         * sign(1.d0, valcrse(ivalc(n,i_coarse,j_coarse+1))    &
                                          - valcrse(ivalc(n,i_coarse,j_coarse-1)))
                    endif

                    ! Set initial values for max/min of current field
                    if (valcrse(ivalc(1,i_coarse,j_coarse)) > dry_tolerance) then
                        vel_max(i_coarse,j_coarse) = valcrse(ivalc(n,i_coarse,j_coarse)) /    &
                                                     valcrse(ivalc(1,i_coarse,j_coarse))
                        vel_min(i_coarse,j_coarse) = valcrse(ivalc(n,i_coarse,j_coarse)) /    &
                                                     valcrse(ivalc(1,i_coarse,j_coarse))
                    else
                        ! refinement has made mass in a fine cell that was in a
                        ! previously dry coarse cell. This should only happen for
                        ! sea level refinement.

                        ! for landslide material, refinement should result in a dry fine cell.
                        ! if this is not the case, there may be something wrong
                        ! with how the refinement is specified.

                        ! the case for which this approach is insufficient is if one
                        ! had a sediment laden lake that we want to have unrefined
                        ! at the beginning of the simulation and want to refine to a
                        ! flat surface and retain the m value when velocity arrives.

                        vel_min(i_coarse,j_coarse) = 0.d0
                        vel_max(i_coarse,j_coarse) = 0.d0

                        !!DIG: verify that these are correct...
                        !! initial values only set here if h<drytol.
                        if (n == 4) then
                            ! solid volume fraction
                            vel_max(i_coarse,j_coarse) = 0.d0 ! if refinement makes new mass in a dry cell, assume it is water.
                            vel_min(i_coarse,j_coarse) = 0.d0 ! KRB and DLG changed this from m0 to 0.d0 4/2/2024

                         elseif (n == 5) then
                            ! pore pressure
                            vel_max(i_coarse,j_coarse) = rho_f*grav
                            vel_min(i_coarse,j_coarse) = rho_f*grav
                         endif
                    endif

                    ! Look for bounds on velocity around each cell
                    ! Necessary since we are interpolating momentum linearly
                    ! but not interpolating depth linearly
                    do i =-1,1,2
                        if (valcrse(ivalc(1,i_coarse + i,j_coarse)) > dry_tolerance) then
                            vel_max(i_coarse,j_coarse) =  max(vel_max(i_coarse,j_coarse),      &
                                                valcrse(ivalc(n,i_coarse + i,j_coarse)) /      &
                                                valcrse(ivalc(1,i_coarse + i,j_coarse)))
                            vel_min(i_coarse,j_coarse) = min(vel_min(i_coarse,j_coarse),       &
                                                valcrse(ivalc(n,i_coarse + i,j_coarse)) /      &
                                                valcrse(ivalc(1,i_coarse + i,j_coarse)))
                        endif

                        if (valcrse(ivalc(1,i_coarse,j_coarse + i)) > dry_tolerance) then
                            vel_max(i_coarse,j_coarse) = max(vel_max(i_coarse,j_coarse),       &
                                                  valcrse(ivalc(n,i_coarse,j_coarse + i)) /    &
                                                  valcrse(ivalc(1,i_coarse,j_coarse + i)))
                            vel_min(i_coarse,j_coarse) = min(vel_min(i_coarse,j_coarse),       &
                                                  valcrse(ivalc(n,i_coarse,j_coarse + i)) /    &
                                                  valcrse(ivalc(1,i_coarse,j_coarse + i)))

                        endif
                    enddo
                enddo
            enddo

            ! Determine momentum in fine cells
            do j_fine = 1, my_patch
                j_coarse     = floor((j_fine + jlo - 1) / ratio_y) - jplo + 1
                ycent_coarse = ylow_coarse + (j_coarse-.5d0)*dy_coarse
                ycent_fine   =  ylower + (j_fine-1+jlo + .5d0)*dy_fine
                eta2         = (ycent_fine-ycent_coarse)/dy_coarse

                do i_fine = 1, mx_patch
                    i_coarse     = floor((i_fine+ilo-1) / ratio_x) - iplo + 1
                    xcent_coarse = xlow_coarse + (i_coarse-.5d0)*dx_coarse
                    xcent_fine   =  xlower + (i_fine-1+ilo + .5d0)*dx_fine
                    eta1         = (xcent_fine-xcent_coarse)/dx_coarse

                    if (flaguse(i_fine,j_fine) == 0) then
                        ! Cell not already set

                        if (.not.(fine_flag(1,i_coarse,j_coarse))) then
                            ! This cell has no coarse cells that are dry
                            hv_fine = valcrse(ivalc(n,i_coarse,j_coarse))         &
                                            + eta1 * slope(1,i_coarse,j_coarse)   &
                                            + eta2 * slope(2,i_coarse,j_coarse)
                            v_fine = hv_fine  / valbig(1,i_fine+nrowst-1, j_fine+ncolst-1)
                            if (v_fine<vel_min(i_coarse,j_coarse) .or.  v_fine>vel_max(i_coarse,j_coarse)) then
                                fine_flag(n,i_coarse,j_coarse) = .true.
                                reloop = .true.
                            else
                                valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = hv_fine
                            endif
                        endif
                    endif
                enddo
            enddo

            ! Reset momentum to conserve momentum in the cases where we may have
            ! gained momentum or if velocity bounds were violated
            if (reloop) then

                do j_fine  = 1, my_patch
                  j_coarse     = floor((j_fine + jlo - 1) / ratio_y) - jplo + 1
                  ycent_coarse = ylow_coarse + (j_coarse-.5d0)*dy_coarse
                  ycent_fine   =  ylower + (j_fine-1+jlo + .5d0)*dy_fine
                  eta2         = (ycent_fine-ycent_coarse)/dy_coarse

                    do i_fine = 1, mx_patch
                       i_coarse     = floor((i_fine+ilo-1) / ratio_x) - iplo + 1
                       xcent_coarse = xlow_coarse + (i_coarse-.5d0)*dx_coarse
                       xcent_fine   =  xlower + (i_fine-1+ilo + .5d0)*dx_fine
                       eta1         = (xcent_fine-xcent_coarse)/dx_coarse

                        if (flaguse(i_fine,j_fine) == 0) then
                            if (fine_flag(1,i_coarse,j_coarse) .or. fine_flag(n,i_coarse,j_coarse)) then
                                if (fine_mass(i_coarse,j_coarse) > dry_tolerance) then
                                    h_coarse = valcrse(ivalc(1,i_coarse,j_coarse))
                                    h_count = real(fine_cell_count(i_coarse,j_coarse),kind=8)
                                    h_fine_average = fine_mass(i_coarse,j_coarse) / h_count
                                    divide_mass = max(h_coarse, h_fine_average)
                                    h_fine = valbig(1, i_fine + nrowst - 1, j_fine + ncolst - 1)
                                    v_new = valcrse(ivalc(n,i_coarse,j_coarse)) / (divide_mass)
                                    if (n > 3) then
                                        !!DIG change
                                        v_new = max(vel_min(i_coarse,j_coarse), v_new)
                                        v_new = min(vel_max(i_coarse,j_coarse), v_new)
                                    endif
                                    valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = &
                                        v_new * valbig(1,i_fine+nrowst-1,j_fine+ncolst-1)
                                else
                                    valbig(n,i_fine+nrowst-1,j_fine+ncolst-1) = 0.d0
                                endif
                            endif
                        endif
                    enddo
                enddo
            endif ! end if reloop
        enddo  ! end loop over nvar
    endif   ! end if patch not set

    ! set bcs, whether or not recursive calls needed. set any part of patch that stuck out
    xhi_fine = xlower + (ihi + 1) * dx_fine
    yhi_fine = ylower + (jhi + 1) * dy_fine
    if (patchOnly) then
       call bc2amr(valbig,aux,mx,my,nvar,naux,dx_fine,dy_fine,level,t,   &
                   xlow_fine,xhi_fine,ylow_fine,yhi_fine)
    endif

contains

    integer pure function ivalc(n,i,j)
        implicit none
        integer, intent(in) :: n,i,j
        ivalc = n + nvar*(i-1) + nvar*mx_coarse*(j-1)
    end function ivalc

    ! Index into first component of aux = topo:
    integer pure function iauxc(i,j)
        implicit none
        integer, intent(in) :: i,j
        iauxc = 1 + naux*(i-1) + naux*mx_coarse*(j-1)
    end function iauxc

    ! Index into vetac:
    integer pure function ivetac(i,j)
        implicit none
        integer, intent(in) :: i,j
        ivetac = 1 + (i-1) + mx_coarse*(j-1)
    end function ivetac

    ! logical for checking if this patch sticks outside of the domain
    logical pure function sticksout(iplo,iphi,jplo,jphi)
        implicit none
        integer, intent(in) :: iplo,iphi,jplo,jphi
        sticksout = (iplo < 0 .or. jplo < 0 .or. &
                     iphi >= iregsz(level - 1) .or. jphi >= jregsz(level - 1))
    end function sticksout

end subroutine filrecur

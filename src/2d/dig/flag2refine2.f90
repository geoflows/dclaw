!! D-Claw specific core file
!! This file is a modified version of
!! clawpack/geoclaw/src/2d/shallow/flag2refine.f90 
!!

! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for D-Claw: This routine uses flowgrades and if specified
! conditions are met, a cell will be flagged for refinement. If
! keep_fine = True, conditions will be evaluated at both this level and
! at finer levels (up to level = mxnest)

! D-Claw refinement does not consider wave or speed tolerances as these
! can be specified by flowgrades.
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
subroutine flag2refine2(mx, my, mbc, mbuff, meqn, maux, xlower, ylower, dx, dy, t, level, &
                        tolsp, q, aux, amrflags)

   use amr_module, only: mxnest, t0, DOFLAG, UNSET
   use amr_module, only: lstart, node, rnode
   use amr_module, only: cornxlo, cornylo, levelptr, mxnest, ndihi, ndilo
   use amr_module, only: ndjhi, ndjlo, store1, storeaux
   use amr_module, only: alloc, hxposs, hyposs

   use refinement_module, only: mflowgrades, keep_fine, fine_level

   use digclaw_module, only: i_h
   use geoclaw_module, only: dry_tolerance

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: mx, my, mbc, meqn, maux, level, mbuff
   real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, tolsp

   real(kind=8), intent(in) :: q(meqn, 1 - mbc:mx + mbc, 1 - mbc:my + mbc)
   real(kind=8), intent(inout) :: aux(maux, 1 - mbc:mx + mbc, 1 - mbc:my + mbc)

   ! Flagging
   real(kind=8), intent(in out) :: amrflags(1 - mbuff:mx + mbuff, 1 - mbuff:my + mbuff)

   ! Generic locals
   integer :: i, j
   real(kind=8) :: x_c, y_c, x_low, y_low, x_hi, y_hi
   logical :: FGFLAG

   real(kind=8) :: qstencil(meqn,5)
   real(kind=8) :: auxstencil(maux,5)

    !! Flowgrade and keep_fine variables
   real(kind=8) :: xlowfg, xhifg, ylowfg, yhifg, xxlow, xxhi, yylow, yyhi
   real(kind=8) :: dxfine, dyfine
   integer :: nx, ny, q_loc, aux_loc, mitot, mjtot, grid_ptr, ii, jj, mm

   ! Loop over interior points on this grid
   ! (i,j) grid cell is [x_low,x_hi] x [y_low,y_hi], cell cen at (x_c,y_c)
   y_loop: do j = 1, my
      y_low = ylower + (j - 1)*dy
      y_c = ylower + (j - 0.5d0)*dy
      y_hi = ylower + j*dy

      x_loop: do i = 1, mx
         x_low = xlower + (i - 1)*dx
         x_c = xlower + (i - 0.5d0)*dx
         x_hi = xlower + i*dx

         if (amrflags(i, j) .ne. UNSET) then
            ! do not consider this cell further
            cycle x_loop
         end if

         ! determine if flowgrades are used
         FGFLAG = .false.
         if ((mflowgrades > 0).and. &
             ((q(i_h, i, j).gt.dry_tolerance)).and. &
             (.not. keep_fine)) then

         ! Need to change this logic so that instead of .not.keep_fine
         ! you can enter this loop if keep fine is True but that
         ! level is .ge. fine_level


         ! only check flowgrades if central cell has thickness
         ! DIG: Discuss

            qstencil(:,:) = 0.d0
            auxstencil(:,:) = 0.d0
            do mm = 1, meqn
                qstencil(mm,1) = q(mm, i, j)
                qstencil(mm,2) = q(mm, i - 1, j)
                qstencil(mm,3) = q(mm, i + 1, j)
                qstencil(mm,4) = q(mm, i, j - 1)
                qstencil(mm,5) = q(mm, i, j + 1)
            end do

            do mm = 1, maux
                auxstencil(mm,1) = aux(mm, i, j)
                auxstencil(mm,2) = aux(mm, i - 1, j)
                auxstencil(mm,3) = aux(mm, i + 1, j)
                auxstencil(mm,4) = aux(mm, i, j - 1)
                auxstencil(mm,5) = aux(mm, i, j + 1)
            end do

            call check_flowgrades(meqn, maux, level, &
                                  qstencil, auxstencil, &
                                  dx, dy, FGFLAG)
            if (FGFLAG) then
               amrflags(i, j) = DOFLAG
               cycle x_loop
            end if

         end if ! if using normal flowgraded and cell is wet.

         ! keep fine added by KRB 2022/12/28
         ! updated by KRB 2024/10/17
         ! if keep_fine is used, refinement must be based on the highest
         ! level of refinement only. checking flowgrades on all levels (as is
         ! done in the code block above) is problematic because if a coarse
         ! grid at any level advects material beyond the extent of a fine grid,
         ! then at the next regrid time point, this material will be placed into
         ! the fine grid.

         ! requirement for use of keep fine is that all flow of concern is on
         ! mxnest only and the sole purpose of keep fine is amr is for
         ! 'tracking' the location of fine grids (e.g. shallow flow problems).

         if (keep_fine .and. mflowgrades .gt. 0) then

            ! if level is lower than fine_level determine whether a grid exists
            ! here on each level between level+1 and fine_level. For each cell
            ! within that grid, test whether refinement should be maintained.
            ! ignore ghost cells on the finer grid

            if (level .lt. fine_level) then

              ! consider flow on fine_level (not highest here, only fine_level)

              ! get pointer to the start of this fine_level
              grid_ptr = lstart(fine_level)

              ! might need to update to know where the end of the fine_level grids are.

              ! get fine dx and dy values
              dxfine = hxposs(fine_level)
              dyfine = hyposs(fine_level)

              do while (grid_ptr /= 0)

                ! calculate extent of fine grid (no ghost cells reflected
                ! in extent)
                nx = node(ndihi, grid_ptr) - node(ndilo, grid_ptr) + 1
                ny = node(ndjhi, grid_ptr) - node(ndjlo, grid_ptr) + 1

                xlowfg = rnode(cornxlo, grid_ptr)
                ylowfg = rnode(cornylo, grid_ptr)
                xhifg = xlowfg + nx*dxfine
                yhifg = ylowfg + ny*dyfine

                ! if there is overlap between the fine grid extent and this
                ! cell on the coarse grid, consider further.
                ! (i,j) coarse grid cell is [x_low,x_hi] x [y_low,y_hi].
                ! fine grid extent is [xlowfg,xhifg] x [ylowfg,yhifg]

                if ((x_hi .gt. xlowfg) .and. &
                    (x_low .lt. xhifg) .and. &
                    (y_hi .gt. ylowfg) .and. &
                    (y_low .lt. yhifg)) then

                   ! If considering further, loop through fine grid cells to
                   ! find the ones that overlap. indexing into q on the fine
                   ! grid includes ghost cells, so update the fine grid
                   ! lower left corner location and the nrows, ncols to
                   ! include ghost cells. other corners no longer needed.

                   xlowfg = rnode(cornxlo, grid_ptr) - mbc*dxfine
                   ylowfg = rnode(cornylo, grid_ptr) - mbc*dyfine

                   ! define mitot and mjtot to include ghost cells.
                   ! iadd and iaddaux modified to reflect this.
                   mitot = nx + 2*mbc
                   mjtot = ny + 2*mbc

                   ! get pointers into q and aux
                   q_loc = node(store1, grid_ptr)
                   aux_loc = node(storeaux, grid_ptr)

                   ! loop through fine grid cells (consider only the non-
                   ! ghost cell)
                   do jj = mbc + 1, mjtot - mbc
                      do ii = mbc + 1, mitot - mbc
                         ! check overlap between fine and coarse cells
                         ! ignore ghost cells
                         xxlow = xlowfg + dxfine*ii
                         xxhi = xxlow + dxfine
                         yylow = ylowfg + dyfine*jj
                         yyhi = yylow + dyfine

                         ! if fine grid cell overlaps with coarse grid cell
                         ! consider further.
                         if ((x_hi .gt. xxlow) .and. &
                             (x_low .lt. xxhi) .and. &
                             (y_hi .gt. yylow) .and. &
                             (y_low .lt. yyhi)) then

                            ! test if central cell has thickness.
                            if (alloc(iadd(i_h, ii, jj)).gt.dry_tolerance) then
                              ! only check flowgrades if central cell in stencil has thickness
                              ! DIG: Discuss

                              ! construct arrays of q and aux for the five cell
                              ! stencil surrounding cell ii,jj.
                              ! stencil may include ghost cell values.
                              qstencil(:,:) = 0.d0
                              auxstencil(:,:) = 0.d0

                              do mm = 1, meqn
                                 qstencil(mm,1) = alloc(iadd(mm, ii, jj))
                                 qstencil(mm,2) = alloc(iadd(mm, ii - 1, jj))
                                 qstencil(mm,3) = alloc(iadd(mm, ii + 1, jj))
                                 qstencil(mm,4) = alloc(iadd(mm, ii, jj - 1))
                                 qstencil(mm,5) = alloc(iadd(mm, ii, jj + 1))
                              end do

                              do mm = 1, maux
                                 auxstencil(mm,1) = alloc(iaddaux(mm, ii, jj))
                                 auxstencil(mm,2) = alloc(iaddaux(mm, ii - 1, jj))
                                 auxstencil(mm,3) = alloc(iaddaux(mm, ii + 1, jj))
                                 auxstencil(mm,4) = alloc(iaddaux(mm, ii, jj - 1))
                                 auxstencil(mm,5) = alloc(iaddaux(mm, ii, jj + 1))

                              end do

                              ! call check_flowgrades
                              call check_flowgrades(meqn, maux, level, &
                                                    qstencil, auxstencil, &
                                                    dxfine, dyfine, &
                                                    FGFLAG)

                              if (FGFLAG) then
                                 amrflags(i, j) = DOFLAG
                                 cycle x_loop
                              end if ! end if FGFLAG

                            end if ! end if central fine cell has thickness
                         end if ! end if coarse and fine cell overlap
                      end do ! endif fine cell loop ii
                   end do ! endif fine cell loop jj
                end if ! endif coarse cell overlaps fine grid

                ! update grid pointer to look at next grid at this level.
                grid_ptr = node(levelptr, grid_ptr)
              end do ! end while (grid_ptr > 0)

            end if ! end level .lt. mxnest
         end if ! end if keep_fine




      end do x_loop
   end do y_loop

contains

    ! Index into q array
    ! KRB modified: assumes mitot includes ghost cells
    pure integer function iadd(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iadd = q_loc + m - 1 + meqn*((j - 1)*(mitot) + i - 1)
    end function iadd

    ! Index into aux array
    ! KRB modified: assumes mitot includes ghost cells
    pure integer function iaddaux(m, i, j)
        implicit none
        integer, intent(in) :: m, i, j
        iaddaux = aux_loc + m - 1 + maux*(i - 1) + maux*(mitot)*(j - 1)
    end function iaddaux

end subroutine flag2refine2

subroutine check_flowgrades(meqn, maux, level, &
                            qstencil, auxstencil, &
                            dx, dy, FGFLAG)

   use refinement_module, only: flowgradevalue, iflowgradevariable, &
                                iflowgradetype, iflowgrademinlevel, mflowgrades
   use geoclaw_module, only: dry_tolerance, sea_level
   use digclaw_module, only: i_h, i_hu, i_hv, i_bdif

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: meqn, maux, level

   real(kind=8), intent(in) :: qstencil(meqn,5)
   real(kind=8), intent(in) :: auxstencil(maux,5)

   real(kind=8), intent(in) :: dx, dy
   logical, intent(inout) :: FGFLAG

    !! Locals
   real(kind=8) :: flowgrademeasure,lef,rig,bot,top
   integer :: iflow

   mfg_loop: do iflow = 1, mflowgrades

      ! flowgradevariable = 1, depth
      if (iflowgradevariable(iflow) .eq. 1) then

         ! flowgradetype = 1, magnitude
         if (iflowgradetype(iflow) .eq. 1) then

            flowgrademeasure = qstencil(i_h,1)

            ! flowgradetype != 1, gradient
         else

            lef = qstencil(i_h,2)
            rig = qstencil(i_h,3)
            bot = qstencil(i_h,4)
            top = qstencil(i_h,5)


            flowgrademeasure = sqrt(((rig - lef)/(2.0d0*dx))**2 + ((top - bot)/(2.0d0*dy))**2)

         end if

         ! flowgradevariable = 2, velocity
      elseif (iflowgradevariable(iflow) .eq. 2) then

         ! flowgradetype = 1, magnitude
         if (iflowgradetype(iflow) .eq. 1) then

            if (qstencil(i_h,1) .gt. dry_tolerance) then ! only consider if center cell is wet.
               flowgrademeasure = sqrt((qstencil(i_hu,1)/qstencil(i_h,1))**2 + (qstencil(i_hv,1)/qstencil(i_h,1))**2)
            else
               flowgrademeasure = 0.d0
            end if

            ! flowgradetype != 1, gradient
         else

            if ((qstencil(i_h,1) .gt. dry_tolerance) .and. &
                (qstencil(i_h,2) .gt. dry_tolerance) .and. &
                (qstencil(i_h,3) .gt. dry_tolerance) .and. &
                (qstencil(i_h,4) .gt. dry_tolerance) .and. &
                (qstencil(i_h,5) .gt. dry_tolerance)) then ! only consider if all stencil cells are wet.

               lef = sqrt((qstencil(i_hu,2)/qstencil(i_h,2))**2 + (qstencil(i_hv,2)/qstencil(i_h,2))**2)
               rig = sqrt((qstencil(i_hu,3)/qstencil(i_h,3))**2 + (qstencil(i_hv,3)/qstencil(i_h,3))**2)
               bot = sqrt((qstencil(i_hu,4)/qstencil(i_h,4))**2 + (qstencil(i_hv,4)/qstencil(i_h,4))**2)
               top = sqrt((qstencil(i_hu,5)/qstencil(i_h,5))**2 + (qstencil(i_hv,5)/qstencil(i_h,5))**2)

               flowgrademeasure = sqrt(((rig - lef)/(2.0d0*dx))**2 + ((top - bot)/(2.0d0*dy))**2)
            else
               flowgrademeasure = 0.d0
            end if

         end if

         ! flowgradevariable = 3, perturbation from sea level
      elseif (iflowgradevariable(iflow) .eq. 3) then

         ! flowgradetype = 1, magnitude
         if (iflowgradetype(iflow) .eq. 1) then

            ! only calculate surface perturbation if center cell is wet.
            if (qstencil(i_h,1) .gt. dry_tolerance) then
               flowgrademeasure = dabs(qstencil(i_h,1) + auxstencil(1,1) - qstencil(i_bdif,1) - sea_level)
            else
               flowgrademeasure = 0.d0
            end if

            ! flowgradetype != 1, gradient
         else
            if ((qstencil(i_h,1) .gt. dry_tolerance) .and. &
                (qstencil(i_h,2) .gt. dry_tolerance) .and. &
                (qstencil(i_h,3) .gt. dry_tolerance) .and. &
                (qstencil(i_h,4) .gt. dry_tolerance) .and. &
                (qstencil(i_h,5) .gt. dry_tolerance)) then ! only consider if all stencil cells are wet.

               lef = qstencil(i_h,2) + auxstencil(1,2) - qstencil(i_bdif,2)
               rig = qstencil(i_h,3) + auxstencil(1,3) - qstencil(i_bdif,3)
               bot = qstencil(i_h,4) + auxstencil(1,4) - qstencil(i_bdif,4)
               top = qstencil(i_h,5) + auxstencil(1,5) - qstencil(i_bdif,5)

               flowgrademeasure = sqrt(((rig - lef)/(2.0d0*dx))**2 + ((top - bot)/(2.0d0*dy))**2)
            else
               flowgrademeasure = 0.d0
            end if
         end if
      else
         write (*, *) '+++++ only flowgradevariable = 1,2,3 supported'
         stop
      end if

      if ((flowgrademeasure .gt. flowgradevalue(iflow)) .and. &
          (level .lt. iflowgrademinlevel(iflow))) then
         FGFLAG = .true.
         exit mfg_loop
      end if

   end do mfg_loop! end mflowgrades loop

end subroutine check_flowgrades

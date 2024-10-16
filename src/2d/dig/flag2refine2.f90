! ::::::::::::::::::::: flag2refine ::::::::::::::::::::::::::::::::::
!
! User routine to control flagging of points for refinement.
!
! Specific for D-Claw: This routine uses flowgrades and if specified
! conditions are met, a cell will be flagged for refinement. If
! keep_fine = True, conditions will be evaluated at both this level and
! the finest level (level = mxnest, not the finest level at this
! location)
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
   use amr_module, only: lfine, lstart, node, rnode
   use amr_module, only: cornxlo, cornylo, levelptr, mxnest, ndihi, ndilo
   use amr_module, only: ndjhi, ndjlo, store1, store2, storeaux
   use amr_module, only: alloc, hxposs, hyposs

   use refinement_module, only: mflowgrades, keep_fine

   implicit none

   ! Subroutine arguments
   integer, intent(in) :: mx, my, mbc, meqn, maux, level, mbuff
   real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, tolsp

   real(kind=8), intent(in) :: q(meqn, 1 - mbc:mx + mbc, 1 - mbc:my + mbc)
   real(kind=8), intent(inout) :: aux(maux, 1 - mbc:mx + mbc, 1 - mbc:my + mbc)

   ! Flagging
   real(kind=8), intent(in out) :: amrflags(1 - mbuff:mx + mbuff, 1 - mbuff:my + mbuff)

   ! Generic locals
   integer :: i, j, m, r
   real(kind=8) :: x_c, y_c, x_low, y_low, x_hi, y_hi
   real(kind=8) :: speed, eta, ds
   logical :: FGFLAG

   real(kind=8) :: qcen(meqn)
   real(kind=8) :: qtop(meqn)
   real(kind=8) :: qbot(meqn)
   real(kind=8) :: qlef(meqn)
   real(kind=8) :: qrig(meqn)
   real(kind=8) :: auxcen(maux)
   real(kind=8) :: auxtop(maux)
   real(kind=8) :: auxbot(maux)
   real(kind=8) :: auxlef(maux)
   real(kind=8) :: auxrig(maux)

    !! Flowgrade and keep_fine variables
   real(kind=8) :: xlow, xhi, ylow, yhi, xxlow, xxhi, yylow, yyhi
   real(kind=8) :: dxfine, dyfine
   integer :: nx, ny, loc, locaux, mitot, mjtot, mptr, iflow, ii, jj, mm

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
         if (mflowgrades > 0) then
            call check_flowgrades(meqn, maux, level, &
                                  q(:, i, j), & ! center
                                  q(:, i - 1, j), & ! left
                                  q(:, i + 1, j), & ! right
                                  q(:, i, j + 1), & ! top
                                  q(:, i, j - 1), & ! bottom
                                  aux(:, i, j), & ! center
                                  aux(:, i - 1, j), & ! left
                                  aux(:, i + 1, j), & ! right
                                  aux(:, i, j + 1), & ! top
                                  aux(:, i, j - 1), & ! bottom
                                  dx, dy, FGFLAG)
            if (FGFLAG) then
               amrflags(i, j) = DOFLAG
               !write(*,*) "++++++++++ flagged"
               cycle x_loop
            end if
         end if

         ! keep fine added by KRB 2022/12/28
         if (keep_fine .and. mflowgrades .gt. 0) then
            ! if level is lower than lfine determine whether a grid exists
            ! here on lfine and test whether refinement should be maintained.
            ! ignore ghost cells on the fine grid

            if (level .lt. lfine) then

               ! loop through lfine grids
               mptr = lstart(lfine)
               do while (mptr > 0)

                  ! calculate extent of fine grid
                  nx = node(ndihi, mptr) - node(ndilo, mptr) + 1
                  ny = node(ndjhi, mptr) - node(ndjlo, mptr) + 1
                  loc = node(store1, mptr)
                  locaux = node(storeaux, mptr)
                  mitot = nx + 2*mbc
                  mjtot = ny + 2*mbc

                  xlow = rnode(cornxlo, mptr)
                  ylow = rnode(cornylo, mptr)
                  xhi = xlow + nx*hxposs(lfine)
                  yhi = ylow + ny*hyposs(lfine)

                  ! if there is overlap between the fine grid and this
                  ! location on the coarse grid, flag.
                  ! (i,j) grid cell is [x1,x2] x [y1,y2].
                  ! fine grid is [xlow,xhi] x [ylow,yhi]
                  if (x_hi .gt. xlow .and. x_low .lt. xhi .and. &
                      y_hi .gt. ylow .and. y_low .lt. yhi) then

                     ! this loop includes ghosts, update extent and
                     mitot = nx + 2*mbc
                     mjtot = ny + 2*mbc
                     xlow = rnode(cornxlo, mptr) - mbc*hxposs(lfine)
                     ylow = rnode(cornylo, mptr) - mbc*hyposs(lfine)
                     xhi = xlow + mitot*hxposs(lfine)
                     yhi = ylow + mjtot*hyposs(lfine)

                     ! get pointers
                     loc = node(store1, mptr)
                     locaux = node(storeaux, mptr)

                     ! loop through fine grid cells.
                     do jj = mbc + 1, mjtot - mbc
                        do ii = mbc + 1, mitot - mbc
                           ! check overlap between fine and coarse cells
                           ! ignore ghost cells
                           xxlow = xlow + hxposs(lfine)*ii
                           xxhi = xxlow + hxposs(lfine)
                           yylow = ylow + hyposs(lfine)*jj
                           yyhi = yylow + hyposs(lfine)

                           ! if fine grid cell is inside of coarse grid cell
                           ! calculate flowgrade values
                           if (x_hi .gt. xxlow .and. x_low .lt. xxhi .and. &
                               y_hi .gt. yylow .and. y_low .lt. yyhi) then

                              ! construct arrays of q and aux for the five cell
                              ! stencil surrounding cell ii,jj.
                              do mm = 1, meqn
                                 qcen(mm) = alloc(iadd(mm, ii, jj))
                                 qlef(mm) = alloc(iadd(mm, ii - 1, jj))
                                 qrig(mm) = alloc(iadd(mm, ii + 1, jj))
                                 qtop(mm) = alloc(iadd(mm, ii, jj + 1))
                                 qbot(mm) = alloc(iadd(mm, ii, jj - 1))
                              end do

                              do mm = i, maux
                                 auxcen(mm) = alloc(iaddaux(mm, ii, jj))
                                 auxlef(mm) = alloc(iaddaux(mm, ii - 1, jj))
                                 auxrig(mm) = alloc(iaddaux(mm, ii + 1, jj))
                                 auxtop(mm) = alloc(iaddaux(mm, ii, jj + 1))
                                 auxbot(mm) = alloc(iaddaux(mm, ii, jj - 1))
                              end do

                              ! get fine dx and dy values
                              dxfine = hxposs(lfine)
                              dyfine = hyposs(lfine)

                              ! call check_flowgrades
                              call check_flowgrades(meqn, maux, level, &
                                                    qcen, qlef, qrig, qtop, qbot, &
                                                    auxcen, auxlef, auxrig, auxtop, auxbot, &
                                                    dxfine, dyfine, &
                                                    FGFLAG)

                              if (FGFLAG) then
                                 amrflags(i, j) = DOFLAG
                                 !write(*,*) "++++++++++ flagged - keep fine"
                                 cycle x_loop
                              end if

                           end if ! end if coarse and fine cell overlap
                        end do ! endif fine cell loop ii
                     end do ! endif fine cell loop jj

                  end if ! endif coarse cell overlaps fine grid
                  mptr = node(levelptr, mptr)
               end do ! end while (mptr > 0)

            end if ! (end level .lt. lfine)
         end if ! (end if keep_fine)
         ! end keep fine

      end do x_loop
   end do y_loop

contains

   ! Index into q array
   pure integer function iadd(m, i, j)
      implicit none
      integer, intent(in) :: m, i, j
      iadd = loc + m - 1 + meqn*((j - 1)*(mitot + 2*mbc) + i - 1)
   end function iadd

   ! Index into aux array
   pure integer function iaddaux(m, i, j)
      implicit none
      integer, intent(in) :: m, i, j
      iaddaux = locaux + m - 1 + maux*(i - 1) + maux*(mitot + 2*mbc)*(j - 1)
   end function iaddaux

end subroutine flag2refine2

subroutine check_flowgrades(meqn, maux, level, &
                            qcen, qlef, qrig, qtop, qbot, &
                            auxcen, auxlef, auxrig, auxtop, auxbot, &
                            dx, dy, &
                            FGFLAG)

   use refinement_module, only: flowgradevalue, iflowgradevariable, &
                                iflowgradetype, iflowgrademinlevel, mflowgrades
   use geoclaw_module, only: dry_tolerance, sea_level
   use digclaw_module, only: i_h, i_hu, i_hv, i_bdif

   ! Subroutine arguments
   integer, intent(in) :: meqn, maux

   real(kind=8), intent(in) :: qcen(meqn)
   real(kind=8), intent(in) :: qlef(meqn)
   real(kind=8), intent(in) :: qrig(meqn)
   real(kind=8), intent(in) :: qtop(meqn)
   real(kind=8), intent(in) :: qbot(meqn)

   real(kind=8), intent(in) :: auxcen(maux)
   real(kind=8), intent(in) :: auxlef(maux)
   real(kind=8), intent(in) :: auxrig(maux)
   real(kind=8), intent(in) :: auxbot(maux)
   real(kind=8), intent(in) :: auxtop(maux)

   real(kind=8), intent(in) :: dy, dx
   logical, intent(inout) :: FGFLAG

    !! Locals
   real(kind=8) :: flowgrademeasure
   integer :: iflow

   FGFLAG = .false.

   mfg_loop: do iflow = 1, mflowgrades

      ! flowgradevariable = 1, depth
      if (iflowgradevariable(iflow) .eq. 1) then

         ! flowgradetype = 1, magnitude
         if (iflowgradetype(iflow) .eq. 1) then

            flowgrademeasure = qcen(i_h)

            ! flowgradetype != 1, gradient
         else

            bot = qbot(i_h)
            top = qtop(i_h)
            lef = qlef(i_h)
            rig = qrig(i_h)

            flowgrademeasure = sqrt(((rig - lef)/(2.0d0*dx))**2 + ((top - bot)/(2.0d0*dy))**2)

         end if

         ! flowgradevariable = 2, velocity
      elseif (iflowgradevariable(iflow) .eq. 2) then

         ! flowgradetype = 1, magnitude
         if (iflowgradetype(iflow) .eq. 1) then

            if (qcen(i_h) .gt. dry_tolerance) then ! only consider if center cell is wet.
               flowgrademeasure = sqrt((qcen(i_hu)/qcen(i_h))**2 + (qcen(i_hv)/qcen(i_h))**2)
            else
               flowgrademeasure = 0.d0
            end if

            ! flowgradetype != 1, gradient
         else

            if ((qcen(i_h) .gt. dry_tolerance) .and. &
                (qtop(i_h) .gt. dry_tolerance) .and. &
                (qbot(i_h) .gt. dry_tolerance) .and. &
                (qlef(i_h) .gt. dry_tolerance) .and. &
                (qrig(i_h) .gt. dry_tolerance)) then ! only consider if all stencil cells are wet.

               bot = sqrt((qbot(i_hu)/qbot(i_h))**2 + (qbot(i_hv)/qbot(i_h))**2)
               top = sqrt((qtop(i_hu)/qtop(i_h))**2 + (qtop(i_hv)/qtop(i_h))**2)
               lef = sqrt((qlef(i_hu)/qlef(i_h))**2 + (qlef(i_hv)/qlef(i_h))**2)
               rig = sqrt((qrig(i_hu)/qrig(i_h))**2 + (qrig(i_hv)/qrig(i_h))**2)

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
            if (qcen(i_h) .gt. dry_tolerance) then
               flowgrademeasure = dabs(qcen(i_h) + auxcen(1) - qcen(i_bdif) - sea_level)
            else
               flowgrademeasure = 0.d0
            end if

            ! flowgradetype != 1, gradient
         else
            if ((qcen(i_h) .gt. dry_tolerance) .and. &
                (qtop(i_h) .gt. dry_tolerance) .and. &
                (qbot(i_h) .gt. dry_tolerance) .and. &
                (qlef(i_h) .gt. dry_tolerance) .and. &
                (qrig(i_h) .gt. dry_tolerance)) then ! only consider if all stencil cells are wet.

               bot = qbot(i_h) + auxbot(1) - qbot(i_bdif)
               top = qtop(i_h) + auxtop(1) - qtop(i_bdif)
               lef = qlef(i_h) + auxlef(1) - qlef(i_bdif)
               rig = qrig(i_h) + auxrig(1) - qrig(i_bdif)

               flowgrademeasure = sqrt(((rig - lef)/(2.0d0*dx))**2 + ((top - bot)/(2.0d0*dy))**2)
            else
               flowgrademeasure = 0.d0
            end if
         end if
      else
         write (*, *) '+++++ only flowgradevariable = 1,2,3 supported'
         stop

      end if

      if (flowgrademeasure .gt. flowgradevalue(iflow) &
          .and. level .lt. iflowgrademinlevel(iflow)) then
         FGFLAG = .true.
         cycle mfg_loop
      end if

   end do mfg_loop! end mflowgrades loop

end subroutine check_flowgrades

! Set auxiliary arrays
!
! Routine that sets the values found in the aux array including topography and
! other fields that are used in the source terms and other things.
!
! In the default routine sets the following fields depending on the input
! parameters present:
!  (1) topography
!  (2) capacity (set if coordinate_system == 2)
!  (3) length ratio of edge (set if coordinate_system == 2)
!  (frition_index). Location of manning's N (variable_friction)
!  (wind_index, wind_index + 1).  Location of x and y wind speeds (wind_forcing)
!  (pressure_index).  Location of pressure field (pressure_forcing)
!
! There are a couple of additional things of note:
!
!  - This routine is called anytime a grid is created during re-gridding to fill
!    in the aux arrays.
!  - Built-in to this routine is an ability to copy aux values from grids from
!    the same level.  If you do not want this behavior you will need to modify
!    this routine.
!  - If you are using periodic BCs then you must handle this here.  Look to the
!    setting of topography for an example of how this may need to be done.
!  - If your aux arrays are time dependent and set somewhere else it may be
!    prudent to set a default value here as is done with the wind and pressure.
!
subroutine setaux(mbc,mx,my,xlow,ylow,dx,dy,maux,aux)

    use amr_module, only: mcapa, xupper, yupper, xlower, ylower, NEEDS_TO_BE_SET
    use amr_module, only: xperdom, yperdom

    use geoclaw_module, only: coordinate_system, earth_radius, deg2rad
    use geoclaw_module, only: sea_level, ambient_pressure

    use storm_module, only: wind_forcing, pressure_forcing
    use storm_module, only: wind_index, pressure_index, set_storm_fields

    use friction_module, only: variable_friction, friction_index
    use friction_module, only: set_friction_field

    use topo_module, only: mtopofiles

    use adjoint_module, only : adjoint_flagging,innerprod_index

    use auxinit_module, only: mauxinitfiles,iauxinit,ylowauxinit
    use auxinit_module, only: yhiauxinit,xlowauxinit,xhiauxinit
    use auxinit_module, only: i0auxinit,auxinitwork,mxauxinit
    use auxinit_module, only: dyauxinit,dxauxinit,myauxinit
    use auxinit_module, only: mauxinit

    use digclaw_module, only: i_dig,i_phi,i_theta,phi_bed,theta_input

    implicit none

    ! Arguments
    integer, intent(in) :: mbc,mx,my,maux
    real(kind=8), intent(in) :: xlow,ylow,dx,dy
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: ii,jj,m, iint,jint
    real(kind=8) :: x,y,xm,ym,xp,yp,topo_integral
    real(kind=8) :: xper,yper,xperm,yperm,xperp,yperp
    character(len=*), parameter :: aux_format = "(2i4,4d15.3)"
    integer :: skipcount,iaux,ilo,jlo

    ! Declare D-Claw variables...
    integer :: mf,istart,iend,jstart,jend,i,j
    real(kind=8) :: xhigher,yhigher,xintlow,xinthi,yintlow,yinthi
    real(kind=8) :: xim,xip,yjm,yjp,xipc,ximc,xc,yjpc,yjmc,yc,daux
    real(kind=8) :: b_x,b_y,gradang,kappa,phi_tread
    logical :: use_phi_bed,use_theta_input,friction_correction
    ! Topography integral function
    real(kind=8) :: topointegral

    ! Lat-Long coordinate system in use, check input variables
    if (coordinate_system == 2) then
        if (mcapa /= 2 .or. maux < 3) then
            print *,'ERROR in setaux:  for coordinate_system==2'
            print *,'     need mcapa == 2 and maux >= 3'
            print *,'     have mcapa = ',mcapa,'  maux = ',maux
            stop
        endif
    endif

    ! Compute integer indices based off same corner of domain to reduce round
    ! off discrepancies
    ilo = floor((xlow - xlower + .05d0*dx)/dx)
    jlo = floor((ylow - ylower + .05d0*dy)/dy)

    ! Set geometry values
    if (coordinate_system == 2) then
        do jj = 1 - mbc, my + mbc
            do ii = 1 - mbc, mx + mbc
                ym = ylower + (jlo+jj-1.d0) * dy
                yp = ylower + (jlo+jj) * dy
                aux(2,ii,jj) = deg2rad * earth_radius**2                      &
                                * (sin(yp * deg2rad) - sin(ym * deg2rad)) / dy
                aux(3,ii,jj) = ym * deg2rad
            end do
        end do
    end if


    ! Set bathymetry based on input files
    if (mtopofiles > 0) then

        do jj=1-mbc,my+mbc
            ym = ylower + (jlo+jj-1.d0) * dy
            yp = ylower + (jlo+jj) * dy
            y = 0.5d0*(ym+yp)


            do ii=1-mbc,mx+mbc
                xm = xlower + (ilo+ii-1.d0) * dx
                xp = xlower + (ilo+ii) * dx
                x = 0.5d0*(xm+xp)

                ! Parameter NEEDS_TO_BE_SET initialized in amr_module.f90
                ! saves time by otherwise copying instead of reinitializing
                if (aux(1,ii,jj) .ne. NEEDS_TO_BE_SET) then
                    cycle
                endif

                topo_integral = 0.d0
                if ((y>yupper).or.(y<ylower).or.(x>xupper).or.(x<xlower)) then
                    if (.not.(xperdom .or. yperdom)) then
                        ! Skip setting as this cell sticks out of the physical
                        ! domain and we are not setting periodic BCs
                        cycle

                    else
                        ! We evaluate the periodic BC topography by computing
                        ! the appropriate periodic grid cell coordinates and
                        ! again calling cellgridintegrate
                        call wrap_coords(x,y,xperm,xper,xperp,yperm,yper, &
                                         yperp,dx,dy)
                        call cellgridintegrate(topo_integral,  &
                                               xperm,xperp,yperm,yperp)
                    endif
                else
                    ! Cell does not extend outside of physical domain
                    call cellgridintegrate(topo_integral,xm,xp,ym,yp)
                endif

                ! Correct for geometry
                if (coordinate_system == 2) then
                    aux(1,ii,jj) = topo_integral / (dx * dy * aux(2,ii,jj))
                else
                    aux(1,ii,jj) = topo_integral / (dx * dy)
                endif
            enddo
        enddo
    else
        print *, "ERROR:  There is no way to set bathymetry!  Either "
        print *, "        provide topography files or request topography "
        print*,  "        defined by a function."
        stop
    end if

    ! Copy topo to ghost cells if outside physical domain and not periodic
    if (.not. yperdom) then
        do jj=1-mbc,my+mbc
            y = ylower + (jlo+jj-.5d0) * dy
            if ((y < ylower) .or. (y>yupper)) then
                do ii=1-mbc,mx+mbc
                    x = xlower + (ilo+ii-.5d0) * dx
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                enddo
            endif
        enddo
    endif
    if (.not. xperdom) then
        do ii=1-mbc,mx+mbc
            x =  xlower + (ilo+ii-.5d0) * dx
            if ((x < xlower) .or. (x > xupper)) then
                do jj=1-mbc,my+mbc
                    y = ylower + (jlo+jj-.5d0) * dy
                    iint = ii + max(0, ceiling((xlower-x)/dx)) &
                         - max(0, ceiling((x-xupper)/dx))
                    jint = jj + max(0, ceiling((ylower-y)/dy)) &
                         - max(0, ceiling((y-yupper)/dy))
                    aux(1,ii,jj) = aux(1,iint,jint)
                enddo
            endif
        enddo
    endif


    !=================================
    !DIG:  Adapted from D-Claw...

    ! ------- zero aux variables that will be set by files
    do mf = 1,mauxinitfiles
       aux(iauxinit(mf),:,:) = 0.d0
    enddo

    !--------zero all d-claw aux variables----
    do mf= i_dig,maux
       aux(mf,:,:) = 0.0
    enddo

!     --------------integrate auxinit files if they exist---------------
      xhigher = xlow + (mx-0.5d0)*dx
      yhigher = ylow + (my-0.5d0)*dy

      do mf =1,mauxinitfiles
         if ((xlow.le.xhiauxinit(mf).and. &
                  xhigher.ge.xlowauxinit(mf)).and. &
                  (ylow.le.yhiauxinit(mf).and. &
                  yhigher.ge.ylowauxinit(mf))) then

            xintlow = dmax1(xlow,xlowauxinit(mf))
            xinthi  = dmin1(xhigher,xhiauxinit(mf))
            istart  = max(1-mbc,int(0.5 + (xintlow-xlow)/dx)-mbc)
            iend    = min(mx+mbc,int(1.0 + (xinthi-xlow)/dx)+mbc)

            yintlow = dmax1(ylow,ylowauxinit(mf))
            yinthi  = dmin1(yhigher,yhiauxinit(mf))
            jstart  = max(1-mbc,int(0.5 + (yintlow-ylow)/dy)-mbc)
            jend    = min(my+mbc,int(1.0 + (yinthi-ylow)/dy)+mbc)

            do i=istart,iend
               x = xlow + (i-0.5d0)*dx
               xim = x - 0.5d0*dx
               xip = x + 0.5d0*dx
               do j=jstart,jend
                  y = ylow + (j-0.5d0)*dy
                  yjm = y - 0.5d0*dy
                  yjp = y + 0.5d0*dy

                  if (xip.gt.xlowauxinit(mf) &
                       .and.xim.lt.xhiauxinit(mf) &
                       .and.yjp.gt.ylowauxinit(mf) &
                       .and.yjm.lt.yhiauxinit(mf)) then

                     xipc=dmin1(xip,xhiauxinit(mf))
                     ximc=dmax1(xim,xlowauxinit(mf))
                     xc=0.5d0*(xipc+ximc)

                     yjpc=dmin1(yjp,yhiauxinit(mf))
                     yjmc=dmax1(yjm,ylowauxinit(mf))
                     yc=0.5d0*(yjmc+yjpc)

                     daux =topointegral(ximc,xipc,yjmc,yjpc, &
                             xlowauxinit(mf),ylowauxinit(mf), &
                             dxauxinit(mf),dyauxinit(mf), &
                             mxauxinit(mf),myauxinit(mf), &
                             auxinitwork(i0auxinit(mf):i0auxinit(mf) &
                             +mauxinit(mf)-1), 1)

                     ! Correct for geometry
                     if (coordinate_system == 2) then
                        daux=daux/((xipc-ximc)*(yjpc-yjmc)*aux(2,i,j))
                     else
                        daux=daux/((xipc-ximc)*(yjpc-yjmc))
                     endif

                     aux(iauxinit(mf),i,j) = aux(iauxinit(mf),i,j)+daux
                  endif
               enddo
            enddo
         endif
      enddo

      use_phi_bed = .true.
      use_theta_input = .true.
      do mf = 1,mauxinitfiles
         if (iauxinit(mf).eq.i_phi) use_phi_bed = .false.
         if (iauxinit(mf).eq.i_theta) use_theta_input = .false.
      enddo

      if (use_phi_bed) then
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               aux(i_phi,i,j) = phi_bed
            enddo
         enddo
      endif

      if (use_theta_input) then
         do j=1-mbc,my+mbc
            do i=1-mbc,mx+mbc
               aux(i_theta,i,j) = theta_input
            enddo
         enddo
      endif

      friction_correction = .false.

      ! Iverson, R. M., & George, D. L. (2019). Basal stress
      ! equations for granular debris masses on smooth or
      ! discretized slopes. Journal of Geophysical Research:
      ! Earth Surface, 124, 1464–1484. https://doi.org/10.1029/2018JF004802

      if (friction_correction) then
        do j=1-mbc+1,my+mbc-1
            do i=1-mbc+1,mx+mbc-1
              b_x = (aux(1,i+1,j)-aux(1,i-1,j))/(2.d0*dx)
              b_y = (aux(1,i,j+1)-aux(1,i,j-1))/(2.d0*dy)
              gradang = atan(sqrt(b_x**2+b_y**2))
              kappa=1.d0
              phi_tread = tan(aux(i_phi,i,j)-gradang) &
                    +2.d0*kappa*tan(gradang)
              aux(i_phi,i,j) = phi_tread
            enddo
         enddo
      endif


contains

    ! Provide wrapper function for providing periodic coordinates
    subroutine wrap_coords(x,y,xperm,xper,xperp,yperm,yper,yperp,dx,dy)

        use amr_module, only: xperdom, yperdom, xupper, yupper, xlower, ylower

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: x,y,dx,dy
        real(kind=8), intent(out) :: xperm,xper,xperp,yperm,yper,yperp
        real(kind=8) :: xdomain, ydomain

        ! need size of domain for wrapping
        xdomain = xupper - xlower
        ydomain = yupper - ylower
        xper = x
        yper = y

        ! test for wrapping of coordinates
        if (x > xupper .and. xperdom) xper = xper - xdomain
        if (x < xlower .and. xperdom) xper = xper + xdomain
        if (y > yupper .and. yperdom) yper = yper - ydomain
        if (y < ylower .and. yperdom) yper = yper + ydomain

        ! adjust rest of variables
        xperm = xper - 0.5d0 * dx
        xperp = xper + 0.5d0 * dx
        yperm = yper - 0.5d0 * dy
        yperp = yper + 0.5d0 * dy

    end subroutine wrap_coords

end subroutine setaux

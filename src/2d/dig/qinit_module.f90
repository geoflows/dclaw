!! D-Claw specific core file
!! This file is a modified version of
!! clawpack/geoclaw/src/2d/shallow/qinit_module.f90
!!

! ============================================================================
!  Module for initialization of the q arrays that might come from files
! ============================================================================


module qinit_module

    ! Updated to include D-Claw files as e.g. set_qinit_dig
    ! which reads data from setqinit_dclaw.data

    use amr_module, only: rinfinity
    use digclaw_module, only: i_h,i_hu,i_hv,i_hm,i_pb,i_hchi

    implicit none
    save

    logical :: module_setup = .false.

    ! Type of q initialization
    integer, public :: qinit_type

    ! Work array
    real(kind=8), private, allocatable :: qinit(:)

    ! Geometry
    real(kind=8) :: x_low_qinit
    real(kind=8) :: y_low_qinit
    real(kind=8) :: t_low_qinit
    real(kind=8) :: x_hi_qinit
    real(kind=8) :: y_hi_qinit
    real(kind=8) :: t_hi_qinit
    real(kind=8) :: dx_qinit
    real(kind=8) :: dy_qinit

    integer, private :: mx_qinit
    integer, private :: my_qinit

    integer, parameter ::  QINIT_PARM_UNIT = 106

    ! for initializing using force_dry to indicate dry regions below sealevel:

    integer :: mx_fdry, my_fdry
    real(kind=8) :: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry, dx_fdry, dy_fdry
    integer(kind=1), allocatable :: force_dry(:,:)
    logical :: use_force_dry
    real(kind=8) :: tend_force_dry  ! always use mask up to this time

    logical :: variable_eta_init

    ! to initialize using different initial eta values in different regions:
    integer :: etain_mx, etain_my
    real(kind=8) :: etain_dx, etain_dy
    real(kind=8), allocatable :: etain_x(:), etain_y(:), etain_eta(:,:)

      !DIG parameters

      ! Work array
      double precision, allocatable :: qinitwork(:)

      ! Topography file data
      character*150, allocatable :: qinitfname(:)
      integer :: mqinitfiles,mqinitsize
      double precision, allocatable :: xlowqinit(:), ylowqinit(:)
      double precision, allocatable :: xhiqinit(:), yhiqinit(:)
      double precision, allocatable :: dxqinit(:), dyqinit(:)

      integer, allocatable ::  mxqinit(:), myqinit(:)
      integer, allocatable :: i0qinit(:), mqinit(:)
      integer, allocatable :: iqinit(:), qinitftype(:)

contains

    subroutine set_qinit(fname)

        implicit none

        ! Subroutine arguments
        character(len=*), optional, intent(in) :: fname

        ! File handling
        integer, parameter :: unit = 7
        character(len=150) :: qinit_fname
        character(len=150) :: fname_force_dry

        integer :: num_force_dry

        if (.not.module_setup) then
            write(QINIT_PARM_UNIT,*) ' '
            write(QINIT_PARM_UNIT,*) '--------------------------------------------'
            write(QINIT_PARM_UNIT,*) 'SETQINIT:'
            write(QINIT_PARM_UNIT,*) '-------------'

            ! Open the data file
            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,"qinit.data")
            endif

            read(unit,"(i1)") qinit_type
            if (qinit_type == 0) then
                ! No perturbation specified
                write(QINIT_PARM_UNIT,*)  '  qinit_type = 0, no perturbation'
                print *,'  qinit_type = 0, no perturbation'
            else
                read(unit,*) qinit_fname
                write(QINIT_PARM_UNIT,*)  qinit_fname

                call read_qinit(qinit_fname)
            endif


            ! If variable_eta_init then function set_eta_init is called
            ! to set initial eta when interpolating onto newly refined patches
            read(unit,*) variable_eta_init


            read(unit,*) num_force_dry
            use_force_dry = (num_force_dry > 0)

            if (num_force_dry > 1) then
                write(6,*) '*** num_force_dry > 1 not yet implemented'
                stop
                endif

            if (use_force_dry) then
                read(unit,*) fname_force_dry
                read(unit,*) tend_force_dry
                call read_force_dry(trim(fname_force_dry))
                endif

            module_setup = .true.
        end if

    end subroutine set_qinit


    subroutine add_perturbation(meqn,mbc,mx,my,xlow_patch,ylow_patch,dx,dy,q,maux,aux)

        use geoclaw_module, only: sea_level, coordinate_system
        use amr_module, only: mcapa

        implicit none

        ! Subroutine arguments
        integer, intent(in) :: meqn,mbc,mx,my,maux
        real(kind=8), intent(in) :: xlow_patch,ylow_patch,dx,dy
        real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local
        integer :: i,j
        real(kind=8) :: ximc,xim,x,xip,xipc,yjmc,yjm,y,yjp,yjpc,dq

        ! Topography integral function
        real(kind=8) :: topointegral

        if (qinit_type > 0) then
            do j=1-mbc,my+mbc
                y = ylow_patch + (j-0.5d0)*dy
                yjm = y - 0.5d0*dy
                yjp = y + 0.5d0*dy

                do i=1-mbc,mx+mbc
                    x = xlow_patch + (i-0.5d0)*dx
                    xim = x - 0.5d0*dx
                    xip = x + 0.5d0*dx

                    ! Check to see if we are in the qinit region at this grid point
                    if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                        (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then

                        xipc=min(xip,x_hi_qinit)
                        ximc=max(xim,x_low_qinit)

                        yjpc=min(yjp,y_hi_qinit)
                        yjmc=max(yjm,y_low_qinit)

                        dq = topointegral(ximc,xipc,yjmc,yjpc,x_low_qinit, &
                                          y_low_qinit,dx_qinit,dy_qinit,mx_qinit, &
                                          my_qinit,qinit,1)
                        if (coordinate_system == 2) then
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(mcapa,i,j))
                        else
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc))
                        endif

                        if (qinit_type < 4) then
                            if (aux(1,i,j) <= sea_level) then
                                q(qinit_type,i,j) = q(qinit_type,i,j) + dq
                            endif
                        else if (qinit_type == 4) then
                            q(i_h,i,j) = max(dq-aux(1,i,j),0.d0)
                        endif
                    endif
                enddo
            enddo
        endif

    end subroutine add_perturbation


    ! currently only supports one file type:
    ! x,y,z values, one per line in standard order from NW corner to SE
    ! z is perturbation from standard depth h,hu,hv set in qinit_geo,
    ! if iqinit = 1,2, or 3 respectively.
    ! if iqinit = 4, the z column corresponds to the definition of the
    ! surface elevation eta. The depth is then set as q(i,j,1)=max(eta-b,0)
    subroutine read_qinit(fname)

        implicit none

        ! Subroutine arguments
        character(len=150) :: fname

        ! Data file opening
        integer, parameter :: unit = 19
        integer :: i,num_points,status
        double precision :: x,y

        print *,'  '
        print *,'Reading qinit data from file  ', fname
        print *,'  '

        write(QINIT_PARM_UNIT,*) '  '
        write(QINIT_PARM_UNIT,*) 'Reading qinit data from'
        write(QINIT_PARM_UNIT,*) fname
        write(QINIT_PARM_UNIT,*) '  '

        open(unit=unit, file=fname, iostat=status, status="unknown", &
             form='formatted',action="read")
        if ( status /= 0 ) then
            print *,"Error opening file", fname
            stop
        endif

        ! Initialize counters
        num_points = 0
        mx_qinit = 0

        ! Read in first values, determines x_low and y_hi
        read(unit,*) x_low_qinit,y_hi_qinit
        num_points = num_points + 1
        mx_qinit = mx_qinit + 1

        ! Sweep through first row figuring out mx
        y = y_hi_qinit
        do while (y_hi_qinit == y)
            read(unit,*) x,y
            num_points = num_points + 1
            mx_qinit = mx_qinit + 1
        enddo
        ! We over count by one in the above loop
        mx_qinit = mx_qinit - 1

        ! Continue to count the rest of the lines
        do
            read(unit,*,iostat=status) x,y
            if (status /= 0) exit
            num_points = num_points + 1
        enddo
        if (status > 0) then
            print *,"ERROR:  Error reading qinit file ",fname
            stop
        endif

        ! Extract rest of geometry
        x_hi_qinit = x
        y_low_qinit = y
        my_qinit = num_points / mx_qinit
        dx_qinit = (x_hi_qinit - x_low_qinit) / (mx_qinit-1)
        dy_qinit = (y_hi_qinit - y_low_qinit) / (my_qinit-1)

        rewind(unit)
        allocate(qinit(num_points))

        ! Read and store the data this time
        do i=1,num_points
            read(unit,*) x,y,qinit(i)
        enddo
        close(unit)

    end subroutine read_qinit

    subroutine read_force_dry(fname)

        use utility_module, only: parse_values
        character(len=*), intent(in) :: fname
        integer :: iunit,i,j,n
        real(kind=8) :: values(16), nodata_value
        character(len=80) :: str

        iunit = 8

        open(unit=iunit,file=fname,status='old',form='formatted')
        !read(iunit,*) tend_force_dry
        !write(6,*) 'tend_force_dry = ',tend_force_dry
        read(iunit,*) mx_fdry
        read(iunit,*) my_fdry
        read(iunit,*) xlow_fdry
        read(iunit,*) ylow_fdry

        read(iunit,'(a)') str
        call parse_values(str, n, values)
        dx_fdry = values(1)
        if (n == 2) then
            dy_fdry = values(2)
          else
            dy_fdry = dx_fdry
          endif

        read(iunit,*) nodata_value
        allocate(force_dry(mx_fdry,my_fdry))

        xhi_fdry = xlow_fdry + mx_fdry*dx_fdry
        yhi_fdry = ylow_fdry + my_fdry*dy_fdry
        write(6,*) '+++ xlow_fdry, xhi_fdry: ',xlow_fdry, xhi_fdry
        write(6,*) '+++ ylow_fdry, yhi_fdry: ',ylow_fdry, yhi_fdry

        do j=1,my_fdry
            read(iunit, *) (force_dry(i,j), i=1,mx_fdry)
            enddo

        close(iunit)
        return
    end subroutine read_force_dry


    subroutine read_eta_init(file_name)
        ! To read in file specifying different eta value in at different
        ! locations, then used in qinit function.
        ! Uses etain module variables.

        implicit none

        ! Input arguments
        character(len=*), intent(in), optional :: file_name

        ! local
        integer, parameter :: iunit = 7
        integer :: i,j
        real(kind=8) :: nodata_value, xllower, yllower

        if (present(file_name)) then
            open(unit=iunit, file=file_name, status='unknown',&
                      form='formatted')
        else
            open(unit=iunit, file='eta_init.data', status='unknown',&
                      form='formatted')
        endif

        read(iunit,*) etain_mx
        !write(6,*) '+++ etain_mx = ',etain_mx
        read(iunit,*) etain_my
        !write(6,*) '+++ etain_my = ',etain_my
        read(iunit,*) xllower
        read(iunit,*) yllower
        read(iunit,*) etain_dx
        etain_dy = etain_dx
        !read(iunit,*) etain_dy
        read(iunit,*) nodata_value

        allocate(etain_x(etain_mx), etain_y(etain_my))
        allocate(etain_eta(etain_mx, etain_my))

        do i=1,etain_mx
            etain_x(i) = xllower + etain_dx*(i-1)
            enddo

        do j=1,etain_my
            etain_y(j) = yllower + etain_dy*(etain_my-j+1)
            read(iunit,*) (etain_eta(i,j),i=1,etain_mx)
            enddo


        close(unit=iunit)
    end subroutine read_eta_init

    ! ========================================================================
    ! Read qinit files as specified in qinit_dclaw.data
    !
    ! Each file has a type stored in qinitftype(i).
    !   qinittype = 1:  standard GIS format: 3 columns: x,y,z(m)
    !   qinittype = 2:  Header as in DEM file, height(m) one value per line
    !   qinittype = 3:  Header as in DEM file, height(m) one row per line
    ! For other formats modify readqinit routine.
    !
    ! advancing northwest to northeast then from north to south. Values should
    ! be uniformly spaced.
    !
    ! Associated with each file is a initialization type, iqinit(file):
    !     as follows:
    !     m: perturbation to q(m) for m=equation_number
    !     m+1:     file defines eta, and h =q(i_h,i,j)= max(eta-b,0)
    ! ========================================================================

   subroutine set_qinit_dig(fname)

      implicit none

      ! Input arguments
      character*25, intent(in), optional :: fname

      ! Locals
      integer, parameter :: iunit = 7
      integer :: i,j,iqinitfile
      character*25 :: file_name
      logical :: found_file


      ! Open and begin parameter file output
      write(QINIT_PARM_UNIT,*) ' '
      write(QINIT_PARM_UNIT,*) '--------------------------------------------'
      write(QINIT_PARM_UNIT,*) 'SETQINIT:'
      write(QINIT_PARM_UNIT,*) '---------'

      if (present(fname)) then
         file_name = fname
      else
         file_name  = 'qinit_dclaw.data'
      endif
      inquire(file=file_name,exist=found_file)
      if (.not. found_file) then
         print *, 'You must provide a file ', file_name
         stop
      endif

      call opendatafile(iunit, file_name)

      read(iunit,*) mqinitfiles

      if (mqinitfiles==0) then
         write(QINIT_PARM_UNIT,*) '   mqinitfiles = 0'
         write(QINIT_PARM_UNIT,*) '   no initial perturbation = 0'
         write(QINIT_PARM_UNIT,*) '   h will be set max(0-b,0)   '
         return
      endif

      write(QINIT_PARM_UNIT,*) '   mqinitfiles = ',mqinitfiles

      ! Read and allocate data parameters for each file
      allocate(mxqinit(mqinitfiles),myqinit(mqinitfiles))
      allocate(xlowqinit(mqinitfiles),ylowqinit(mqinitfiles))
      allocate(xhiqinit(mqinitfiles),yhiqinit(mqinitfiles))
      allocate(dxqinit(mqinitfiles),dyqinit(mqinitfiles))
      allocate(qinitfname(mqinitfiles),qinitftype(mqinitfiles))
      allocate(iqinit(mqinitfiles))
      allocate(i0qinit(mqinitfiles),mqinit(mqinitfiles))

      do i=1,mqinitfiles
         read(iunit,*) qinitfname(i)
         read(iunit,*) qinitftype(i),iqinit(i)

         write(QINIT_PARM_UNIT,*) '   '
         write(QINIT_PARM_UNIT,*) '   ',qinitfname(i)
         write(QINIT_PARM_UNIT,*) '  qinitftype = ', qinitftype(i)
         write(QINIT_PARM_UNIT,*) '  iqinit = ', iqinit(i)

         call read_qinit_dig_header(qinitfname(i),qinitftype(i),mxqinit(i), &
                myqinit(i),xlowqinit(i),ylowqinit(i),xhiqinit(i),yhiqinit(i), &
                dxqinit(i),dyqinit(i))
            mqinit(i) = mxqinit(i)*myqinit(i)
      enddo

      ! Indexing into work array
      i0qinit(1)=1
      if (mqinitfiles > 1) then
         do i=2,mqinitfiles
            i0qinit(i)=i0qinit(i-1) + mqinit(i-1)
         enddo
      endif

      ! Read and allocate space in work array for each file
      mqinitsize = sum(mqinit)
      allocate(qinitwork(mqinitsize))

      do i=1,mqinitfiles
            call read_qinit_dig(mxqinit(i),myqinit(i),qinitftype(i),qinitfname(i), &
                qinitwork(i0qinit(i):i0qinit(i)+mqinit(i)-1))
      enddo

   end subroutine set_qinit_dig


    ! ========================================================================
    !
    !  Read qinit_dig file.
    !
    ! ========================================================================
    subroutine read_qinit_dig(mx,my,filetype,fname,qinit)

        use utility_module, only: parse_values, to_lower

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,filetype
        character*150, intent(in) :: fname
        double precision, intent(inout) :: qinit(1:mx*my)

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        double precision, parameter :: qinit_missing = -150.d0
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,qinit_start,n
        double precision :: no_data_value,x,y,z
        real(kind=8) :: values(16)
        character(len=80) :: str

        print *, ' '
        print *, 'Reading qinit file  ', fname

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(filetype))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                i = 0
                status = 0
                do while (status == 0)
                    i = i + 1
                    read(iunit,fmt=*,iostat=status) x,y,qinit(i)
                enddo

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if filetype=2 or
            ! mx values per line if filetype=3
            ! ================================================================
            case(2:3)
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                no_data_value = values(1)

                ! Read in data
                missing = 0
                select case(abs(filetype))
                    case(2)
                        do i=1,mx*my
                            read(iunit,*) qinit(i)
                            if (qinit(i) == no_data_value) then
                                missing = missing + 1
                                qinit(i) = qinit_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (qinit((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (qinit((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    qinit((j-1)*mx + i) = qinit_missing
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print *, '   WARNING...some missing data values this file'
                    print *, '       ',missing,' missing data values'
                    print *, '              (see fort.missing)'
                    print *, '   These values have arbitrarily been set to ',&
                        qinit_missing
                endif
        end select

        close(unit=iunit)

   end subroutine read_qinit_dig


    ! ========================================================================
    ! subroutine read_qinit_dig_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read qinit file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - qinit_type - (int) Type of file format (1 < qinit_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_qinit_dig_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)

        use utility_module, only: parse_values, to_lower

        implicit none

        ! Input and Output
        character(len=150), intent(in) :: fname
        integer, intent(in) :: qinit_type
        integer, intent(out) :: mx,my
        real(kind=8), intent(out) :: xll,yll,xhi,yhi,dx,dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: qinit_size, status, n
        real(kind=8) :: x,y,z,nodata_value
        logical :: found_file
        real(kind=8) :: values(16)
        character(len=80) :: str
        logical :: verbose
        logical :: xll_registered, yll_registered

        verbose = .false.

        inquire(file=fname,exist=found_file)
        if (.not. found_file) then
            print *, 'Missing qinit file:'
            print *, '   ', fname
            stop
        endif

        select case(abs(qinit_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                open(unit=iunit, file=fname, status='unknown',form='formatted')

                ! Initial size variables
                qinit_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                qinit_size = qinit_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    qinit_size = qinit_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) x,y,z
                    qinit_size = qinit_size + 1
                enddo
                if (status > 0) then
                    print *, "ERROR:  Error reading header of qinit file ",fname
                    stop
                endif

                ! Calculate remaining values
                my = qinit_size / mx
                xhi = x
                yll = y
                dx = (xhi-xll) / (mx-1)
                dy = (yhi-yll) / (my-1)

            ! ASCII file with header followed by z data
            case(2:3)
                open(unit=iunit, file=fname, status='unknown',form='formatted')
                read(iunit,'(a)') str
                call parse_values(str, n, values)
                mx = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                my = nint(values(1))

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                xll = values(1)
                str = to_lower(str)  ! convert to lower case
                xll_registered = (index(str, 'xllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                yll = values(1)
                str = to_lower(str)  ! convert to lower case
                yll_registered = (index(str, 'yllcorner') > 0)

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                dx = values(1)
                if (n == 2) then
                    dy = values(2)
                  else
                    dy = dx
                  endif

                read(iunit,'(a)') str
                call parse_values(str, n, values)
                nodata_value = values(1)

                if (xll_registered) then
                    xll = xll + 0.5d0*dx
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting xllcorner by 0.5*dx to cell center'
                    endif

                if (yll_registered) then
                    yll = yll + 0.5d0*dy
                    write(6,*) '*** in file: ',trim(fname)
                    write(6,*) '    Shifting yllcorner by 0.5*dy to cell center'
                    endif

                xhi = xll + (mx-1)*dx
                yhi = yll + (my-1)*dy

            case default
                print *, 'ERROR:  Unrecognized qinit_type'
                print *, '    qinit_file_type = ',qinit_type
                print *, '  for qinit file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)

        write(QINIT_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(QINIT_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(QINIT_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_qinit_dig_header

end module qinit_module

! Module containing refinement flagging criteria
module refinement_module

    use geoclaw_module, only: GEO_PARM_UNIT

    implicit none
    save

    logical, private :: module_setup = .false.

    ! ========================================================================
    !  Refinement Criteria
    ! ========================================================================
    real(kind=8) :: wave_tolerance
    real(kind=8), allocatable :: speed_tolerance(:)
    logical :: varRefTime = .FALSE. ! Choose dt refinement automatically


    ! ========================================================================
    !  Flow grades flagging support adapted from D-Claw
    ! ========================================================================
    real(kind=8), allocatable :: flowgradevalue(:)

    integer, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, allocatable :: iflowgrademinlevel(:)
    integer :: mflowgrades
    logical :: keep_fine

contains

    ! =========================================================================
    !  Reads in the refinement control parameters
    ! =========================================================================
    subroutine set_refinement(file_name)

        use amr_module, only: mxnest
        use utility_module, only: get_value_count

        implicit none

        ! Arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: unit = 127
        integer :: i
        character(len=128) :: line

        if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Refinement Control Parameters:'
            write(GEO_PARM_UNIT,*) '------------------------------'

            if (present(file_name)) then
                call opendatafile(unit, file_name)
            else
                call opendatafile(unit, 'refinement.data')
            endif

            ! Basic criteria
            read(unit,*) wave_tolerance
            read(unit,'(a)') line
            allocate(speed_tolerance(get_value_count(line)))
            read(line,*) speed_tolerance
            read(unit,*)
            read(unit,*) varRefTime
            close(unit)

            ! Write out data to parameter file
            write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
            write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
            write(GEO_PARM_UNIT,*) '   Variable dt Refinement Ratios:',varRefTime
            write(GEO_PARM_UNIT,*) ''

            module_setup = .true.
        end if

    end subroutine set_refinement


    subroutine set_flow_grades(fname)

        implicit none

        ! Input arguments
        character*25, intent(in), optional :: fname

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i
        character*25 :: file_name
        logical :: found_file


        write(GEO_PARM_UNIT,*) ' '
        write(GEO_PARM_UNIT,*) '--------------------------------------------'
        write(GEO_PARM_UNIT,*) 'SET FLOW GRADES:'
        write(GEO_PARM_UNIT,*) '------------'

        ! Read user parameters from setflowgrades.data

        if (present(fname)) then
            file_name = fname
        else
            file_name = 'setflowgrades.data'
        endif
        inquire(file=file_name,exist=found_file)
        if (.not. found_file) then
            print *, 'You must provide a file ', file_name
            print *, 'Or comment out call set_flow_grades in setprob'
            stop
        endif

        call opendatafile(iunit, file_name)

        read(iunit,*) mflowgrades

        if (mflowgrades == 0) then
            write(GEO_PARM_UNIT,*) '  No flow grades specified'
            return
        endif

        ! Allocate arrays
        allocate(flowgradevalue(mflowgrades),iflowgradevariable(mflowgrades))
        allocate(iflowgradetype(mflowgrades),iflowgrademinlevel(mflowgrades))

        do i=1,mflowgrades
            read(iunit,*) flowgradevalue(i),iflowgradevariable(i), &
                iflowgradetype(i),iflowgrademinlevel(i)
        enddo

        read(iunit,*) keep_fine
        write(*,*) 'keep_fine', keep_fine

        close(iunit)

        write(GEO_PARM_UNIT,*) '   mflowgrades:',  mflowgrades

        do i=1,mflowgrades
            write(GEO_PARM_UNIT,"(d12.3,3i4)") flowgradevalue(i), &
                iflowgradevariable(i),iflowgradetype(i),iflowgrademinlevel(i)

        enddo

    end subroutine set_flow_grades


end module refinement_module

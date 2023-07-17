program compute_rms_density
  
  ! compute_rms_density.x [qc file] [fgh file] [npoints]

  implicit none

  double precision, dimension(:), allocatable :: qcden, fghden
  integer :: npts
  
  character(255)     :: qcden_file, fghden_file, cmdstr
  integer, parameter :: qcden_unit = 10, fghden_unit = 11
  
  double precision :: rms_error

  integer :: ios
  
  ! Process command line arguments
  call get_command_argument(number=1, value=qcden_file, status=ios)
  if (ios .ne. 0) then
      print *, "WARNING: No file name provided. Assuming fhf_dens.data"
      write(qcden_file, "(A)") "fhf_dens.data"
  end if
  call get_command_argument(number=2, value=fghden_file, status=ios)
  if (ios .ne. 0) then
      print *, "WARNING: No file name provided. Assuming fhf_dens.fgh.data"
      write(fghden_file, "(A)") "fhf_dens.fgh.data"
  end if
  call get_command_argument(number=3, value=cmdstr, status=ios)
  if (ios .ne. 0) then
      print *, "WARNING: No point number provided. Assuming 64 points."
      npts = 64
  else
      read(cmdstr, "(I3)") npts
  end if
  
  ! Get density information
  allocate(qcden(npts))
  allocate(fghden(npts))
  call read_density_file(qcden_file, qcden_unit, qcden, npts)
  call read_density_file(fghden_file, fghden_unit, fghden, npts)
  
  call compute_rms_error(qcden, fghden, npts, rms_error)

  print "(F20.8)", rms_error

contains
  
  ! Compute the RMS error between two densities
  ! Assumed den2 is reference
  subroutine compute_rms_error(den1, den2, npts, rmse)
     
    implicit none
    
    integer, intent(in) :: npts
    double precision, dimension(npts), intent(in) :: den1, den2
    
    double precision, intent(out) :: rmse

    double precision, dimension(npts) :: diff
    double precision, dimension(npts) :: wts
    
    integer :: i

    ! Create weight vector
    call create_weight_vector(wts, den2, npts)

    diff = den1 - den2
    rmse = 0d0
    do i = 1, npts
        rmse = rmse + wts(i) * diff(i) * diff(i)
    end do
    rmse = rmse / dble(npts)
    rmse = sqrt(rmse)

  end subroutine compute_rms_error

  ! Create weight vector based on value
  subroutine create_weight_vector(wts, ref, npts)
    
    implicit none
    
    integer, intent(in) :: npts
    double precision, dimension(npts), intent(in) :: ref
    double precision, dimension(npts), intent(inout) :: wts
    
    ! Fit all points at wt=1.0
    wts = 1.0

  end subroutine create_weight_vector

  ! Read in density file
  subroutine read_density_file(flname, flunit, dens, npts)
    
    implicit none
    
    character(255), intent(in) :: flname
    integer, intent(in) :: flunit
    integer, intent(in) :: npts
    double precision, dimension(npts), intent(inout) :: dens
    integer :: i

    open(unit=flunit, file=trim(adjustl(flname)), action="read", &
        status="unknown", position="rewind", iostat=ios)
    if (ios .ne. 0) stop "*** Can not open density file. ***"
    do i = 1, npts
        read(flunit,fmt="(20X,F20.8)") dens(i)
    end do
    close(unit=flunit)
    
  end subroutine read_density_file

end program compute_rms_density

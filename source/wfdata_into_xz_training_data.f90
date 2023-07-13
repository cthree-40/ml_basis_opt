program wfdata_into_xz_training_data

  implicit none

  double precision, dimension(:,:), allocatable :: xdens
  double precision, dimension(:,:), allocatable :: ydens
  
  double precision, dimension(:,:), allocatable :: xyzwf
  
  integer, parameter :: data_dim = 32768 ! 32 * 32 * 32 entries
  integer, parameter :: cube_dim = 32
  
  character(255) :: wf_file
  integer :: ios
  integer :: i
  
  ! Process command line argument
  call get_command_argument(number=1, value=wf_file, status=ios)
  if (ios /= 0) then
      print *, "WARNING: No file name provided. Assuming WF_1_1.dat"
      write(wf_file, "(A)") "WF_1_1.dat"
  end if
  
  allocate(xyzwf(4, data_dim))
  call read_xyz_wf(wf_file, xyzwf, data_dim)

  allocate(xdens(2, cube_dim))
  call get_xdens(xyzwf, data_dim, xdens, cube_dim)
  allocate(ydens(2, cube_dim))
  call get_ydens(xyzwf, data_dim, ydens, cube_dim)
  
  ! Print x density
  !print *, "=== x-axis Density (bohr / a.u.) ==="
  do i = 1, cube_dim
      print "(F20.6,F20.8)", xdens(1,i), xdens(2,i)
  end do

  ! Print z density
  !print *, "=== z-axis Density (bohr / a.u.) ==="
  do i = 1, cube_dim
      print "(F20.6,F20.8)", ydens(1,i), ydens(2,i)
  end do

contains
  
  ! Get density along x axis
  subroutine get_xdens(xyz_wf, xyz_dim, x_dens, x_dim)
    
    implicit none
    
    integer, intent(in) :: xyz_dim, x_dim
    double precision, dimension(4, xyz_dim), intent(in) :: xyz_wf
    double precision, dimension(2, x_dim), intent(inout) :: x_dens
    
    double precision :: vol

    integer :: xcnt
    integer :: i
    
    ! Compute volume
    vol = (xyz_wf(3,2) - xyz_wf(3,1)) / 0.529177
    vol = vol * vol * vol

    xcnt = 0
    do i = 1, xyz_dim
        
        if (abs(xyz_wf(2,i)) < 1d-9 .and. abs(xyz_wf(3,i)) < 1d-9) then
            xcnt = xcnt + 1
            x_dens(1,xcnt) = xyz_wf(1,i) / 0.529177 ! x coordinate (in bohr)
            x_dens(2,xcnt) = xyz_wf(4,i) * xyz_wf(4,i) / vol ! density
        end if
    end do
    
  end subroutine get_xdens

  ! Get density along y axis
  subroutine get_ydens(xyz_wf, xyz_dim, y_dens, y_dim)
    
    implicit none
    
    integer, intent(in) :: xyz_dim, y_dim
    double precision, dimension(4, xyz_dim), intent(in) :: xyz_wf
    double precision, dimension(2, y_dim), intent(inout) :: y_dens
    
    double precision :: vol

    integer :: ycnt
    integer :: i
    
    ! Compute volume
    vol = (xyz_wf(3,2) - xyz_wf(3,1)) / 0.529177
    vol = vol * vol * vol

    ycnt = 0
    do i = 1, xyz_dim
        
        if (abs(xyz_wf(2,i)) < 1d-9 .and. abs(xyz_wf(3,i)) < 1d-9) then
            ycnt = ycnt + 1
            y_dens(1,ycnt) = xyz_wf(1,i) / 0.529177 ! x coordinate (in bohr)
            y_dens(2,ycnt) = xyz_wf(4,i) * xyz_wf(4,i) / vol ! density
        end if
    end do
    
  end subroutine get_ydens

  ! Read in xyz WF data
  subroutine read_xyz_wf(wf_file, xyz_wf, dim)
    
    implicit none
    
    character(255), intent(in) :: wf_file
    integer, intent(in) :: dim
    double precision, dimension(4, dim), intent(inout) :: xyz_wf
    
    integer, parameter :: wf_unit = 11
    integer :: i, ios

    open(file=trim(adjustl(wf_file)),unit=wf_unit,status="unknown", &
        action="read",iostat=ios)
    if (ios /= 0) stop "*** Could not open WF.dat file ***"
    
    ! Read each entry: x, y, z, WF
    do i = 1, dim
        read(unit=wf_unit,fmt="(4F20.10)") xyz_wf(1,i), xyz_wf(2,i), xyz_wf(3,i), xyz_wf(4,i)
    end do
    
    close(wf_unit)
  end subroutine read_xyz_wf

end program wfdata_into_xz_training_data

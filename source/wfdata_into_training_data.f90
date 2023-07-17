program wfdata_into_xz_training_data

  implicit none

  double precision, dimension(:,:), allocatable :: xdens
  double precision, dimension(:,:), allocatable :: ydens
  double precision, dimension(:,:), allocatable :: zdens

  double precision, dimension(:,:), allocatable :: xyzwf
  
  integer, parameter :: data_dim = 32768 ! 32 * 32 * 32 entries
  integer, parameter :: cube_dim = 32
  
  character(255) :: wf_file, jobtype_name
  integer :: ios, jobtype
  integer :: i
  
  ! Process command line argument
  call get_command_argument(number=1, value=wf_file, status=ios)
  if (ios .ne. 0) then
      print *, "WARNING: No file name provided. Assuming WF_1_1.dat"
      write(wf_file, "(A)") "WF_1_1.dat"
  end if
  call get_command_argument(number=2, value=jobtype_name, status=ios)
  if (ios .ne. 0) then
      print *, "WARNING: No jobtype selected. Assuming print of x, y, and z."
      jobtype = 3
  else
      read(jobtype_name, "(I3)") jobtype
  end if
  
  allocate(xyzwf(4, data_dim))
  call read_xyz_wf(wf_file, xyzwf, data_dim)

  allocate(xdens(2, cube_dim))
  allocate(ydens(2, cube_dim))
  allocate(zdens(2, cube_dim))
  call get_axis_dens(xyzwf, data_dim, xdens, cube_dim, 1)
  call get_axis_dens(xyzwf, data_dim, ydens, cube_dim, 2)
  call get_axis_dens(xyzwf, data_dim, zdens, cube_dim, 3)
  
  if (jobtype .eq. 1) then
      ! Print x density only
      do i = 1, cube_dim
          print "(F20.5,F20.8)", xdens(1,i), xdens(2,i)
      end do
  else if (jobtype .eq. 2) then
      ! Print x and z densities
      do i = 1, cube_dim
          print "(F20.5,F20.8)", xdens(1,i), xdens(2,i)
      end do
      do i = 1, cube_dim
          print "(F20.5,F20.8)", zdens(1,i), zdens(2,i)
      end do

  else if (jobtype .eq. 3) then
      ! Print x and y and z densities
      do i = 1, cube_dim
          print "(F20.5,F20.8)", xdens(1,i), xdens(2,i)
      end do
      do i = 1, cube_dim
          print "(F20.5,F20.8)", ydens(1,i), ydens(2,i)
      end do
      do i = 1, cube_dim
          print "(F20.5,F20.8)", zdens(1,i), zdens(2,i)
      end do

  else
      stop "*** Error: Invalid Job Type! ***"
  end if

contains
  
  ! Get density along specified axis
  subroutine get_axis_dens(xyz_wf, xyz_dim, y_dens, y_dim, axis)
    
    implicit none
    
    integer, intent(in) :: xyz_dim, y_dim, axis
    double precision, dimension(4, xyz_dim), intent(in) :: xyz_wf
    double precision, dimension(2, y_dim), intent(inout) :: y_dens
    
    double precision :: vol
    
    integer :: a1, a2
    integer :: ycnt
    integer :: i
    
    ! Compute axes we are not interested in
    if (axis .eq. 1) then
        a1 = 2
        a2 = 3
    else if (axis .eq. 2) then
        a1 = 1
        a2 = 3
    else
        a1 = 1
        a2 = 2
    end if


    ! Compute volume
    ! z axis is fast index, so difference is edge length of cube
    vol = (xyz_wf(3,2) - xyz_wf(3,1)) / 0.529177
    vol = vol * vol * vol

    ycnt = 0
    do i = 1, xyz_dim
        
        if (abs(xyz_wf(a1,i)) < 1d-9 .and. abs(xyz_wf(a2,i)) < 1d-9) then
            ycnt = ycnt + 1
            y_dens(1,ycnt) = xyz_wf(axis,i) / 0.529177 ! y coordinate (in bohr)
            y_dens(2,ycnt) = xyz_wf(4,i) * xyz_wf(4,i) / vol ! density
        end if
    end do
    
  end subroutine get_axis_dens

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

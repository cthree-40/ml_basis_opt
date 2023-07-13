program cubedata_into_xz_training_data
  
  implicit none
  
  double precision, dimension(:,:), allocatable :: xdens
  double precision, dimension(:,:), allocatable :: zdens
  
  double precision, dimension(:,:,:), allocatable :: cubedata

  
  integer :: cube_dim
  integer, dimension(3) :: ndisps
  double precision, dimension(3) :: disp_size, start 
  
  character(255) :: cubefile_name
  integer :: ios
  integer :: i
  
  ! Process command line argument
  call get_command_argument(number=1, value=cubefile_name, status=ios)
  if (ios /= 0) then
      print *, "WARNING: No file name provided. Assuming den_p_0.cube"
      write(cubefile_name, "(A)") "den_p_0.cube"
  end if
  
  call read_cubefile_header(cubefile_name, cube_dim, ndisps, disp_size, start)
  
  allocate(cubedata(ndisps(1), ndisps(2), ndisps(3)))
  if (ndisps(1) .ne. ndisps(2)) stop "*** NDISPS(1) != NDISPS(2) ***"
  if (ndisps(1) .ne. ndisps(3)) stop "*** NDISPS(1) != NDISPS(3) ***"
  if (ndisps(2) .ne. ndisps(3)) stop "*** NDISPS(2) != NDISPS(3) ***"
  call read_cubefile_cubedata(cubefile_name, cubedata, ndisps(1))
  
  allocate(xdens(2, ndisps(1)))
  allocate(zdens(2, ndisps(3)))
  call generate_densdata(cubedata, 1, start(1), ndisps(1), disp_size(1), xdens)
  call generate_densdata(cubedata, 3, start(3), ndisps(3), disp_size(3), zdens)

  ! Print x density
  !print *, "=== x-axis Density (bohr / a.u.) ==="
  do i = 1, ndisps(1)
      print "(F20.5,F20.8)", xdens(1,i), xdens(2,i)
  end do

  ! Print z density
  !print *, "=== z-axis Density (bohr / a.u.) ==="
  do i = 1, ndisps(3)
      print "(F20.5,F20.8)", zdens(1,i), zdens(2,i)
  end do
  
  !print *, ""

contains
  
  ! Generate density data along some specified axis (1, 2, or 3)
  subroutine generate_densdata(cubedata, axis, start, ndisps, disp_size, dens)
    
    implicit none
    
    integer, intent(in) :: axis
    integer, intent(in) :: ndisps
    double precision, dimension(ndisps, ndisps, ndisps), intent(in) :: cubedata
    double precision, intent(in) :: start, disp_size

    double precision, dimension(2, ndisps), intent(inout) :: dens
    
    integer :: i, mid

    mid = ndisps / 2
    
    do i = 1, ndisps
        dens(1, i) = start + disp_size * (i - 1)
        select case (axis)
        case (1)
            dens(2, i) = cubedata(i, mid, mid)
        case (2)
            dens(2, i) = cubedata(mid, i, mid)
        case (3)
            dens(2, i) = cubedata(mid, mid, i)
        case default
            stop "*** Error selecting axis for cube data ***"
        end select
    end do

            
  end subroutine generate_densdata

  ! Read cube data from cubefile
  subroutine read_cubefile_cubedata(flnm, cubedata, ndisps)
    
    implicit none
    
    integer, intent(in) :: ndisps
    character(255), intent(in) :: flnm
    double precision, dimension(ndisps * ndisps * ndisps), intent(inout) :: cubedata
    integer, parameter :: flun = 11    
    integer :: i, ios

    open(unit=flun, file=trim(adjustl(flnm)), action="read", status="unknown", &
        position="rewind", iostat=ios)
    if (ios /= 0) stop "*** Could not open den_p.cube file. ***"
    
    ! Skip first nine lines
    do i = 1, 9
        read(unit=flun,fmt=*) 
    end do
    
    ! Read the cube data
    read(unit=flun,fmt="(4e20.8)") cubedata
    
    close(unit=flun)

  end subroutine read_cubefile_cubedata

  ! Read cube header from cubefile
  subroutine read_cubefile_header(flnm, dim, ndisps, disp_size, start)

    implicit none
    
    character(255), intent(in) :: flnm
    integer, intent(out) :: dim
    integer, dimension(3), intent(inout) :: ndisps
    double precision, dimension(3), intent(inout) :: disp_size
    double precision, dimension(3), intent(inout) :: start
    
    integer, parameter :: flun = 11
    integer :: i, ios

    open(unit=flun, file=trim(adjustl(flnm)), action="read", &
        status="unknown", position="rewind", iostat=ios)
    if (ios /= 0) stop "*** Could not open den_p.cube file. ***"

    read(unit=flun,fmt=*) ! Skip first line
    read(unit=flun,fmt=*) ! Skip second line
    
    ! Read dimension, and then starting positions for grid
    read(unit=flun,fmt="(I6,3F13.8)") dim, start(1), start(2), start(3)
    
    ! Read displacement data
    do i = 1, 3
        select case (i)
        case (1)
            read(unit=flun,fmt="(I6,F13.8,26X)") ndisps(i), disp_size(i)
        case (2)
            read(unit=flun,fmt="(I6,13X,F13.8,13X)") ndisps(i), disp_size(i)
        case (3)
            read(unit=flun,fmt="(I6,26X,F13.8)") ndisps(i), disp_size(i)
        case default
            stop "*** Should not be here! Error reading displacements. ***"
        end select
    end do
    
    ! Close file
    close(unit=flun)
    
  end subroutine read_cubefile_header

end program

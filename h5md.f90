module h5md
  use hdf5
  implicit none

  integer :: h5_error

  type h5md_obs
     integer(HID_T) :: obs_id
     double precision, allocatable :: d_buffer(:)
     integer, allocatable :: i_buffer(:)
     integer :: buffer_len
     integer :: buffer_i
     integer :: last_val
  end type h5md_obs

contains

  ! opens a h5md file
  ! creates the '/h5md' group and gives it the attributes 'creator' and 
  ! 'version'. currently, only supports creating a new file.
  ! file_id is the returned hdf5 location of the file
  ! filename is the name of the file
  ! prog_name is the name that appears in the 'creator' global attribute
  subroutine h5md_open_file(file_id, filename, prog_name)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename, prog_name

    character(len=5) :: h5md_version
    integer(HID_T) :: h5_g_id
    integer(HID_T) :: a_type, a_space, a_id
    integer(HSIZE_T) :: a_size(1)

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5_error)

    call h5gcreate_f(file_id, 'h5md', h5_g_id, h5_error)

    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    a_size(1) = len(prog_name)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    
    call h5acreate_f(h5_g_id, 'creator', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, prog_name, a_size, h5_error)

    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    
    h5md_version = '0.1.0'
    a_size(1) = len(h5md_version)

    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    
    call h5acreate_f(h5_g_id, 'version', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, h5md_version, a_size, h5_error)

    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)

    call h5gclose_f(h5_g_id, h5_error)

  end subroutine h5md_open_file

  ! adds a trajectory group in a h5md file
  ! group_name, if present, is the name of a subgroup of 'trajectory'. else,
  ! the group is directly 'trajectory'
  ! also, adds 'time' and 'step' information
  ! if 'time_is_group' is .true., two groups are created, in which individual datasets
  ! will be placed. else, two datasets are created in the given group
  subroutine h5md_create_trajectory_group(file_id, group_name, time_is_group)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in), optional :: group_name
    logical, intent(in), optional :: time_is_group
  end subroutine h5md_create_trajectory_group

  ! add a trajectory dataset to a trajectory group of name trajectory_name
  ! traj_group_id is the id of the trajectory troup
  ! traj_id is the id of the resulting dataset
  ! trajectory_name can be 'position', 'velocity', 'force' or 'species'
  ! N is the number of atoms, D the spatial dimension
  ! species_react is an optional argument. if set to .true., 'species' will be
  ! time dependent, if set to .false., 'species' will not possess the time
  ! dimension
  subroutine h5md_add_trajectory_data(file_id, group_name, trajectory_name, N, D, species_react)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in), optional :: group_name
    character(len=*), intent(in) :: trajectory_name
    integer, intent(in) :: N, D
    logical, intent(in), optional :: species_react
    
  end subroutine h5md_add_trajectory_data

  ! beware: should be doubled for integer and real values, with explicit interfacing
  ! adds a "time frame" to a trajectory
  ! traj_id is the trajectory dataset
  ! data is the actual data of dim (D,N)
  ! step is the integer step of simulation. if the linked 'step' dataset's latest value
  ! is present_step, 'step' and 'time' are not updated. Else, they are.
  subroutine h5md_write_trajectory_data(traj_id, d_data, i_data, present_step, time)
    integer(HID_T), intent(in) :: traj_id
    double precision, intent(in), optional :: d_data(:,:)
    integer, intent(in), optional :: i_data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time
    
  end subroutine h5md_write_trajectory_data

  ! takes a single value and appends it to the appropriate buffer.
  ! if the buffer size is reached, the buffer is dumped to the file.
  subroutine h5md_append_observable_value
  end subroutine h5md_append_observable_value
  
end module h5md

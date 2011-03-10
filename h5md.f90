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
  ! also creates 'trajectory' and 'observables' groups.
  ! file_id is the returned hdf5 location of the file
  ! filename is the name of the file
  ! prog_name is the name that appears in the 'creator' global attribute
  subroutine h5md_open_file(file_id, filename, prog_name)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename, prog_name

    character(len=5) :: h5md_version
    integer(HID_T) :: h5_g_id, g_id
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
    call h5tclose_f(a_type, h5_error)
    
    h5md_version = '0.1.0'
    a_size(1) = len(h5md_version)

    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    
    call h5acreate_f(h5_g_id, 'version', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, h5md_version, a_size, h5_error)

    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    call h5tclose_f(a_type, h5_error)

    call h5gclose_f(h5_g_id, h5_error)

    call h5gcreate_f(file_id, 'trajectory', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

    call h5gcreate_f(file_id, 'observables', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_open_file

  ! adds a trajectory group in a h5md file
  ! group_name is the name of a subgroup of 'trajectory'.
  subroutine h5md_create_trajectory_group(file_id, group_name)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: group_name

    integer(HID_T) :: traj_g_id
    
    call h5gcreate_f(file_id, 'trajectory/'//group_name , traj_g_id, h5_error)

    call h5gclose_f(traj_g_id, h5_error)

  end subroutine h5md_create_trajectory_group

  ! creates the step and time datasets in the specified group
  subroutine h5md_create_step_time(group_id)
    integer(HID_T), intent(inout) :: group_id

    integer :: rank
    integer(HSIZE_T) :: dims(1), max_dims(1), chunk_dims(1)
    integer(HID_T) :: s_id, d_id, plist

    rank = 1
    dims = (/ 0 /)
    max_dims = (/ H5S_UNLIMITED_F /)
    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    chunk_dims = (/ 1024 /)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(group_id, 'step', H5T_NATIVE_INTEGER, s_id, d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5dclose_f(d_id, h5_error)
    call h5sclose_f(s_id, h5_error)
    
    rank = 1
    dims = (/ 0 /)
    max_dims = (/ H5S_UNLIMITED_F /)
    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    chunk_dims = (/ 1024 /)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(group_id, 'time', H5T_NATIVE_DOUBLE, s_id, d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5dclose_f(d_id, h5_error)
    call h5sclose_f(s_id, h5_error)
    
  end subroutine h5md_create_step_time

  ! add a trajectory dataset to a trajectory group of name trajectory_name
  ! group_name is an optional group name for the trajectory group
  ! trajectory_name can be 'position', 'velocity', 'force' or 'species'
  ! N is the number of atoms, D the spatial dimension
  ! species_react is an optional argument. if set to .true., 'species' will be
  ! time dependent, if set to .false., 'species' will not possess the time
  ! dimension
  ! link_from is the name of another trajectory from which the time can be copied
  subroutine h5md_add_trajectory_data(file_id, trajectory_name, N, D, d_id, step_id, time_id, group_name, species_react, link_from)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: trajectory_name
    integer, intent(in) :: N, D
    integer(HID_T), intent(out) :: d_id, time_id, step_id
    character(len=*), intent(in), optional :: group_name
    logical, intent(in), optional :: species_react
    character(len=*), intent(in), optional :: link_from
    
    character(len=128) :: path
    integer :: rank
    integer(HSIZE_T) :: dims(3), max_dims(3), chunk_dims(3)
    integer(HID_T) :: traj_g_id, g_id, s_id, plist

    if ( (trajectory_name .ne. 'position') .and. (trajectory_name .ne. 'velocity') .and. (trajectory_name .ne. 'force') .and. (trajectory_name .ne. 'species') ) then
       write(*,*) 'non conforming trajectory name in h5md_add_trajectory_data'
       stop
    end if

    if (present(group_name)) then
       call h5gopen_f(file_id, 'trajectory/'//group_name, traj_g_id, h5_error)
    else
       call h5gopen_f(file_id, 'trajectory', traj_g_id, h5_error)
    end if
    path = trajectory_name

    ! g_id is opened as the container of trajectory_name
    call h5gcreate_f(traj_g_id, path, g_id, h5_error)
       
    if (trajectory_name .eq. 'species') then
       if (present(species_react) .and. (species_react) ) then
          rank = 2
          dims = (/ N, 0, 0 /)
          max_dims = (/ N, H5S_UNLIMITED_F, 0 /)
          chunk_dims = (/ N, 1, 0 /)
       else
          rank = 1
          dims = (/ N, 0, 0 /)
          max_dims = (/ N, 0, 0 /)
          chunk_dims = (/ N, 0, 0 /)
       end if
    else
       rank = 3
       dims = (/ D, N, 0 /)
       max_dims = (/ D, N, H5S_UNLIMITED_F /)
       chunk_dims = (/ D, N, 1 /)
    end if

    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    if (trajectory_name .ne. 'species') then
       call h5dcreate_f(g_id, 'coordinates', H5T_NATIVE_DOUBLE, s_id, d_id, h5_error, plist)
    else
       call h5dcreate_f(g_id, 'coordinates', H5T_NATIVE_INTEGER, s_id, d_id, h5_error, plist)
    end if
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(s_id, h5_error)

    if (present(link_from)) then
       call h5lcreate_hard_f(traj_g_id, link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(traj_g_id, link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', step_id, h5_error)
    call h5dopen_f(g_id, 'time', time_id, h5_error)

    call h5gclose_f(g_id, h5_error)
    call h5gclose_f(traj_g_id, h5_error)
    
  end subroutine h5md_add_trajectory_data

  ! beware: should be doubled for integer and real values, with explicit interfacing
  ! adds a "time frame" to a trajectory
  ! traj_id is the trajectory dataset
  ! data is the actual data of dim (D,N)
  ! step is the integer step of simulation. if the linked 'step' dataset's latest value
  ! is present_step, 'step' and 'time' are not updated. Else, they are.
  subroutine h5md_write_trajectory_data_d(traj_id, data, present_step, time)
    integer(HID_T), intent(in) :: traj_id
    double precision, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    
    
  end subroutine h5md_write_trajectory_data_d

  ! takes a single value and appends it to the appropriate buffer.
  ! if the buffer size is reached, the buffer is dumped to the file.
  subroutine h5md_append_observable_value
  end subroutine h5md_append_observable_value
  
end module h5md

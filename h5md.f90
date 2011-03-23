
!> This module allows to handle H5MD files.

module h5md
  use hdf5
  implicit none

  !> Global variable to keep the error from HDF5 instructions.
  integer :: h5_error

  !> A type to hold an observable.
  !! Provides a buffering facility.
  type h5md_obs
     !> The ID of the observable in the file
     integer(HID_T) :: id
     !> The double precision buffer.
     double precision, allocatable :: d_buffer(:)
     !> The integer buffer.
     integer, allocatable :: i_buffer(:)
     !> The length of the buffer.
     integer :: buffer_len
     !> The index in the buffer.
     integer :: buffer_i
  end type h5md_obs

  !> A type to hold a reference to a dataset alongside with the step and time information.
  type h5md_t
     !> The ID of the dataset in the file.
     integer(HID_T) :: d_id
     !> The ID of the related step dataset.
     integer(HID_T) :: s_id
     !> The ID of the related time dataset.
     integer(HID_T) :: t_id
  end type h5md_t

  interface h5md_create_obs
     module procedure h5md_create_obs_is
     module procedure h5md_create_obs_i1
     module procedure h5md_create_obs_i2
     module procedure h5md_create_obs_ds
     module procedure h5md_create_obs_d1
     module procedure h5md_create_obs_d2
  end interface

  interface h5md_write_obs
     module procedure h5md_write_obs_is
     module procedure h5md_write_obs_i1
     module procedure h5md_write_obs_i2
     module procedure h5md_write_obs_ds
     module procedure h5md_write_obs_d1
     module procedure h5md_write_obs_d2
  end interface

  
contains

  !> Creates a h5md file
  !!
  !! Creates the '/h5md' group and gives it the attributes 'creator' and 
  !! 'version'. currently, only supports creating a new file.
  !! also creates 'trajectory' and 'observables' groups.
  !! @param file_id the returned hdf5 location of the file.
  !! @param filename name of the file.
  !! @prog_name name that appears in the 'creator' global attribute.
  subroutine h5md_create_file(file_id, filename, prog_name)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename, prog_name

    integer :: h5md_version(2), creation_time, val(8)
    integer(HID_T) :: h5_g_id, g_id
    integer(HID_T) :: a_type, a_space, a_id
    integer(HSIZE_T) :: a_size(1)
    integer :: months(12), i

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5_error)

    call h5gcreate_f(file_id, 'h5md', h5_g_id, h5_error)

    ! write creator attribute
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    a_size(1) = len(prog_name)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    call h5acreate_f(h5_g_id, 'creator', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, prog_name, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    call h5tclose_f(a_type, h5_error)
    
    ! write h5md/version
    h5md_version = (/ 0, 1 /)
    a_size(1) = 2
    call h5screate_simple_f(1, a_size, a_space, h5_error)
    call h5acreate_f(h5_g_id, 'version', H5T_NATIVE_INTEGER, a_space, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_INTEGER, h5md_version, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)

    ! write h5md/creation_time as seconds since Epoch
    months = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    call date_and_time(values=val)
    creation_time = sum(months) * 24*60*60 * (val(1)-1970)
    if (val(2).gt.1) then
       creation_time = creation_time + sum(months(1:val(2)-1)) * 24*60*60
    end if
    do i=1970,val(1)
       if (mod(i,4).eq.0) then
          if (mod(i,100).eq.0) then
             if (mod(i,400).eq.0) creation_time = creation_time + 24*60*60
          else
             creation_time = creation_time + 24*60*60
          end if
       end if
    end do
    creation_time = creation_time + (val(3)-1)*24*60*60
    creation_time = creation_time - val(4)*60  
    creation_time = creation_time + val(5)*60*60
    creation_time = creation_time + val(6)*60
    creation_time = creation_time + val(7)

    a_size(1) = 1
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    call h5acreate_f(h5_g_id, 'creation_time', H5T_NATIVE_INTEGER, a_space, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_INTEGER, creation_time, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)

    call h5gclose_f(h5_g_id, h5_error)

    call h5gcreate_f(file_id, 'trajectory', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

    call h5gcreate_f(file_id, 'observables', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_file

  !> Opens a h5md file
  !!
  !! Returns file_id as the location of the opened file.
  !! @param file_id the returned hdf5 location of the file.
  !! @param filename name of the file.
  !! @param rw flag that allows to open the file in read/write mode.
  subroutine h5md_open_file(file_id, filename, rw)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: rw
    
    if (present(rw) .and. rw) then
       call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5_error)
    else
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5_error)
    end if

  end subroutine h5md_open_file

  !> Adds a trajectory group in a h5md file
  !! @param group_name name of a subgroup of 'trajectory'.
  subroutine h5md_create_trajectory_group(file_id, group_name)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: group_name

    integer(HID_T) :: traj_g_id
    
    call h5gcreate_f(file_id, 'trajectory/'//group_name , traj_g_id, h5_error)

    call h5gclose_f(traj_g_id, h5_error)

  end subroutine h5md_create_trajectory_group

  !> Creates the step and time datasets in the specified group
  !! @param group_id location to place the datasets.
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

  !> Adds a trajectory dataset to a trajectory group of name trajectory_name
  !! @param group_name optional group name for the trajectory group
  !! @param trajectory_name can be 'position', 'velocity', 'force' or 'species'
  !! @param N the number of atoms
  !! @param D the spatial dimension
  !! @param species_react optional argument. if set to .true., 'species' will be
  !! time dependent, if set to .false., 'species' will not possess the time
  !! dimension
  !! @param link_from is the name of another trajectory from which the time can be copied
  subroutine h5md_add_trajectory_data(file_id, trajectory_name, N, D, ID, group_name, species_react, link_from)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: trajectory_name
    integer, intent(in) :: N, D
    type(h5md_t), intent(out) :: ID
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
       call h5dcreate_f(g_id, 'coordinates', H5T_NATIVE_DOUBLE, s_id, ID% d_id, h5_error, plist)
    else
       call h5dcreate_f(g_id, 'coordinates', H5T_NATIVE_INTEGER, s_id, ID% d_id, h5_error, plist)
    end if
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(s_id, h5_error)

    if (present(link_from)) then
       call h5lcreate_hard_f(traj_g_id, link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(traj_g_id, link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)
    call h5gclose_f(traj_g_id, h5_error)
    
  end subroutine h5md_add_trajectory_data

  !> Sets the minimum and maximum attributes to a position trajectory dataset.
  !! @param ID h5md_t of the dataset
  !! @param xmin lower coordinates of the box
  !! @param xmax upper coordinates of the box
  subroutine h5md_set_box_size(ID, xmin, xmax)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: xmin(:)
    double precision, intent(in) :: xmax(:)

    integer(HID_T) :: att_id, att_s
    integer(HSIZE_T) :: att_size(1)

    ! Checking that xmin and xmax have the same size
    if (size(xmin) /= size(xmax)) stop 'non equal dimensions for xmin and xmax in h5md_set_box_size'
    
    ! Setting the size for the attribute
    att_size(1) = size(xmin)

    ! Writing the minimum attribute
    call h5screate_simple_f(1, att_size, att_s, h5_error)
    call h5acreate_f(ID%d_id, 'minimum', H5T_NATIVE_DOUBLE, att_s, att_id, h5_error)
    call h5awrite_f(att_id, H5T_NATIVE_DOUBLE, xmin, att_size, h5_error)
    call h5aclose_f(att_id, h5_error)
    call h5sclose_f(att_s, h5_error)

    ! Writing the maximum attribute
    call h5screate_simple_f(1, att_size, att_s, h5_error)
    call h5acreate_f(ID%d_id, 'maximum', H5T_NATIVE_DOUBLE, att_s, att_id, h5_error)
    call h5awrite_f(att_id, H5T_NATIVE_DOUBLE, xmax, att_size, h5_error)
    call h5aclose_f(att_id, h5_error)
    call h5sclose_f(att_s, h5_error)

  end subroutine h5md_set_box_size

  !> Opens a trajectory dataset and its associated step and time datasets.
  !! @param file_id the HDF5 ID of the file.
  !! @param trajectory_name The name of the trajectory: position, velocity, force or species.
  !! @param ID The h5md_t variable for the dataset.
  !! @param group_name The optional group name.
  subroutine h5md_open_trajectory(file_id, trajectory_name, ID, group_name)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: trajectory_name
    type(h5md_t), intent(out) :: ID
    character(len=*), intent(in), optional :: group_name

    if (present(group_name)) then
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/coordinates', ID%d_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/step', ID%s_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/time', ID%t_id, h5_error)
    else
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/coordinates', ID%d_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/step', ID%s_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/time', ID%t_id, h5_error)
    end if

  end subroutine h5md_open_trajectory

  !> adds a "time frame" to a trajectory
  !! @param traj_id is the trajectory dataset
  !! @param data is the actual data of dim (D,N)
  !! @param step is the integer step of simulation. if the linked 'step' dataset's latest value
  !! is present_step, 'step' and 'time' are not updated. Else, they are.
  !! @todo should be doubled for integer and real values, with explicit interfacing
  subroutine h5md_write_trajectory_data_d(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time
    
    integer(HID_T) :: d_file_space, mem_s, step_id, step_s
    integer(HSIZE_T) :: dims(3), max_dims(3), start(3), num(3)
    integer :: last_step(1)
    double precision :: last_time(1)

    dims(1:2) = shape(data)
    call h5screate_simple_f(2, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id, d_file_space, h5_error)
    call h5sget_simple_extent_dims_f(d_file_space, dims, max_dims, h5_error)
    dims(3) = dims(3) + 1
    call h5sclose_f(d_file_space, h5_error)
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id, d_file_space, h5_error)

    start(1) = 0 ; start(2) = 0 ; start(3) = dims(3)-1
    num(1) = dims(1) ; num(2) = dims(2) ; num(3) = 1

    call h5sselect_hyperslab_f(d_file_space, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, dims, h5_error, mem_space_id=mem_s, file_space_id=d_file_space)
    call h5sclose_f(d_file_space, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)

  end subroutine h5md_write_trajectory_data_d

  !> Load a time frame from a time frame in a H5MD file.
  !>
  !> @param ID h5md_t variable of the dataset
  !> @param data the coordinates array
  !> @param step the integer step from which to choose the data. Currently the only supported value is -1 for the last step.
  subroutine h5md_load_trajectory_data_d(ID, data, step)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(inout) :: data(:,:)
    integer, intent(in) :: step

    integer(HID_T) :: file_space, mem_space
    integer :: ndims, data_shape(2)
    integer(HSIZE_T) :: dims(3), max_dims(3), start(3)

    if (step /= -1) then
       stop 'unknown step in h5md_load_trajectory_data'
    end if

    call h5dget_space_f(ID%d_id , file_space, h5_error)
    call h5sget_simple_extent_ndims_f(file_space, ndims, h5_error)
    if (ndims /= 3) stop 'wrong number of dimensions in h5md_load_trajectory_data'
    call h5sget_simple_extent_dims_f(file_space, dims, max_dims, h5_error)
    data_shape = shape(data)
    if (data_shape(1) /= dims(1) .or. data_shape(2) /= dims(2)) then
       stop 'wrong dimensions for data in h5md_load_trajectory_data'
    end if
    if (dims(3) < 1) stop 'no data in h5md_load_trajectory_data'
    start = (/ 0, 0, 0 /)
    start(3) = dims(3)-1
    dims(3) = 1
    call h5sselect_hyperslab_f(file_space, H5S_SELECT_SET_F, start, dims, h5_error)
    call h5dread_f(ID%d_id, H5T_NATIVE_DOUBLE, data, dims, h5_error, file_space_id=file_space)
    call h5sclose_f(file_space, h5_error)

  end subroutine h5md_load_trajectory_data_d

  !> Adds a "time frame" to a trajectory
  !! @param traj_id is the trajectory dataset
  !! @param data is the actual data of dim (D,N)
  !! @param step is the integer step of simulation. if the linked 'step' dataset's latest value
  !! is present_step, 'step' and 'time' are not updated. Else, they are.
  !! @todo should be doubled for integer and real values, with explicit interfacing
  subroutine h5md_write_trajectory_data_d1d(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time
    
    integer(HID_T) :: d_file_space, mem_s, step_id, step_s
    integer(HSIZE_T) :: dims(3), max_dims(3), start(3), num(3)
    integer :: last_step(1)
    double precision :: last_time(1)

    dims(1:1) = shape(data)
    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id, d_file_space, h5_error)
    call h5sget_simple_extent_dims_f(d_file_space, dims, max_dims, h5_error)
    dims(3) = dims(3) + 1
    call h5sclose_f(d_file_space, h5_error)
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id, d_file_space, h5_error)

    start(1) = 0 ; start(2) = 0 ; start(3) = dims(3)-1
    num(1) = dims(1) ; num(2) = dims(2) ; num(3) = 1

    call h5sselect_hyperslab_f(d_file_space, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, dims, h5_error, mem_space_id=mem_s, file_space_id=d_file_space)
    call h5sclose_f(d_file_space, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)

  end subroutine h5md_write_trajectory_data_d1d

  !> Appends step and time information.
  !! If the last step is the present step, does nothing. Else, append the data.
  !! @param s_id ID of the step dataset
  !! @param t_id ID of the time dataset
  !! @param step The present step
  !! @param time The present time
  subroutine h5md_append_step_time(s_id, t_id, step, time)
    integer(HID_T), intent(inout) :: s_id, t_id
    integer, intent(in) :: step
    double precision, intent(in) :: time

    integer(HID_T) :: file_s, mem_s
    integer(HSIZE_T) :: dims(1), max_dims(1), start(1), num(1)
    integer :: last_step(1)

    ! open step
    call h5dget_space_f(s_id, file_s, h5_error)
    call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)

    if (dims(1) .le. 0) then
       last_step(1) = -1
    else
       start(1) = dims(1)-1
       num(1) = 1
       call h5screate_simple_f(1, num, mem_s, h5_error)
       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dread_f(s_id, H5T_NATIVE_INTEGER, last_step, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(mem_s, h5_error)
    end if

    ! check last
    if (last_step(1) .gt. step) then ! if last > present_step -> error
       write(*,*) 'error, last step is bigger than present step'
    else if (last_step(1) .lt. step) then ! else if last < present_step -> extend step and append present_step, same for time
       ! add step value to the end of the step dataset

       dims(1) = 1
       call h5screate_simple_f(1, dims, mem_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       call h5sclose_f(file_s, h5_error)
       start(1) = dims(1)
       num(1) = 1
       dims(1) = dims(1) + 1
       call h5dset_extent_f(s_id, dims, h5_error)
       call h5dget_space_f(s_id, file_s, h5_error)

       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dwrite_f(s_id, H5T_NATIVE_INTEGER, step, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(file_s, h5_error)
       call h5sclose_f(mem_s, h5_error)

       ! add time value to the end of the time dataset
       dims(1) = 1
       call h5screate_simple_f(1, dims, mem_s, h5_error)

       call h5dget_space_f(t_id, file_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       call h5sclose_f(file_s, h5_error)
       dims(1) = dims(1) + 1
       call h5dset_extent_f(t_id, dims, h5_error)
       call h5dget_space_f(t_id, file_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       start(1) = dims(1) - 1
       num(1) = 1
       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dwrite_f(t_id, H5T_NATIVE_DOUBLE, time, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(file_s, h5_error)
       call h5sclose_f(mem_s, h5_error)

    end if ! else if last = present_step, do nothing

  end subroutine h5md_append_step_time


  !> Close a h5md_obs variable.
  !! @param obs h5md_obs variable.
  subroutine h5md_close_obs(obs)
    type(h5md_obs), intent(inout) :: obs

    if (obs% buffer_i .gt. 0) then
       !call h5md_append_obs_value_d(obs, force_dump = .true.)
    end if

    call h5dclose_f(obs% id, h5_error)

  end subroutine h5md_close_obs
  
  !> Close a h5md_t variable.
  !! @param ID h5md_t variable.
  subroutine h5md_close_ID(ID)
    type(h5md_t), intent(inout) :: ID

    call h5dclose_f(ID% d_id, h5_error)
    call h5dclose_f(ID% s_id, h5_error)
    call h5dclose_f(ID% t_id, h5_error)

  end subroutine h5md_close_ID

  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_i1(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:1) = shape(data)
    max_dims(1:1) = shape(data)
    chunk_dims(1:1) = shape(data)
    

    dims(2)     = 0
    max_dims(2) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(2) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i1


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_is(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1)     = 0
    max_dims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(1) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_is


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_i2(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:2) = shape(data)
    max_dims(1:2) = shape(data)
    chunk_dims(1:2) = shape(data)
    

    dims(3)     = 0
    max_dims(3) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(3) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i2


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_d1(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:1) = shape(data)
    max_dims(1:1) = shape(data)
    chunk_dims(1:1) = shape(data)
    

    dims(2)     = 0
    max_dims(2) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(2) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d1


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_ds(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1)     = 0
    max_dims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(1) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_ds


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  subroutine h5md_create_obs_d2(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:2) = shape(data)
    max_dims(1:2) = shape(data)
    chunk_dims(1:2) = shape(data)
    

    dims(3)     = 0
    max_dims(3) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(3) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'samples', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d2

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_i1(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = dims(2)
    num(2) = 1
    dims(2) = dims(2) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i1

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_is(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1) = 1

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)
    start(1) = dims(1)
    num(1) = 1
    dims(1) = dims(1) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_is

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_i2(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:2) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = dims(3)
    num(3) = 1
    dims(3) = dims(3) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i2

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_d1(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = dims(2)
    num(2) = 1
    dims(2) = dims(2) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d1

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_ds(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1) = 1

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)
    start(1) = dims(1)
    num(1) = 1
    dims(1) = dims(1) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_ds

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  subroutine h5md_write_obs_d2(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:2) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = dims(3)
    num(3) = 1
    dims(3) = dims(3) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d2


end module h5md

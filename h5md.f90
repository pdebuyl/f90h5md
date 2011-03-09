module h5md
  use hdf5
  implicit none

  type h5md_file
     character(len=128) :: file_name
     integer(HID_T) :: file_id
     integer(HID_T) :: traj_group_id
     integer(HID_T) :: obs_group_id
     integer :: error
  end type h5md_file

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
  ! prog_name is the name that appears in the 'creator' global attribute
  subroutine h5md_open_file(file, filename, prog_name)
    type(h5md_file), intent(out) :: file
    character(len=*), intent(in) :: filename, prog_name
  end subroutine h5md_open_file

  ! adds a trajectory group in a h5md file
  ! group_name, if present, is the name of a subgroup of 'trajectory'. else,
  ! the group is directly 'trajectory'
  subroutine h5md_create_trajectory_group(file, group_name)
    type(h5md_file), intent(inout) :: file
    character(len=*), intent(in), optional :: group_name
  end subroutine h5md_create_trajectory_group

  ! adds 'time' and 'step' information
  ! if 'is_group' is .true., two groups are created, in which individual datasets
  ! will be placed. else, two datasets are created in the given group
  subroutine h5md_create_time_step_dataset(traj_group_id, is_group)
    integer(HID_T), intent(in) :: traj_group_id
    logical, intent(in), optional :: is_group
  end subroutine h5md_create_time_step_dataset

  ! add a trajectory dataset to a trajectory group of name trajectory_name
  ! traj_group_id is the id of the trajectory troup
  ! traj_id is the id of the resulting dataset
  ! trajectory_name can be 'position', 'velocity', 'force' or 'species'
  ! N is the number of atoms, D the spatial dimension
  ! species_react is an optional argument. if set to .true., 'species' will be
  ! time dependent, if set to .false., 'species' will not possess the time
  ! dimension
  subroutine h5md_add_trajectory_data(traj_group_id, traj_id, trajectory_name, N, D, species_react)
    integer(HID_T), intent(in) :: traj_group_id
    integer(HID_T), intent(out) :: traj_id
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

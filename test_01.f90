!! Copyright 2011-2013 Pierre de Buyl
!!
!! This file is part of f90h5md
!!
!! f90h5md is free software and is licensed under the modified BSD license (see
!! LICENSE file).

program h5md_test_01
  use h5md
  implicit none

  integer(HID_T) :: file_id
  type(h5md_t) :: pos_ID, vel_ID, kin_ID, extra_ID, fourd_ID
  integer(HSIZE_T), parameter :: D = 3, N = 5
  double precision :: r(D,N), v(D,N), t, r_read(D,N)
  double precision :: fourd(3,2,32,32), fourd_read(3,2,32,32)
  integer :: i, j, i_t, i_t_read
  integer :: mylist(10), mylist_read(10)
  logical :: sw1, sw2(5), sw1_read, sw2_read(5)

  integer :: broken_count
  logical :: test_var

  double precision :: kin, kin_read

  call h5open_f(h5_error)

  call h5md_create_file(file_id, 'data.h5md', 'Pierre de Buyl <pdebuyl@ulb.ac.be>', 'h5md_test_01', 'no version information')

  call h5md_write_par(file_id, 'reason', 'testing')
  mylist = (/ ( i**2, i=1, 10 ) /)
  call h5md_write_par(file_id, 'mylist', mylist)
  sw1 = .true.
  call h5md_write_par(file_id, 'sw1', sw1)
  sw2 = (/ .true., .false., .false., .true., .true. /)
  call h5md_write_par(file_id, 'sw2', sw2)

  call h5md_create_trajectory_group(file_id, 'solvent')

  call h5md_add_trajectory_data(file_id, 'position', N, D, pos_ID, group_name='solvent', compress=.true.)
  call h5md_add_trajectory_data(file_id, 'velocity', N, D, vel_ID, group_name='solvent', link_from='position')
  call h5md_add_trajectory_data(file_id, 'so_do_reac', N, D, extra_ID, group_name='solvent',force_kind='integer',force_rank=2)
  
  call h5md_create_obs(file_id, 'kinetic', kin_ID, kin)
  call h5md_create_obs(file_id, 'fourd', fourd_ID, fourd)
  
  fourd = 0.d0

  r(1,:) = 1.d0
  r(2,:) = 2.d0
  r(3,:) = 3.d0
  r = 100.5d0
  v(1,:) = 10.d0
  v(2,:) = 20.d0
  v(3,:) = 30.d0

  i_t = 0
  t = 0.d0

  do i=0,3
     do j=1,10
        i_t = i_t + 1
        t = t + 0.01d0
        kin = 0.5d0*sum(v**2)
        call h5md_write_obs(kin_ID, kin, i_t, t)
     end do

     r = r + v*0.1d0

     fourd = fourd + 1.d0

     call h5md_write_obs(pos_ID, r, i_t, t)
     call h5md_write_obs(vel_ID, v, i_t, t)
     call h5md_write_obs(fourd_ID, fourd, i_t, t)
  
  end do
  write(*,*) 'i_t = ', i_t

  call h5md_write_par(file_id, 'last_step', i_t)

  call h5md_close_ID(pos_ID)
  call h5md_close_ID(vel_ID)
  call h5md_close_ID(kin_ID)
  call h5md_close_ID(extra_ID)
  call h5md_close_ID(fourd_ID)

  call h5fclose_f(file_id, h5_error)

  call h5md_open_file(file_id, 'data.h5md')
  broken_count = 0

  call h5md_read_par(file_id, 'mylist', mylist_read)
  write(*,*) 'mylist = ', mylist
  write(*,*) 'mylist = ', mylist_read
  call check( (minval(abs(mylist-mylist_read)).eq.0) )

  call h5md_read_par(file_id, 'sw1', sw1_read)
  write(*,*) 'sw1 = ', sw1
  write(*,*) 'sw1 = ', sw1_read
  call check(sw1.eqv.sw1_read)
  call h5md_read_par(file_id, 'sw2', sw2_read)
  write(*,*) 'sw2 = ', sw2
  write(*,*) 'sw2 = ', sw2_read
  test_var = .true.
  do i=1,size(sw2)
     if (sw2(i).neqv.sw2_read(i)) test_var=.false.
  end do
  call check( test_var )

  call h5md_read_par(file_id, 'last_step', i_t_read)
  write(*,*) 'i_t = ', i_t
  write(*,*) 'i_t = ', i_t_read
  call check( i_t .eq. i_t_read )

  call h5md_open_ID(file_ID, pos_ID, 'trajectory', 'solvent/position')
  call h5md_read_obs(pos_ID, r_read, i_t, t)
  write(*,*) 't = ', t, ' , r1 = ', r(:,1)
  call h5md_close_ID(pos_ID)
  call check( minval(abs(r-r_read)).eq.0)

  call h5md_open_ID(file_ID, kin_ID, 'observables', 'kinetic')
  call h5md_read_obs(kin_ID, kin_read, i_t, t)
  write(*,*) 't = ', t, ' , kin = ', kin, kin_read
  call h5md_close_ID(kin_ID)

  call h5md_open_ID(file_ID, fourd_ID, 'observables', 'fourd')
  call h5md_read_obs(fourd_ID, fourd_read, i_t, t)
  write(*,*) 'checking fourd'
  call check( minval(abs(fourd-fourd_read)) < 1d-14 )
  call h5md_close_ID(fourd_ID)

  call h5fclose_f(file_id, h5_error)

  call h5close_f(h5_error)
  
  write(*,*) broken_count, "test(s) failed"

contains
  
  subroutine check(test)
    implicit none
    logical, intent(in) :: test

    if (test) then
       write(*,*) "OK"
    else
       write(*,*) "FAILED"
       broken_count = broken_count + 1
    end if
    write(*,*) '---------------'

  end subroutine check

end program h5md_test_01

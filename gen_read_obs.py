#!/usr/bin/env python

types = dict()
types['i'] = 'integer'
types['d'] = 'double precision'
H5T = dict()
H5T['i'] = 'H5T_NATIVE_INTEGER'
H5T['d'] = 'H5T_NATIVE_DOUBLE'
dims = dict()
dims['s'] = ''
dims['1'] = '(:)'
dims['2'] = '(:,:)'

for t_k,t_v in types.iteritems():
    for d_k,d_v in dims.iteritems():
        if (d_k == 's'):
            rank = 1
        else:
            rank = int(d_k)+1
        s=''
        s+="""  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private"""
        s+="""
  subroutine h5md_read_obs_%s%s(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    %s, intent(out) :: data%s
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = %i
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))
""" % (t_k,d_k,t_v,d_v,rank)
        if (d_k!='s'):
            s+="""
    dims(1:%i) = shape(data)
""" % (rank-1,)
        else:
            s+="""
    dims(1) = 1
"""

        s+="""
    call h5screate_simple_f(%i, dims, mem_s, h5_error)

    call h5dget_space_f(ID%% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
""" % (rank-1, )
        if (d_k!='s'):
            s+="""
    start(1:%i) = 0
    num(1:%i) = dims(1:%i)
""" % (rank-1, rank-1, rank-1)
        s+="""    start(%i) = idx
    num(%i) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID%% d_id, %s, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_%s%s
""" % (rank,rank,H5T[t_k],t_k,d_k)
        print s


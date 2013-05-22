#!/usr/bin/env python

# Copyright 2011-2013 Pierre de Buyl
#
# This file is part of f90h5md
#
# f90h5md is free software and is licensed under the modified BSD license (see
# LICENSE file).

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
dims['3'] = '(:,:,:)'
dims['4'] = '(:,:,:,:)'

for t_k,t_v in types.iteritems():
    for d_k,d_v in dims.iteritems():
        if (d_k == 's'):
            rank = 1
        else:
            rank = int(d_k)+1
        s=''
        s+="""  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private""" 
        s+="""
  subroutine h5md_create_obs_%s%s(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    %s, intent(in) :: data%s
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = %i
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))
""" % (t_k,d_k, t_v, d_v, rank )

        if (d_k!='s'):
            s+="""
    dims(1:%i) = shape(data)
    max_dims(1:%i) = shape(data)
    chunk_dims(1:%i) = shape(data)
    
""" % (rank-1, rank-1, rank-1)


        s+="""
    dims(%i)     = 0
    max_dims(%i) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(%i) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', %s, file_s, ID%% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID%% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID%% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_%s%s

""" % (rank,rank,rank, H5T[t_k],t_k,d_k)
        print s

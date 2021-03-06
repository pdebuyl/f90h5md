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
types['l'] = 'logical'
#types['c'] = 'character'
H5T = dict()
H5T['i'] = 'H5T_NATIVE_INTEGER'
H5T['d'] = 'H5T_NATIVE_DOUBLE'
#H5T['c'] = 'a_type'
H5T['l'] = 'H5T_NATIVE_INTEGER'
dims = dict()
dims['s'] = ''
dims['1'] = '(:)'
dims['2'] = '(:,:)'

for t_k,t_v in types.iteritems():
    if (t_v=='character'):
        dims_var = dict()
        dims_var['s'] = ''
        t_v = 'character(len=*)'
    else:
        dims_var = dims
    for d_k,d_v in dims_var.iteritems():
        if (d_k == 's'):
            rank = 1
        else:
            rank = int(d_k)
        s=''
        s+="""  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private"""
        s+="""
  subroutine h5md_read_par_%s%s(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    %s, intent(out) :: data%s

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(%i)"""  % (t_k,d_k,t_v,d_v,rank)
        if (t_k=='c'):
            s+="""
    integer(HSIZE_T) :: a_size(1)
    integer(HID_T) :: a_type"""
        if (t_k=='l'):
            if (d_k=='s'):
                s+="""
    integer :: data_int
    if (data) then
        data_int = 1
    else
        data_int = 0
    end if
"""
            elif(d_k=='1'):
                s+="""
    integer, allocatable :: data_int%s
    allocate(data_int(size(data)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere
""" % (d_v,)
            elif(d_k=='2'):
                s+="""
    integer, allocatable :: data_int%s
    allocate(data_int(size(data,dim=1),size(data,dim=2)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere
""" % (d_v,)

        if (t_k=='c'):
            s+="""
    a_size(1) = len(data)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)"""

        if (d_k=='s'):
            s+="""
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)
"""
        else:
            s+="""
    dims = shape(data)
    call h5screate_simple_f(%i, dims, par_s, h5_error)
""" % (rank, )


        s+="""
    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)
""" 
        
        if (t_k=='l'):
            s+="""
    call h5dread_f(par_d, %s, data_int, dims, h5_error)
""" % (H5T[t_k],)
        else:
            s+="""
    call h5dread_f(par_d, %s, data, dims, h5_error)
""" % (H5T[t_k],)

        if (t_k=='l'):
            if (d_k=='s'):
                s+="""
    if (data_int.eq.0) then
        data = .false.
    else
        data = .true.
    end if
"""
            else: #(d_k=='1' or d_k=='2'):
                s+="""
    where (data_int .eq. 0)
        data = .false.
    elsewhere
        data = .true.
    endwhere
    deallocate(data_int)
""" 

        s+="""
    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_%s%s
""" % (t_k, d_k)
        print s


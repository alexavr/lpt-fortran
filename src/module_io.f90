module module_io
!===============================================================================
!
! module_io: The module provides some procedures for 
!               input|output operations. Here we have everything 
!               related to NetCDF, reading point.dat and the namelist.
!               
!===============================================================================

use,intrinsic :: iso_fortran_env,only:real32,real64
use netcdf
use module_globals
use omp_lib

private

public :: read_namelist 
public :: save_results 
public :: get_scale
public :: get_offset
public :: get_var_xyz
public :: get_var_xyzt
public :: get_var_xyt
public :: get_dims
public :: get_ndims
public :: get_attr_str
public :: check
public :: open_nc
public :: close_nc
public :: get_pointtxt_num
public :: get_pointtxt
public :: get_pointgrid_num
public :: get_pointgrid

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_namelist
! 
! Reads the namelist. All variables, including name of the namelist
! Are set in module_globals.f90
! 
implicit none
  integer :: unit
  pt_height = fFillValue
  
  open(newunit=unit,file=trim(namelist_name))
  read(unit,nml=data)
  read(unit,nml=scheme) 
  read(unit,nml=particles)
  close(unit)

  pt_nlevels = count(pt_height.ne.fFillValue)

end subroutine read_namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine save_results(file_out,file_in,ncid_in,result,newtime)
implicit none
    character(len=*),intent(in)  :: file_out,file_in
    integer,intent(in) :: ncid_in
    real(kind=real32),intent(in) :: result(:,:,:)
    real(kind=real64),intent(in) :: newtime(:)

    integer :: tdim, pdim, vdim
    integer :: unit, status
    integer :: ncid_out, ia, natts
    integer :: tdim_id, pdim_id, vdim_id, tdim_in_id
    integer :: tvar_id, var_id
    character(len=70) :: tunits, tcalendar, tlong_name, attname
    character(len=120) :: dirname, str_datetime
    character(len=8)  :: str_date
    character(len=10) :: str_time

    call getcwd( dirname )
    write(dirname,'(a,"/")') trim(dirname)
    call date_and_time(DATE=str_date,ZONE=str_time)
    write(str_datetime,'("Created ",a,a)') str_date,str_time

    tdim  = UBOUND(result,1)
    pdim  = UBOUND(result,2)
    vdim  = UBOUND(result,3)

    ! print*,tdim,pdim,vdim

    open(newunit=unit, file=trim(file_out), status="unknown")
    close(unit,status='delete')

    call check( nf90_create( path = file_out , cmode = or(nf90_clobber,nf90_64bit_offset), ncid=ncid_out) )

    ! Dimension
    call check( nf90_def_dim(ncid_out, "lon_lat_hgt",vdim, vdim_id) )  
    call check( nf90_def_dim(ncid_out, "n",          pdim, pdim_id) )     
    call check( nf90_def_dim(ncid_out, "time",       tdim, tdim_id) )

    tunits = get_attr_str(ncid_in,"time","units")
    tcalendar = get_attr_str(ncid_in,"time","calendar")
    ! tlong_name = get_attr_str(ncid_in,"time","long_name")
    
    ! TIME
    call check( nf90_def_var(ncid_out, name="time", xtype=NF90_DOUBLE, dimids = (/tdim_id/), varid = tvar_id ) )
    call check( nf90_inq_varid(ncid_in, "time", tdim_in_id) )
    call check( nf90_inquire_variable(ncid_in, tdim_in_id, natts = natts) )
     do ia = 1, natts
        call check( nf90_inq_attname(ncid_in, tdim_in_id, ia, attname) )
        call check( nf90_copy_att(ncid_in, tdim_in_id, attname, ncid_out, tvar_id) )
    end do

    ! GLOBALS
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "description", "Particle Tracing") )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "codename"   , "pt.f90") )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "history"    , str_datetime) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "pwd"        , dirname) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "src_file"   , file_in) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "out_file"   , file_out) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "start_time" , stime) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "end_time"   , etime) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "timestep"   , timestep) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "cell_detector", cell_detector) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "accuracy"   , merge(1, 0, accuracy  )) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "global"   , merge(1, 0, global  )) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "zoutput"    , merge(1, 0, zoutput   )) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "pt_grid"    , merge(1, 0, pt_grid   )) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "pt_step"    , pt_step) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "pt_height"  , pt_height(:pt_nlevels)) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "horizontal" , merge(1, 0, horizontal)) )
    if (horizontal) call check( nf90_put_att(ncid_out, NF90_GLOBAL, "horizontal_level"   , horizontal_level) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "unstructured_grid" , merge(1, 0, unstructured_grid)) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "cartesian_grid" , merge(1, 0, cartesian_grid)) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "ideal_case" , merge(1, 0, ideal_case)) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "dx"    , dx) )
    call check( nf90_put_att(ncid_out, NF90_GLOBAL, "dy"    , dy) )

    ! RESULT
    call check( nf90_def_var(ncid_out, name="points", xtype=NF90_FLOAT, dimids = (/ tdim_id, pdim_id, vdim_id /), varid = var_id ) )
        call check( nf90_put_att(ncid_out, var_id, "long_name", "Tracing points") )
        call check( nf90_put_att(ncid_out, var_id, "_FillValue", fFillValue) )
        call check( nf90_put_att(ncid_out, var_id, "missing_value", fFillValue) )


    call check( nf90_enddef(ncid_out) )
    call check( nf90_put_var(ncid_out, tvar_id , newtime  ) )
    call check( nf90_put_var(ncid_out, var_id  , result   ) )

    call close_nc(ncid_out)

end subroutine save_results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks if scaling was applied to the data 
! returns 0 - scale factor is NOT applied
function get_scale(ncid,var_id)
implicit none
    integer,intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=70) :: scale_names(7) = (/"SCALE", "Scale", "_scale", "scale_factor", "Scale_factor", "Slope" , "slope"/)
    ! character(len=70) :: offset_name(6) = (/"add_offset", "OFFSET", "Offset", "_offset", "Intercept", "intercept"/)
    logical :: check_scale, scale, offset
    integer ii, xdim,status
    real(kind=real64) :: get_scale

    get_scale = 0

    xdim  = UBOUND(scale_names,1)
    ! if found scale -> stop loop and return the value
    do ii = 1, xdim
        status = nf90_get_att(ncid, var_id, scale_names(ii), get_scale )
        if( status .eq. nf90_noerr ) exit
    end do

end function get_scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks if offset was applied to the data 
! returns 0 - offset is NOT applied
function get_offset(ncid,var_id)
implicit none
    integer,intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=70) :: offset_name(6) = (/"add_offset", "OFFSET", "Offset", "_offset", "Intercept", "intercept"/)
    logical :: check_scale, scale, offset
    integer ii, xdim,status
    real(kind=real64) :: get_offset

    get_offset = 0

    xdim  = UBOUND(offset_name,1)
    ! if found offset -> stop loop and return the value
    do ii = 1, xdim
        status = nf90_get_att(ncid, var_id, offset_name(ii), get_offset )
        if( status .eq. nf90_noerr ) exit
    end do

end function get_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_var_xyz(ncid,var_name,itime,array)
implicit none
    integer,intent(in)  :: ncid, itime
    character(len=*),intent(in)  :: var_name
    real(kind=real32),intent(out) :: array(:,:,:)
    integer :: status, var_id, xdim, ydim, zdim
    real(kind=real64) :: scale_factor,add_offset
    real(kind=real64),dimension(:,:,:),allocatable :: sarray

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)
    zdim  = UBOUND(array,3)

    call check( nf90_inq_varid(ncid, trim(var_name), var_id) )

    add_offset = get_offset(ncid,var_id)
    scale_factor = get_scale(ncid,var_id)

    if(add_offset .ne. 0) then
        allocate( sarray(xdim,ydim,zdim) )

        call check( nf90_get_var(ncid, var_id, sarray , start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )
        array = sarray*scale_factor+add_offset

        deallocate( sarray )
    else
        call check( nf90_get_var(ncid, var_id, array , start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )
    end if

    ! print*,trim(var_name),var_id,minval(array),maxval(array)

end subroutine get_var_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Прочитать данные одной горизонтальной плоскости
subroutine get_var_xyt(ncid,var_name,level,st,tduration,array)
implicit none
    integer,intent(in)  :: ncid,st,tduration,level
    character(len=*), intent(in)  :: var_name
    real(kind=real32),intent(out) :: array(:,:,:)
    integer :: status, var_id, dim_id, xdim, ydim, zdim, tdim
    real(kind=real64) :: scale_factor,add_offset
    real(kind=real32),dimension(:,:,:),allocatable :: sarray
    ! real(kind=real32),dimension(:),    allocatable :: levels

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)
    tdim  = UBOUND(array,3)

    ! ! понять какой индекс у этого уровня
    ! call check( nf90_inq_varid(ncid, "level", var_id) )
    ! call check( nf90_inq_dimid(ncid, "level", dim_id) )
    ! call check( nf90_inquire_dimension(ncid = ncid, dimid = dim_id, len = zdim) )
    ! allocate( levels(zdim) )
    ! call check( nf90_get_var(ncid, var_id, levels ) )

    ! zloc=minloc(abs(levels-level),1)

    call check( nf90_inq_varid(ncid, trim(var_name), var_id) )

    add_offset = get_offset(ncid,var_id)
    scale_factor = get_scale(ncid,var_id)

    if(add_offset .eq. 0) add_offset = 1
    if(scale_factor .eq. 0) scale_factor = 1

    allocate( sarray(xdim,ydim,tdim) )

    call check( nf90_get_var(ncid, var_id, sarray, &
                    start = (/ 1, 1, level, st /), &
                    count = (/ xdim, ydim, 1, (tduration) /)   ) )
    
    array = sarray*scale_factor+add_offset  ! Unscaling если нужно

    deallocate( sarray )
    ! deallocate( levels )

    ! print*,trim(var_name),var_id,minval(array),maxval(array)

end subroutine get_var_xyt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_var_xyzt(ncid,var_name,st,et,array)
! array has 4D: xdim,ydim,zdim,tdim
implicit none
    integer,intent(in)  :: ncid,st,et
    character(len=*),intent(in)  :: var_name
    real(kind=real32),intent(out) :: array(:,:,:,:)
    integer :: status, var_id, xdim, ydim, zdim, tdim
    real(kind=real64) :: scale_factor,add_offset
    real(kind=real64),dimension(:,:,:,:),allocatable :: sarray

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)
    zdim  = UBOUND(array,3)
    tdim  = UBOUND(array,4)

    call check( nf90_inq_varid(ncid, trim(var_name), var_id) )

    add_offset = get_offset(ncid,var_id)
    scale_factor = get_scale(ncid,var_id)

    if(add_offset .eq. 0) add_offset = 1
    if(scale_factor .eq. 0) scale_factor = 1

    allocate( sarray(xdim,ydim,zdim,tdim) )

    call check( nf90_get_var(ncid, var_id, sarray, start = (/ 1, 1, 1, st /), count = (/ xdim, ydim, zdim, et /) ) )
    array = sarray*scale_factor+add_offset

    deallocate( sarray )

    ! print*,trim(var_name),var_id,minval(array),maxval(array)

end subroutine get_var_xyzt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_dims(ncid, lon2d, lat2d, z, time)
implicit none
    integer,intent(in)  :: ncid
    real(kind=real32),intent(out) :: lon2d(:,:), lat2d(:,:), z(:)
    real(kind=real64),intent(out) :: time(:)

    character(len=20) :: lon_names(4) = (/"longitude","lon","west_east","x"/)
    character(len=20) :: lat_names(4) = (/"latitude","lat","south_north","y"/)
    character(len=20) :: z_names(4) = (/"level","bottom_top","z","hight"/)
    character(len=20) :: time_names(1) = (/"time"/)

    integer :: ii, jj, xdim, ydim, zdim, tdim
    integer :: status, var_id

    lon2d = fFillValue 
    lat2d = fFillValue
    z = fFillValue 
    time = fFillValue

    xdim = ubound(lon2d,1)
    ydim = ubound(lon2d,2) 
    zdim = ubound(z,1) 
    tdim = ubound(time,1)


    if (cartesian_grid) then

        do ii = 1, xdim 
            lon2d(ii,:) = ii
        enddo

        do jj = 1, ydim 
            lat2d(:,jj) = jj
        enddo

    else

        do ii = 1, ubound(lon_names,1)
            status = nf90_inq_varid(ncid, trim(lon_names(ii)), var_id)
            if(status == nf90_noerr) then

                if( lon_names(ii).EQ."XLONG") then
                    call check( nf90_get_var(ncid, varid=var_id, values=lon2d ) )
                else
                    ! It is 1D coordinates. So read vector and spread in 2D array
                    call check( nf90_get_var(ncid, varid=var_id, values=lon2d(:,1)) )
                    ! spread in 2d
                    do jj = 2, ydim
                        lon2d(:,jj) = lon2d(:,1)
                    end do
                end if

                exit

            end if
        end do

        do ii = 1, ubound(lat_names,1)
            status = nf90_inq_varid(ncid, trim(lat_names(ii)), var_id)
            if(status == nf90_noerr) then

                if( lat_names(ii).EQ."XLAT") then
                    call check( nf90_get_var(ncid, varid=var_id, values=lat2d ) )
                else
                    ! It is 1D coordinates. So read vector and spread in 2D array
                    call check( nf90_get_var(ncid, varid=var_id, values=lat2d(1,:)) )
                    ! spread in 2d
                    do jj = 2, xdim
                        lat2d(jj,:) = lat2d(1,:)
                    end do
                end if

                exit
                
            end if
        
        end do

    endif

    ! no levels in WRF data, so leaving it blank
    ! also no coordinate spreading needed
    do ii = 1, ubound(z_names,1)
        status = nf90_inq_varid(ncid, trim(z_names(ii)), var_id)
        if(status == nf90_noerr) then
            call check( nf90_get_var(ncid, varid=var_id, values=z) )
            exit 
        end if
    end do

    do ii = 1, ubound(time_names,1)
        status = nf90_inq_varid(ncid, trim(time_names(ii)), var_id)
        if(status == nf90_noerr) then
            call check( nf90_get_var(ncid, varid=var_id, values=time) )
            exit 
        end if
    end do

end subroutine get_dims
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_ndims(ncid, nlon2d, nlat2d, nz, ntime)
implicit none
    integer,intent(in)  :: ncid
    integer,intent(out) :: nlon2d, nlat2d, nz, ntime

    character(len=20) :: lon_names(4) = (/"longitude","lon","west_east","x"/)
    character(len=20) :: lat_names(4) = (/"latitude","lat","south_north","y"/)
    character(len=20) :: z_names(4) = (/"level","bottom_top","z","hight"/)
    character(len=20) :: time_names(1) = (/"time"/)
    
    integer :: ii
    integer :: status, dim_id

    nlon2d = -1 
    nlat2d = -1
    nz = -1 
    ntime = -1

    do ii = 1, ubound(lon_names,1)
        status = nf90_inq_dimid(ncid, trim(lon_names(ii)), dim_id)
        if(status == nf90_noerr) then
            call check( nf90_inquire_dimension(ncid, dimid=dim_id, len=nlon2d) )
            exit
        end if
    end do

    do ii = 1, ubound(lat_names,1)
        status = nf90_inq_dimid(ncid, trim(lat_names(ii)), dim_id)
        if(status == nf90_noerr) then
            call check( nf90_inquire_dimension(ncid, dimid=dim_id, len=nlat2d) )
            exit
        end if
    end do

    do ii = 1, ubound(z_names,1)
        status = nf90_inq_dimid(ncid, trim(z_names(ii)), dim_id)
        if(status == nf90_noerr) then
            call check( nf90_inquire_dimension(ncid, dimid=dim_id, len=nz) )
            exit
        end if
    end do

    do ii = 1, ubound(time_names,1)
        status = nf90_inq_dimid(ncid, trim(time_names(ii)), dim_id)
        if(status == nf90_noerr) then
            call check( nf90_inquire_dimension(ncid, dimid=dim_id, len=ntime) )
            exit
        end if
    end do


end subroutine get_ndims
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_attr_str(ncid,var_name,att_name)
implicit none
    integer,intent(in) :: ncid
    character(len=*), intent(in) :: var_name,att_name
    character(len=120) :: get_attr_str
    integer :: VarId
    
    call check( nf90_inq_varid(ncid, trim(var_name), VarId) )
    call check( nf90_get_att(ncid, VarId, att_name, get_attr_str ))


end function get_attr_str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check(status)
implicit none
    integer, intent (in) :: status
      
    if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
    end if
end subroutine check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine open_nc(filename,ncid)
implicit none
    character(len=*), intent(in) :: filename
    integer,intent(out) :: ncid
    integer :: status

        call check( nf90_open(trim(filename),  NF90_NOWRITE, ncid ) )

end subroutine open_nc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine close_nc(ncid)
implicit none
    integer,intent(in) :: ncid
    call check( nf90_close(ncid) )
end subroutine close_nc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_pointtxt_num(filename)
  implicit none
  character(len=*),INTENT(IN)  :: filename
  integer :: get_pointtxt_num, line, unit
  integer :: ios = 0
  character(len=70)   :: str_tmp

  line = 0
  open(newunit=unit, file=filename, status="old")
  do while(ios == 0)
    read(unit, '(a)', iostat=ios) str_tmp
    if (ios == 0) line = line + 1
  end do
  close(unit)

  get_pointtxt_num = line - 1

end function get_pointtxt_num
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_pointtxt(filename,points)
  implicit none
  character(len=*),INTENT(IN)  :: filename
  real(kind=real32),intent(out) :: points(:,:)

  integer :: get_pointtxt_num, unit, line
  integer :: ios = 0
  character(len=70)   :: buffer
  ! character(len=1)   :: a

  line = 0
  open(newunit=unit, file=filename, status="old")
  do while(ios == 0)
    line = line + 1
    if (line.eq.1) then
      read(unit, *, iostat=ios) buffer
    else 
      read(unit, *, iostat=ios) points(line-1,:)
    end if
  end do
  close(unit)


end subroutine get_pointtxt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_pointgrid_num(istart, nlon, nlat, pt_nlevels, step)
  implicit none
  integer,intent(in)  :: istart, nlon, nlat, pt_nlevels, step
  integer :: get_pointgrid_num, ii, jj

  get_pointgrid_num = 0
  do ii = istart, nlon, step
  do jj = istart, nlat, step
    get_pointgrid_num = get_pointgrid_num + 1
  end do
  end do

  get_pointgrid_num = get_pointgrid_num*pt_nlevels

end function get_pointgrid_num
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_pointgrid(istart,lon2d,lat2d,height,step,points)
  implicit none
  real(kind=real32),intent(in)  :: lon2d(:,:), lat2d(:,:), height(:)
  integer,intent(in)  :: istart, step
  real(kind=real32),intent(out) :: points(:,:)

  integer :: xdim, ydim, iii, ii, jj, kk

  xdim = ubound(lon2d,1)
  ydim = ubound(lon2d,2) 

  iii = 0
  do kk = 1, pt_nlevels
  do ii = istart, xdim, step
  do jj = istart, ydim, step
    iii = iii + 1
    points(iii,1) = lon2d(ii,jj)
    points(iii,2) = lat2d(ii,jj)
    points(iii,3) = height(kk)
  end do
  end do
  end do

end subroutine get_pointgrid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module module_io
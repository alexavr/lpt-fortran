module module_utils

  real(8), parameter :: dFillValue = -(HUGE(1.d0))
  real(4), parameter :: fFillValue = 1E-35
  integer, parameter :: iFillValue = -(HUGE(1))
  real(8), parameter :: g = 9.80665D0, rg = 1.D0/g
  integer, parameter :: maxheightlevels = 100
  character(len=80), parameter  :: namelist_name = "namelist.pt"
  integer :: pt_nlevels

! NAMELIST DATA
  character(len=19)  :: stime,etime
  logical :: regional=.false.
! NAMELIST SCHEME
  logical :: accuracy=.true., horizontal=.false., &
             perfomance=.false., zoutput=.false.
  real(4) :: timestep=12
! NAMELIST PACTICLES
  logical :: pt_grid=.false.
  real(4) :: pt_height(maxheightlevels)
  integer :: pt_step

  namelist /data/   stime,etime,regional
  namelist /scheme/ accuracy,timestep,horizontal,perfomance,zoutput
  namelist /pacticles/ pt_grid,pt_height,pt_step

contains 
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
  real(4),intent(out) :: points(:,:)

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
  real(4),intent(in)  :: lon2d(:,:), lat2d(:,:), height(:)
  integer,intent(in)  :: istart, step
  real(4),intent(out) :: points(:,:)

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
subroutine read_namelist
implicit none
  pt_height = fFillValue
  
  open(555,file=trim(namelist_name))
  read(555,nml=data)
  read(555,nml=scheme) 
  read(555,nml=pacticles)
  close(555)

  pt_nlevels = count(pt_height.ne.fFillValue)

end subroutine read_namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! https://stackoverflow.com/questions/17354100/date-time-difference-in-fortran
! https://github.com/wavebitscientific/datetime-fortran
! 
function strtime2units(time_str,tunits)
use datetime_module
! USE iso_c_binding
! USE iso_fortran_env,ONLY:REAL64
implicit none

  character(len=*),INTENT(IN)  :: time_str
  character(len=*),INTENT(IN)  :: tunits
  character(len=19)   :: tunits_timestamp
  real(8) :: strtime2units, dtime
  TYPE(datetime)  :: date1,date2
  TYPE(timedelta) :: timediff
  integer :: ii, iii
  integer :: icheck

  if (index(tunits,"seconds").NE.0) then
    dtime = 1.D0
  else if (index(tunits,"hours").NE.0) then
    dtime = 3600.D0
  else if (index(tunits,"days").NE.0) then
    dtime = 86400.D0
  else
    print*, "strtime2units"
    print*, "========================================================"
    print*, "Could not find seconds/hours/days in the time attribute "
    print*, "Was the src data prepared? It should have relative time "
    print*, "with 'seconds/hours/days since YYYY-MM-DD HH:mm:SS' attribute."
    print*, "========================================================"
    print*, "STOP."
    STOP
  end if

  icheck = scan(tunits,"1234567890")
  tunits_timestamp = tunits(icheck:icheck+18)

  date1 = strptime(time_str,"%Y-%m-%d %H:%M:%S")
  date2 = strptime(tunits_timestamp,"%Y-%m-%d %H:%M:%S")
  timediff = date1-date2

  ! print*,trim(tunits)
  ! print*,trim(time_str)," ",trim(tunits_timestamp)

  strtime2units = timediff%total_seconds()/dtime

end function strtime2units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! https://stackoverflow.com/questions/17354100/date-time-difference-in-fortran
! https://github.com/wavebitscientific/datetime-fortran
! 
subroutine convert_utime2str(utime,tunits,strtime)
use datetime_module
! USE iso_c_binding
! USE iso_fortran_env,ONLY:REAL64
implicit none

  real(8),intent(in)  :: utime(:)
  character(len=*),intent(in)  :: tunits
  character(len=19),intent(out)  :: strtime(:)

  character(len=19)   :: tunits_timestamp
  TYPE(datetime) :: date1
  TYPE(datetime),dimension(:),allocatable  :: date2(:)
  TYPE(timedelta),dimension(:),allocatable :: timediff(:)
  integer :: xdim
  integer :: icheck, nspace

  xdim = ubound(utime,1)
  allocate( timediff(xdim)  )
  allocate( date2(xdim)  )

  ! everything into minutes
  if (index(tunits,"seconds").NE.0) then
    timediff = timedelta(minutes=int(utime/60.))
  else if (index(tunits,"minutes").NE.0) then
    timediff = timedelta(minutes=int(utime))
  else if (index(tunits,"hours").NE.0) then
    timediff = timedelta(minutes=int(utime*60))
  else if (index(tunits,"days").NE.0) then
    timediff = timedelta(minutes=int(utime*24*60))
  else
    print*, "convert_utime2str:"
    print*, "========================================================"
    print*, "Could not find seconds/hours/days in the time attribute "
    print*, "Was the src data prepared? It should have relative time "
    print*, "with 'seconds/hours/days since YYYY-MM-DD HH:mm:SS' attribute."
    print*, "========================================================"
    print*, "STOP."
    STOP
  end if

  icheck = scan(tunits,"1234567890")
  tunits_timestamp = tunits(icheck:icheck+18)

  date1 = strptime(tunits_timestamp,"%Y-%m-%d %H:%M:%S")
  date2 = date1 + timediff
  strtime = date2 % isoformat()

  deallocate (timediff,date2)

end subroutine convert_utime2str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function get_tunit(tunits)
implicit none

  character(len=*),INTENT(IN)  :: tunits
  character(len=19)   :: get_tunit
  integer :: ii, iii

  if (index(tunits,"seconds").NE.0) then
    get_tunit = "seconds"
  else if (index(tunits,"hours").NE.0) then
    get_tunit = "hours"
  else if (index(tunits,"days").NE.0) then
    get_tunit = "days"
  else
    print*, "get_tunit:"
    print*, "========================================================"
    print*, "Could not find seconds/hours/days in the time attribute "
    print*, "Was the src data prepared? It should have relative time "
    print*, "with 'seconds/hours/days since YYYY-MM-DD HH:mm:SS' attribute."
    print*, "========================================================"
    print*, "STOP."
    STOP
  end if

end function get_tunit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module module_utils
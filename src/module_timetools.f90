module module_timetools
!===============================================================================
!
! module_timetools: The module provides procedures related to time operations.  
!               
!===============================================================================

use,intrinsic :: iso_fortran_env,only:real32,real64
use module_globals

private
public :: strtime2units  
public :: convert_utime2str
public :: get_tunit

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function strtime2units(time_str,tunits)
! https://stackoverflow.com/questions/17354100/date-time-difference-in-fortran
! https://github.com/wavebitscientific/datetime-fortran
! 
! Takes time as string and converts it into relative time based on description 
! in tunits. The time_str has to be in %Y-%m-%d %H:%M:%S format.
! The tunits has to be in standart "seconds/hours/days since YYYY-MM-DD HH:mm:SS"
! 
use datetime_module
implicit none

  character(len=*),INTENT(IN)  :: time_str
  character(len=*),INTENT(IN)  :: tunits
  character(len=19)   :: tunits_timestamp
  real(kind=real64) :: strtime2units, dtime
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

  strtime2units = timediff%total_seconds()/dtime

end function strtime2units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine convert_utime2str(utime,tunits,strtime)
! https://stackoverflow.com/questions/17354100/date-time-difference-in-fortran
! https://github.com/wavebitscientific/datetime-fortran
! 
! Converts relative time into string time based on 
! tunits description. 
! The tunits of utime has to be in standart 
! "seconds/hours/days since YYYY-MM-DD HH:mm:SS". 
! 
use datetime_module
implicit none

  real(kind=real64),intent(in)  :: utime(:)
  character(len=*),intent(in)   :: tunits
  character(len=19),intent(out) :: strtime(:)

  character(len=19)   :: tunits_timestamp
  TYPE(datetime) :: date1
  TYPE(datetime),dimension(:),allocatable  :: date2(:)
  TYPE(timedelta),dimension(:),allocatable :: timediff(:)
  integer :: xdim
  integer :: icheck, nspace

  xdim = ubound(utime,1)
  allocate( timediff(xdim)  )
  allocate( date2(xdim)  )

  ! Turning everything into minutes
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
! 
! Parses unit of relative time
! 
function get_tunit(tunits)
implicit none

  character(len=*),INTENT(IN)  :: tunits
  character(len=19)   :: get_tunit
  integer :: ii, iii

  if (index(tunits,"seconds").NE.0) then
    get_tunit = "seconds"
  else if (index(tunits,"minutes").NE.0) then
    get_tunit = "minutes"
  else if (index(tunits,"hours").NE.0) then
    get_tunit = "hours"
  else if (index(tunits,"days").NE.0) then
    get_tunit = "days"
  else
    print*, "get_tunit:"
    print*, "========================================================"
    print*, "Could not find seconds/minutes/hours/days in the time attribute "
    print*, "Was the src data prepared? It should have relative time "
    print*, "with 'seconds/minutes/hours/days since YYYY-MM-DD HH:mm:SS' attribute."
    print*, "========================================================"
    print*, "STOP."
    STOP
  end if

end function get_tunit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module module_timetools
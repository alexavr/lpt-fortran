module module_globals
use,intrinsic :: iso_fortran_env,only:real32,real64
implicit none

  real(kind=real64), parameter :: dFillValue = -(HUGE(1.d0))
  real(kind=real32), parameter :: fFillValue = 1E-35
  integer, parameter :: iFillValue = -(HUGE(1))
  real(kind=real64), parameter :: g = 9.80665D0, rg = 1.D0/g
  integer, parameter :: maxheightlevels = 100
  character(len=80), parameter  :: namelist_name = "namelist.pt"
  character(len=80), parameter  :: file_start_particles = "init_particles.dat"

  integer :: pt_nlevels

! NAMELIST DATA
  character(len=19)  :: stime,etime
  logical :: regional=.false.
! NAMELIST WRF_PREP
  character(len=80)  :: type = "model_levels"
  real(kind=real32)  :: levels(maxheightlevels)
! NAMELIST SCHEME
  logical :: accuracy=.true., horizontal=.false., &
             perfomance=.false., zoutput=.false.
  real(kind=real32) :: timestep=12
  real(kind=real32) :: horizontal_level=0
! NAMELIST PACTICLES
  logical :: pt_grid=.false.
  real(kind=real32) :: pt_height(maxheightlevels)
  integer :: pt_step,cell_detector=0

  namelist /data/   stime,etime,regional
  namelist /wrf_prep/ type,levels
  namelist /scheme/ accuracy,timestep,horizontal,horizontal_level,perfomance,zoutput,cell_detector
  namelist /particles/ pt_grid,pt_height,pt_step

end module module_globals

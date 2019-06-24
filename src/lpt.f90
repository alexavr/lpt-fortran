! erai_2015-01_NA.nc : no point interpolation and output: 0m23.317s
! Plans:
! 1. Add point interpolation
! 2. Add output 
! 3. Add perfomance measurement
! 
program lpt
  use,intrinsic :: iso_fortran_env,only:real32,real64
  use module_globals
  use module_io
  use module_grid
  use module_timetools
  implicit none

  character(len=120)  :: file_in, file_out
  integer  :: ncid_in, ncid_out

  character(len=120)  :: tunits, tunit
  real(kind=real64)  :: stimeu,etimeu
  integer  :: stimeui,etimeui

  integer  :: nlon, nlat, nz, ntime_in
  real(kind=real32),dimension(:,:),allocatable :: lon2d, lat2d
  real(kind=real32),dimension(:),allocatable :: levels
  real(kind=real64),dimension(:),allocatable :: time_in

  real(kind=real64)  :: dt,dtm,dtt
  integer  :: ntimesteps, duration, npoints

  real(kind=real32),dimension(:,:),allocatable :: points

  integer  :: ii, itime, it, ip

  real(kind=real32),dimension(:,:,:)      ,allocatable :: result
  real(kind=real64),dimension(:)          ,allocatable :: newtime_ind,newtime
  character(len=19),dimension(:),allocatable :: newtime_datatime
  
  logical :: update

  real(kind=real32),dimension(:,:,:),allocatable :: u,v,w,z
  real(kind=real32),dimension(:,:,:),allocatable :: u1,v1,w1,z1
  real(kind=real32),dimension(:,:,:),allocatable :: u2,v2,w2,z2
  real(kind=real64) :: ditime
  real(kind=real32) :: pu, pv, pw

  ! For interpolation used rodust scheme:
  ! horisontal: min distance to 4 points (x,y) [ij(4,2)]
  ! vertical: takes ij points and finds levels below and above [4,2]
  integer :: status, ij(4,2), kk(4,2) 

  call get_command_argument(1,file_in)
  if(trim(file_in).EQ."") then
    stop "Usage: ./pt.exe input_file.nc"
  end if
  file_out = "_pt_"//trim(file_in)

  print*,"-> Files:"
  print*,"   Input file:      ",trim(file_in)
  print*,"   Output file: ",trim(file_out)
  
  print*,"-> Reading namelist..."
  call read_namelist

  print*,"-> Particles..."
  if(pt_grid) then
    print*,"   Particles are spread all over the domain"
    print*,"           with step  : ",pt_step
    print*,"           at level(s): ",pt_height(:pt_nlevels)
  else 
    print*,"   Input stat points file: ",trim(file_start_particles)
  end if

  print*, "-> Opening src file..."
  call open_nc(file_in,ncid_in)
  print*, "-> Converting stime and etime to src time..."
  tunits = get_attr_str(ncid_in,"time","units")
  print*, "   ",trim(tunits)
  stimeu = strtime2units(stime,tunits)
  etimeu = strtime2units(etime,tunits)

  print*, "-> Getting coordinates..."
  call get_ndims(ncid_in, nlon, nlat, nz, ntime_in)
  print*, "   Dimensions are "
  write(*,'( "    nlon = ",i3,"; nlat = ",i3,"; nz = ",i3,"; ntime = ",i3 )') nlon, nlat, nz, ntime_in 

  allocate( lon2d   ( nlon, nlat ) )
  allocate( lat2d   ( nlon, nlat ) )
  allocate( time_in ( ntime_in   ) )
  allocate( levels  ( nz         ) )
  call get_dims(ncid_in, lon2d, lat2d, levels, time_in)
  print*, "   Coordinates ... Ok"

  print*, "-> Finding stime and etime indexes in src..."
  stimeui = minloc( abs(time_in-stimeu),1 )
  etimeui = minloc( abs(time_in-etimeu),1 )
  ! print*,stimeu,etimeu
  ! print*,time_in(stimeui),time_in(etimeui)

  print*, "-> Read or create star points..."
  if(pt_grid) then
    npoints =  get_pointgrid_num(3, nlon, nlat, pt_nlevels, pt_step)
    allocate( points(npoints,3) )
    call get_pointgrid(3,lon2d,lat2d,pt_height(:pt_nlevels),pt_step,points)
  else
    npoints =  get_pointtxt_num(file_start_particles)
    allocate( points(npoints,3) )
    call get_pointtxt(file_start_particles,points)
  end if

 print*, "-> Time driver..."
 ! If accuracy: calculate number in timesteps between src timesteps
 ! Since timestep is in mins, than we need to convert in min (dtm) 
  dt = time_in(2) - time_in(1)
  tunit = get_tunit(tunits)

  if(trim(tunit).eq."seconds") then
    dtm = dt/60
  else if(trim(tunit).eq."minutes") then
    dtm = dt
  else if(trim(tunit).eq."hours") then 
    dtm = dt*60
  else if(trim(tunit).eq."days") then
    dtm = dt*24*60
  end if

  ntimesteps = 1
  if(accuracy) ntimesteps = int(dtm)/int(timestep)

  if(.not.accuracy) timestep = dtm

  duration = etimeui - stimeui
  dtt = dt/ntimesteps

  allocate( result           (3,npoints,duration*ntimesteps+1)  )
  allocate( newtime_ind      (          duration*ntimesteps+1)  )
  allocate( newtime          (          duration*ntimesteps+1)  )
  allocate( newtime_datatime (          duration*ntimesteps+1)  )

  result = fFillValue

  newtime_ind = [ (stimeui+(ii-1)*1./ntimesteps, ii=1,(duration*ntimesteps+1) ) ]
  newtime = [ (time_in(stimeui)+(ii-1)*dtt, ii=1,(duration*ntimesteps+1) ) ]

  call convert_utime2str(newtime,tunits,newtime_datatime)

  allocate ( u(nlon, nlat, nz) )
  allocate ( v(nlon, nlat, nz) )
  allocate ( w(nlon, nlat, nz) )
  allocate ( z(nlon, nlat, nz) )
  if(accuracy) then
    allocate ( u1(nlon, nlat, nz) )
    allocate ( v1(nlon, nlat, nz) )
    allocate ( w1(nlon, nlat, nz) )
    allocate ( z1(nlon, nlat, nz) )
    allocate ( u2(nlon, nlat, nz) )
    allocate ( v2(nlon, nlat, nz) )
    allocate ( w2(nlon, nlat, nz) )
    allocate ( z2(nlon, nlat, nz) )
  end if

  print*, "-> Tracing..."
  update = .true.
  do itime = 1, (duration*ntimesteps+1)

    write(*,'(a," ",f5.1,"%")') newtime_datatime(itime), float(count(points(:,1).NE.fFillValue))/float(npoints)*100.

    ! If all points are out -- no need to waste our computation time
    if(all(points.eq.fFillValue)) exit

    it = int(floor(newtime_ind(itime)))

    if(itime.ne.1) update = ( it .gt. int(floor(newtime_ind(itime-1))) )

    ! print*, newtime_datatime(itime), newtime_ind(itime), it

    if(accuracy) then

      if (update) then

        ! print*,"   reading..."

        call get_var_xyz(ncid_in,"u",it  ,u1) 
        call get_var_xyz(ncid_in,"v",it  ,v1) 
        call get_var_xyz(ncid_in,"w",it  ,w1) 
        call get_var_xyz(ncid_in,"z",it  ,z1) 
        call get_var_xyz(ncid_in,"u",it+1,u2) 
        call get_var_xyz(ncid_in,"v",it+1,v2) 
        call get_var_xyz(ncid_in,"w",it+1,w2) 
        call get_var_xyz(ncid_in,"z",it+1,z2) 
        update = .false.

      end if
        
        ditime = newtime_ind(itime)-it   ! 
        z = (1.-ditime)*z1 + (ditime)*z2 ! Linear inperpolation between 
        u = (1.-ditime)*u1 + (ditime)*u2 !    src timesteps 
        v = (1.-ditime)*v1 + (ditime)*v2 !    weights are 0..1 =>
        w = (1.-ditime)*w1 + (ditime)*w2 !    no need to normalize

    else
        call get_var_xyz(ncid_in,"u",it,u) 
        call get_var_xyz(ncid_in,"v",it,v) 
        call get_var_xyz(ncid_in,"w",it,w) 
        call get_var_xyz(ncid_in,"z",it,z) 
        ditime = 0
    end if

    ! calculating heights (for z interpolation)
    z = z * rg ! geopotential / g = h

    if(itime.eq.1) then
        do ip = 1, npoints
          result(:,ip,itime) = points(ip,:)
          ! print'(a,24x, 2f10.4, f10.2)', newtime_datatime(itime), points(ip,:)
        end do
    else

      do ip = 1, npoints

        call get_cell_hor ( points(ip,:), lon2d, lat2d, ij) 

        call get_cell_vert( points(ip,:), ij, z, kk )

        call interpolate  ( points(ip,:), u, v, w, z, ij, kk, &
                            lon2d, lat2d, levels, &
                            pu, pv, pw) 
        
        call locate       ( points(ip,:), pu, pv, pw, timestep )
        
        if( any(points(ip,:).eq.fFillValue) ) cycle

        ! print'(a,2f7.2,f10.4, 2f10.4, f10.2)', newtime_datatime(itime), pu, pv, pw, points(ip,:)

        result(:,ip,itime) = points(ip,:)

      end do

    end if

  end do

  deallocate( u,v,w,z )
  if(accuracy) then
    deallocate( u1,v1,w1,z1 )
    deallocate( u2,v2,w2,z2 )
  end if

  ! OUTPUT
  call save_results(file_out,file_in,ncid_in,result,newtime)

  deallocate( points,lon2d,lat2d,time_in,levels           )
  deallocate( result,newtime_ind,newtime,newtime_datatime )

  call close_nc(ncid_in)

! contains 

end program lpt
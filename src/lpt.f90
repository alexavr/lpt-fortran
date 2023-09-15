! RUN:
! export OMP_NUM_THREADS=16 
! ./lpt.exe src_file.nc


program lpt
  use,intrinsic :: iso_fortran_env,only:real32,real64
  use module_globals
  use module_io
  use module_grid
  use module_timetools
  use omp_lib
  implicit none

  character(len=120)  :: file_in, file_out
  integer  :: ncid_in, ncid_out

  character(len=120)  :: tunits, tunit
  real(kind=real64)  :: stimeu,etimeu
  integer  :: stimeui,etimeui, t1, t2

  integer  :: nlon, nlat, nz, ntime_in
  real(kind=real32),dimension(:,:),allocatable :: lon2d, lat2d
  ! real(kind=real32),dimension(:),allocatable :: levels
  real(kind=real64),dimension(:),allocatable :: time_in

  real(kind=real64)  :: dt,dtm,dtt, timestep_r, ditime
  integer  :: ntimesteps, duration, npoints
  integer  :: iz

  real(kind=real32),dimension(:,:),allocatable :: points

  integer  :: ii, jj, itime, it, ip, ipoint, ps, pe

  real(kind=real32),dimension(:,:,:)      ,allocatable :: result
  real(kind=real64),dimension(:)          ,allocatable :: newtime_ind,newtime
  character(len=19),dimension(:),allocatable :: newtime_datatime
  
  logical :: first, good

  real(kind=real32),dimension(:,:,:),allocatable :: u2d,v2d
  logical,dimension(:,:),allocatable :: mask
  real(kind=real32),dimension(:,:),allocatable :: u2d_tmp,v2d_tmp
  ! real(kind=real32),dimension(:,:,:,:),allocatable :: u,v,w,z
  ! real(kind=real32),dimension(:,:,:,:),allocatable :: u1,v1,w1,z1
  ! real(kind=real32),dimension(:,:,:,:),allocatable :: u2,v2,w2,z2
  real(kind=real32) :: pu, pv, pw

  ! For interpolation used rodust scheme:
  ! horisontal: min distance to 4 points (x,y) [ij(4,2)]
  ! vertical: takes ij points and finds levels below and above [4,2]
  integer :: status, ij(4,2), kk(4,2) 

  integer :: tid, nthreads, nmaxthreads ! OpenMP OMP

  call get_command_argument(1,file_in)
  if(trim(file_in).EQ."") then
    stop "Usage: ./lpt.exe input_file.nc"
  endif
  
  ! file_out = trim(file_in)
  file_out = file_in(1:len(trim(file_in))-3)//"_tracks.nc"

  print*,"-> Files:"
  print*,"   Input file:      ",trim(file_in)
  print*,"   Output file: "    ,trim(file_out)

  print*,"-> Reading namelist..."
  call read_namelist

  if (horizontal) then
    print*,"-> DOING HORIZONTAL SIMULATION"
  else
    print*,"-> DOING 3D SIMULATION"
  endif

  print*,"-> Particles..."
  if(pt_grid) then
    print*,"    * Particles are spread all over the domain"
    print*,"           with step  : ",pt_step
    print*,"           at level(s): ",pt_height(:pt_nlevels)
  else 
    print*,"    * Input stat points file: ",trim(file_start_particles)
  endif

  print*, "-> Opening src file..."
  call open_nc(file_in,ncid_in)
  print*, "#####################################################################"
  print*, "-> Converting stime and etime to src time..."
  tunits = get_attr_str(ncid_in,"time","units")
  print*, "   ",trim(tunits)
  stimeu = strtime2units(stime,tunits)
  etimeu = strtime2units(etime,tunits)
  print*, "#####################################################################"

  print*, "-> Getting coordinates..."
  call get_ndims(ncid_in, nlon, nlat, nz, ntime_in)
  print*,    "   * Dimensions are "
  write(*,'( "   * nlon = ",i3,"; nlat = ",i3,"; nz = ",i3,"; ntime = ",i3 )') nlon, nlat, nz, ntime_in

  allocate( lon2d   ( nlon, nlat ) )
  allocate( lat2d   ( nlon, nlat ) )
  allocate( time_in ( ntime_in   ) )
  ! allocate( levels  ( nz         ) )
  call get_dims(ncid_in, lon2d, lat2d, levels(:nz), time_in)
  print*, "   * Coordinates ... Ok"
  print*, "#####################################################################"

  print*, "-> Finding stime and etime indexes in src..."
  stimeui = minloc( abs(time_in-stimeu),1 )
  etimeui = minloc( abs(time_in-etimeu),1 )
  duration = etimeui - stimeui
  write(*,'( "   * Working with ",i0," time steps ")') duration
  print*, "#####################################################################"

  iz = MINLOC(abs(levels(:nz)-horizontal_level), 1)

  print*, "-> Read or create start points..."
  if(pt_grid) then
    npoints =  get_pointgrid_num(3, nlon, nlat, pt_nlevels, pt_step)
    allocate( points(npoints,3) )
    call get_pointgrid(3,lon2d,lat2d,pt_height(:pt_nlevels),pt_step,points)
  else
    npoints =  get_pointtxt_num(file_start_particles)
    allocate( points(npoints,3) )
    call get_pointtxt(file_start_particles,points)
  endif
  ! points = fFillValue
  ! npoints = 96 ! npoints - 1

  print*, "   * Number of points: ", npoints
  print*, "#####################################################################"

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
  endif

  ! Шаги по времени (если accuracy, то учитываем шаги между шагами в исходном файле)
  ntimesteps = 1
  if (accuracy) then
    ntimesteps = int(dtm)/int(timestep)
  else
    timestep = dtm
  endif

  dtt = dt/ntimesteps   ! минимальный шаг по времени

  allocate( result           (3,npoints,duration*ntimesteps+1)  )
  allocate( newtime_ind      (          duration*ntimesteps+1)  )
  allocate( newtime          (          duration*ntimesteps+1)  )
  allocate( newtime_datatime (          duration*ntimesteps+1)  )

  result = fFillValue

  timestep_r = 1./ntimesteps
  newtime = [ (time_in(stimeui)+(ii-1)*dtt, ii=1,(duration*ntimesteps+1) ) ]
  newtime_ind = [ (stimeui+(ii-1)*timestep_r, ii=1,(duration*ntimesteps+1) ) ]

  call convert_utime2str(newtime,tunits,newtime_datatime)

  ! allocate ( w(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  ! allocate ( z(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  ! if(accuracy) then
  !   allocate ( u1(nlon, nlat, nz, ntime_in) )
  !   allocate ( v1(nlon, nlat, nz, ntime_in) )
  !   allocate ( w1(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  !   allocate ( z1(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  !   allocate ( u2(nlon, nlat, nz, ntime_in) )
  !   allocate ( v2(nlon, nlat, nz, ntime_in) )
  !   allocate ( w2(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  !   allocate ( z2(nlon, nlat, nz, ntime_in) )  ! if (.not.horizontal)
  ! endif
  print*, "#####################################################################"

  print*, "-> Tracing..."


  write(*,'("   * Reading data......")')


  first = .true.

  allocate ( u2d(nlon, nlat, duration) )
  allocate ( v2d(nlon, nlat, duration) )

  call get_var_xyt(ncid_in, "u", iz, stimeui, duration, u2d)
  call get_var_xyt(ncid_in, "v", iz, stimeui, duration, v2d)

  allocate ( u2d_tmp(nlon, nlat) )
  allocate ( v2d_tmp(nlon, nlat) )
  allocate ( mask(nlon, nlat) )


  write(*,'("   * Entering the parallel block......")')
!$omp parallel private(tid,mask,u2d_tmp,v2d_tmp,itime,ij,pu, pv,t1,t2,ipoint,ditime,ii,jj,good,ps, pe) ! shared(nthreads,npoints,stimeui,nlon,nlat,ntimesteps,duration,u2d,v2d,newtime_ind,horizontal,result,points)
! !$omp parallel default(none) private(tid,mask,u2d_tmp,v2d_tmp,itime,ij,pu, pv,t1,t2,ipoint,ditime,ii,jj,good,ps, pe) shared(nthreads,npoints,stimeui,nlon,nlat,ntimesteps,duration,u2d,v2d,newtime_ind,horizontal,result,points,lon2d,lat2d,timestep)
! !$omp parallel default(none) private(tid,ps,pe,t1,t2,ditime,mask,good,u2d_tmp,              &
                                    !  v2d_tmp,ij,pu,pv,ipoint,itime,ii,jj) shared(nthreads,  &
                                    ! npoints,stimeui,nlon,nlat,ntimesteps,duration,          &
                                    ! u2d,v2d,newtime_ind,horizontal,result,points,           &
                                    ! lon2d,lat2d,timestep)

  tid = OMP_GET_THREAD_NUM()
  nthreads = OMP_GET_NUM_THREADS()

  ps = 1+(tid*npoints/nthreads)
  pe = (tid+1)*npoints/nthreads

  write(*,'("   * OpenMP calculations on ",i3," core of ",i3," | from  ",i6," to ",i6 " / ",i3)') tid, nthreads, ps, pe, (pe - ps + 1)

  outer: do ipoint = ps, pe

    if (tid .EQ. 0) then
      write(*,'("   * * Progress appoximatley ",f5.1," % ")') float(ipoint)/float(pe-ps+1)*100. ! minval((/float(ipoint)/float(pe-ps+1)*100.,100./))
    endif

! print*, ipoint
    if (horizontal) then ! Если считаем только на проскости, то

! print*, "horizontal"
      inner: do itime = 1, duration*ntimesteps+1

! print*, itime
        ! Считаем веса для интерполяции по времени
        t1 = floor(newtime_ind(itime))
        t2 = ceiling(newtime_ind(itime))
        ditime = newtime_ind(itime)-t1   !
        ! А теперь считаем индексы для обрезанного массива (с 1 до duration)
        t1 = t1 - stimeui + 1
        t2 = t2 - stimeui + 1

        ! write(*,'(i0" | ",i0," | ",i0, " " ,f15.10, " ",i0, " : ",f10.4)') tid, itime, t1, newtime_ind(itime), t2, ditime

! print*,"HERE 1"
        mask = .false.
! print*,"HERE 1.1"
! print*,"HERE 2"
        do ii=1,nlon
          do jj=1,nlat
! print*,ii,jj
            good = u2d(ii,jj,t1).NE.fFillValue .OR. v2d(ii,jj,t1).NE.fFillValue .OR. &
                   u2d(ii,jj,t2).NE.fFillValue .OR. v2d(ii,jj,t2).NE.fFillValue

            if ( good ) then
              u2d_tmp(ii,jj) = (1.-ditime)*u2d(ii,jj,t1) + (ditime)*u2d(ii,jj,t2)
              v2d_tmp(ii,jj) = (1.-ditime)*v2d(ii,jj,t1) + (ditime)*v2d(ii,jj,t2)
            else
              u2d_tmp(ii,jj) = fFillValue
              v2d_tmp(ii,jj) = fFillValue
              mask(ii,jj) = .true.
            endif

          enddo
        enddo

! print*,"HERE XXXX"
! !$omp critical
        call get_cell_hor ( points(ipoint,:), lon2d, lat2d, ij, mask                    )
        call interpolate2d( points(ipoint,:), u2d_tmp, v2d_tmp, ij, lon2d, lat2d, pu, pv)
        call locate       ( points(ipoint,:), pu, pv, 0., timestep                      )

! print*,tid,ipoint,points(ipoint,:)

! !$omp end single
! !$omp end critical

! print*,tid,ipoint

! !$omp critical
        ! write(*,'(i0," ",i0," ", 8i12," | ", 3f10.2," | ", 2f7.2)') ipoint, itime, ij, points(ipoint,:), pu, pv
        if( any(points(ipoint,:).eq.fFillValue) ) points(ipoint,:) = fFillValue
        result(:,ipoint,itime) = points(ipoint,:)

! !$omp barrier

! !$omp end critical

      enddo inner

! else
!     !   call get_var_xyzt(ncid_in,"u",u)
!     !   call get_var_xyzt(ncid_in,"v",v)
!     !   call get_var_xyzt(ncid_in,"w",w)
!     !   call get_var_xyzt(ncid_in,"z",z)

    endif
  enddo outer

! print*,tid,"HERE!"
!$omp end parallel

  print*, "#####################################################################"

  ! do ii=1, npoints
  !   print*,ii,minval(result(1,ii,:),mask=result(1,ii,:).NE.fFillValue), &
  !             maxval(result(1,ii,:),mask=result(1,ii,:).NE.fFillValue), &
  !             minval(result(2,ii,:),mask=result(2,ii,:).NE.fFillValue), &
  !             maxval(result(2,ii,:),mask=result(2,ii,:).NE.fFillValue)
  ! enddo

  ! print*,result
  call save_results(file_out,file_in,ncid_in,result,newtime)

  deallocate( result )
  deallocate( points )
  deallocate( time_in )
  deallocate( newtime_ind      )
  deallocate( newtime          )
  deallocate( newtime_datatime )
  deallocate( lon2d   )
  deallocate( lat2d   )
  deallocate( u2d,v2d )
  deallocate( u2d_tmp,v2d_tmp,mask )

!       if (update) then

!         if (horizontal) then
!           w1 = 0; w2 = 0 ! if it is 2D interpolation
!           z1 = 0; z2 = 0 ! if it is 2D interpolation
!           call get_var_xyzt(ncid_in,"u",u1)
!           call get_var_xyzt(ncid_in,"v",v1)
!           call get_var_xyzt(ncid_in,"u",u2)
!           call get_var_xyzt(ncid_in,"v",v2)
!         else
!           call get_var_xyzt(ncid_in,"u",u1)
!           call get_var_xyzt(ncid_in,"v",v1)
!           call get_var_xyzt(ncid_in,"w",w1)
!           call get_var_xyzt(ncid_in,"z",z1)
!           call get_var_xyzt(ncid_in,"u",u2)
!           call get_var_xyzt(ncid_in,"v",v2)
!           call get_var_xyzt(ncid_in,"w",w2)
!           call get_var_xyzt(ncid_in,"z",z2)
!         endif

!         update = .false.

!       endif
        
!       ditime = newtime_ind(itime)-it   ! 

!       ! Linear inperpolation between src timesteps 
!       ! weights are 0..1 => no need to normalize  (!!!!!!!)
!       if (horizontal) then
!         u = (1.-ditime)*u1 + (ditime)*u2
!         v = (1.-ditime)*v1 + (ditime)*v2
!       else
!         u = (1.-ditime)*u1 + (ditime)*u2
!         v = (1.-ditime)*v1 + (ditime)*v2
!         z = (1.-ditime)*z1 + (ditime)*z2 !
!         w = (1.-ditime)*w1 + (ditime)*w2 !
!       endif

!     else

!       if (horizontal) then
!         call get_var_xyzt(ncid_in,"u",it,u)
!         call get_var_xyzt(ncid_in,"v",it,v)
!       else
!         call get_var_xyzt(ncid_in,"u",it,u)
!         call get_var_xyzt(ncid_in,"v",it,v)
!         call get_var_xyzt(ncid_in,"w",it,w)
!         call get_var_xyzt(ncid_in,"z",it,z)
!       endif

!       ditime = 0

!     endif

! ! print*,' *** 1'

!     ! calculating heights (for z interpolation)
!     ! z = z * rg ! geopotential / g = h

! ! print*,' *** 2'

!     if(itime.eq.1) then
!         do ip = 1, npoints
!           result(:,ip,itime) = points(ip,:)
!           ! print'(a,24x, 2f10.4, f10.2)', newtime_datatime(itime), points(ip,:)
!         end do
!     else

!       if (horizontal) then ! if it's horizontal plane
        
! ! print*,' *** 3'

!         iz = MINLOC(abs(levels-horizontal_level), 1)

!         do ip = 1, npoints

!           call get_cell_hor ( points(ip,:), lon2d, lat2d, ij)
!           call interpolate2d( points(ip,:), u(:,:,iz), v(:,:,iz), ij, &
!                               lon2d, lat2d, pu, pv, pw)
!           call locate       ( points(ip,:), pu, pv, pw, timestep )

!           if( any(points(ip,:).eq.fFillValue) ) cycle

!           result(:,ip,itime) = points(ip,:)

!         end do

!       else ! it is 3D plane

!         do ip = 1, npoints

!           call get_cell_hor ( points(ip,:), lon2d, lat2d, ij)
!           call get_cell_vert( points(ip,:), ij, z, kk )
!           call interpolate3d( points(ip,:), u, v, w, z, ij, kk, &
!                               lon2d, lat2d, levels, pu, pv, pw)
!           call locate       ( points(ip,:), pu, pv, pw, timestep )

!           if( any(points(ip,:).eq.fFillValue) ) cycle

!           ! print'(a,2f7.2,f10.4, 2f10.4, f10.2)', newtime_datatime(itime), pu, pv, pw, points(ip,:)
!           ! print'(a, f7.2," |",2f10.4, f10.2)', newtime_datatime(itime), sqrt(pu**2+pv**2+pw**2), points(ip,:)

!           result(:,ip,itime) = points(ip,:)

!         end do



!       endif

    ! endif


  ! end do

  ! deallocate( u,v,w,z )
  ! if(accuracy) then
  !   deallocate( u1,v1,w1,z1 )
  !   deallocate( u2,v2,w2,z2 )
  ! endif

  ! OUTPUT
  
  ! deallocate( points,lon2d,lat2d,time_in,levels            )
  ! deallocate( points,lon2d,lat2d,time_in            )
  ! deallocate( result,newtime_ind,newtime,newtime_datatime  )

  call close_nc(ncid_in)

contains 

end program lpt
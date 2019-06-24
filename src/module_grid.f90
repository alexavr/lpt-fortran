module module_grid
!===============================================================================
!
! module_grid: Module that provides the grid-related procedures 
!               like finging indexes and interpolation
!               
!               
!===============================================================================

use,intrinsic :: iso_fortran_env,only:real32,real64
use module_globals
implicit none ! religion first

private

public :: haversine 
public :: get_cell_hor 
public :: get_cell_vert
public :: interpolate
public :: locate

    real(kind=real64), parameter :: pi = 3.14159265359, er = 6371.D0

    interface haversine
        procedure dhaversine1d, dhaversine2d
    end interface haversine
    
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine get_cell_hor ( point, lon2d, lat2d, ij)
implicit none ! religion first
real(kind=real32), intent(inout) :: point(3)
real(kind=real32), intent(in) :: lon2d(:,:), lat2d(:,:)
integer, intent(out) :: ij(4,2)
real(kind=real64),dimension(:,:),allocatable :: dist

integer :: xdim, ydim, npoints
integer :: ii
real(kind=real32) :: plon, plat, dist_min
logical :: south, north, east, west
logical :: skip = .false.

    skip = ( count( point.eq.fFillValue ) .NE. 0 ) 

    ij = iFillValue

    if(.not.skip) then

        plon = point(1)
        plat = point(2)

        xdim = ubound(lon2d,1)
        ydim = ubound(lon2d,2)
        npoints = ubound(ij,1)

        allocate( dist(xdim, ydim) ) 

        dist = haversine(lon2d, lat2d, plon, plat)

        do ii = 1, npoints
            ij(ii,:) = minloc(dist)
            ! print*,"ij ", ii, ij(ii,:),dist( ij(ii,1),ij(ii,2) )

            dist( ij(ii,1),ij(ii,2) ) = huge(dist)
        end do

        ! If point reached the border, that loosing it
        if(regional) then
            west  = (count( ij(:,1).EQ.1 ) .NE. 0)
            south = (count( ij(:,2).EQ.1 ) .NE. 0)
            east  = (count( ij(:,1).EQ.xdim ) .NE. 0)
            north = (count( ij(:,2).EQ.ydim ) .NE. 0)
            if( south .or. north .or. east .or. west ) point = fFillValue
        end if

        deallocate (dist)

    end if

end subroutine get_cell_hor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_cell_vert ( point, ij, z, kk)
implicit none ! religion first
    real(kind=real32), intent(inout) :: point(3)
    real(kind=real32), intent(in)    :: z(:,:,:)
    integer, intent(in)    :: ij(:,:)
    integer, intent(out)   :: kk(:,:)
    
    real(kind=real32) :: phgt
    integer :: nkk, nij, nz, k1, k2
    integer :: ii, ik
    logical :: skip = .false.
    logical :: ascending
    
    real(kind=real64),dimension(:,:),allocatable :: zz
  
    skip = ( count( point.eq.fFillValue ) .NE. 0 ) 

    kk = iFillValue

    ! if it was regional gris and we at the border -- skip all to the end
    if(.not.skip) then
        nij = ubound(ij,1)
        nz = ubound(z,3)
        phgt = point(3)
        
        allocate( zz(nij, nz) ) 

        ascending = ( z(1,1,nz) .gt. z(1,1,1) )

        k1 = 1
        k2 = 0
        if(ascending) then
            k1 = 0
            k2 = 1
        end if

        do ii = 1, nij
            zz(ii,:) = z( ij(ii,1),ij(ii,2),: )
            do ik = 1, nz-1
                ! print*, zz(ii,ik+k1), phgt, zz(ii,ik+k2)
                if( phgt .ge. zz(ii,ik+k1) .and. phgt .lt. zz(ii,ik+k2) ) then
                    kk(ii,1) = ik
                    kk(ii,2) = ik+1
                    exit
                end if
            end do

            ! print*,"kk ", ii, kk(ii,1), kk(ii,2)

            if(ascending) then
                if( (phgt .gt. zz(ii,nz)) .OR. (phgt .lt. zz(ii,1))) point = fFillValue
            else
                if( (phgt .gt. zz(ii,1)) .OR. (phgt .lt. zz(ii,nz))) point = fFillValue
            end if

            ! print*,"kk ", ii, kk(ii,1), kk(ii,2)

        end do

    end if



end subroutine get_cell_vert
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interpolate(point, u, v, w, z, ij, kk,       &
                            lon2d, lat2d, levs,         &
                            pu, pv, pw)
implicit none ! religion first
    real(kind=real32), intent(inout) :: point(3)
    real(kind=real32), intent(in)    :: u(:,:,:),v(:,:,:),w(:,:,:),z(:,:,:)
    real(kind=real32), intent(in)    :: lon2d(:,:), lat2d(:,:),levs(:)
    integer, intent(in)    :: ij(:,:),kk(:,:)
    real(kind=real32), intent(out)   :: pu, pv, pw
    logical :: skip = .false., levels
    real(kind=real32) :: plon, plat, phgt
    real(kind=real64) :: dlon, dlat, dhgt
    real(kind=real64),dimension(:),allocatable :: dist, wgt
    integer :: nn, nk, ndist, ik, in, ii
    real(kind=real32) :: dp, dz
    real(kind=real64) :: sum_dist

    skip = ( count( point.eq.fFillValue ) .NE. 0 ) 

    if(.not.skip) then

        plon = point(1)
        plat = point(2)
        phgt = point(3)

        nn = ubound(ij,1)
        nk = ubound(kk,2)
        ndist = nn*nk

        allocate( dist(ndist) )
        allocate(  wgt(ndist) )

        ii = 1
        do ik = 1, nk
        do in = 1, nn
            dlon = haversine(lon2d(ij(in,1),ij(in,2)), plat, plon, plat)*1000.d0 ! in m
            dlat = haversine(plon, lat2d(ij(in,1),ij(in,2)), plon, plat)*1000.d0 ! in m
            dhgt = z(ij(in,1),ij(in,2), kk(in,ik)) - phgt
            dist(ii) = sqrt( dlon**2 + dlat**2 + dhgt**2  )
            ii = ii + 1
        end do
        end do

        sum_dist = sum(dist)
        wgt = (1.d0 - dist/sum_dist)/dble(ndist-1)

        levels = all( levs.ne.fFillValue ) 

        pu = 0.
        pv = 0.
        pw = 0.
        ii = 1
        do ik = 1, nk
        do in = 1, nn
            dp = 1.
            dz = 1.
            if(levels) then ! if we have levels => it is reanalisis and have w[Pa/m]
                dp = abs( levs(kk(in,1)) - levs(kk(in,2)) )*100.
                dz = abs( z(ij(in,1),ij(in,2),kk(in,1)) - z(ij(in,1),ij(in,2),kk(in,2)) )
            end if

            pu = pu + wgt(ii)*u(  ij(in,1),ij(in,2),kk(in,ik)  )
            pv = pv + wgt(ii)*v(  ij(in,1),ij(in,2),kk(in,ik)  )
            pw = pw + wgt(ii)*w(  ij(in,1),ij(in,2),kk(in,ik)  )*dz/dp

            ii = ii + 1
        end do
        end do

        deallocate( dist )
        deallocate( wgt  )

    end if

end subroutine interpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dt - in min
subroutine locate(point, pu, pv, pw, dt)
implicit none ! religion first
    real(kind=real32), intent(inout) :: point(3)
    real(kind=real32), intent(in)    :: pu, pv, pw 
    real(kind=real32), intent(in)    :: dt 
    logical :: skip = .false.
    real(kind=real32) :: plon, plat, phgt
    real(kind=real64) :: R

    skip = ( count( point.eq.fFillValue ) .NE. 0 ) 

    if(.not.skip) then

        plon = fFillValue
        plat = fFillValue
        phgt = fFillValue

        R = er*1000.d0

        plon = point(1) + pu*(dt*60.d0)/R*180.d0/pi*cos(point(2)*pi/180.d0)
        plat = point(2) + pv*(dt*60.d0)/R*180.d0/pi
        phgt = point(3) + pw*(dt*60.d0)


        ! FIX THIS!!!!
        if (plat.ge.90) then 
            plat = 180. - plat
            plon = plon - 180
        else if (plat .le. -90) then
            plat = -180 - plat
            plon = plon - 180
        end if

        if (plon .ge. 180) then 
            plon = plon - 360
        else if (plon .le. -180) then
            plon = 360 - plon
        end if

        point(1) = plon
        point(2) = plat
        point(3) = phgt

    end if

end subroutine locate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the great circle distance between two points 
! on the earth (specified in decimal degrees)
! in km
function dhaversine1d(lon1, lat1, lon2, lat2)
implicit none ! religion first
real(kind=real32), intent(in) :: lon1, lat1, lon2, lat2
real(kind=real64)             :: rlon1, rlat1, rlon2, rlat2, dlon, dlat, a, c
real(kind=real64)             :: dhaversine1d

    rlon1 = dble(lon1)*pi/180.D0 
    rlat1 = dble(lat1)*pi/180.D0 
    rlon2 = dble(lon2)*pi/180.D0 
    rlat2 = dble(lat2)*pi/180.D0 

    dlon = rlon2 - rlon1 
    dlat = rlat2 - rlat1 

    a = sin(dlat/2.D0)**2 + cos(rlat1) * cos(rlat2) * sin(dlon/2.D0)**2
    c = 2.D0 * asin(sqrt(a))

    dhaversine1d = c * er 

end function dhaversine1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the great circle distance on the earth (specified in decimal degrees)
! form a point with reapect to 2d array of coordinates 
! in km
function dhaversine2d(lon2d, lat2d, plon, plat)
implicit none ! religion first
real(kind=real32), intent(in)                 :: lon2d(:,:), lat2d(:,:), plon, plat
real(kind=real64)                             :: rplon, rplat
real(kind=real64),dimension(:,:),allocatable  :: rlon2d, rlat2d, dlon, dlat, a, c
real(kind=real64),dimension(:,:),allocatable  :: dhaversine2d

integer :: xdim, ydim
    
    xdim = UBOUND(lon2d,1)
    ydim = UBOUND(lon2d,2)

    allocate( rlon2d(xdim,xdim),rlat2d(xdim,xdim) ) 
    allocate( dlon(xdim,xdim), dlat(xdim,xdim) ) 
    allocate( a(xdim,xdim), c(xdim,xdim) )
    allocate( dhaversine2d(xdim,xdim) )

    rlon2d = dble(lon2d)*pi/180.D0 
    rlat2d = dble(lat2d)*pi/180.D0 
    rplon  = dble(plon)*pi/180.D0 
    rplat  = dble(plat)*pi/180.D0 

    dlon = rlon2d - rplon 
    dlat = rlat2d - rplat 

    a = sin(dlat/2.D0)**2 + cos(rlat2d) * cos(rplat) * sin(dlon/2.D0)**2
    c = 2.D0 * asin(sqrt(a))

    dhaversine2d = c * er 
    deallocate( rlon2d, rlat2d, dlon, dlat, a, c )

end function dhaversine2d

end module module_grid
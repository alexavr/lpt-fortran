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

public :: create_grid 
public :: haversine 
public :: get_cell_hor 
public :: get_cell_vert
public :: interpolate3d
public :: interpolate2d
public :: locate

    real(kind=real64), parameter :: pi = 3.14159265359d0, er = 6371.D0

    interface haversine
        procedure dhaversine1d, dhaversine2d
    end interface haversine
    
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! создаем коордтнатную сетку индексов, если нет своей
subroutine create_grid ( x2d, y2d)
implicit none ! religion first
real(kind=real32), intent(inout)    :: x2d(:,:), y2d(:,:)
integer :: ii,jj, xdim, ydim

xdim = ubound(x2d,1)
ydim = ubound(x2d,2)

do ii = 1, xdim
do jj = 1, ydim

    x2d(ii,jj) = ii
    y2d(ii,jj) = jj

enddo
enddo

end subroutine create_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! получаем координаты 4-х окружающих точек
subroutine get_cell_hor ( point1d, x2d, y2d, ij, mask)
implicit none ! religion first
real(kind=real32), intent(inout) :: point1d(3)
real(kind=real32), intent(in)    :: x2d(:,:), y2d(:,:)
logical,           intent(in)    :: mask(:,:)
integer,           intent(out)   :: ij(4,2)
real(kind=real64),dimension(:,:),allocatable :: dist
real(kind=real64)                :: grid(3,3,5) = dFillValue

integer :: xdim, ydim, npoints
integer :: ii
real(kind=real32) :: pointx, pointy, dist_min
logical :: south, north, east, west
logical :: skip = .false.

! !$omp critical

    skip = ( count( point1d.eq.fFillValue ) .NE. 0 )


    ij = iFillValue

! print*,' *** 4.1', skip

    if(.not.skip) then


        pointx = point1d(1)
        pointy = point1d(2)

        xdim = ubound(x2d,1)
        ydim = ubound(x2d,2)
        npoints = ubound(ij,1)

        allocate( dist(xdim, ydim) ) 

        if (cartesian_grid) then ! Если сетка равных расстояний, то считаем обычным Пифагором

            ! ПИФАГОР
            dist = sqrt( (x2d-pointx)**2 + (y2d - pointy)**2)
            call get_grid(dist, x2d, y2d, pointx, pointy, grid)

        else ! иначе считаем на сфере честно

            dist = haversine(x2d, y2d, pointx, pointy)    
            call get_grid(dist, x2d, y2d, pointx, pointy, grid)

        endif

        if(cell_detector.eq.0) then
            ! print*,"Simple scheme"
            call closest_distance(grid, pointx, pointy, ij)
        else if(cell_detector.eq.1) then
            call triangle_method(grid, pointx, pointy, ij)
        else
            write(*,'("get_cell_hor: Can not recognize the cell_detector option eq ",i)') cell_detector
            stop
        end if

        skip = ( count( ij.eq.iFillValue ) .NE. 0 ) 


        if(.not.global .OR. .not.periodic) then ! If point reached the border, that loosing it
            west  = (count( ij(:,1).EQ.1 ) .NE. 0)
            south = (count( ij(:,2).EQ.1 ) .NE. 0)
            east  = (count( ij(:,1).EQ.xdim ) .NE. 0)
            north = (count( ij(:,2).EQ.ydim ) .NE. 0)
            if( south .or. north .or. east .or. west ) point1d = fFillValue
        end if

        skip = ( mask( ij(1,1),ij(1,2)) .OR. mask( ij(2,1),ij(2,2)) .OR. &
                 mask( ij(3,1),ij(3,2)) .OR. mask( ij(4,1),ij(4,2))  )


            ! print*,skip
        if(skip) point1d = fFillValue ! If point reached the missing area
        
        deallocate (dist)


    end if

! !$omp end critical

! print*,' *** 4.6', point

end subroutine get_cell_hor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Searching for grid indexes around the point 
! by using the triangle method: the sum of trangles areas  
! (2 grid to the point) will be equal the gridcell area only if 
! the point is inside the cell
! Issues: 
! - projection uncertanties make it usless to find 
!           the actual equlity (even with certain epsilon).
!           So for now I'm using min(sum(trianlges)) amoung 
!           all 4 cells. Tested numerious times - works very well.
! - Only 4 cells and 4 points acounted yet.
! 
subroutine triangle_method(grid, plon, plat, ij)
implicit none ! religion first
    real(kind=real64), intent(in) :: grid(:,:,:)
    real(kind=real32), intent(in) :: plon, plat
    integer          , intent(out):: ij(:,:)
    logical :: skip = .false.
    integer,parameter :: ncells = 4, npoints = 4 
    real(kind=real32) :: cell(ncells,npoints,5)
    real(kind=real32) :: segm(ncells,npoints,2)
    real(kind=real64) :: darea(ncells)

    integer :: ic, ip, ic_min(1)

    ij = iFillValue

    skip = ( count( grid.eq.dFillValue ) .NE. 0 ) 
    
    if(.not.skip) then
        ! initializing the traversal order (+1 comes from python indexing)
        segm(1,1,:) = (/1,1/)+1
        segm(1,2,:) = (/1,0/)+1
        segm(1,3,:) = (/2,0/)+1
        segm(1,4,:) = (/2,1/)+1
        segm(2,1,:) = (/1,1/)+1
        segm(2,2,:) = (/2,1/)+1
        segm(2,3,:) = (/2,2/)+1
        segm(2,4,:) = (/1,2/)+1
        segm(3,1,:) = (/1,1/)+1
        segm(3,2,:) = (/1,0/)+1
        segm(3,3,:) = (/0,0/)+1
        segm(3,4,:) = (/0,1/)+1
        segm(4,1,:) = (/1,1/)+1
        segm(4,2,:) = (/1,2/)+1
        segm(4,3,:) = (/0,2/)+1
        segm(4,4,:) = (/0,1/)+1

        do ic = 1, ncells
        do ip = 1, npoints
            cell(ic,ip,:) = grid( segm(ic,ip,1 ), segm(ic,ip,2 ), :)
        end do
        end do

        do ic = 1, ncells

            if (cartesian_grid) then
                darea(ic) = darea_triangle_cartesian( plon, plat, cell(ic,:,3), cell(ic,:,4) )
            else
                darea(ic) = darea_triangle( plon, plat, cell(ic,:,3), cell(ic,:,4) )
            endif

        end do

        ic_min = minloc(darea)
        ij = cell(ic_min(1),:,1:2)

    end if

end subroutine triangle_method
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=real64) function darea_triangle( plon, plat, lons1d, lats1d)
implicit none ! religion first
    real(kind=real32), intent(in) :: plon, plat
    real(kind=real32), intent(in) :: lons1d(:), lats1d(:)

    real(kind=real32) :: x0, x1, x2, x3, x4 
    real(kind=real32) :: y0, y1, y2, y3, y4 
    real(kind=real64) :: a_x, a_y, b_x, b_y, c_x, c_y, d_x, d_y
    real(kind=real64) :: ao_x, ao_y, bo_x, bo_y, co_x, co_y, do_x, do_y
    real(kind=real64) :: Sabc,Scda,Saob,Sboc,Scod,Sdoa,Sabcd,Ssum

    x0 = plon
    x1 = lons1d(1)
    x2 = lons1d(2)
    x3 = lons1d(3)
    x4 = lons1d(4)
    y0 = plat
    y1 = lats1d(1)
    y2 = lats1d(2)
    y3 = lats1d(3)
    y4 = lats1d(4)

    a_x = haversine(x1, y1, x2, y1) 
    a_y = haversine(x1, y1, x1, y2) 
    b_x = haversine(x2, y2, x3, y2) 
    b_y = haversine(x2, y2, x2, y3) 
    c_x = haversine(x3, y3, x4, y3) 
    c_y = haversine(x3, y3, x3, y4) 
    d_x = haversine(x4, y4, x1, y4) 
    d_y = haversine(x4, y4, x4, y1) 

    ao_x = haversine(x1, y1, x0, y1)
    ao_y = haversine(x1, y1, x1, y0)
    bo_x = haversine(x2, y2, x0, y2)
    bo_y = haversine(x2, y2, x2, y0)
    co_x = haversine(x3, y3, x0, y3)
    co_y = haversine(x3, y3, x3, y0)
    do_x = haversine(x4, y4, x0, y4)
    do_y = haversine(x4, y4, x4, y0)

    Sabc = abs(0.5d0*(a_x*b_y-a_y*b_x))
    Scda = abs(0.5d0*(c_x*d_y-c_y*d_x))
    Sabcd = Sabc + Scda

    Saob = abs(0.5d0*(a_x*ao_y-a_y*ao_x))
    Sboc = abs(0.5d0*(b_x*bo_y-b_y*bo_x))
    Scod = abs(0.5d0*(c_x*co_y-c_y*co_x))
    Sdoa = abs(0.5d0*(d_x*do_y-d_y*do_x))
    Ssum = Saob+Sboc+Scod+Sdoa

    darea_triangle = abs( Sabcd-Ssum )

end function darea_triangle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! для равномерной по расстоянию сетки
real(kind=real64) function darea_triangle_cartesian( plon, plat, lons1d, lats1d)
implicit none ! religion first
    real(kind=real32), intent(in) :: plon, plat
    real(kind=real32), intent(in) :: lons1d(:), lats1d(:)

    real(kind=real32) :: x0, x1, x2, x3, x4 
    real(kind=real32) :: y0, y1, y2, y3, y4 
    real(kind=real64) :: a_x, a_y, b_x, b_y, c_x, c_y, d_x, d_y
    real(kind=real64) :: ao_x, ao_y, bo_x, bo_y, co_x, co_y, do_x, do_y
    real(kind=real64) :: Sabc,Scda,Saob,Sboc,Scod,Sdoa,Sabcd,Ssum

    x0 = plon
    x1 = lons1d(1)
    x2 = lons1d(2)
    x3 = lons1d(3)
    x4 = lons1d(4)
    y0 = plat
    y1 = lats1d(1)
    y2 = lats1d(2)
    y3 = lats1d(3)
    y4 = lats1d(4)

    a_x = abs(x1 - x2) 
    a_y = abs(y1 - y2) 
    b_x = abs(x2 - x3) 
    b_y = abs(y2 - y3) 
    c_x = abs(x3 - x4) 
    c_y = abs(y3 - y4) 
    d_x = abs(x4 - x1) 
    d_y = abs(y4 - y1) 

    ao_x = abs(x1 - x0)
    ao_y = abs(y1 - y0)
    bo_x = abs(x2 - x0)
    bo_y = abs(y2 - y0)
    co_x = abs(x3 - x0)
    co_y = abs(y3 - y0)
    do_x = abs(x4 - x0)
    do_y = abs(y4 - y0)

    Sabc = abs(0.5d0*(a_x*b_y-a_y*b_x))
    Scda = abs(0.5d0*(c_x*d_y-c_y*d_x))
    Sabcd = Sabc + Scda

    Saob = abs(0.5d0*(a_x*ao_y-a_y*ao_x))
    Sboc = abs(0.5d0*(b_x*bo_y-b_y*bo_x))
    Scod = abs(0.5d0*(c_x*co_y-c_y*co_x))
    Sdoa = abs(0.5d0*(d_x*do_y-d_y*do_x))
    Ssum = Saob+Sboc+Scod+Sdoa

    darea_triangle_cartesian = abs( Sabcd-Ssum )

end function darea_triangle_cartesian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gets GRID data as output and finds npoints 
! closest points. The npoints determines as 
! the dimention of the ij array. 
! 
! grid(:,:,1) - ii (global)
! grid(:,:,2) - jj (global)
! grid(:,:,3) - lon 
! grid(:,:,4) - lat 
! grid(:,:,5) - dist 
! 
subroutine closest_distance(grid, plon, plat, ij)
implicit none ! religion first
    real(kind=real64), intent(inout) :: grid(:,:,:)
    real(kind=real32), intent(in) :: plon, plat
    integer          , intent(out):: ij(:,:)

    integer :: npoints
    logical :: skip = .false.

    integer :: ii, ij_tmp(2)

    skip = ( count( grid.eq.dFillValue ) .NE. 0 ) 
    
    if(.not.skip) then
        npoints = ubound(ij,1)

        do ii = 1, npoints
            ij_tmp = minloc(grid(:,:,5))
            ij(ii,1) = nint( grid(ij_tmp(1),ij_tmp(2),1) )
            ij(ii,2) = nint( grid(ij_tmp(1),ij_tmp(2),2) )
            ! print*,"ij ", ii, ij(ii,:), grid( ij_tmp(1),ij_tmp(2), 5 )
            grid( ij_tmp(1),ij_tmp(2), 5 ) = huge(grid)
        end do


    end if


end subroutine closest_distance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! grid(:,:,1) - ii (global)
! grid(:,:,2) - jj (global)
! grid(:,:,3) - lon 
! grid(:,:,4) - lat 
! grid(:,:,5) - dist 
subroutine get_grid(dist, lon2d, lat2d, plon, plat, grid)
implicit none ! religion first
    real(kind=real64), intent(in) :: dist(:,:)
    real(kind=real32), intent(in) :: lon2d(:,:),lat2d(:,:)
    real(kind=real32), intent(in) :: plon, plat
    real(kind=real64), intent(out) :: grid(3,3,5)

    integer :: xdim, ydim
    integer :: ii, jj
    integer :: iii, iiip1, iiim1, jjj, jjjp1, jjjm1
    integer :: iii_jjjm1,iii_jjjp1,iiip1_jjjm1,iiim1_jjjm1,iiip1_jjjp1,iiim1_jjjp1
    integer :: tmp2(2), tmp1(1)
    real(kind=real32) :: lon2

    grid = dFillValue

    xdim = ubound(dist,1)
    ydim = ubound(dist,2)

    tmp2 = minloc(dist)
    ii = tmp2(1)
    jj = tmp2(2)

! print*,' *** 4.3.1', xdim, ydim, ii, jj
! print*,' *** 4.3.1.1', minval(dist), maxval(dist)


    ! if its the N/S Pole - we need to point the y-direction
    if(global) then
        if(lat2d(ii, jj) == 90 .or. lat2d(ii, jj) == -90) then
            tmp1 = minloc( abs(lon2d(:, jj)-plon) )
            ii = tmp1(1)
        end if
    end if

    iii   = ii
    iiip1 = ii+1
    iiim1 = ii-1
    jjj   = jj
    jjjp1 = jj+1
    jjjm1 = jj-1

! print*,' *** 4.3.2', iii, iiip1, iiim1, jjj, jjjp1, jjjm1


    if ( global ) then
        ! глобальная сетка в координатах щирота/долгота
        ! подразумевают периодические условия по умолчанию

        if(ii == xdim) iiip1 = 1
        if(ii == 1   ) iiim1 = xdim

        iii_jjjm1 = iii
        iii_jjjp1 = iii
        iiip1_jjjm1 = iiip1
        iiim1_jjjm1 = iiim1
        iiip1_jjjp1 = iiip1
        iiim1_jjjp1 = iiim1

        ! handling the N/S Poles
        if(jj == 1) then
            if (lon2d(ii,jj+1) > 0) then
                lon2 = lon2d(ii,jj+1)-180.
            else
                lon2 = lon2d(ii,jj+1)+180.
            end if

            jjjm1 = jj+1

            tmp1 = minloc( abs(lon2d(:, jjjm1)-lon2) )
            iii_jjjm1 = tmp1(1)
            iiip1_jjjm1 = iii_jjjm1 - 1
            iiim1_jjjm1 = iii_jjjm1 + 1

        end if

        if(jj == ydim) then

            if (lon2d(ii,jj-1) >  0) then
                lon2 = lon2d(ii,jj-1)-180.
            else
                lon2 = lon2d(ii,jj-1)+180.
            end if

            jjjp1 = jj-1

            tmp1 = minloc( abs(lon2d(:, jjjp1)-lon2) )
            iii_jjjp1 = tmp1(1)
            iiip1_jjjp1 = iii_jjjp1 - 1
            iiim1_jjjp1 = iii_jjjp1 + 1
        end if

        grid(3,1,:) = (/dble(iiim1_jjjp1) , dble(jjjp1), dble(lon2d(iiim1_jjjp1,jjjp1)), dble(lat2d(iiim1_jjjp1,jjjp1)), dist(iiim1_jjjp1,jjjp1)/)
        grid(3,2,:) = (/dble(iii_jjjp1  ) , dble(jjjp1), dble(lon2d(iii_jjjp1  ,jjjp1)), dble(lat2d(iii_jjjp1  ,jjjp1)), dist(iii_jjjp1  ,jjjp1)/)
        grid(3,3,:) = (/dble(iiip1_jjjp1) , dble(jjjp1), dble(lon2d(iiip1_jjjp1,jjjp1)), dble(lat2d(iiip1_jjjp1,jjjp1)), dist(iiip1_jjjp1,jjjp1)/)
        grid(2,1,:) = (/dble(iiim1      ) , dble(jjj  ), dble(lon2d(iiim1      ,jjj  )), dble(lat2d(iiim1      ,jjj  )), dist(iiim1      ,jjj  )/)
        grid(2,2,:) = (/dble(iii        ) , dble(jjj  ), dble(lon2d(iii        ,jjj  )), dble(lat2d(iii        ,jjj  )), dist(iii        ,jjj  )/)
        grid(2,3,:) = (/dble(iiip1      ) , dble(jjj  ), dble(lon2d(iiip1      ,jjj  )), dble(lat2d(iiip1      ,jjj  )), dist(iiip1      ,jjj  )/)
        grid(1,1,:) = (/dble(iiim1_jjjm1) , dble(jjjm1), dble(lon2d(iiim1_jjjm1,jjjm1)), dble(lat2d(iiim1_jjjm1,jjjm1)), dist(iiim1_jjjm1,jjjm1)/)
        grid(1,2,:) = (/dble(iii_jjjm1  ) , dble(jjjm1), dble(lon2d(iii_jjjm1  ,jjjm1)), dble(lat2d(iii_jjjm1  ,jjjm1)), dist(iii_jjjm1  ,jjjm1)/)
        grid(1,3,:) = (/dble(iiip1_jjjm1) , dble(jjjm1), dble(lon2d(iiip1_jjjm1,jjjm1)), dble(lat2d(iiip1_jjjm1,jjjm1)), dist(iiip1_jjjm1,jjjm1)/)

    else if( .not.global .AND. .not.periodic ) then
        ! региональная модель со свободными границами (все из домена утекает вникуда)

       ! If closest grid is at the border, than loose the particle else do:
        if (.not. (((ii == 1) .or. (ii == xdim) .or. (jj == 1) .or. (jj == ydim)))) then
            grid(3,1,:) = (/ dble(iiim1), dble(jjj+1), dble(lon2d(iiim1,jjj+1)), dble(lat2d(iiim1,jjj+1)), dist(iiim1,jjj+1) /)
            grid(3,2,:) = (/ dble(iii  ), dble(jjj+1), dble(lon2d(iii  ,jjj+1)), dble(lat2d(iii  ,jjj+1)), dist(iii  ,jjj+1) /)
            grid(3,3,:) = (/ dble(iiip1), dble(jjj+1), dble(lon2d(iiip1,jjj+1)), dble(lat2d(iiip1,jjj+1)), dist(iiip1,jjj+1) /)
            grid(2,1,:) = (/ dble(iiim1), dble(jjj  ), dble(lon2d(iiim1,jjj  )), dble(lat2d(iiim1,jjj  )), dist(iiim1,jjj  ) /)
            grid(2,2,:) = (/ dble(iii  ), dble(jjj  ), dble(lon2d(iii  ,jjj  )), dble(lat2d(iii  ,jjj  )), dist(iii  ,jjj  ) /)
            grid(2,3,:) = (/ dble(iiip1), dble(jjj  ), dble(lon2d(iiip1,jjj  )), dble(lat2d(iiip1,jjj  )), dist(iiip1,jjj  ) /)
            grid(1,1,:) = (/ dble(iiim1), dble(jjj-1), dble(lon2d(iiim1,jjj-1)), dble(lat2d(iiim1,jjj-1)), dist(iiim1,jjj-1) /)
            grid(1,2,:) = (/ dble(iii  ), dble(jjj-1), dble(lon2d(iii  ,jjj-1)), dble(lat2d(iii  ,jjj-1)), dist(iii  ,jjj-1) /)
            grid(1,3,:) = (/ dble(iiip1), dble(jjj-1), dble(lon2d(iiip1,jjj-1)), dble(lat2d(iiip1,jjj-1)), dist(iiip1,jjj-1) /)
        end if

    else if( .not.global .AND. periodic ) then

        if(ii == xdim) iiim1 = 1; iii = 2; iiip1 = 3
        if(ii == 1   ) iiim1 = xdim-1; iii = xdim-1; iiip1 = xdim

        iii_jjjm1 = iii
        iii_jjjp1 = iii
        iiip1_jjjm1 = iiip1
        iiim1_jjjm1 = iiim1
        iiip1_jjjp1 = iiip1
        iiim1_jjjp1 = iiim1

        ! handling the N/S Poles
        if(jj == ydim) jjjm1 = 1; jjj = 2; jjjp1 = 3
        if(jj == 1   ) jjjm1 = ydim-1; jjj = ydim-1; jjjp1 = ydim
        ! if(jj == 1) then
        !     if (lon2d(ii,jj+1) > 0) then
        !         lon2 = lon2d(ii,jj+1)-180.
        !     else
        !         lon2 = lon2d(ii,jj+1)+180.
        !     end if

        !     jjjm1 = jj+1

        !     tmp1 = minloc( abs(lon2d(:, jjjm1)-lon2) )
        !     iii_jjjm1 = tmp1(1)
        !     iiip1_jjjm1 = iii_jjjm1 - 1
        !     iiim1_jjjm1 = iii_jjjm1 + 1

        ! end if

        ! if(jj == ydim) then

        !     if (lon2d(ii,jj-1) >  0) then
        !         lon2 = lon2d(ii,jj-1)-180.
        !     else
        !         lon2 = lon2d(ii,jj-1)+180.
        !     end if

        !     jjjp1 = jj-1

        !     tmp1 = minloc( abs(lon2d(:, jjjp1)-lon2) )
        !     iii_jjjp1 = tmp1(1)
        !     iiip1_jjjp1 = iii_jjjp1 - 1
        !     iiim1_jjjp1 = iii_jjjp1 + 1
        ! end if

        grid(3,1,:) = (/dble(iiim1_jjjp1) , dble(jjjp1), dble(lon2d(iiim1_jjjp1,jjjp1)), dble(lat2d(iiim1_jjjp1,jjjp1)), dist(iiim1_jjjp1,jjjp1)/)
        grid(3,2,:) = (/dble(iii_jjjp1  ) , dble(jjjp1), dble(lon2d(iii_jjjp1  ,jjjp1)), dble(lat2d(iii_jjjp1  ,jjjp1)), dist(iii_jjjp1  ,jjjp1)/)
        grid(3,3,:) = (/dble(iiip1_jjjp1) , dble(jjjp1), dble(lon2d(iiip1_jjjp1,jjjp1)), dble(lat2d(iiip1_jjjp1,jjjp1)), dist(iiip1_jjjp1,jjjp1)/)
        grid(2,1,:) = (/dble(iiim1      ) , dble(jjj  ), dble(lon2d(iiim1      ,jjj  )), dble(lat2d(iiim1      ,jjj  )), dist(iiim1      ,jjj  )/)
        grid(2,2,:) = (/dble(iii        ) , dble(jjj  ), dble(lon2d(iii        ,jjj  )), dble(lat2d(iii        ,jjj  )), dist(iii        ,jjj  )/)
        grid(2,3,:) = (/dble(iiip1      ) , dble(jjj  ), dble(lon2d(iiip1      ,jjj  )), dble(lat2d(iiip1      ,jjj  )), dist(iiip1      ,jjj  )/)
        grid(1,1,:) = (/dble(iiim1_jjjm1) , dble(jjjm1), dble(lon2d(iiim1_jjjm1,jjjm1)), dble(lat2d(iiim1_jjjm1,jjjm1)), dist(iiim1_jjjm1,jjjm1)/)
        grid(1,2,:) = (/dble(iii_jjjm1  ) , dble(jjjm1), dble(lon2d(iii_jjjm1  ,jjjm1)), dble(lat2d(iii_jjjm1  ,jjjm1)), dist(iii_jjjm1  ,jjjm1)/)
        grid(1,3,:) = (/dble(iiip1_jjjm1) , dble(jjjm1), dble(lon2d(iiip1_jjjm1,jjjm1)), dble(lat2d(iiip1_jjjm1,jjjm1)), dist(iiip1_jjjm1,jjjm1)/)


    else
        print*,"get_cell_hor -> get_grid: ERROR!"
        print*,"                          Configuration in namelist global=",global," AND periodic="
        print*,"                          global=",global," AND periodic=",periodic
        print*,"                          is unacceptable"
        print*,"                          STOP"
        stop

    end if

end subroutine get_grid
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
!   ADD DESCRIPTION 
subroutine interpolate3d(point, u, v, w, z, ij, kk,       &
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
    real(kind=real32) :: dp, dz, dx, dy
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
            
            if (cartesian_grid) then
            
                dlon = abs(plon - lon2d(ij(in,1),ij(in,2)))*dx ! haversine(lon2d(ij(in,1),ij(in,2)), plat, plon, plat)*1000.d0 ! in m
                dlat = abs(plat - lat2d(ij(in,1),ij(in,2)))*dy ! haversine(plon, lat2d(ij(in,1),ij(in,2)), plon, plat)*1000.d0 ! in m

            else

                dlon = haversine(lon2d(ij(in,1),ij(in,2)), plat, plon, plat)*1000.d0 ! in m
                dlat = haversine(plon, lat2d(ij(in,1),ij(in,2)), plon, plat)*1000.d0 ! in m

            endif
            
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
            if(levels) then ! if we have levels => it is pressure|hight levels with w[Pa/m] or w[m/s]
                dp = abs( levs(kk(in,1)) - levs(kk(in,2)) )
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

end subroutine interpolate3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Имея координаты окружающих точек, интеполируем скорость искомой точки
subroutine interpolate2d(point1d, u, v, ij, &
                            lon2d, lat2d, &
                            pu, pv)
implicit none ! religion first
    real(kind=real32), intent(in)   :: point1d(:)
    real(kind=real32), intent(in)   :: u(:,:),v(:,:)
    real(kind=real32), intent(in)   :: lon2d(:,:), lat2d(:,:)
    integer, intent(in)    :: ij(:,:)
    real(kind=real32), intent(out)  :: pu, pv
    logical :: skip = .false.
    real(kind=real32) :: plon, plat, phgt
    real(kind=real64) :: dlon, dlat, dhgt
    real(kind=real64),dimension(:),allocatable :: dist, wgt
    integer :: nn, nk, ndist, ik, ind 
    real(kind=real32) :: dp, dz
    real(kind=real64) :: sum_dist

    skip = ( count( point1d.eq.fFillValue ) .NE. 0 )
 
    if(.not.skip) then

        plon = point1d(1)
        plat = point1d(2)
        phgt = point1d(3)

        nn = ubound(ij,1)
        ndist = nn

        allocate( dist(ndist) )
        allocate(  wgt(ndist) )

        do ind = 1, nn

            if (cartesian_grid) then
            
                dlon = plon - lon2d(ij(ind,1),ij(ind,2)) ! безразмерные
                dlat = plat - lat2d(ij(ind,1),ij(ind,2)) ! безразмерные
            
            else
            
                dlon = haversine(lon2d(ij(ind,1),ij(ind,2)), plat, plon, plat)*1000.d0 ! in m
                dlat = haversine(plon, lat2d(ij(ind,1),ij(ind,2)), plon, plat)*1000.d0 ! in m
            
            endif

            dist(ind) = sqrt( dlon**2 + dlat**2  )

        end do

        sum_dist = sum(dist)
        wgt = (1.d0 - dist/sum_dist)/dble(ndist-1)

        pu = 0.
        pv = 0.

        do ind = 1, nn
            if ( u( ij(ind,1),ij(ind,2) ).EQ.fFillValue .OR. v( ij(ind,1),ij(ind,2) ).EQ.fFillValue ) then
                pu = fFillValue
                pv = fFillValue
            else
                pu = pu + wgt(ind)*u( ij(ind,1),ij(ind,2) )
                pv = pv + wgt(ind)*v( ij(ind,1),ij(ind,2) )
            endif

           ! print*, count( u==0 ) ! write(*,'( i4, " | pv:",f4.2," | wgt:",f4.2," |  dist:",f4.2," | u:" ,f5.2," | xy:",2i4 )') ii,pv,wgt(ii),dist(ii),u( ij(in,1),ij(in,2) ), ij(in,1), ij(in,2)

        end do

        deallocate( dist )
        deallocate( wgt  )

    end if

end subroutine interpolate2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gets Earth radius in km 
pure elemental real(kind=real32) function earth_radius(lat)
implicit none
    real(kind=real32), intent(in)    :: lat 
    real(kind=real32) :: er_a,er_b

    er_a = 6378.1370 ! 6370.000 # The Earth's equatorial radius
    er_b = 6356.7523 ! 6370.000 # The Earth's polar radius
        
    earth_radius = sqrt( ( (er_a**2*cos(lat))**2 + &
        (er_b**2*sin(lat))**2 ) / ( (er_a*cos(lat))**2 + &
        (er_b*sin(lat))**2 ))

end function earth_radius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Поиск координат частицы исходя из прошых координат и скорости точки
! dt - in min
! nlon, nlat - для ideal_case (чтобы выкидывать точки при достижении границ домена)
!              ПЕРИОДИЧЕСКИЕ ГРАНИЦЫ В ИДЕАЛИЗИРОВАННЫХ СЛУЧЯХ ПОКА НЕ РАССМАТРИВАЕМ!
subroutine locate(point1d, pu, pv, pw, dt, dx, dy, dz, nlon, nlat, nlevel)
implicit none ! religion first
    real(kind=real32), intent(inout) :: point1d(3)
    real(kind=real32), intent(in)    :: pu, pv, pw 
    real(kind=real32), intent(in)    :: dt, dx, dy, dz
    integer,           intent(in)    :: nlon, nlat, nlevel
    logical :: skip = .false.
    real(kind=real32) :: plon, plat, phgt
    real(kind=real64) :: R

    skip = ( count( point1d.eq.fFillValue ) .NE. 0 )

    if(.not.skip) then


        plon = fFillValue
        plat = fFillValue
        phgt = fFillValue

        if (cartesian_grid) then


            ! знаем прошлые координаты, получаем скорости в метрах -> переводим в расстояние
            !       и переводим в градусы (поскольку исходящие данные в градусах)
            plon = point1d(1) + pu*dt*60/dx
            plat = point1d(2) + pv*dt*60/dy
            phgt = 0 !point1d(3) + pw*dt*60/dz ! *dt*60

            ! if (periodic) then
            !     ! if ( plon <= 1 )    plon = nlon - plon
            !     if ( plon <= 1 )    print*,plon, (nlon - plon)
            !     if ( plat <= 1 )    plat = nlat - plat
            !     ! if ( plon >= nlon ) plon = 1 + plon
            !     if ( plon >= nlon ) print*,plon, (nlon + plon)
            !     if ( plat >= nlat ) plat = 1 + plat
            ! else
            !     if ( plon <= 1 .OR. plon >= nlon ) plon = fFillValue
            !     if ( plat <= 1 .OR. plat >= nlat ) plat = fFillValue
            ! endif

        else

            ! R = 6371.0088 (simple)
            R = earth_radius(point1d(2))
            R = R*1000.d0 ! в метры

            plon = point1d(1) + pu*60.d0*dt/(pi*R/180.d0*cos(point1d(2)*pi/180.d0))
            plat = point1d(2) + pv*60.d0*dt/(pi*R/180.d0)
            phgt = 0 !point1d(3) + pw*60.d0*dt

            ! Forbits to go over abs(360)
            if (abs(plon) > 360) plon = mod(plon,360.)

            ! Если вылетели выше полюса, то переходим в другое полушарие
            if (plat.gt.90) then 
                plat = 180. - plat
                plon = plon - 180.
            else if (plat .lt. -90) then
                plat = -180. - plat
                plon = plon - 180.
            end if

            ! придерживаемся -180..+180 координатной сетки
            if (plon .gt. 180) then 
                plon = plon - 360.
            else if (plon .lt. -180) then
                plon = plon + 360.
            end if

        endif

        ! if (phgt <= 1 .OR. phgt >= nlevel )  phgt = fFillValue !!!!!!!!!!!!!!!!!!!!!!!

        ! write(*,'(6f7.2)') pu,pv,pw,plon,plat,phgt
        point1d(1) = plon
        point1d(2) = plat
        point1d(3) = phgt

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

    dhaversine1d = c * earth_radius( lat2 ) 

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

    allocate( rlon2d(xdim,ydim),rlat2d(xdim,ydim) ) 
    allocate( dlon(xdim,ydim), dlat(xdim,ydim) ) 
    allocate( a(xdim,ydim), c(xdim,ydim) )
    allocate( dhaversine2d(xdim,ydim) )

    rlon2d = dble(lon2d)*pi/180.D0 
    rlat2d = dble(lat2d)*pi/180.D0 
    rplon  = dble(plon)*pi/180.D0 
    rplat  = dble(plat)*pi/180.D0 

    dlon = rlon2d - rplon 
    dlat = rlat2d - rplat 

    a = sin(dlat/2.D0)**2 + cos(rlat2d) * cos(rplat) * sin(dlon/2.D0)**2
    c = 2.D0 * asin(sqrt(a))

    dhaversine2d = c * earth_radius( plat ) 
    deallocate( rlon2d, rlat2d, dlon, dlat, a, c )

end function dhaversine2d

end module module_grid
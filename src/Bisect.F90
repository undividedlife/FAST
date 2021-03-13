!-------------------------------------------------------------------------------
!
!   MODULE Bisect
!
!>  @brief
!>   - Module for index searching
!
!   REVISION HISTORY
!
!   16AUG2015 : In-Sun Song : First complete version
!   09JUN2017 : In-Sun Song : Out-of-domain index at top and bottom boundaries
!                             are separately defined in Locver (version 1.0)
!   08JUL2017 : In-Sun Song : Reviewed and comments added
!   01APR2018 : In-Sun Song : Double checked and reviewed.
!   22FEB2021 : In-Sun Song : Symbol Locgrd is made using Locver
!
!-------------------------------------------------------------------------------
!
MODULE Bisect
!
    USE Base, ONLY: i4, r4, r8
!
    IMPLICIT NONE
!
    PRIVATE
!
    REAL(r8), PARAMETER :: EPS = 1.e-10_r8
!
    INTERFACE Bisect_Chkgrd
        MODULE PROCEDURE Bisect_Chkgrd_r4
        MODULE PROCEDURE Bisect_Chkgrd_r8
    END INTERFACE
!
    INTERFACE Bisect_Locgrd
        MODULE PROCEDURE Bisect_Locver_r4
        MODULE PROCEDURE Bisect_Locver_r8
    END INTERFACE
!
    INTERFACE Bisect_Locver
        MODULE PROCEDURE Bisect_Locver_r4
        MODULE PROCEDURE Bisect_Locver_r8
    END INTERFACE
!
    INTERFACE Bisect_Lochor
        MODULE PROCEDURE Bisect_Lochor_r4
        MODULE PROCEDURE Bisect_Lochor_r8
    END INTERFACE
!
    PUBLIC :: Bisect_Chkgrd
    PUBLIC :: Bisect_Locgrd
    PUBLIC :: Bisect_Locver
    PUBLIC :: Bisect_Lochor
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Bisect_Chkgrd
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Bisect_Chkgrd_r4(n, grd)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN) :: n
        REAL(r4),    DIMENSION(n), INTENT(IN) :: grd 
!
        INTEGER(i4) :: i
        LOGICAL :: flag
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Chkgrd_r4)'
!!!
!
        flag = .TRUE.
        DO i = 1, n-1
            IF (grd(i) >= grd(i+1)) flag = .FALSE.
        END DO
!
        IF (flag) THEN
            WRITE(6,'(A)') HEADER//': Grid monotonically increases'
        ELSE 
            WRITE(6,'(A)') HEADER//': Grid is not monotonic'
            WRITE(6,'(A)') HEADER//': Program stops'
            STOP
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Chkgrd_r4
!
    SUBROUTINE Bisect_Chkgrd_r8(n, grd)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN) :: n
        REAL(r8),    DIMENSION(n), INTENT(IN) :: grd 
!
        INTEGER(i4) :: i
        LOGICAL :: flag
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Chkgrd_r8)'
!!!
!
        flag = .TRUE.
        DO i = 1, n-1
            IF (grd(i) >= grd(i+1)) flag = .FALSE.
        END DO
!
        IF (flag) THEN
            WRITE(6,'(A)') HEADER//': Grid monotonically increases'
        ELSE 
            WRITE(6,'(A)') HEADER//': Grid is not monotonic'
            WRITE(6,'(A)') HEADER//': Program stops'
            STOP
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Chkgrd_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Bisect_Locver
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Bisect_Locver_r4(n, zv, z, k0, k1)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN)  :: n
        REAL(r4),    DIMENSION(n), INTENT(IN)  :: zv
        REAL(r4),                  INTENT(IN)  :: z
        INTEGER(i4),               INTENT(OUT) :: k0, k1
!
        INTEGER(i4) :: kl, km, ku
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Locver_r4)'
!!!
!
!       This routine works for monotonically increasing grid zv
!
        kl = 0
        ku = n+1
        IF (z == zv(1)) THEN
            k0 = 1
        ELSE IF (z == zv(n)) THEN
            k0 = n-1
        ELSE
            DO WHILE (ku-kl > 1) 
                km = (ku+kl)/2
                IF (z >= zv(km)) THEN
                    kl = km
                ELSE
                    ku = km
                END IF
            END DO
            k0 = kl
        END IF
!
!       Comments below shows what value k0 has when z < zv(1) or z > zv(n).
!       k0 = 0 when z < zv(1), and k0 = n when z > zv(n)
!
!       For n = 4,
!
!       If z < zv(1),  : below the lowest vertical grid: k0 = 0
!           ku=5, kl=0, ku-kl=5>1, km=2
!           z<zv(km)(yes), ku=km=2, ku-kl=2>1(yes), km=1
!           z<zv(km)(yes), ku=km=1, ku-kl=1>1(no),  k0=kl=0
!
!       If z > zv(4),  : above the highest vertical grid: k0 = n
!           ku=5, kl=0, ku-kl=5>1, km=2
!           z>zv(km)(yes), kl=km=2, ku-kl=3>1(yes), km=3
!           z>zv(km)(yes), kl=km=3, ku-kl=2>1(yes), km=4
!           z>zv(km)(yes), kl=km=4, ku-kl=1>1(no) , k0=kl=4
!
!       Determines k1. Basically k1 = k0+1
!
        IF (k0 == 0) THEN
            k0 = -99999
            k1 = -99999
        ELSE IF (k0 == n) THEN
            k0 =  99999
            k1 =  99999
        ELSE
            k1 = k0+1
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Locver_r4
!
    SUBROUTINE Bisect_Locver_r8(n, zv, z, k0, k1)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN)  :: n
        REAL(r8),    DIMENSION(n), INTENT(IN)  :: zv
        REAL(r8),                  INTENT(IN)  :: z
        INTEGER(i4),               INTENT(OUT) :: k0, k1
!
        INTEGER(i4) :: kl, km, ku
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Locver_r8)'
!!!
!
!       This routine works for monotonically increasing grid zv
!
        kl = 0
        ku = n+1
        IF (z == zv(1)) THEN
            k0 = 1
        ELSE IF (z == zv(n)) THEN
            k0 = n-1
        ELSE
            DO WHILE (ku-kl > 1) 
                km = (ku+kl)/2
                IF (z >= zv(km)) THEN
                    kl = km
                ELSE
                    ku = km
                END IF
            END DO
            k0 = kl
        END IF
!
!       Comments below shows what value k0 has when z < zv(1) or z > zv(n).
!       k0 = 0 when z < zv(1), and k0 = n when z > zv(n)
!
!       For n = 4,
!
!       If z < zv(1),  : below the lowest vertical grid: k0 = 0
!           ku=5, kl=0, ku-kl=5>1, km=2
!           z<zv(km)(yes), ku=km=2, ku-kl=2>1(yes), km=1
!           z<zv(km)(yes), ku=km=1, ku-kl=1>1(no),  k0=kl=0
!
!       If z > zv(4),  : above the highest vertical grid: k0 = n
!           ku=5, kl=0, ku-kl=5>1, km=2
!           z>zv(km)(yes), kl=km=2, ku-kl=3>1(yes), km=3
!           z>zv(km)(yes), kl=km=3, ku-kl=2>1(yes), km=4
!           z>zv(km)(yes), kl=km=4, ku-kl=1>1(no) , k0=kl=4
!
!       Determines k1. Basically k1 = k0+1
!
        IF (k0 == 0) THEN
            k0 = -99999
            k1 = -99999
        ELSE IF (k0 == n) THEN
            k0 =  99999
            k1 =  99999
        ELSE
            k1 = k0+1
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Locver_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Bisect_Lochor
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Bisect_Lochor_r4(n, xv, x, m, yv, y,  &
                               i0, i1, j0, j1, xorg, yorg)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN)  :: n, m
        REAL(r4),    DIMENSION(n), INTENT(IN)  :: xv
        REAL(r4),    DIMENSION(m), INTENT(IN)  :: yv
        REAL(r4),                  INTENT(INOUT)  :: x, y
        INTEGER(i4),               INTENT(OUT) :: i0, i1
        INTEGER(i4),               INTENT(OUT) :: j0, j1
        LOGICAL,                   INTENT(OUT), OPTIONAL :: xorg, yorg
!
        REAL(r4), DIMENSION(n+1) :: xvw   ! wrapped version of xv
        REAL(r4) :: x0, y0, sgn
        INTEGER(i4) :: nn
        INTEGER(i4) :: il, im, iu
        INTEGER(i4) :: jl, jm, ju
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Lochor_r4)'
!!!
!
!       This routine works for monotonically increasing grids xv and yv
!
        x0 = x     ! original value is saved for future use
        y0 = y    
!
!       Meridional index
!
!       Note that meridional index should always be found. A given latitude is
!       converted to a value within the range of -90 and 90.
!  
!       y0 = -450,-405,-390,-360,-315,-270,-225,-180,-135, -90, -45,   0
!       y  =  -90, -45, -30,   0,  45,  90,  45,   0, -45, -90, -45,   0
!
!       y0 =   0,  45,  90, 135, 180, 225, 270, 315, 360, 390, 405, 450
!       y  =   0,  45,  90,  45,   0, -45, -90, -45,   0,  30,  45,  90
!
        IF (ABS(y) > 90._r4) THEN
!
!           Makes initial |y| smaller than 360 because it is unnecessary to
!           deal with |y| > 360 without loss of generality. For example,
!           starting from y = 719 is equivalent to starting from 359.
!
            y = MOD(ABS(y), 360._r4)   ! Initially, set y > 0
            DO
                IF (y > 90._r4) THEN
                    y = 180._r4 - y
                    x = x + 180._r4      ! x always increases by 180 (IMPORTANT)
                ELSE IF (y < -90._r4) THEN
                    y = -180._r4 - y
                    x = x + 180._r4      ! x always increases by 180 (IMPORTANT)
                ELSE
                    EXIT
                END IF
            END DO
            y = SIGN(1._r4, y0) * y
!
        END IF
!
!       yorg = .TRUE.  if the computed y is the same as the given y (i.e., y0).
!       yorg = .FALSE. Otherwise
!
        IF (PRESENT(yorg)) THEN
            IF (y0 == y) THEN
                yorg = .TRUE.
            ELSE
                yorg = .FALSE.
            END IF
        END IF 
!
!       Searches meridional index using the (re)computed value of y
!       based on the bisection method
!
        jl = 0
        ju = m+1
        IF (y == yv(1)) THEN
            j0 = 1
        ELSE IF (y == yv(m)) THEN
            j0 = m-1
        ELSE
            DO WHILE (ju-jl > 1) 
                jm = (ju+jl)/2
                IF (y >= yv(jm)) THEN
                    jl = jm
                ELSE
                    ju = jm
                END IF
            END DO
            j0 = jl
        END IF
!
!       Determines j1. Basically j1 = j0+1 except for some special cases
!
        IF (j0 == 0) THEN  ! -90 - delta lat <= lat < -90
            j1 = 2           ! is equivalent to -90 <= lat < -90 + delta lat
            j0 = 1
        ELSE IF (j0 == m) THEN  ! 90 <= lat < 90 + delta lat
            j1 = m                ! is equivalent to 90 - delta lat <= lat < 90
            j0 = m-1
        ELSE
            j1 = j0+1
        END IF
!
!       Makes wrapped zonal grid
!
        nn = n+1
        xvw(1:n) = xv(1:n)
        xvw(nn) = xv(1) + 360._r4  ! xv(1) = -180 => xvw(nn) = 180
!                                  ! xv(1) =    0 => xvw(nn) = 360
!       Zonal index
!
        DO WHILE (x >= xvw(nn))
            x = x - 360._r4
        END DO
        DO WHILE (x < xvw( 1))
            x = x + 360._r4
        END DO
!
!       For the two cases below, x is set equal to xvw(1), which is usually
!       either 0 or -180, ignoring small error of EPS (10**-1).
!
!              170 --- 180 --- 190  or
!              350 --- 360 --- 370
!       CASE 1:            x
!       CASE 2:          x
!
        IF (x > xvw(nn) .AND. ABS(x-xvw(nn)) < EPS) x = xvw(1)   ! CASE 1
        IF (x < xvw(nn) .AND. ABS(x-xvw(nn)) < EPS) x = xvw(1)   ! CASE 2
!
!       xorg = .TRUE.  if the computed x is the same as the given x (i.e., x0).
!       xorg = .FALSE. Otherwise
!
        IF (PRESENT(xorg)) THEN
            IF (x0 == x) THEN
                xorg = .TRUE.
            ELSE
                xorg = .FALSE.
            END IF
        END IF
!
!       Searches for zonal index using the (re)computed value of x
!       based on the bisection method
!
        il = 0
        iu = nn+1
        IF (x == xvw(1)) THEN
            i0 = 1
        ELSE IF (x == xvw(nn)) THEN
            i0 = nn-1
        ELSE
            DO WHILE (iu-il > 1) 
                im = (iu+il)/2
                IF (x >= xvw(im)) THEN
                    il = im
                ELSE
                    iu = im
                END IF
            END DO
            i0 = il
        END IF
!
!       Determines i1. Basically i1 = i0+1
!       No special cases are considered in the longitudinal grid
!       because the longitudinal grid is cyclic.
!
        IF (i0 == n) THEN
            i1 = 1
        ELSE
            i1 = i0+1
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Lochor_r4
!
    SUBROUTINE Bisect_Lochor_r8(n, xv, x, m, yv, y,  &
                               i0, i1, j0, j1, xorg, yorg)
!
        IMPLICIT NONE
!
        INTEGER(i4),               INTENT(IN)  :: n, m
        REAL(r8),    DIMENSION(n), INTENT(IN)  :: xv
        REAL(r8),    DIMENSION(m), INTENT(IN)  :: yv
        REAL(r8),                  INTENT(INOUT)  :: x, y
        INTEGER(i4),               INTENT(OUT) :: i0, i1
        INTEGER(i4),               INTENT(OUT) :: j0, j1
        LOGICAL,                   INTENT(OUT), OPTIONAL :: xorg, yorg
!
        REAL(r8), DIMENSION(n+1) :: xvw   ! wrapped version of xv
        REAL(r8) :: x0, y0, sgn
        INTEGER(i4) :: nn
        INTEGER(i4) :: il, im, iu
        INTEGER(i4) :: jl, jm, ju
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Bisect_Lochor_r8)'
!!!
!
!       This routine works for monotonically increasing grids xv and yv
!
        x0 = x     ! original value is saved for future use
        y0 = y    
!
!       Meridional index
!
!       Note that meridional index should always be found. A given latitude is
!       converted to a value within the range of -90 and 90.
!  
!       y0 = -450,-405,-390,-360,-315,-270,-225,-180,-135, -90, -45,   0
!       y  =  -90, -45, -30,   0,  45,  90,  45,   0, -45, -90, -45,   0
!
!       y0 =   0,  45,  90, 135, 180, 225, 270, 315, 360, 390, 405, 450
!       y  =   0,  45,  90,  45,   0, -45, -90, -45,   0,  30,  45,  90
!
        IF (ABS(y) > 90._r8) THEN
!
!           Makes initial |y| smaller than 360 because it is unnecessary to
!           deal with |y| > 360 without loss of generality. For example,
!           starting from y = 719 is equivalent to starting from 359.
!
            y = MOD(ABS(y), 360._r8)   ! Initially, set y > 0
            DO
                IF (y > 90._r8) THEN
                    y = 180._r8 - y
                    x = x + 180._r8      ! x always increases by 180 (IMPORTANT)
                ELSE IF (y < -90._r8) THEN
                    y = -180._r8 - y
                    x = x + 180._r8      ! x always increases by 180 (IMPORTANT)
                ELSE
                    EXIT
                END IF
            END DO
            y = SIGN(1._r8, y0) * y
!
        END IF
!
!       yorg = .TRUE.  if the computed y is the same as the given y (i.e., y0).
!       yorg = .FALSE. Otherwise
!
        IF (PRESENT(yorg)) THEN
            IF (y0 == y) THEN
                yorg = .TRUE.
            ELSE
                yorg = .FALSE.
            END IF
        END IF 
!
!       Searches meridional index using the (re)computed value of y
!       based on the bisection method
!
        jl = 0
        ju = m+1
        IF (y == yv(1)) THEN
            j0 = 1
        ELSE IF (y == yv(m)) THEN
            j0 = m-1
        ELSE
            DO WHILE (ju-jl > 1) 
                jm = (ju+jl)/2
                IF (y >= yv(jm)) THEN
                    jl = jm
                ELSE
                    ju = jm
                END IF
            END DO
            j0 = jl
        END IF
!
!       Determines j1. Basically j1 = j0+1 except for some special cases
!
        IF (j0 == 0) THEN  ! -90 - delta lat <= lat < -90
            j1 = 2           ! is equivalent to -90 <= lat < -90 + delta lat
            j0 = 1
        ELSE IF (j0 == m) THEN  ! 90 <= lat < 90 + delta lat
            j1 = m                ! is equivalent to 90 - delta lat <= lat < 90
            j0 = m-1
        ELSE
            j1 = j0+1
        END IF
!
!       Makes wrapped zonal grid
!
        nn = n+1
        xvw(1:n) = xv(1:n)
        xvw(nn) = xv(1) + 360._r8  ! xv(1) = -180 => xvw(nn) = 180
!                                  ! xv(1) =    0 => xvw(nn) = 360
!       Zonal index
!
        DO WHILE (x >= xvw(nn))
            x = x - 360._r8
        END DO
        DO WHILE (x < xvw( 1))
            x = x + 360._r8
        END DO
!
!       For the two cases below, x is set equal to xvw(1), which is usually
!       either 0 or -180, ignoring small error of EPS (10**-1).
!
!              170 --- 180 --- 190  or
!              350 --- 360 --- 370
!       CASE 1:            x
!       CASE 2:          x
!
        IF (x > xvw(nn) .AND. ABS(x-xvw(nn)) < EPS) x = xvw(1)   ! CASE 1
        IF (x < xvw(nn) .AND. ABS(x-xvw(nn)) < EPS) x = xvw(1)   ! CASE 2
!
!       xorg = .TRUE.  if the computed x is the same as the given x (i.e., x0).
!       xorg = .FALSE. Otherwise
!
        IF (PRESENT(xorg)) THEN
            IF (x0 == x) THEN
                xorg = .TRUE.
            ELSE
                xorg = .FALSE.
            END IF
        END IF
!
!       Searches for zonal index using the (re)computed value of x
!       based on the bisection method
!
        il = 0
        iu = nn+1
        IF (x == xvw(1)) THEN
            i0 = 1
        ELSE IF (x == xvw(nn)) THEN
            i0 = nn-1
        ELSE
            DO WHILE (iu-il > 1) 
                im = (iu+il)/2
                IF (x >= xvw(im)) THEN
                    il = im
                ELSE
                    iu = im
                END IF
            END DO
            i0 = il
        END IF
!
!       Determines i1. Basically i1 = i0+1
!       No special cases are considered in the longitudinal grid
!       because the longitudinal grid is cyclic.
!
        IF (i0 == n) THEN
            i1 = 1
        ELSE
            i1 = i0+1
        END IF
!
        RETURN
    END SUBROUTINE Bisect_Lochor_r8
!
END

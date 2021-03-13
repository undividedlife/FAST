PROGRAM Example2
!
   USE Base,      ONLY: i4, r4, r8
   USE Bisect,    ONLY: Bisect_Lochor, Bisect_Locver
   USE PhysConst, ONLY: PI
!
   IMPLICIT NONE
!
   INTEGER(i4), PARAMETER :: NLON = 128
   INTEGER(i4), PARAMETER :: NLAT = 65
   INTEGER(i4), PARAMETER :: NLEV = 51
!
   REAL(r8), DIMENSION(NLON) :: rlon, rlam
   REAL(r8), DIMENSION(NLAT) :: rlat, rphi
   REAL(r8), DIMENSION(NLEV) :: rlev
!
   REAL(r8) :: xlon, xlat, xlev
   INTEGER(i4) :: i, i0, j0, k0, i1, j1, k1
   LOGICAL :: flagx, flagy
!
!!!
!
   CALL SetGrid(NLON, NLAT, NLEV, rlon, rlam, rlat, rphi, rlev)
!
   DO i = 1, 100000
      xlon = 30.0_r8
      xlat = 0.0_r8 + i*0.2_r8
      xlev = 0.0_r8 + i*10.0_r8
      CALL Bisect_Lochor(NLON, rlon, xlon, NLAT, rlat, xlat, i0, i1, j0, j1, xorg=flagx, yorg=flagy)
      CALL Bisect_Locver(NLEV, rlev, xlev, k0, k1)
      IF (k0 /= -99999 .AND. k0 /= 99999) THEN
         WRITE(10,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'O ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
         IF ((i0 < 1 .OR. i0 > NLON) .OR. (i1 < 1 .OR. i1 > NLON) .OR.  &
             (j0 < 1 .OR. j0 > NLAT) .OR. (j1 < 1 .OR. j1 > NLAT) .OR.  &
             (k0 < 1 .OR. k0 > NLEV) .OR. (k1 < 1 .OR. k1 > NLEV)) THEN
            WRITE(6,*) 'Invalid operation 0'
            STOP
         END IF
      ELSE
         WRITE(11,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'X ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
      END IF
   END DO
!
   DO i = 1, 100000
      xlon = 0.0_r8 + i*0.5_r8
      xlat = 0.0_r8 + i*0.2_r8
      xlev = 0.0_r8 + i*10.0_r8
      CALL Bisect_Lochor(NLON, rlon, xlon, NLAT, rlat, xlat, i0, i1, j0, j1, xorg=flagx, yorg=flagy)
      CALL Bisect_Locver(NLEV, rlev, xlev, k0, k1)
      IF (k0 /= -99999 .AND. k0 /= 99999) THEN
         WRITE(20,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'O ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
         IF ((i0 < 1 .OR. i0 > NLON) .OR. (i1 < 1 .OR. i1 > NLON) .OR.  &
             (j0 < 1 .OR. j0 > NLAT) .OR. (j1 < 1 .OR. j1 > NLAT) .OR.  &
             (k0 < 1 .OR. k0 > NLEV) .OR. (k1 < 1 .OR. k1 > NLEV)) THEN
            WRITE(6,*) 'Invalid operation 1'
            STOP
         END IF
      ELSE
         WRITE(21,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'X ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
      END IF
   END DO
!
   DO i = 1, 100000
      xlon = 0.0_r8 - i*0.5_r8
      xlat = 0.0_r8 - i*0.2_r8
      xlev = 0.0_r8 + i*10.0_r8
      CALL Bisect_Lochor(NLON, rlon, xlon, NLAT, rlat, xlat, i0, i1, j0, j1, xorg=flagx, yorg=flagy)
      CALL Bisect_Locver(NLEV, rlev, xlev, k0, k1)
      IF (k0 /= -99999 .AND. k0 /= 99999) THEN
         WRITE(30,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'O ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
         IF ((i0 < 1 .OR. i0 > NLON) .OR. (i1 < 1 .OR. i1 > NLON) .OR.  &
             (j0 < 1 .OR. j0 > NLAT) .OR. (j1 < 1 .OR. j1 > NLAT) .OR.  &
             (k0 < 1 .OR. k0 > NLEV) .OR. (k1 < 1 .OR. k1 > NLEV)) THEN
            WRITE(6,*) 'Invalid operation 2'
            STOP
         END IF
      ELSE
         WRITE(31,'(A,I9,1X,2(L2,1X),3(F15.4,1X),6(I4,1X))')  &
            'X ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
      END IF
   END DO
!
CONTAINS
!
!-------------------------------------------------------------------------------
!  SUBROUTINE SetGrid
!-------------------------------------------------------------------------------
!
   SUBROUTINE SetGrid(mlon, mlat, mlev, alon, alam, alat, aphi, alev)
!
      IMPLICIT NONE
!
      INTEGER(i4),                  INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon), INTENT(INOUT) :: alon, alam
      REAL(r8),    DIMENSION(mlat), INTENT(INOUT) :: alat, aphi
      REAL(r8),    DIMENSION(mlev), INTENT(INOUT) :: alev
!
      INTEGER(i4) :: i, j, k, ip
!!!
!
      DO i = 1, mlon
         alon(i) = DBLE(i-1)*360._r8/DBLE(mlon)
      END DO
      DO i = 1, mlon
         alam(i) = alon(i)*PI/180._r8
      END DO
      DO j = 1, mlat
         alat(j) = DBLE(j-1)*180._r8/DBLE(mlat-1) - 90._r8
      END DO
      DO i = 1, mlat
         aphi(i) = alat(i)*PI/180._r8
      END DO
      DO k = 1, mlev
         alev(k) = DBLE(k-1)*2000._r8
      END DO
!
      WRITE(6,'(A,3(F15.7,1X))') 'lon : ', MINVAL(alon), MAXVAL(alon), SUM(alon)/DBLE(mlon)
!     DO i = 1, mlon
!        ip = NLON/2+i
!        IF (ip > NLON) ip = ip - NLON
!        WRITE(6,'(2(I4,1X),2(F15.7,1X))') i, ip, alon(i), alon(ip)
!     END DO
      WRITE(6,'(A,3(F15.7,1X))') 'lam : ', MINVAL(alam), MAXVAL(alam), SUM(alam)/DBLE(mlon)
      WRITE(6,'(A,3(F15.7,1X))') 'lat : ', MINVAL(alat), MAXVAL(alat), SUM(alat)/DBLE(mlat)
      WRITE(6,'(A,3(F15.7,1X))') 'phi : ', MINVAL(aphi), MAXVAL(aphi), SUM(aphi)/DBLE(mlat)
      WRITE(6,'(A,3(F15.7,1X))') 'lev : ', MINVAL(alev), MAXVAL(alev), SUM(alev)/DBLE(mlev)
!
      RETURN
   END SUBROUTINE SetGrid

END

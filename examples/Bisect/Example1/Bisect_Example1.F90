PROGRAM Example1
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
   REAL(r8) :: xlon, xlat, xlev, dz
   INTEGER(i4) :: i, k, i0, j0, k0, i1, j1, k1
   LOGICAL :: flagx, flagy
!
!!!
!
   CALL SetGrid(NLON, NLAT, NLEV, rlon, rlam, rlat, rphi, rlev)
!
   CALL RANDOM_SEED()
!
   CALL Bisect_Locver(NLEV, rlev, -10.0_r8, k0, k1)
   PRINT *, k0, k1
   CALL Bisect_Locver(NLEV, rlev, 1.1e+5_r8, k0, k1)
   PRINT *, k0, k1
!  STOP
!
   DO i = 1, 10
      IF (i == 1) THEN
         dz = 1000._r8
      ELSE IF (i == 2) THEN
         dz = 2500._r8
      ELSE
         CALL RANDOM_NUMBER(dz); dz = dz*100._r8
      END IF
      xlev = 0.0_r8
      WRITE(9,'(A,F15.7,A)') '-----', dz, '-----'
      DO
         CALL Bisect_Locver(NLEV, rlev, xlev, k0, k1)
         IF (k0 == -99999 .OR. k0 == 99999) EXIT
         WRITE(9,'(3(F15.7,1X),2(I4,1X))') xlev, rlev(k0), rlev(k1), k0, k1
         xlev = xlev + dz
      END DO
   END DO
!
   DO i = 1, 100000
      CALL RANDOM_NUMBER(xlon); xlon = xlon*720._r8 - 180.0_r8
      CALL RANDOM_NUMBER(xlat); xlat = xlat*360._r8 - 180.0_r8
      CALL RANDOM_NUMBER(xlev); xlev = xlev*100000._r8 - 10000._r8
      CALL Bisect_Lochor(NLON, rlon, xlon, NLAT, rlat, xlat, i0, i1, j0, j1, xorg=flagx, yorg=flagy)
      CALL Bisect_Locver(NLEV, rlev, xlev, k0, k1)
      IF (k0 /= -99999 .AND. k0 /= 99999) THEN
         IF ((i0 < 1 .OR. i0 > NLON) .OR. (i1 < 1 .OR. i1 > NLON) .OR.  &
             (j0 < 1 .OR. j0 > NLAT) .OR. (j1 < 1 .OR. j1 > NLAT) .OR.  &
             (k0 < 1 .OR. k0 > NLEV) .OR. (k1 < 1 .OR. k1 > NLEV)) THEN
            WRITE(6,*) 'Invalid operation'
            STOP
         END IF
         WRITE(10,'(A,I9,2(L6,1X),3(F15.7,1X),6(I4,1X))')  &
            'O ', i, flagx, flagy, xlon, xlat, xlev, i0, i1, j0, j1, k0, k1
      ELSE
         WRITE(11,'(A,I9,2(L6,1X),3(F15.7,1X),6(I4,1X))')  &
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

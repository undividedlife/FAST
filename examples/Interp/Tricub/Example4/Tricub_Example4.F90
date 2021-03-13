PROGRAM Tricub_Example4
!
    USE Base,   ONLY: i4, r4, r8
    USE TriCub, ONLY: Tricub_Init, Tricub_Coef, Tricub_Eval
!
    IMPLICIT NONE
!
    INTEGER(i4), PARAMETER :: NLON = 128
    INTEGER(i4), PARAMETER :: NLAT = 65
    INTEGER(i4), PARAMETER :: NLEV = 51
!
    REAL(r8), PARAMETER :: PI = 3.1415926535897932385_r8
    REAL(r8), PARAMETER :: PI2 = PI*2._r8
    REAL(r8), PARAMETER :: GRAV = 9.806_r8
    REAL(r8), PARAMETER :: RADE = 6371.0_r8
!
    REAL(r8), DIMENSION(NLON) :: rlon, rlam
    REAL(r8), DIMENSION(NLAT) :: rlat, rphi
    REAL(r8), DIMENSION(NLEV) :: rlev
    REAL(r8), DIMENSION(NLON,NLAT,NLEV)   :: fld
    REAL(r8), DIMENSION(NLON,NLAT,NLEV,7) :: drv, drvanal
!
    REAL(r8), DIMENSION(64) :: a
    REAL(r8), DIMENSION(8) :: f
    REAL(r8), DIMENSION(8,7) :: deriv 
    REAL(r8), DIMENSION(7) :: dlvec
    REAL(r8) :: xr, yr, zr, dx, dy, dz
    REAL(r8) :: x, y, z
    INTEGER(i4) :: i0, j0, k0, i1, j1, k1
    INTEGER(i4) :: l
!!!!
!
    CALL SetGrid(NLON, NLAT, NLEV, rlon, rlam, rlat, rphi, rlev)
    CALL ComputField(NLON, NLAT, NLEV, rlam, rphi, rlev, fld)
    CALL ComputDeriv(NLON, NLAT, NLEV, rlam, rphi, rlev, fld, drv, drvanal)
!
    CALL Tricub_Init
!
    i0 = 10  ; j0 = 60  ; k0 = 40
    i1 = i0+1; j1 = j0+1; k1 = k0+1
!
    IF (i1 == NLON+1) THEN 
        i1 = 1
        dx = 360._r8-rlam(NLON)
    ELSE
        dx = rlam(i1)-rlam(i0)
    END IF
    dy = rphi(j1)-rphi(j0)
    dz = rlev(k1)-rlev(k0)
!
    xr = rlam(i0)*0.3_r8 + rlam(i1)*0.7_r8
    yr = rphi(j0)*0.3_r8 + rphi(j1)*0.7_r8
    zr = rlev(k0)*0.3_r8 + rlev(k1)*0.7_r8
!
    f(:) = (/fld(i0,j0,k0),fld(i1,j0,k0),  &
             fld(i0,j1,k0),fld(i1,j1,k0),  &
             fld(i0,j0,k1),fld(i1,j0,k1),  &
             fld(i0,j1,k1),fld(i1,j1,k1)/)
!
    dlvec = (/dx, dy, dz, dx*dy, dx*dz, dy*dz, dx*dy*dz/)
    DO l = 1, 7
        deriv(:,l) = (/drv(i0,j0,k0,l),drv(i1,j0,k0,l),  &
                       drv(i0,j1,k0,l),drv(i1,j1,k0,l),  &
                       drv(i0,j0,k1,l),drv(i1,j0,k1,l),  &
                       drv(i0,j1,k1,l),drv(i1,j1,k1,l)/) * dlvec(l)
    END DO
!
    CALL Tricub_Coef(a, f, deriv)
!
!   WRITE(6,'(A)') '---'
!   WRITE(6,'(4(F15.7,1X))') f
!   WRITE(6,'(A)') '---'
!   WRITE(6,'(4(F15.7,1X))') deriv
!   WRITE(6,'(A)') '---'
!   WRITE(6,'(4(F15.7,1X))') a
!
    x = (xr-rlam(i0))/dx
    y = (yr-rphi(j0))/dy
    z = (zr-rlev(k0))/dz
!
    WRITE(6,'(A)') '---'
    WRITE(6,'(5(F15.7,1X))') xr, yr, zr, Tricub_Eval(a, x, y, z),  &
        AnalEval(xr, yr, zr)
    WRITE(6,'(6(F15.7,1X))') xr, yr, zr, Tricub_Eval(a, x, y, z, 1, 0, 0),  &
        Tricub_Eval(a, x, y, z, 1, 0, 0)/dx,  &
        AnalEval(xr, yr, zr, derx=1, dery=0, derz=0)
    WRITE(6,'(6(F15.7,1X))') xr, yr, zr, Tricub_Eval(a, x, y, z, 0, 1, 0),  &
        Tricub_Eval(a, x, y, z, 0, 1, 0)/dy,  &
        AnalEval(xr, yr, zr, derx=0, dery=1, derz=0)
    WRITE(6,'(6(F15.7,1X))') xr, yr, zr, Tricub_Eval(a, x, y, z, 0, 0, 1),  &
        Tricub_Eval(a, x, y, z, 0, 0, 1)/dz,  &
        AnalEval(xr, yr, zr, derx=0, dery=0, derz=1)
!
CONTAINS
!
!-------------------------------------------------------------------------------
!   SUBROUTINE SetGrid
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
!!!!
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
!       DO i = 1, mlon
!           ip = NLON/2+i
!           IF (ip > NLON) ip = ip - NLON
!           WRITE(6,'(2(I4,1X),2(F15.7,1X))') i, ip, alon(i), alon(ip)
!       END DO
        WRITE(6,'(A,3(F15.7,1X))') 'lam : ', MINVAL(alam), MAXVAL(alam), SUM(alam)/DBLE(mlon)
        WRITE(6,'(A,3(F15.7,1X))') 'lat : ', MINVAL(alat), MAXVAL(alat), SUM(alat)/DBLE(mlat)
        WRITE(6,'(A,3(F15.7,1X))') 'phi : ', MINVAL(aphi), MAXVAL(aphi), SUM(aphi)/DBLE(mlat)
        WRITE(6,'(A,3(F15.7,1X))') 'lev : ', MINVAL(alev), MAXVAL(alev), SUM(alev)/DBLE(mlev)
!
        RETURN
    END SUBROUTINE SetGrid
!
!-------------------------------------------------------------------------------
!   SUBROUTINE ComputField
!-------------------------------------------------------------------------------
!
    SUBROUTINE ComputField(mlon, mlat, mlev, alam, aphi, alev, fld)
!
        IMPLICIT NONE
!
        INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
        REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
        REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
        REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
        REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: fld
!
        REAL(r8), DIMENSION(mlon) :: funclon
        REAL(r8), DIMENSION(mlat) :: funclat
        REAL(r8), DIMENSION(mlev) :: funclev
!
        INTEGER(i4) :: i, j, k
!!!
!
!       funclon = SIN[3*lambda]
!
        DO i = 1, mlon
            funclon(i) = SIN(3._r8 * alam(i))
        END DO
!
!     funclat = LegrendreP[n, m, Cos[phi]] = LegendreP[7, 3, Cos[phi]]
!
      DO j = 1, mlat
         funclat(j) = - (315._r8/8._r8) *  &
                        (3._r8 - 66._r8 * COS(aphi(j))**2 + 143._r8 * COS(aphi(j))**4) *  &
                        ABS(SIN(aphi(j)))**3
      END DO
!
!     funclev = Exp[-(z-50)^2/20^2]
!
      DO k = 1, mlev
         funclev(k) = 0.5_r8 * EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2)
      END DO
!
!     Total field
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               fld(i,j,k) = 30._r8 + funclon(i) * funclat(j) * funclev(k)
            END DO
         END DO
      END DO
!
      WRITE(6,'(A,3(F15.7,1X))') 'funclat : ', MINVAL(funclat), MAXVAL(funclat), SUM(funclat)/DBLE(mlat)
      WRITE(6,'(A,3(F15.7,1X))') 'funclon : ', MINVAL(funclon), MAXVAL(funclon), SUM(funclon)/DBLE(mlon)
      WRITE(6,'(A,3(F15.7,1X))') 'funclev : ', MINVAL(funclev), MAXVAL(funclev), SUM(funclev)/DBLE(mlev)
      WRITE(6,'(A,3(F15.7,1X))') 'fld  : ', MINVAL(fld), MAXVAL(fld), SUM(fld)/DBLE(mlon*mlat*mlev)
!
      RETURN
   END SUBROUTINE ComputField
!
!-------------------------------------------------------------------------------
!  FUNCTION AnalEval
!-------------------------------------------------------------------------------
!
   FUNCTION AnalEval(x, y, z, derx, dery, derz) RESULT(ret)
!
      IMPLICIT NONE
!
      REAL(r8), INTENT(IN) :: x, y, z
      INTEGER(i4), INTENT(IN), OPTIONAL :: derx, dery, derz
!
      REAL(r8) :: funclat, funclon, funclev, ret
!!!
!
!     See Derivative.nb for details
!
!     r[lambda, phi, z] = 30 + f[phi] g[lambda] h[z]
!
!     where f[phi] = LegendreP[7, 3, Cos[phi]] = 
!     -(315/8) * (3 - 66 Cos[phi]^2 + 143 Cos[phi]^4) * (Sin[phi]^2)^{3/2};
!     g[lambda] = Sin[3 lambda]; and h[z] = 0.5 Exp[-(z/1000 - 50)^2 / 20^2].
!
      IF ((.NOT. PRESENT(derx)) .AND. (.NOT. PRESENT(dery)) .AND.  &
          (.NOT. PRESENT(derz))) THEN
         funclon = SIN(3._r8 * x)
         funclat = - (315._r8/8._r8) *  &
                     (3._r8 - 66._r8 * COS(y)**2 + 143._r8 * COS(y)**4) *  &
                     ABS(SIN(y))**3
         funclev = 0.5_r8 * EXP(-(z/1000._r8 - 50._r8)**2 / (20._r8)**2)
         ret = 30._r8 + funclon * funclat * funclev
         RETURN
      END IF
!
!     Zonal derivative of r[lambda, phi, z]
!
      IF (PRESENT(derx)) THEN
         IF (derx == 1) THEN 
            ret = EXP(-(z/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  COS(3._r8 * x) *  &
                  (-177.1875_r8 + 4075.3125_r8 * COS(y)**2 -  &
                   12344.0625_r8 * COS(y)**4 + 8445.9375_r8 * COS(y)**6) *  &
                  ABS(SIN(y))
            RETURN
         END IF
      END IF
!
!     Meridional derivative of r[lambda, phi, z]
!
      IF (PRESENT(dery)) THEN
          IF (dery == 1) THEN 
            ret = EXP(-(z/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (-2775.9375_r8 + 20534.0625_r8 * COS(y)**2  &
                   -37465.3125_r8 * COS(y)**4 + 19707.1875_r8 * COS(y)**6) *  &
                  COS(y) * SIGN(1._r8, SIN(y)) * SIN(3._r8 * x)
            RETURN
         END IF
      END IF
!
!     Vertical derivative of r[lambda, phi, z]
!
      IF (PRESENT(derz)) THEN
         IF (derz == 1) THEN 
            ret = 0.0000140765626_r8 * (z - 50000._r8) *  &
                  EXP(-(z/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  SIN(y)**2 * ABS(SIN(y)) *  &
                  (0.020979020979_r8 - 0.461538461538_r8 * COS(y)**2 + COS(y)**4 ) *  &
                  SIN(3._r8 * x)
         END IF
      END IF
!

      RETURN
   END FUNCTION AnalEval
!
!-------------------------------------------------------------------------------
!  SUBROUTINE ComputDeriv
!-------------------------------------------------------------------------------
!
   SUBROUTINE ComputDeriv(mlon, mlat, mlev, alam, aphi, alev, fld, drv, drvanal)
!
      IMPLICIT NONE
!
      INTEGER(i4),                              INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),             INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),             INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),             INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev),   INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev,7), INTENT(INOUT) :: drv, drvanal
!
      REAL(r8), PARAMETER :: RADE = 6371.0_r8
      REAL(r8) :: drvtmp,drverr,fac
      INTEGER(i4) :: i, j, k, ip
!!!
!
!     (1) dF/dl (\frac{\partial F}{\partial \lambda})
!
      CALL Deriv100    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,1))
      CALL Deriv100Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,1))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dl    ', drv(:,:,:,1), drvanal(:,:,:,1))
!
!     (2) dF/dp (\frac{\partial F}{\partial \phi})
!
      CALL Deriv010    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,2))
      CALL Deriv010Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,2))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dp    ', drv(:,:,:,2), drvanal(:,:,:,2))
!
!     (3) dF/dz (\frac{\partial F}{\partial z})
!
      CALL Deriv001    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,3))
      CALL Deriv001Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,3))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dz    ', drv(:,:,:,3), drvanal(:,:,:,3))
!
!     (4) d^2F/dldp
!
      CALL Deriv110    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,4))
      CALL Deriv110Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,4))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dldp  ', drv(:,:,:,4), drvanal(:,:,:,4))
!
!     (5) d^2F/dldz
!
      CALL Deriv101    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,5))
      CALL Deriv101Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,5))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dldz  ', drv(:,:,:,5), drvanal(:,:,:,5))
!
!     (6) d^2F/dpdz
!
      CALL Deriv011    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,6))
      CALL Deriv011Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,6))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dpdz  ', drv(:,:,:,6), drvanal(:,:,:,6))
!
!     (7) d^2F/dldpdz
!
      CALL Deriv111    (mlon, mlat, mlev, alam, aphi, alev, fld, drv    (:,:,:,7))
      CALL Deriv111Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drvanal(:,:,:,7))
      CALL DerivChk    (mlon, mlat, mlev, 'dF/dldpdz', drv(:,:,:,7), drvanal(:,:,:,7))
!
      RETURN
   END SUBROUTINE ComputDeriv
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DerivX
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv100(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8), PARAMETER :: PI = 3.1415926535897932385_r8
      REAL(r8), PARAMETER :: PI2 = PI*2._r8
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 2, mlon-1
               drv(i,j,k) = (fld(i+1,j,k)-fld(i-1,j,k))/(alam(i+1)-alam(i-1))
            END DO
            drv(   1,j,k) = (fld(2,j,k)-fld(mlon,  j,k))/(alam(2)+PI2-alam(mlon  ))
            drv(mlon,j,k) = (fld(1,j,k)-fld(mlon-1,j,k))/(alam(1)+PI2-alam(mlon-1))
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv100
! 
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv100Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv100Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  COS(3._r8 * alam(i)) *  &
                  (-177.1875_r8 + 4075.3125_r8 * COS(aphi(j))**2  &
                   -12344.0625_r8 * COS(aphi(j))**4 + 8445.9375_r8 * COS(aphi(j))**6) *  &
                  ABS(SIN(aphi(j)))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv100Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv010
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv010(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 2, mlat-1
            DO i = 1, mlon
               drv(i,j,k) = (fld(i,j+1,k)-fld(i,j-1,k))/(aphi(j+1)-aphi(j-1))
            END DO
         END DO
      END DO
      DO k = 1, mlev
         DO i = 1, mlon   ! South pole (southward extrapolation)
            drvtmp = (fld(i,2,k)-fld(i,1,k))/(aphi(2)-aphi(1))
            drv(i,1,k) = 2._r8 * drvtmp - drv(i,2,k)
         END DO
      END DO
      DO k = 1, mlev
         DO i = 1, mlon   ! North pole (northward extrapolation)
            drvtmp = (fld(i,mlat,k)-fld(i,mlat-1,k))/(aphi(mlat)-aphi(mlat-1))
            drv(i,mlat,k) = 2._r8 * drvtmp - drv(i,mlat-1,k)
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv010
! 
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv010Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv010Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (-2775.9375_r8 + 20534.0625_r8 * COS(aphi(j))**2  &
                   -37465.3125_r8 * COS(aphi(j))**4 + 19707.1875_r8 * COS(aphi(j))**6) *  &
                  COS(aphi(j)) * SIGN(1._r8, SIN(aphi(j))) * SIN(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv010Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv001
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv001(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 2, mlev-1
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) = (fld(i,j,k+1)-fld(i,j,k-1))/(alev(k+1)-alev(k-1))
            END DO
         END DO
      END DO
      DO j = 1, mlat
         DO i = 1, mlon  ! Bottom extrapolation
            drvtmp = (fld(i,j,2)-fld(i,j,1))/(alev(2)-alev(1))
            drv(i,j,1) = 2._r8 * drvtmp - drv(i,j,2)
         END DO
      END DO
      DO j = 1, mlat
         DO i = 1, mlon  ! Bottom extrapolation
            drvtmp = (fld(i,j,mlev)-fld(i,j,mlev-1))/(alev(mlev)-alev(mlev-1))
            drv(i,j,mlev) = 2._r8 * drvtmp - drv(i,j,mlev-1)
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv001
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv001Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv001Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  0.0000140765625_r8 * (alev(k) - 50000._r8) *  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  SIN(aphi(j))**2 * ABS(SIN(aphi(j))) *  &
                  (0.020979020979_r8 - 0.461538461538_r8 * COS(aphi(j))**2 + COS(aphi(j))**4 ) *  &
                  SIN(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv001Anal
! 
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv110
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv110(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8), DIMENSION(mlon,mlat,mlev) :: drv100
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      CALL Deriv100(mlon, mlat, mlev, alam, aphi, alev, fld, drv100)
      CALL Deriv010(mlon, mlat, mlev, alam, aphi, alev, drv100, drv)
!
      RETURN
   END SUBROUTINE Deriv110
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv110Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv110Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (-8327.8125_r8 + 61602.1875_r8 * COS(aphi(j))**2   &
                   -112395.9375_r8 * COS(aphi(j))**4 + 59121.5625_r8 * COS(aphi(j))**6) *  &
                  COS(aphi(j)) * SIGN(1._r8,SIN(aphi(j))) *  &
                  COS(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv110Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv101
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv101(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8), DIMENSION(mlon,mlat,mlev) :: drv100
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      CALL Deriv100(mlon, mlat, mlev, alam, aphi, alev, fld, drv100)
      CALL Deriv001(mlon, mlat, mlev, alam, aphi, alev, drv100, drv)
!
      RETURN
   END SUBROUTINE Deriv101
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv101Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv101Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) = 0.0000422296875_r8 *  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (alev(k) - 50000._r8) * SIN(aphi(j))**2 *  &
                  (0.020979020979_r8 - 0.461538461538_r8 * COS(aphi(j))**2 + COS(aphi(j))**4) *  &
                  ABS(SIN(aphi(j))) * COS(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv101Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv011
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv011(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8), DIMENSION(mlon,mlat,mlev) :: drv010
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      CALL Deriv010(mlon, mlat, mlev, alam, aphi, alev, fld, drv010)
      CALL Deriv001(mlon, mlat, mlev, alam, aphi, alev, drv010, drv)
!
      RETURN
   END SUBROUTINE Deriv011
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv011Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv011Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (-0.693984375_r8 + 0.0000138796875_r8 * alev(k) +  &
                   (5.133515625_r8 - 0.0001026703125_r8 * alev(k)) * COS(aphi(j))**2 +   &
                   (-9.366328125_r8 + 0.0001873265625_r8 * alev(k)) * COS(aphi(j))**4 +  &
                   (4.926796875_r8 - 0.0000985359375_r8 * alev(k)) * COS(aphi(j))**6) *  &
                  COS(aphi(j)) * SIGN(1._r8,SIN(aphi(j))) *  &
                  SIN(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv011Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv111
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv111(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8), DIMENSION(mlon,mlat,mlev) :: drv100
      REAL(r8), DIMENSION(mlon,mlat,mlev) :: drv110
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      CALL Deriv100(mlon, mlat, mlev, alam, aphi, alev, fld, drv100)
      CALL Deriv010(mlon, mlat, mlev, alam, aphi, alev, drv100, drv110)
      CALL Deriv001(mlon, mlat, mlev, alam, aphi, alev, drv110, drv)
!
      RETURN
   END SUBROUTINE Deriv111
!
!-------------------------------------------------------------------------------
!  SUBROUTINE Deriv111Anal
!-------------------------------------------------------------------------------
!
   SUBROUTINE Deriv111Anal(mlon, mlat, mlev, alam, aphi, alev, fld, drv)
!
      IMPLICIT NONE
!
      INTEGER(i4),                            INTENT(IN)    :: mlon, mlat, mlev
      REAL(r8),    DIMENSION(mlon),           INTENT(IN)    :: alam
      REAL(r8),    DIMENSION(mlat),           INTENT(IN)    :: aphi
      REAL(r8),    DIMENSION(mlev),           INTENT(IN)    :: alev
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(IN)    :: fld
      REAL(r8),    DIMENSION(mlon,mlat,mlev), INTENT(INOUT) :: drv
!
      REAL(r8) :: drvtmp
      INTEGER(i4) :: i, j, k
!!!
!
      DO k = 1, mlev
         DO j = 1, mlat
            DO i = 1, mlon
               drv(i,j,k) =  &
                  EXP(-(alev(k)/1000._r8 - 50._r8)**2 / (20._r8)**2) *  &
                  (-2.081953125_r8 + 0.0000416390625_r8 * alev(k) +  &
                   (15.400546875_r8 - 0.0003080109375_r8 * alev(k)) * COS(aphi(j))**2 +   &
                   (-28.098984375_r8 + 0.0005619796875_r8 * alev(k)) * COS(aphi(j))**4 +  &
                   (14.7803906250_r8 - 0.0002956078125_r8 * alev(k)) * COS(aphi(j))**6) *  &
                  COS(aphi(j)) * SIGN(1._r8, SIN(aphi(j))) *  &
                  COS(3._r8 * alam(i))
            END DO
         END DO
      END DO
!
      RETURN
   END SUBROUTINE Deriv111Anal
!
!-------------------------------------------------------------------------------
!  SUBROUTINE DerivChk
!-------------------------------------------------------------------------------
!
   SUBROUTINE DerivChk(mlon, mlat, mlev, vname, drv, drvanal)
!
      IMPLICIT NONE
!
      INTEGER(i4),                                 INTENT(IN) :: mlon, mlat, mlev
      CHARACTER(LEN=*),                            INTENT(IN) :: vname
      REAL(r8),         DIMENSION(mlon,mlat,mlev), INTENT(IN) :: drv, drvanal
!
      REAL(r8) :: drverr
!!!
!
      drverr = SUM(ABS(drv-drvanal))/DBLE(mlon*mlat*mlev)
      WRITE(6,'(A,3(F15.7,1X))') 'Discrete   '//vname//' : ',  &
         MINVAL(drv    ), MAXVAL(drv    ), SUM(drv    )/DBLE(mlon*mlat*mlev)
      WRITE(6,'(A,3(F15.7,1X))') 'Continuous '//vname//' : ',  &
         MINVAL(drvanal), MAXVAL(drvanal), SUM(drvanal)/DBLE(mlon*mlat*mlev)
      WRITE(6,'(A,F15.7)')       'Error      '//vname//' : ', drverr
!
      RETURN
   END SUBROUTINE DerivChk
!
END

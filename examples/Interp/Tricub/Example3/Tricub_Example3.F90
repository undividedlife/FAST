PROGRAM Tricub_Example3
!
   USE Base,   ONLY: i4, r4, r8
   USE TriCub
!
   REAL(r8),    PARAMETER :: MX = 10.0_r8
   INTEGER(i4), PARAMETER :: NT = 101
!
   CALL Tricub_Init
!
   WRITE(6,'(A)') 'LIBTRICUBIC TEST PROGRAM'
   WRITE(6,'(A)') '************************'
   WRITE(6,'(A)') 'libtricubic v'//TRIM(Tricub_Vers())
   WRITE(6,'(A)') '************************'
   CALL Test1()
   WRITE(6,'(A)') '************************'
   CALL Test2()
   WRITE(6,'(A)') '************************'
   WRITE(6,'(A)') 'LIBTRICUBIC END OF TESTS'
!
CONTAINS
!
   SUBROUTINE Test1()
!
      IMPLICIT NONE
      INTEGER(i4) :: i
      REAL(r8) :: v1, v2
      REAL(r8) :: x, y, z
      REAL(r8), DIMENSION(8) :: f, dfdx, dfdy, dfdz
      REAL(r8), DIMENSION(8) :: d2fdxdy, d2fdxdz
      REAL(r8), DIMENSION(8) :: d2fdydz, d3fdxdydz 
      REAL(r8), DIMENSION(8,7) :: deriv
      REAL(r8), DIMENSION(64) :: a
      REAL(r8)               :: randomval
!
      EXTERNAL randomval
!
!     Input data are generated using cpp rand function
!
      DO i = 1, 8
         f(i) = -MX+2*MX*randomval()
         dfdx(i) = -MX+2*MX*randomval();      deriv(i,1) = dfdx(i)
         dfdy(i) = -MX+2*MX*randomval();      deriv(i,2) = dfdy(i)
         dfdz(i) = -MX+2*MX*randomval();      deriv(i,3) = dfdz(i)
         d2fdxdy(i) = -MX+2*MX*randomval();   deriv(i,4) = d2fdxdy(i)
         d2fdxdz(i) = -MX+2*MX*randomval();   deriv(i,5) = d2fdxdz(i)
         d2fdydz(i) = -MX+2*MX*randomval();   deriv(i,6) = d2fdydz(i)
         d3fdxdydz(i) = -MX+2*MX*randomval(); deriv(i,7) = d3fdxdydz(i)
      END DO
!
      CALL Tricub_Coef(a,f,deriv)
!
      WRITE(6,'(A)') 'TESTING F VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = f(i)
         v2 = Tricub_Eval(a,x,y,z)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO 
!
      WRITE(6,'(A)') 'TESTING DFDX VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdx(i)
         v2 = Tricub_Eval(a,x,y,z,1,0,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO 
!
      WRITE(6,'(A)') 'TESTING DFDY VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdy(i)
         v2 = Tricub_Eval(a,x,y,z,0,1,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      WRITE(6,'(A)') 'TESTING DFDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdz(i)
         v2 = Tricub_Eval(a,x,y,z,0,0,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDXDY VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdxdy(i)
         v2 = Tricub_Eval(a,x,y,z,1,1,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDXDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdxdz(i)
         v2 = Tricub_Eval(a,x,y,z,1,0,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDYDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdydz(i)
         v2 = Tricub_Eval(a,x,y,z,0,1,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      WRITE(6,'(A)') 'TESTING D3FDXDYDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d3fdxdydz(i)
         v2 = Tricub_Eval(a,x,y,z,1,1,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
         IF (ABS(v1-v2) > 1.0e-10_r8) RETURN
      END DO
!
      RETURN
   END SUBROUTINE Test1
!
   SUBROUTINE Test2()
!
      IMPLICIT NONE
!
      REAL(r8), DIMENSION(8) :: f1, df1dx, df1dy, df1dz
      REAL(r8), DIMENSION(8) :: f2, df2dx, df2dy, df2dz
      REAL(r8), DIMENSION(8) :: d2f1dxdy, d2f1dxdz, d2f1dydz, d3f1dxdydz
      REAL(r8), DIMENSION(8) :: d2f2dxdy, d2f2dxdz, d2f2dydz, d3f2dxdydz
      REAL(r8), DIMENSION(8,7) :: deriv1, deriv2
      REAL(r8), DIMENSION(64) :: a1, a2
      REAL(r8) :: rho,v1,v2,x,y,z
      INTEGER(i4), DIMENSION(4) :: iarr1,iarr2
      INTEGER(i4) :: i,j
!
      REAL(r8) :: randomval
      EXTERNAL randomval
!!!
!
      DO i = 1, 8
         f1(i) = -MX+2*MX*randomval()
         df1dx(i) = -MX+2*MX*randomval()
         df1dy(i) = -MX+2*MX*randomval()
         df1dz(i) = -MX+2*MX*randomval()
         d2f1dxdy(i) = -MX+2*MX*randomval()
         d2f1dxdz(i) = -MX+2*MX*randomval()
         d2f1dydz(i) = -MX+2*MX*randomval()
         d3f1dxdydz(i) = -MX+2*MX*randomval()
         f2(i) = -MX+2*MX*randomval()
         df2dx(i) = -MX+2*MX*randomval()
         df2dy(i) = -MX+2*MX*randomval()
         df2dz(i) = -MX+2*MX*randomval()
         d2f2dxdy(i) = -MX+2*MX*randomval()
         d2f2dxdz(i) = -MX+2*MX*randomval()
         d2f2dydz(i) = -MX+2*MX*randomval()
         d3f2dxdydz(i) = -MX+2*MX*randomval()
      END DO
!
      iarr1(1:4) = (/4,5,6,7/) + 1
      iarr2(1:4) = (/0,1,2,3/) + 1
!
      DO i = 1, 4
         f1(iarr1(i)) = f2(iarr2(i))
         df1dx(iarr1(i)) = df2dx(iarr2(i))
         df1dy(iarr1(i)) = df2dy(iarr2(i))
         df1dz(iarr1(i)) = df2dz(iarr2(i))
         d2f1dxdy(iarr1(i)) = d2f2dxdy(iarr2(i))
         d2f1dxdz(iarr1(i)) = d2f2dxdz(iarr2(i))
         d2f1dydz(iarr1(i)) = d2f2dydz(iarr2(i))
         d3f1dxdydz(iarr1(i)) = d3f2dxdydz(iarr2(i))
      END DO
!
      deriv1(:,1) = df1dx(:)
      deriv1(:,2) = df1dy(:)
      deriv1(:,3) = df1dz(:)
      deriv1(:,4) = d2f1dxdy(:)
      deriv1(:,5) = d2f1dxdz(:)
      deriv1(:,6) = d2f1dydz(:)
      deriv1(:,7) = d3f1dxdydz(:)
!
      deriv2(:,1) = df2dx(:)
      deriv2(:,2) = df2dy(:)
      deriv2(:,3) = df2dz(:)
      deriv2(:,4) = d2f2dxdy(:)
      deriv2(:,5) = d2f2dxdz(:)
      deriv2(:,6) = d2f2dydz(:)
      deriv2(:,7) = d3f2dxdydz(:)
!
      CALL Tricub_Coef(a1,f1,deriv1)
      CALL Tricub_Coef(a2,f2,deriv2)
!
      WRITE(6,'(A)') 'CONTINUITY CHECK...'
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on F:        ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)  
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,0,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,0,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on DFDX:     ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,1,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,1,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on DFDY:     ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,0,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,0,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on DFDZ:     ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,1,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,1,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDXDY:  ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,0,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,0,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDXDZ:  ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,1,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,1,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDYDZ:  ', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,1,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,1,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDXDYDZ:', SQRT(rho)
      IF (SQRT(rho) > 1.0e-8_r8) RETURN
!
      WRITE(6,'(A)') 'THESE ARE NOT NECESSARILY CONTINUOUS...'
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,2,0,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,2,0,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDX2   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
         z = 1.0_r8
         v1 = Tricub_Eval(a1,x,y,z,0,2,0)
         z = 0.0_r8
         v2 = Tricub_Eval(a2,x,y,z,0,2,0)
         IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDY2   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,0,2)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,0,2)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDZ2   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,3,0,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,3,0,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDX3   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,3,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,3,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDY3   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,0,3)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,0,3)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDZ3   :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,2,1,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,2,1,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDX2DY :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,2,0,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,2,0,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDX2DZ :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,2,0)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,2,0)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDXDY2 :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,1,0,2)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,1,0,2)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDXDZ2 :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT 
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,2,1)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,2,1)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDY2DZ :', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = DBLE(i-1)/DBLE(NT-1)
            y = DBLE(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = Tricub_Eval(a1,x,y,z,0,1,2)
            z = 0.0_r8
            v2 = Tricub_Eval(a2,x,y,z,0,1,2)
            IF (rho < (v1-v2)**2) rho = (v1-v2)**2
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D3FDYDZ2 :', SQRT(rho)
!
      RETURN
   END SUBROUTINE Test2
!
END

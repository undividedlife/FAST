PROGRAM Tricub_Example2
!
   USE Base,   ONLY: i4, r4, r8
   USE TriCub
!
   REAL(r8),    PARAMETER :: MX = 10.0_r8
   INTEGER(i4), PARAMETER :: NT = 101
!
   REAL(r8), DIMENSION(0:3,0:3) :: monoa = &
                         RESHAPE( (/ 1._r8, 0._r8, 0._r8, 0._r8, &
                                     0._r8, 0._r8, 1._r8, 0._r8, &
                                    -3._r8, 3._r8,-2._r8,-1._r8, &
                                     2._r8,-2._r8, 1._r8, 1._r8/), (/4,4/))
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
      REAL(r8), dimension(8) :: f, dfdx, dfdy, dfdz
      REAL(r8), dimension(8) :: d2fdxdy, d2fdxdz
      REAL(r8), dimension(8) :: d2fdydz, d3fdxdydz
      REAL(r8), dimension(8) :: a
      REAL(r8)               :: randomval
!
      EXTERNAL randomval
!
!     Input data are generated using cpp rand function
!
      DO i = 1, 8
         f(i) = -MX+2*MX*randomval()
         dfdx(i) = -MX+2*MX*randomval()
         dfdy(i) = -MX+2*MX*randomval()
         dfdz(i) = -MX+2*MX*randomval()
         d2fdxdy(i) = -MX+2*MX*randomval()
         d2fdxdz(i) = -MX+2*MX*randomval()
         d2fdydz(i) = -MX+2*MX*randomval()
         d3fdxdydz(i) = -MX+2*MX*randomval()
      END DO
!
      WRITE(6,'(A)') 'TESTING F VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = f(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,0,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO 
!
      WRITE(6,'(A)') 'TESTING DFDX VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdx(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,0,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO 
!
      WRITE(6,'(A)') 'TESTING DFDY VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdy(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,1,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO
!
      WRITE(6,'(A)') 'TESTING DFDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = dfdz(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,0,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDXDY VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdxdy(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,1,0)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDXDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdxdz(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,0,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO
!
      WRITE(6,'(A)') 'TESTING D2FDYDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d2fdydz(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,0,1,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
      END DO
!
      WRITE(6,'(A)') 'TESTING D3FDXDYDZ VALUES...'
      DO i = 1, 8
         CALL Tricub_PXYZ(i-1, x, y, z)
         v1 = d3fdxdydz(i)
         v2 = MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z,1,1,1)
         WRITE(6,'(I1,6X,F10.6,6X,F10.6,6X,A,E15.7)') i, v1, v2, 'Error=', ABS(v1-v2)
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
      REAL(r8) :: rho
      REAL(r8) :: v1,v2
      REAL(r8) :: x,y,z
      INTEGER(i4), DIMENSION(4) :: iarr1,iarr2
      INTEGER(i4)  :: i,j
!
      REAL(r8) :: randomval
!
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
      WRITE(6,'(A)') 'CONTINUITY CHECK...'
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,0,0)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,0,0)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on F:        ', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)  
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,0,0)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,0,0)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on DFDX:     ', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,1,0)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,1,0)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on DFDY:     ', SQRT(rho)
!
      rho = 0.0_r8
      OPEN(10, FILE='testmonoc.dat', STATUS='unknown')
      WRITE(10,'(A)') 'VARIABLE= x y v1 v2 d'
      WRITE(10,'(A,I3,A,I3)') 'ZONE I=', NT, ' J=', NT
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,0,1)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,0,1)
            WRITE(10,'(2(F5.2,1X),2(F8.5,1X),E13.5)') x, y, v1,v2, v1-v2
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      CLOSE(10)
      WRITE(6,'(A,E15.7)') '  C1-error on DFDZ:     ', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,1,0)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,1,0)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDXDY:  ', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,1,0,1)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,1,0,1)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDXDZ:  ', SQRT(rho)
!
      rho = 0.0_r8
      DO j = 1, NT
         DO i = 1, NT
            x = 0._r8+1._r8*(i-1)/DBLE(NT-1)
            y = 0._r8+1._r8*(j-1)/DBLE(NT-1)
            z = 1.0_r8
            v1 = MonoCubic(f1,df1dx,df1dy,df1dz,d2f1dxdy,d2f1dxdz,d2f1dydz,d3f1dxdydz,x,y,z,0,1,1)
            z = 0.0_r8
            v2 = MonoCubic(f2,df2dx,df2dy,df2dz,d2f2dxdy,d2f2dxdz,d2f2dydz,d3f2dxdydz,x,y,z,0,1,1)
            IF (rho < (v1-v2)*(v1-v2)) rho = (v1-v2)*(v1-v2)
         END DO
      END DO
      WRITE(6,'(A,E15.7)') '  C1-error on D2FDYDZ:  ', SQRT(rho)
!
      RETURN
   END SUBROUTINE Test2
!
   FUNCTION MonoCubic_Run(a, x) RESULT(ret)
!
      IMPLICIT NONE
!
      REAL(r8), DIMENSION(0:3), INTENT(IN) :: a
      REAL(r8),                 INTENT(IN) :: x
      REAL(r8) :: ret
      INTEGER(i4) :: i
!!!
!
      ret = 0._r8
      DO i = 0, 3
         IF (i /= 0) then
            ret = ret + a(i)*(x**DBLE(i))
         ELSE
            ret = ret + a(i)
         END IF
      END DO
!
      RETURN
   END FUNCTION MonoCubic_Run
!
   SUBROUTINE MonoCubic_Mat(a, x)
!
      IMPLICIT NONE
!
      REAL(r8), DIMENSION(0:3), INTENT(INOUT) :: a
      REAL(r8), DIMENSION(0:3), INTENT(IN)    :: x
      INTEGER(i4) :: i, j
!!!
!
      DO i = 0, 3
         a(i) = 0.0_r8
         DO j = 0, 3
            a(i) = a(i) + monoa(j,i)*x(j)
         END DO
      end do
!
      RETURN
   END SUBROUTINE MonoCubic_Mat
!
   RECURSIVE FUNCTION MonoCubic(f, dfdx, dfdy, dfdz,  &
                                d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz,  &
                                x, y, z, dx, dy, dz) RESULT(ret)
!
      IMPLICIT NONE
!
      REAL(r8),    DIMENSION(0:7), INTENT(IN) :: f, dfdx, dfdy, dfdz
      REAL(r8),    DIMENSION(0:7), INTENT(IN) :: d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz
      REAL(r8),                    INTENT(IN) :: x, y, z
      INTEGER(i4),                 INTENT(IN) :: dx, dy, dz
      REAL(r8), DIMENSION(0:3) :: fx, dfdyx, dfdzx, coeff, val
      REAL(r8), DIMENSION(0:2) :: fy, dfdzy, d2fdydzx
      REAL(r8) :: dq, ret
      INTEGER(i4) :: i
!!!
!
      dq = 1.0e-5_r8
      IF (dx > 0) THEN
         ret = (MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x+dq,y,z,dx-1,dy,dz)- &
                MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x-dq,y,z,dx-1,dy,dz))/(2._r8*dq)
         RETURN
      END IF
      IF (dy > 0) THEN
         ret = (MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y+dq,z,dx,dy-1,dz)- &
                MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y-dq,z,dx,dy-1,dz))/(2._r8*dq)
         RETURN
      END IF
      IF (dz > 0) THEN
         ret = (MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z+dq,dx,dy,dz-1)- &
                MonoCubic(f,dfdx,dfdy,dfdz,d2fdxdy,d2fdxdz,d2fdydz,d3fdxdydz,x,y,z-dq,dx,dy,dz-1))/(2._r8*dq)
         RETURN
      END IF
!
!     First we interpolate 4 values wrt x
!
      DO i = 0, 3
         val(0) = f(2*i)
         val(1) = f(2*i+1)
         val(2) = dfdx(2*i)
         val(3) = dfdx(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         fx(i) = MonoCubic_Run(coeff,x)
      END DO
!
!     Next we interpolate 4 values dfdy wrt x
!
      DO i = 0, 3
         val(0) = dfdy(2*i)
         val(1) = dfdy(2*i+1)
         val(2) = d2fdxdy(2*i)
         val(3) = d2fdxdy(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         dfdyx(i) = MonoCubic_Run(coeff,x)
      END DO
!
!     Next we interpolate 4 values dfdz wrt x
!
      DO i = 0, 3
         val(0) = dfdz(2*i)
         val(1) = dfdz(2*i+1)
         val(2) = d2fdxdz(2*i)
         val(3) = d2fdxdz(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         dfdzx(i) = MonoCubic_Run(coeff,x)
      END DO
!
!     Next we interpolate 4 values d2fdydz wrt x
!
      DO i = 0, 3
         val(0) = d2fdydz(2*i)
         val(1) = d2fdydz(2*i+1)
         val(2) = d3fdxdydz(2*i)
         val(3) = d3fdxdydz(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         d2fdydzx(i) = MonoCubic_Run(coeff,x)
      END DO
!
!     Next we interpolate 2 values of f wrt y
!
      DO i = 0, 1
         val(0) = fx(2*i)
         val(1) = fx(2*i+1)
         val(2) = dfdyx(2*i)
         val(3) = dfdyx(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         fy(i) = MonoCubic_Run(coeff,y)
      END DO
!
!     Next we interpolate 2 values of dfdz wrt y
!
      DO i = 0, 1
         val(0) = dfdzx(2*i)
         val(1) = dfdzx(2*i+1)
         val(2) = d2fdydzx(2*i)
         val(3) = d2fdydzx(2*i+1)
         CALL MonoCubic_Mat(coeff,val)
         dfdzy(i) = MonoCubic_Run(coeff,y)
      END DO
!
!     Finally interpolation of f wrt z
!
      val(0) = fy(0)
      val(1) = fy(1)
      val(2) = dfdzy(0)
      val(3) = dfdzy(1)
      CALL MonoCubic_Mat(coeff,val) 
      ret = MonoCubic_Run(coeff,z)
!
      RETURN
   END FUNCTION MonoCubic  
!
END

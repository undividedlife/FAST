PROGRAM Tricub_Example1
!
!-------------------------------------------------------------------------------
!  example1.cpp : illustrates the use of libtricubic
!  Francois Lekien <lekien@mit.edu> 2004-01-20
!  For FORTRAN, In-Sun Song <insun.song@gmail.com>
!-------------------------------------------------------------------------------
!
   USE Base,   ONLY: r8
!
!-------------------------------------------------------------------------------
!  Required include file: tricubic.h
!  If tricubic.h has not been installed in a directory accessible by
!  the compiler, use -I/path/to/tricubic
!  For FORTRAN, use -I/path/to/tricubic/module
!-------------------------------------------------------------------------------
!
   USE Tricub, ONLY: Tricub_Init, Tricub_Coef, Tricub_Eval
!
   IMPLICIT NONE
!
!  Define the box
!  Tricubic is written for cubes of side 1. Multiplications and divisions
!  are needed along the way for rectangular box with arbitrary sides.
!  See below for details.
!
   REAL(r8), PARAMETER :: DX = 0.2_r8
   REAL(r8), PARAMETER :: DY = 2.0_r8
   REAL(r8), PARAMETER :: DZ = 0.3_r8
!
   REAL(r8), DIMENSION(64) :: a
   REAL(r8) :: f1, f2, dfdx, d2fdxdy
   REAL(r8) :: x, y, z
!
   REAL(r8), DIMENSION(8) :: fval
   REAL(r8), DIMENSION(8) :: dfdxval,dfdyval,dfdzval
   REAL(r8), DIMENSION(8) :: d2fdxdyval,d2fdxdzval,d2fdydzval
   REAL(r8), DIMENSION(8) :: d3fdxdydzval
   REAL(r8), DIMENSION(8,7) :: deriv
!!!
!
!  These are the 8 functions that need to be known at the 8 corners.
!  The functions are f, the three first derivatives dfdx, dfdy, dfdz,
!  the 3 mixed 2nd order derivatives and the mixed 3rd order derivative.
!
!  The derivatives can be obtained by numerical differentiation of f at the
!  other corners of the grid. See Numerical Recipes for examples in 2D.
!
!  The order of the points is as follow:
!  0: x=0, y=0, z=0
!  1: x=1, y=0, z=0
!  2: x=0, y=1, z=0
!  3: x=1, y=1, z=0
!  4: x=0, y=0, z=1
!  5: x=1, y=0, z=1
!  6: x=0, y=1, z=1
!  7: x=1, y=1, z=1
!
!  For convenience, the ordering of the points is available at run time
!  using tricubic_pointID2xyz() (for FORTRAN, using TCPID2XYZ).
!
   fval         = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   dfdxval      = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   dfdyval      = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   dfdzval      = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   d2fdxdyval   = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   d2fdxdzval   = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   d2fdydzval   = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
   d3fdxdydzval = (/1.2_r8, 2.3_r8, 3.4_r8, 4.5_r8, 5.6_r8, 6.7_r8, 7.8_r8, 8.9_r8/)
!
!  First we compute the 64 coefficients for the cube
!
   CALL Tricub_Init
!
!  The first step is to scale the derivatives that have been computed
!  in a rectangular box instead of the cube of side 1.
!
!  Notice that this step can be avoided by computing the numerical derivatives
!  without dividing by the length of the boxes.
!
   fval    = fval*1._r8
   dfdxval = dfdxval*DX
   dfdyval = dfdyval*DY
   dfdzval = dfdzval*DZ
   d2fdxdyval = d2fdxdyval*DX*DY
   d2fdxdzval = d2fdxdzval*DX*DZ
   d2fdydzval = d2fdydzval*DY*DZ
   d3fdxdydzval = d3fdxdydzval*DX*DY*DZ
!
   deriv(:,1) = dfdxval
   deriv(:,2) = dfdyval
   deriv(:,3) = dfdzval
   deriv(:,4) = d2fdxdyval
   deriv(:,5) = d2fdxdzval
   deriv(:,6) = d2fdydzval
   deriv(:,7) = d3fdxdydzval
!
!  Next we get the set of coefficients for this cube
!
   CALL Tricub_Coef(a, fval, deriv)
!
!  To get the value in the middle of the cube, we always use (.5,.5,.5)
!  (i.e., relative coordinate)
!
   f1 = Tricub_Eval(a, 0.5_r8, 0.5_r8, 0.5_r8)
   WRITE(6,'(A,E20.12)') ' f1 = ', f1
!
   x = DX*0.5
   y = DY*0.5
   z = DZ*0.5  
!
!  To get the value at a point x,y,z (with reference to corner ID 0),
!  we divide by each legnth
!
   f2 = Tricub_Eval(a, x/dx, y/dy, z/dz)
   WRITE(6,'(A,E20.12)') ' f2 = ', f2
!
!  Derivatives can be computed similarly but need to be scaled by dx,dy,dz
!
   dfdx = Tricub_Eval(a, x/dx, y/dy, z/dz, 1, 0, 0)/dx
   WRITE(6,'(A,E20.12)') ' dfdx = ', dfdx
!
   d2fdxdy = Tricub_Eval(a, x/dx, y/dy, z/dz, 1, 1, 0)/(dx*dy)
   WRITE(6,'(A,E20.12)') ' d2fdxdy = ', d2fdxdy
!
END

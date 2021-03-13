!-------------------------------------------------------------------------------
!
! MODULE Interp
!
!-------------------------------------------------------------------------------
!
MODULE Interp
!
    USE Base,   ONLY: r4, r8, i4
    USE Bisect, ONLY: Bisect_Locgrd
    USE Tricub, ONLY: Tricub_Init, Tricub_Coef, Tricub_Eval
!
    IMPLICIT NONE
!
    PRIVATE
!
!   Interp_Lagr, Interp_Poly
!
!       Lagrange and polynomial interpolations are adapted from
!       J. Caleb Wherry's F90-Numerical-Lib on github
!       https://github.com/calebwherry/F90-Numerical-Lib
!
!   Interp_Rgrd
!
!       A suite (REGRIDPACK) of Fortran routines of interpolation for
!       one-, two-, three-, and four-dimensional arrays on uniform or
!       nonuniform orthogonal grids. This operation is commonlly referred to as
!       "regridding". Linear or cubic interpolation can be selected
!       independently in each dimension. Extrapolation is not allowed.
!
!       John C. Adams at NCAR : 1997: Original REGRIDPACK
!       Jacob Williams        : 2019: Modernized and refactored
!                               https://github.com/jacobwilliams/regridpack
!
!   Interp_Cspl
!
!       de Boor's cubic spline algorithm (https://www.netlib.org/pppack)
!       From a practical guide to splines by C. de Boor
!
!   Interp_Tcub 
!
!       Tricubic interpolation routine refactored from C++ routine written
!       by Lekien and Marsen (2005) in Internaltional Journal for Numerical 
!       Methods in Engineering
!       
!
    INTERFACE Interp_Lagr
        MODULE PROCEDURE Interp_Lagr_r8
    END INTERFACE
!
    INTERFACE Interp_Poly
        MODULE PROCEDURE Interp_Poly_r8
    END INTERFACE
!
    INTERFACE Interp_Rgrd
        MODULE PROCEDURE Interp_Rgrd1_r8
!!!!    MODULE PROCEDURE Interp_Rgrd2_r8
!!!!    MODULE PROCEDURE Interp_Rgrd3_r8
!!!!    MODULE PROCEDURE Interp_Rgrd4_r8
    END INTERFACE
!
    INTERFACE Interp_Cspl
        MODULE PROCEDURE Interp_Cspl_r8
    END INTERFACE
!
    INTERFACE Interp_Tcub
        MODULE PROCEDURE Interp_Tcub_r8
    END INTERFACE
!
    PUBLIC :: Interp_Lagr
    PUBLIC :: Interp_Poly
    PUBLIC :: Interp_Rgrd
    PUBLIC :: Interp_Cspl
    PUBLIC :: Interp_Tcub
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Interp_Lagr
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Interp_Lagr_r8(n, xgrd, ygrd, x, y)
!
        IMPLICIT NONE
!
        INTEGER(i4),                 INTENT(IN) :: n
        REAL(r8),    DIMENSION(0:n), INTENT(IN) :: xgrd, ygrd
        REAL(r8),    DIMENSION(:),   INTENT(IN) :: x
        REAL(r8),    DIMENSION(:),   INTENT(OUT) :: y
!
        REAL(r8) :: lnk
        INTEGER(i4) :: nnew
        INTEGER(i4) :: j, k, l
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Interp_Lagr_r8)'
!!!
!
        nnew = SIZE(x)
!
        DO l = 1, nnew
            y(l) = 0._r8
            DO k = 0, n
                lnk = 1._r8
                DO j = 0, n
                    IF (j /= k) THEN
                        lnk = lnk * (x(l) - xgrd(j))/(xgrd(k) - xgrd(j))
                    END IF
                END DO
                y(l) = y(l) + ygrd(k) * lnk
            END DO
        END DO
!
        RETURN
    END SUBROUTINE Interp_Lagr_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Interp_Poly
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Interp_Poly_r8(n, xgrd, ygrd, x, y)
!
        IMPLICIT NONE
!
        INTEGER(i4),                 INTENT(IN) :: n
        REAL(r8),    DIMENSION(0:n), INTENT(IN) :: xgrd, ygrd
        REAL(r8),    DIMENSION(:),   INTENT(IN) :: x
        REAL(r8),    DIMENSION(:),   INTENT(OUT) :: y
!
        REAL(r8), DIMENSION(0:n,0:n) :: a
        INTEGER(i4) :: nnew
        INTEGER(i4) :: i, k, l
        REAL(r8) :: lnk
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Interp_Poly_r8)'
!!!
!
        nnew = SIZE(x)
!
        DO l = 1, nnew
            a(0:n,0:n) = 1._r8
            a(0:n,0)   = ygrd(0:n)
            DO k = 1, n
                DO i = k, n
                    a(i,k) = ((x(l)-xgrd(i-k))*a(i  ,k-1)-  &
                              (x(l)-xgrd(i  ))*a(i-1,k-1))/(xgrd(i)-xgrd(i-k))
                END DO
            END DO
            y(l) = a(n,n)
        END DO
!
        RETURN
    END SUBROUTINE Interp_Poly_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Interp_Rgrd1
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Interp_Rgrd1_r8(x, p, xx, q, intpol, ier)
!
        IMPLICIT NONE
!
        REAL(r8),    DIMENSION(:), INTENT(IN)  :: x, p  ! original x and values
        REAL(r8),    DIMENSION(:), INTENT(IN)  :: xx    ! target grid (xx)
        REAL(r8),    DIMENSION(:), INTENT(OUT) :: q     ! regridded values at xx
        INTEGER(i4),               INTENT(IN)  :: intpol
        INTEGER(i4),               INTENT(OUT) :: ier   ! status code
                                                        ! 0    : no error
                                                        ! 1-6  : error
                                                        ! 10   : input vectors
                                                        !        are wrong size
                                                        ! 100  : output memory
!
        INTEGER(i4) :: lw, liw
        INTEGER(i4) :: nx, mx
        INTEGER(i4) :: np, nq
        REAL(r8),    DIMENSION(:), ALLOCATABLE :: w
        INTEGER(i4), DIMENSION(:), ALLOCATABLE :: iw
        INTEGER(i4) :: ierr1, ierr2
!!!
!
!       Get array sizes:
!
        nx = SIZE(x)
        np = SIZE(p)
!
        mx = size(xx)
        nq = size(q)
!
        IF (nx /= np .OR. mx /= nq) THEN ! Error: vectors are the wrong size
            ier = 10
            RETURN
        END IF
!
!       Allocate work matrices:
!
        SELECT CASE(intpol)
        CASE(1)
            lw = mx
        CASE(3)
            lw = 4*mx
        CASE DEFAULT
            ier = 6     ! Error: invalid intpol value
            RETURN
        END SELECT
!
        liw = mx
!
        ALLOCATE(w (lw),  STAT=ierr1)
        ALLOCATE(iw(liw), STAT=ierr2)
!
        IF (ierr1 == 0 .AND. ierr2 == 0) THEN ! call the main routine:
            CALL Rgrd1_r8(nx,x,p,mx,xx,q,intpol,w,lw,iw,liw,ier)
        ELSE !error: out of memory
            ier = 100
        END IF
!
!       Clean up:
!
        IF (ALLOCATED(w))  DEALLOCATE(w)
        IF (ALLOCATED(iw)) DEALLOCATE(iw)
!
        RETURN
    END SUBROUTINE Interp_Rgrd1_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Rgrd1
!
!   SUBROUTINE Rgrd1 interpolates the values p(i) on the grid x(i)
!   for i=1,...,nx onto q(ii) on the grid xx(ii),ii=1,...,mx.
!
!   Requirements
!
!   x must be a strictly increasing grid and xx must be an increasing
!   grid (see ier = 4).  in addition the interval
!
!   [xx(1),xx(mx)]
!
!   must lie within the interval
!
!   [x(1),x(nx)].
!
!   extrapolation is not allowed (see ier=3).  if these intervals
!   are identical and the x and xx grids are UNIFORM then subroutine
!   rgrd1u should be used in place of rgrd1.
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Rgrd1_r8(nx,x,p,mx,xx,q,intpol,w,lw,iw,liw,ier)
!
        IMPLICIT NONE
!
        INTEGER(i4),                 INTENT(IN)    :: nx
        REAL(r8),    DIMENSION(nx),  INTENT(IN)    :: x
        REAL(r8),    DIMENSION(nx),  INTENT(IN)    :: p
        INTEGER(i4),                 INTENT(IN)    :: mx
        REAL(r8),    DIMENSION(mx),  INTENT(IN)    :: xx
        INTEGER(i4),                 INTENT(IN)    :: intpol  ! 1 (lin) or 3 (cub)
        REAL(r8),                    INTENT(OUT)   :: q(mx)
        INTEGER(i4),                 INTENT(IN)    :: lw 
        REAL(r8),    DIMENSION(lw),  INTENT(INOUT) :: w
        INTEGER(i4),                 INTENT(IN)    :: liw
        INTEGER(i4), DIMENSION(liw), INTENT(INOUT) :: iw
        INTEGER(i4),                 INTENT(OUT)   :: ier
!
        INTEGER(i4) :: i,ii,i1,i2,i3,i4
!!!
!
!       Check arguments for errors
!
!       check xx grid resolution
        ier = 1
        IF (mx < 1) RETURN
!
!       Check intpol
!
        ier = 6
        IF (intpol /= 1 .AND. intpol /= 3) RETURN
!
!       Check x grid resolution
!
        ier = 2
        IF (intpol == 1 .AND. nx < 2) RETURN
        IF (intpol == 3 .AND. nx < 4) RETURN
!
!       Check xx grid contained in x grid
!
        ier = 3
        IF (xx(1) < x(1) .OR. xx(mx) > x(nx)) RETURN
!
!       Check montonicity of grids
!
        DO i = 2, nx
            IF (x(i-1) >= x(i)) THEN
                ier = 4
                RETURN
            END IF
        END DO
        DO ii = 2, mx
            IF (xx(ii-1) > xx(ii)) THEN
                 ier = 4
                 RETURN
            END IF
        END DO
!
!       Check minimum work space lengths
!
        ier = 5
        IF (intpol == 1) THEN
            IF (lw < mx) RETURN
        ELSE
            IF (lw < 4*mx) RETURN
        END IF
        IF (liw < mx) RETURN
!
!       Arguments o.k.
!
        ier = 0
!
        IF (intpol == 1) THEN
!       ! linear interpolation in x
            CALL Linmx_r8(nx,x,mx,xx,iw,w)
            CALL Lint1_r8(nx,p,mx,q,iw,w)
        ELSE
!       ! cubic interpolation in x
            i1 = 1
            i2 = i1+mx
            i3 = i2+mx
            i4 = i3+mx
            CALL Cubnmx_r8(nx,x,mx,xx,iw,w(i1),w(i2),w(i3),w(i4))
            CALL Cubt1_r8(nx,p,mx,q,iw,w(i1),w(i2),w(i3),w(i4))
        END IF
!
    END SUBROUTINE Rgrd1_r8
!
!-------------------------------------------------------------------------------
!   SUBROUTINE Lint1
!-------------------------------------------------------------------------------
!
    SUBROUTINE Lint1_r8(nx,p,mx,q,ix,dx)
!
        IMPLICIT NONE
!
        INTEGER(i4) :: mx,ix(mx),nx,ii,i
        REAL(r8) :: p(nx),q(mx),dx(mx)
!!!
!
        DO ii = 1, mx
            i = ix(ii)
            q(ii) = p(i)+dx(ii)*(p(i+1)-p(i))
        END DO
!
        RETURN
    END SUBROUTINE Lint1_r8
!
!-------------------------------------------------------------------------------
!   SUBROUTINE Cubt1
!-------------------------------------------------------------------------------
!
   SUBROUTINE Cubt1_r8(nx,p,mx,q,ix,dxm,dx,dxp,dxpp)
!
        IMPLICIT NONE
!
        INTEGER(i4) :: mx,ix(mx),nx,i,ii
        REAL(r8) :: p(nx),q(mx),dxm(mx),dx(mx),dxp(mx),dxpp(mx)
!
        DO ii = 1, mx
            i = ix(ii)
            q(ii) = dxm(ii)*p(i-1)+dx(ii)*p(i)+dxp(ii)*p(i+1)+dxpp(ii)*p(i+2)
        END DO
!
        RETURN
    END SUBROUTINE Cubt1_r8
!
!-------------------------------------------------------------------------------
!   SUBROUTINE Linmx
!-------------------------------------------------------------------------------
!
    SUBROUTINE Linmx_r8(nx,x,mx,xx,ix,dx)
!
        IMPLICIT NONE
!
        REAL(r8) :: x(*),xx(*),dx(*)
        INTEGER(i4) :: ix(*),isrt,ii,i,nx,mx
!
        isrt = 1
!
        DO ii = 1, mx
!           Find x(i) s.t. x(i) < xx(ii) <= x(i+1)
            DO i = isrt, nx-1
                IF (x(i+1) >= xx(ii)) THEN
                    isrt = i
                    ix(ii) = i
                    EXIT
                END IF
            END DO
        END DO
!
!       Set linear scale term
!
        DO ii = 1, mx
            i = ix(ii)
            dx(ii) = (xx(ii)-x(i))/(x(i+1)-x(i))
        END DO
!
        RETURN
    END SUBROUTINE Linmx_r8
!
!-------------------------------------------------------------------------------
!   SUBROUTINE Cubnmx
!-------------------------------------------------------------------------------
!
    SUBROUTINE Cubnmx_r8(nx,x,mx,xx,ix,dxm,dx,dxp,dxpp)
!
        IMPLICIT NONE
!
        REAL(r8) :: x(*),xx(*),dxm(*),dx(*),dxp(*),dxpp(*)
        INTEGER(i4) :: ix(*),mx,nx,i,ii,isrt
!
        isrt = 1
!
        DO ii = 1, mx
!           Set i in [2,nx-2] closest s.t.
!           x(i-1),x(i),x(i+1),x(i+2) can interpolate xx(ii)
            DO i=isrt,nx-1
                IF (x(i+1) >= xx(ii)) THEN
                    ix(ii) = MIN(nx-2,MAX(2,i))
                    isrt = ix(ii)
                    EXIT
                END IF
            END DO
        END DO
!
!       Set cubic scale terms
!
        DO ii=1,mx
            i = ix(ii)
            dxm (ii) = (xx(ii)-x(i))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/  &
                       ((x(i-1)-x(i))*(x(i-1)-x(i+1))*(x(i-1)-x(i+2)))
            dx  (ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i+1))*(xx(ii)-x(i+2))/  &
                       ((x(i)-x(i-1))*(x(i)-x(i+1))*(x(i)-x(i+2)))
            dxp (ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+2))/  &
                       ((x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i+1)-x(i+2)))
            dxpp(ii) = (xx(ii)-x(i-1))*(xx(ii)-x(i))*(xx(ii)-x(i+1))/  &
                       ((x(i+2)-x(i-1))*(x(i+2)-x(i))*(x(i+2)-x(i+1)))
        END DO
!
        RETURN
    END SUBROUTINE Cubnmx_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Interp_Cspl
!
!   This routine is obtained from PPPACK on netlib.org
!   From a practical guide to splines by C. de Boor
!    
!   n      : number of data points. assumed to be >= 2.
!   x      : Absciaase of data points (strictly increasing)
!   y      : Ordinates of data points
!   ybc1   : Boundary condition at x(1)
!   ybcn   : Boundary condition at x(n)
!   ibcbeg : Boundary condition indicator at x(1)
!   ibcend : Boundary condition indicator at x(n)
!
!   ibcbeg = 0 : No boundary condition at x(1) is given.
!                In this case, the not-a-knot condition is used, i.e. the
!                jump in the third derivative across x(2) is forced to
!                zero, thus the first and the second cubic polynomial pieces
!                are made to coincide.
!   ibcbeg = 1 : The slope at x(1) is made to equal ybc1, supplied by input.
!   ibcbeg = 2 : The second derivative at x(1) is made to equal ybc1
!                supplied by input.
!   ibcend = 0, 1, 2 : Analogous to ibcbeg = 0, 1, 2, respectively
!
!   On output,
!   
!   c(1,:) = y(:), 
!   c(2,1) = ybc1 and c(2,n) = ybcn
!   c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coeffs of the cubic
!   interpolating spline with interior knots (or joints) x(2), ..., x(n-1),
!   precisely, in the interval x(i), x(i+1)), the spline f is given by
!
!        f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!
!   where h = x - x(i). the function program *ppvalu* may be
!   used to evaluate f or its derivatives from x, c, l = n-1, and k=4.
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Interp_Cspl_r8(n, x, y, ybc1, ybcn, ibcbeg, ibcend,  &
                              nn, xval, yval, jderiv)
!
        IMPLICIT NONE
!
        INTEGER(i4),                INTENT(IN)  :: n
        REAL(r8),    DIMENSION(n),  INTENT(IN)  :: x, y
        REAL(r8),                   INTENT(IN)  :: ybc1, ybcn 
        INTEGER(i4),                INTENT(IN)  :: ibcbeg, ibcend
        INTEGER(i4),                INTENT(IN)  :: nn
        REAL(r8),    DIMENSION(nn), INTENT(IN)  :: xval
        REAL(r8),    DIMENSION(nn), INTENT(OUT) :: yval
        INTEGER(i4),                INTENT(IN)  :: jderiv  ! jderiv-th deriv
!
        REAL(r8), DIMENSION(4,n) :: c
        REAL(r8) :: divdf1, divdf3, dtau, g
        REAL(r8) :: h, fmmjdr, fmmjdr0
        INTEGER(i4) :: i, j, l, m, iflag, ip1
        CHARACTER(LEN=*), PARAMETER :: HEADER = '(Interp_CSpl_r8)'
!!!
!
        IF (n < 2) THEN
            WRITE(6,'(A)') HEADER//': Error: n < 2.'
            STOP
        END IF
!
!       A tridiagonal linear system for the unknown slopes s(i) of f
!       at x(i), i = 1, ..., n, is generated and then solved by gauss
!       elimination, with s(i) ending up in c(2,i), all i.
!       c(3,:) and c(4,:) are used initially for temporary storage.
!
        c(1:4,1:n) = 0._r4
        c(1,1:n) = y(1:n)
        c(2,1) = ybc1
        c(2,n) = ybcn
!
!       Compute first differences of x sequence and store in c(3,:).
!       Also, compute first divided difference of data and store in c(4,:).
!
#if 1
        l = n-1
        DO m = 2, n
            c(3,m) = x(m) - x(m-1)
            c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
        END DO
!
!       Construct first equation from the boundary condition of the form
!       c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
!
        SELECT CASE (ibcbeg-1)
        CASE (:-1)             ! Label 11 ! no condition at left end
            IF (n > 2) THEN    ! not-a-knot condition at left end and n > 2
                c(4,1) = c(3,3)
                c(3,1) = c(3,2) + c(3,3)
                c(2,1) = ((c(3,2)+2._r8*c(3,1))*c(4,2)*c(3,3)+  &
                           c(3,2)**2*c(4,3))/c(3,1)
                DO m = 2, l
                    g = -c(3,m+1)/c(4,m-1)
                    c(2,m) = g*c(2,m-1) + 3._r8*(c(3,m  )*c(4,m+1)+  &
                                                 c(3,m+1)*c(4,m  ))
                    c(4,m) = g*c(3,m-1) + 2._r8*(c(3,m  )+c(3,m+1))
                END DO
!               IF (ibcend-1) 21,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = -1  ! GO TO 21
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            ELSE  ! no condition at left and n = 2
                c(4,1) = 1._r8
                c(3,1) = 1._r8
                c(2,1) = 2._r8*c(4,2)
!               IF (ibcend-1) 26,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = 0   ! GO TO 26
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            END IF
        CASE (0)   ! slope prescribed at left end
            c(4,1) = 1._r8
            c(3,1) = 0._r8
            IF (n /= 2) THEN
                DO m = 2, l
                    g = -c(3,m+1)/c(4,m-1)
                    c(2,m) = g*c(2,m-1) + 3._r8*(c(3,m)*c(4,m+1)+  &
                                                 c(3,m+1)*c(4,m))
                    c(4,m) = g*c(3,m-1) + 2._r8*(c(3,m) + c(3,m+1))
                END DO
!               IF (ibcend-1) 21,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = -1  ! GO TO 21
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            ELSE IF (n == 2) THEN
!               IF (ibcend-1) 26,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = 0   ! GO TO 26
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            END IF
        CASE (1:)  ! second derivative prescribed at left end
            c(4,1) = 2._r8
            c(3,1) = 1._r8
            c(2,1) = 3._r8*c(4,2) - c(3,2)/2._r8*c(2,1)
            IF (n /= 2) THEN
                DO m = 2, l
                    g = -c(3,m+1)/c(4,m-1)
                    c(2,m) = g*c(2,m-1) + 3._r8*(c(3,m  )*c(4,m+1)+  &
                                                 c(3,m+1)*c(4,m  ))
                    c(4,m) = g*c(3,m-1) + 2._r8*(c(3,m) + c(3,m+1))
                END DO
!               IF (ibcend-1) 21,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = -1  ! GO TO 21
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            ELSE IF (n == 2) THEN
!               IF (ibcend-1) 26,30,24
                SELECT CASE (ibcend-1)
                CASE(:-1); iflag = 0   ! GO TO 26
                CASE(0);   iflag = 1   ! GO TO 30
                CASE(1:);  iflag = 2   ! GO TO 24
                END SELECT
            END IF
        END SELECT
!
        SELECT CASE (iflag)
        CASE(-1)   ! label 21
            IF (n == 3 .AND. ibcbeg == 0) THEN
                c(2,n) = 2._r8*c(4,n)
                c(4,n) = 1._r8
                g = -1._r8/c(4,n-1)
                c(4,n) = g*c(3,n-1) + c(4,n)
                c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
                j = l 
                DO 
                    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                    j = j - 1
                    IF (j <= 0) EXIT
                END DO
                DO i = 2, n
                    dtau = c(3,i)
                    divdf1 = (c(1,i) - c(1,i-1))/dtau
                    divdf3 = c(2,i-1) + c(2,i) - 2._r8*divdf1
                    c(3,i-1) = 2._r8*(divdf1 - c(2,i-1) - divdf3)/dtau
                    c(4,i-1) = (divdf3/dtau)*(6._r8/dtau)
                END DO
            ELSE
                g = c(3,n-1) + c(3,n)
                c(2,n) = ((c(3,n)+2._r8*g)*c(4,n)*c(3,n-1)  &
                         + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
                g = -g/c(4,n-1)
                c(4,n) = g*c(3,n-1) + c(4,n)
                c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
                j = l
                DO
                    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                    j = j - 1
                    IF (j <= 0) EXIT         
                END DO
                DO i = 2, n
                    dtau = c(3,i)
                    divdf1 = (c(1,i) - c(1,i-1))/dtau
                    divdf3 = c(2,i-1) + c(2,i) - 2._r4*divdf1
                    c(3,i-1) = 2._r4*(divdf1 - c(2,i-1) - divdf3)/dtau
                    c(4,i-1) = (divdf3/dtau)*(6._r4/dtau)
                END DO
            END IF
        CASE(0)  ! label 26
            IF (ibcbeg > 0) THEN
                c(2,n) = 2._r8*c(4,n)
                c(4,n) = 1._r8
                g = -1._r8/c(4,n-1)
                c(4,n) = g*c(3,n-1) + c(4,n)
                c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
                j = l
                DO
                    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                    j = j - 1
                    IF (j <= 0) EXIT
                END DO
                DO i = 2, n
                    dtau = c(3,i)
                    divdf1 = (c(1,i) - c(1,i-1))/dtau
                    divdf3 = c(2,i-1) + c(2,i) - 2._r4*divdf1
                    c(3,i-1) = 2._r4*(divdf1 - c(2,i-1) - divdf3)/dtau
                    c(4,i-1) = (divdf3/dtau)*(6._r8/dtau)
                END DO
            ELSE
                c(2,n) = c(4,n)
                j = l 
                DO
                    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                    j = j - 1
                    IF (j <= 0) EXIT
                END DO
                DO i = 2, n
                    dtau = c(3,i)
                    divdf1 = (c(1,i) - c(1,i-1))/dtau
                    divdf3 = c(2,i-1) + c(2,i) - 2._r4*divdf1
                    c(3,i-1) = 2._r4*(divdf1 - c(2,i-1) - divdf3)/dtau
                    c(4,i-1) = (divdf3/dtau)*(6._r8/dtau)
                END DO
            END IF
        CASE(1)   ! label 30
            j = l
            DO
                c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                j = j - 1
                IF (j <= 0) EXIT
            END DO
            DO i = 2, n
                dtau = c(3,i)
                divdf1 = (c(1,i) - c(1,i-1))/dtau
                divdf3 = c(2,i-1) + c(2,i) - 2._r8*divdf1
                c(3,i-1) = 2._r8*(divdf1 - c(2,i-1) - divdf3)/dtau
                c(4,i-1) = (divdf3/dtau)*(6._r8/dtau)
            END DO
        CASE(2)   ! label 24
            c(2,n) = 3._r8*c(4,n) + c(3,n)/2._r8*c(2,n)
            c(4,n) = 2._r8
            g = -1._r8/c(4,n-1)
            c(4,n) = g*c(3,n-1) + c(4,n)
            c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
            j = l
            DO  
                c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
                j = j - 1 
                IF (j <= 0) EXIT
            END DO
            DO i = 2, n
                dtau = c(3,i)
                divdf1 = (c(1,i) - c(1,i-1))/dtau
                divdf3 = c(2,i-1) + c(2,i) - 2._r8*divdf1
                c(3,i-1) = 2._r8*(divdf1 - c(2,i-1) - divdf3)/dtau
                c(4,i-1) = (divdf3/dtau)*(6._r8/dtau)
            END DO
        END SELECT
#endif
#if 0
      l = n-1
!compute first differences of tau sequence and store in c(3,.). also,
!compute first divided difference of data and store in c(4,.).
      do 10 m=2,n
         c(3,m) = x(m) - x(m-1)
   10    c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
!construct first equation from the boundary condition, of the form
!c             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      if (ibcbeg-1)                     11,15,16
   11 if (n .gt. 2)                     go to 12
!c     no condition at left end and n = 2.
      c(4,1) = 1.
      c(3,1) = 1.
      c(2,1) = 2.*c(4,2)
                                        go to 25
!c     not-a-knot condition at left end and n .gt. 2.
   12 c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
                                        go to 19
!c     slope prescribed at left end.
   15 c(4,1) = 1.
      c(3,1) = 0.
                                        go to 18
!c     second derivative prescribed at left end.
   16 c(4,1) = 2.
      c(3,1) = 1.
      c(2,1) = 3.*c(4,2) - c(3,2)/2.*c(2,1)
   18 if(n .eq. 2)                      go to 25
!c  if there are interior knots, generate the corresp. equations and car-
!c  ry out the forward pass of gauss elimination, after which the m-th
!c  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   19 do 20 m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1) + 3.*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
   20    c(4,m) = g*c(3,m-1) + 2.*(c(3,m) + c(3,m+1))
!construct last equation from the second boundary condition, of the form
!c           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
!c     if slope is prescribed at right end, one can go directly to back-
!c     substitution, since c array happens to be set up just right for it
!c     at this point.
      if (ibcend-1)                     21,30,24
   21 if (n .eq. 3 .and. ibcbeg .eq. 0) go to 22
!c     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
!c     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.*g)*c(4,n)*c(3,n-1)  &
                  + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
                                        go to 29
!c     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
!c     knot at left end point).
   22 c(2,n) = 2.*c(4,n)
      c(4,n) = 1.
                                        go to 28
!c     second derivative prescribed at right endpoint.
   24 c(2,n) = 3.*c(4,n) + c(3,n)/2.*c(2,n)
      c(4,n) = 2.
                                        go to 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg .gt. 0)                go to 22
!c     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n) = c(4,n)
                                        go to 30
   28 g = -1./c(4,n-1)
!complete forward pass of gauss elimination.
   29 c(4,n) = g*c(3,n-1) + c(4,n)
      c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
!carry out back substitution
   30 j = l 
   40    c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
         j = j - 1
         if (j .gt. 0)                  go to 40
!c****** generate cubic coefficients in each interval, i.e., the deriv.s
!c  at its left endpoint, from value and slope at its endpoints.
      do 50 i=2,n
         dtau = c(3,i)
         divdf1 = (c(1,i) - c(1,i-1))/dtau
         divdf3 = c(2,i-1) + c(2,i) - 2.*divdf1
         c(3,i-1) = 2.*(divdf1 - c(2,i-1) - divdf3)/dtau
   50    c(4,i-1) = (divdf3/dtau)*(6./dtau)
#endif
!
!!!     DO i = 1, n
!!!         WRITE(6,'(4(F15.7,1X))') (c(j,i), j = 1, 4)
!!!     END DO
!
!       End of construction of c matrix
!       Evaluate interpolation using the c matrix
!       Below is obtained from PPVALU.F in PPPACK
!
!       jderiv = 0: Function value
!       jderiv = 1: First-order derivative
!       jderiv = 2: Second-order derivative
!       jderiv = 3: Third-order derivative
!
        yval(1:nn) = 0._r8
!
        fmmjdr0 = DBLE(4 - jderiv)
        IF (fmmjdr0 <= 0._r8) RETURN
!
        DO j = 1, nn
            fmmjdr = fmmjdr0
            CALL Bisect_Locgrd(n, x, xval(j), i, ip1)
            h = xval(j) - x(i)
            !!!PRINT *, j, xval(j), x(i)
            m = 4
            DO 
                yval(j) = (yval(j)/fmmjdr)*h + c(m,i)
                !!!PRINT *, m, yval(j)
                m = m - 1
                fmmjdr = fmmjdr - 1._r8
                IF (fmmjdr <= 0._r8) EXIT
            END DO
        END DO
!
        RETURN
    END SUBROUTINE Interp_Cspl_r8
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Interp_Tcub
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Interp_Tcub_r8
!
        IMPLICIT NONE
!
        RETURN
    END SUBROUTINE Interp_Tcub_r8
!
END

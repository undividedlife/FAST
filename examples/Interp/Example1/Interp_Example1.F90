PROGRAM Test_Interp
!
    USE Base,   ONLY: i4, r4, r8
    USE Interp, ONLY: Interp_Lagr, Interp_Poly,  &
                      Interp_Rgrd, Interp_Cspl
!
    IMPLICIT NONE
!
    REAL(r8), DIMENSION(:), ALLOCATABLE :: xgrd, ygrd
    REAL(r8), DIMENSION(:), ALLOCATABLE :: xnew, ynew_lagr, ynew_poly
    REAL(r8), DIMENSION(:), ALLOCATABLE ::       ynew_line, ynew_cube
    REAL(r8), DIMENSION(:), ALLOCATABLE ::       ynew_cspl
    REAL(r8) :: tmp1, tmp2
    INTEGER(i4) :: ndat, nnew
    INTEGER(i4) :: istat, ierr
    INTEGER(i4) :: i, j
    CHARACTER(LEN=*), PARAMETER :: HEADER = '(Test_Interp)'
!!!
!
    OPEN(UNIT=10, FILE='./points.dat', STATUS='OLD')
    ndat = 0
    DO
        READ(UNIT=10, FMT='(F3.1)', IOSTAT=istat) tmp1
        IF (istat /= 0) EXIT 
        ndat = ndat + 1
    END DO
    CLOSE(UNIT=10)
!
    WRITE(6,'(A,I5)') HEADER//': ndat = ', ndat
!
    ALLOCATE(xgrd(ndat))
    ALLOCATE(ygrd(ndat))
!
    OPEN(UNIT=10, FILE='./points.dat', STATUS='OLD')
    DO i = 1, ndat
        READ(UNIT=10, FMT='(F3.1,1X,F6.4)', IOSTAT=istat) xgrd(i), ygrd(i)
        WRITE(6,'(F3.1,1X,F6.4)') xgrd(i), ygrd(i)
        IF (istat /= 0) EXIT
    END DO
    CLOSE(UNIT=10)
!
!   Set a new grid
!
    nnew = 31
!
    ALLOCATE(xnew(nnew))
    ALLOCATE(ynew_lagr(nnew))
    ALLOCATE(ynew_poly(nnew))
    ALLOCATE(ynew_line(nnew))
    ALLOCATE(ynew_cube(nnew))
    ALLOCATE(ynew_cspl(nnew))
!
    DO j = 1, nnew
        xnew(j) = xgrd(1) + (xgrd(ndat)-xgrd(1)) * REAL(j-1) / REAL(nnew-1)
    END DO
!
!   Interpolation
!
    CALL Interp_Lagr(ndat-1, xgrd, ygrd, xnew, ynew_lagr)
    CALL Interp_Poly(ndat-1, xgrd, ygrd, xnew, ynew_poly)
    CALL Interp_Rgrd(xgrd, ygrd, xnew, ynew_line, 1, ierr)
    CALL Interp_Rgrd(xgrd, ygrd, xnew, ynew_cube, 3, ierr)
    CALL Interp_Cspl(ndat, xgrd, ygrd, 0._r8, 0._r8, 2, 2,  &
                     nnew, xnew, ynew_cspl, 0)
!
!   Output
!
    DO j = 1, nnew
        WRITE(6,'(6(F7.4,1X))') xnew(j), ynew_lagr(j), ynew_poly(j),  &
                                         ynew_line(j), ynew_cube(j),  &
                                         ynew_cspl(j)
    END DO
!
END

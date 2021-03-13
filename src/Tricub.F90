!-------------------------------------------------------------------------------
!
!   MODULE TriCub
!
!>  @brief
!>
!>  - Module for tricubic interpolation
!
!>  @authors
!>
!>  - Francois Lekien
!>    lekien@mit.edu
!>    Department of Mechanical Engineering
!>    Massachusetts Institute of Technology
!>
!>  - Chad Coulliette
!>    chad@griffonlabs.com
!>    Griffon Lab's
!>
!>  - Jerry Marsden
!>    marsden@cds.caltech.edu
!>    Control and Dynamical Systems
!>    California Institute of Technology
!
!   Original sources for tricubic-1.0 can be downloaded at
!   http://gyre.cds.caltech.edu/pub/software/tricubic, but
!   this link is broken as of 2016, but some files were found
!   at https://github.com/nbigaouette/libtricubic 
!
!   Reference
!
!   [1] Lekien, F., and J. Marsden, 2005: Tricubic interpolation
!       in three dimensions. Int. J. Numer. Meth. Engng, 63, 455-471.
!
!   [2] Lekien, F., C. Coulliette, and J. Marsden, 2004: Tricubic engine.
!       Technical notes and full matrix. California Institute of Technology
!       (this document can be downloaded from the github repository).
!
!   Fortran 90/95 module and single precision version are refactored by
!
!   In-Sun Song
!   in-sun.song@nasa.gov
!   Global Modeling and Assimilation Office
!   NASA Goddard Space Flight Center
!
!   09 MAR 2021: Name of module procedures are changed for FAST library
!
!-------------------------------------------------------------------------------
!
MODULE Tricub
!
    USE Base, ONLY: i4, r4, r8
!
    IMPLICIT NONE
!
    PRIVATE
!
    REAL(r8), DIMENSION(64,64) :: aa
    LOGICAL :: aa_initialized = .FALSE.
!
    PUBLIC :: Tricub_Init    ! Fortran only
    PUBLIC :: Tricub_Vers
    PUBLIC :: Tricub_PXYZ
    PUBLIC :: Tricub_Coef
    PUBLIC :: Tricub_Eval
    PUBLIC :: Tricub_IJK2N
    PUBLIC :: Tricub_POINT2XYZ
!
    CHARACTER(LEN=*), PARAMETER :: Tricub_VERSION_STORED = '0.2'
!
    INTERFACE Tricub_PXYZ
        MODULE PROCEDURE Tricub_PXYZ_i4
        MODULE PROCEDURE Tricub_PXYZ_r4
        MODULE PROCEDURE Tricub_PXYZ_r8
    END INTERFACE
!
    INTERFACE Tricub_Coef_Stacked
        MODULE PROCEDURE Tricub_Coef_Stacked_r4
        MODULE PROCEDURE Tricub_Coef_Stacked_r8
    END INTERFACE
!
    INTERFACE Tricub_Coef
        MODULE PROCEDURE Tricub_Coef_r4
        MODULE PROCEDURE Tricub_Coef_r8
    END INTERFACE
!
    INTERFACE Tricub_Eval
        MODULE PROCEDURE Tricub_Eval_Smpl_r4
        MODULE PROCEDURE Tricub_Eval_Smpl_r8
        MODULE PROCEDURE Tricub_Eval_Full_r4
        MODULE PROCEDURE Tricub_Eval_Full_r8
    END INTERFACE
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Tricub_Init
!
!   Initializes A matrix
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_Init
!
        IMPLICIT NONE
!
        INTEGER(i4), DIMENSION(64,64) :: ia
!
        ia(:, 1) = (/  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 2) = (/  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 3) = (/ -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 4) = (/  2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 5) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 6) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 7) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 8) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:, 9) = (/ -3,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,10) = (/  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,11) = (/  9, -9, -9,  9,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  &
                       6, -6,  3, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       4,  2,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,12) = (/ -6,  6,  6, -6,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0,  &
                      -4,  4, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2, -2, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,13) = (/  2,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,14) = (/  0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,15) = (/ -6,  6,  6, -6,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0,  &
                      -3,  3, -3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2, -1, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,16) = (/  4, -4, -4,  4,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  &
                       2, -2,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,17) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,18) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,19) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,20) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,21) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,22) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,23) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  3,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0/)
        ia(:,24) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2, -2,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0/)
        ia(:,25) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,26) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  3,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0,  0,  0,  0,  0/)
        ia(:,27) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  9, -9, -9,  9,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  6,  3, -6, -3,  0,  0,  0,  0,  &
                       6, -6,  3, -3,  0,  0,  0,  0,  4,  2,  2,  1,  0,  0,  0,  0/)
        ia(:,28) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -3, -3,  3,  3,  0,  0,  0,  0,  &
                      -4,  4, -2,  2,  0,  0,  0,  0, -2, -2, -1, -1,  0,  0,  0,  0/)
        ia(:,29) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,30) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0/)
        ia(:,31) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -6,  6,  6, -6,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -4, -2,  4,  2,  0,  0,  0,  0,  &
                      -3,  3, -3,  3,  0,  0,  0,  0, -2, -1, -2, -1,  0,  0,  0,  0/)
        ia(:,32) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  4, -4, -4,  4,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  2,  2, -2, -2,  0,  0,  0,  0,  &
                       2, -2,  2, -2,  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0/)
        ia(:,33) = (/ -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,34) = (/  0,  0,  0,  0,  0,  0,  0,  0, -3,  0,  0,  0,  3,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,35) = (/  9, -9,  0,  0, -9,  9,  0,  0,  6,  3,  0,  0, -6, -3,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  3, -3,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,36) = (/ -6,  6,  0,  0,  6, -6,  0,  0, -3, -3,  0,  0,  3,  3,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -4,  4,  0,  0, -2,  2,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,37) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,38) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0,  0,  0, -1,  0,  0,  0/)
        ia(:,39) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       9, -9,  0,  0, -9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       6,  3,  0,  0, -6, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       6, -6,  0,  0,  3, -3,  0,  0,  4,  2,  0,  0,  2,  1,  0,  0/)
        ia(:,40) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3, -3,  0,  0,  3,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -4,  4,  0,  0, -2,  2,  0,  0, -2, -2,  0,  0, -1, -1,  0,  0/)
        ia(:,41) = (/  9,  0, -9,  0, -9,  0,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       4,  0,  2,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,42) = (/  0,  0,  0,  0,  0,  0,  0,  0,  9,  0, -9,  0, -9,  0,  9,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       6,  0,  3,  0, -6,  0, -3,  0,  6,  0, -6,  0,  3,  0, -3,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  2,  0,  2,  0,  1,  0/)
        ia(:,43) = (/-27, 27, 27,-27, 27,-27,-27, 27,-18, -9, 18,  9, 18,  9,-18, -9,  &
                     -18, 18, -9,  9, 18,-18,  9, -9,-18, 18, 18,-18, -9,  9,  9, -9,  &
                     -12, -6, -6, -3, 12,  6,  6,  3,-12, -6, 12,  6, -6, -3,  6,  3,  &
                     -12, 12, -6,  6, -6,  6, -3,  3, -8, -4, -4, -2, -4, -2, -2, -1/)
        ia(:,44) = (/ 18,-18,-18, 18,-18, 18, 18,-18,  9,  9, -9, -9, -9, -9,  9,  9,  &
                      12,-12,  6, -6,-12, 12, -6,  6, 12,-12,-12, 12,  6, -6, -6,  6,  &
                       6,  6,  3,  3, -6, -6, -3, -3,  6,  6, -6, -6,  3,  3, -3, -3,  &
                       8, -8,  4, -4,  4, -4,  2, -2,  4,  4,  2,  2,  2,  2,  1,  1/)
        ia(:,45) = (/ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0, -2,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,46) = (/  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  0, -3,  0,  3,  0,  3,  0, -4,  0,  4,  0, -2,  0,  2,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -2,  0, -1,  0, -1,  0/)
        ia(:,47) = (/ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6,  &
                       9, -9,  9, -9, -9,  9, -9,  9, 12,-12,-12, 12,  6, -6, -6,  6,  &
                       6,  3,  6,  3, -6, -3, -6, -3,  8,  4, -8, -4,  4,  2, -4, -2,  &
                       6, -6,  6, -6,  3, -3,  3, -3,  4,  2,  4,  2,  2,  1,  2,  1/)
        ia(:,48) = (/-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6,  &
                      -6,  6, -6,  6,  6, -6,  6, -6, -8,  8,  8, -8, -4,  4,  4, -4,  &
                      -3, -3, -3, -3,  3,  3,  3,  3, -4, -4,  4,  4, -2, -2,  2,  2,  &
                      -4,  4, -4,  4, -2,  2, -2,  2, -2, -2, -2, -2, -1, -1, -1, -1/)
        ia(:,49) = (/  2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,50) = (/  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0, -2,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,51) = (/ -6,  6,  0,  0,  6, -6,  0,  0, -4, -2,  0,  0,  4,  2,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -3,  3,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,52) = (/  4, -4,  0,  0, -4,  4,  0,  0,  2,  2,  0,  0, -2, -2,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  2, -2,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,53) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,54) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2,  0,  0,  0, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,  0,  0,  0/)
        ia(:,55) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -6,  6,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -4, -2,  0,  0,  4,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -3,  3,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0, -2, -1,  0,  0/)
        ia(:,56) = (/  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       4, -4,  0,  0, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2,  2,  0,  0, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2, -2,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0/)
        ia(:,57) = (/ -6,  0,  6,  0,  6,  0, -6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -2,  0, -1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,58) = (/  0,  0,  0,  0,  0,  0,  0,  0, -6,  0,  6,  0,  6,  0, -6,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                      -4,  0, -2,  0,  4,  0,  2,  0, -3,  0,  3,  0, -3,  0,  3,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0, -2,  0, -1,  0, -2,  0, -1,  0/)
        ia(:,59) = (/ 18,-18,-18, 18,-18, 18, 18,-18, 12,  6,-12, -6,-12, -6, 12,  6,  &
                      12,-12,  6, -6,-12, 12, -6,  6,  9, -9, -9,  9,  9, -9, -9,  9,  &
                       8,  4,  4,  2, -8, -4, -4, -2,  6,  3, -6, -3,  6,  3, -6, -3,  &
                       6, -6,  3, -3,  6, -6,  3, -3,  4,  2,  2,  1,  4,  2,  2,  1/)
        ia(:,60) = (/-12, 12, 12,-12, 12,-12,-12, 12, -6, -6,  6,  6,  6,  6, -6, -6,  &
                      -8,  8, -4,  4,  8, -8,  4, -4, -6,  6,  6, -6, -6,  6,  6, -6,  &
                      -4, -4, -2, -2,  4,  4,  2,  2, -3, -3,  3,  3, -3, -3,  3,  3,  &
                      -4,  4, -2,  2, -4,  4, -2,  2, -2, -2, -1, -1, -2, -2, -1, -1/)
        ia(:,61) = (/  4,  0, -4,  0, -4,  0,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0/)
        ia(:,62) = (/  0,  0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0, -4,  0,  4,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  &
                       2,  0,  2,  0, -2,  0, -2,  0,  2,  0, -2,  0,  2,  0, -2,  0,  &
                       0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0/)
        ia(:,63) = (/-12, 12, 12,-12, 12,-12,-12, 12, -8, -4,  8,  4,  8,  4, -8, -4,  &
                      -6,  6, -6,  6,  6, -6,  6, -6, -6,  6,  6, -6, -6,  6,  6, -6,  &
                      -4, -2, -4, -2,  4,  2,  4,  2, -4, -2,  4,  2, -4, -2,  4,  2,  &
                      -3,  3, -3,  3, -3,  3, -3,  3, -2, -1, -2, -1, -2, -1, -2, -1/)
        ia(:,64) = (/  8, -8, -8,  8, -8,  8,  8, -8,  4,  4, -4, -4, -4, -4,  4,  4,  &
                       4, -4,  4, -4, -4,  4, -4,  4,  4, -4, -4,  4,  4, -4, -4,  4,  &
                       2,  2,  2,  2, -2, -2, -2, -2,  2,  2, -2, -2,  2,  2, -2, -2,  &
                       2, -2,  2, -2,  2, -2,  2, -2,  1,  1,  1,  1,  1,  1,  1,  1/)
!
        aa = DBLE(ia)
        aa_initialized = .TRUE.
!
        RETURN
    END SUBROUTINE Tricub_Init
!
!-------------------------------------------------------------------------------
!
!   Tricub_Vers
!
!   Fortran implementation of the function
!   char *tricubic_version(void) in libtricubic.cpp
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_Vers() RESULT(ver)
!
        IMPLICIT NONE
        CHARACTER(LEN=2048) :: ver
!!!
!
        ver = Tricub_VERSION_STORED
!
        RETURN
    END FUNCTION Tricub_Vers
!
!-------------------------------------------------------------------------------
!
!   Tricub_PXYZ_i4
!
!   Fortran implementation of the function
!   void tricubic_pointID2xyz(int id, int *x, int *y, int *z)
!   in libtricubic.cpp
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_PXYZ_i4(id, x, y, z)
!
        IMPLICIT NONE
!
        INTEGER(i4), INTENT(IN)    :: id
        INTEGER(i4), INTENT(INOUT) :: x, y, z
!!!
!
        CALL Tricub_POINT2XYZ(id, x, y, z)
!
        RETURN
    END SUBROUTINE Tricub_PXYZ_i4
!
!-------------------------------------------------------------------------------
!
!   Tricub_PXYZ_r8
!
!   Fortran implementation of the function
!   void tricubic_pointID2xyz(int id, double *x, double *y, double *z)
!   in libtricubic.cpp
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_PXYZ_r8(id, x, y, z)
!
        IMPLICIT NONE
!
        INTEGER(i4), INTENT(IN)    :: id
        REAL(r8),    INTENT(INOUT) :: x, y, z
!
        INTEGER(i4) :: x2, y2, z2
!!!
!
        CALL Tricub_POINT2XYZ(id, x2, y2, z2)
!
        x = DBLE(x2)
        y = DBLE(y2)
        z = DBLE(z2)
!
        RETURN
    END SUBROUTINE Tricub_PXYZ_r8
!
!-------------------------------------------------------------------------------
!
!   Tricub_PXYZ_r4
!
!   Single precision version of tricubic_pointID2xyz
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_PXYZ_r4(id, x, y, z)
!
        IMPLICIT NONE
!
        INTEGER(i4), INTENT(IN)    :: id
        REAL(r4),    INTENT(INOUT) :: x, y, z
!
        INTEGER(i4) :: x2, y2, z2
!!!
!
        CALL Tricub_POINT2XYZ(id, x2, y2, z2)
!
        x = FLOAT(x2)
        y = FLOAT(y2)
        z = FLOAT(z2)
!
        RETURN
    END SUBROUTINE Tricub_PXYZ_r4
!
!-------------------------------------------------------------------------------
!
!    Tricub_Coef_Stacked_r8
!
!    Fortran implementation of the function
!    void tricubic_get_coeff_stacked(double a[64], double x[64]) in tricubic.cpp
!
!-------------------------------------------------------------------------------
!
     SUBROUTINE Tricub_Coef_Stacked_r8(a, x)
!
        IMPLICIT NONE
!
        REAL(r8), DIMENSION(64), INTENT(OUT) :: a
        REAL(r8), DIMENSION(64), INTENT(IN)  :: x
!
        INTEGER(i4) :: i, j
!!!
!
        DO i = 1, 64
            a(i) = 0.0_r8
            DO j = 1, 64
                a(i) = a(i) + aa(j,i) * x(j)
            END DO
        END DO
!
        RETURN
    END SUBROUTINE Tricub_Coef_Stacked_r8
!
!-------------------------------------------------------------------------------
!
!   Tricub_Coef_Stacked_r4
!
!   Single precision version of tricubic_get_coeff_stacked
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_Coef_Stacked_r4(a, x)
!
        IMPLICIT NONE
!
        REAL(r4), DIMENSION(64), INTENT(OUT) :: a
        REAL(r4), DIMENSION(64), INTENT(IN)  :: x
!
        INTEGER(i4) :: i, j
!!!
!
        DO i = 1, 64
            a(i) = 0.0_r4
            DO j = 1, 64
                a(i) = a(i) + aa(j,i) * x(j)
            END DO
        END DO
!
        RETURN
    END SUBROUTINE Tricub_Coef_Stacked_r4
!
!-------------------------------------------------------------------------------
!
!   Tricub_Coef_r8
!
!   Fortran implementation of the function
!   void tricubic_get_coeff(double a[64], double f[8],
!   double dfdx[8], double dfdy[8], double dfdz[8], 
!   double d2fdxdy[8], double d2fdxdz[8], double d2fdydz[8], double d3fdxdydz[8])
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Tricub_Coef_r8(a, f, dfdx, dfdy, dfdz,  &
!                             d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz)
    SUBROUTINE Tricub_Coef_r8(a, f, deriv)
!
        IMPLICIT NONE
!
        REAL(r8), DIMENSION(64),  INTENT(OUT) :: a
        REAL(r8), DIMENSION(8),   INTENT(IN)  :: f
!       REAL(r8), DIMENSION(8),   INTENT(IN)  :: dfdx, dfdy, dfdz
!       REAL(r8), DIMENSION(8),   INTENT(IN)  :: d2fdxdy, d2fdxdz, d2fdydz
!       REAL(r8), DIMENSION(8),   INTENT(IN)  :: d3fdxdydz
        REAL(r8), DIMENSION(8,7), INTENT(IN)  :: deriv
!
        REAL(r8), DIMENSION(64) :: x
        INTEGER(i4) :: i
!!!
!
        DO i = 1, 8
            x(   i) = f(i)
            x( 8+i) = deriv(i,1) !dfdx(i)
            x(16+i) = deriv(i,2) !dfdy(i)
            x(24+i) = deriv(i,3) !dfdz(i)
            x(32+i) = deriv(i,4) !d2fdxdy(i)
            x(40+i) = deriv(i,5) !d2fdxdz(i)
            x(48+i) = deriv(i,6) !d2fdydz(i)
            x(56+i) = deriv(i,7) !d3fdxdydz(i)
        END DO
!
        CALL Tricub_Coef_Stacked_r8(a, x)
!
        RETURN
    END SUBROUTINE Tricub_Coef_r8
!
!-------------------------------------------------------------------------------
!
!   Tricub_Coef_r4
!
!   Single precision version of tricubic_get_coeff_dbl
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Tricub_Coef_r4(a, f, dfdx, dfdy, dfdz,  &
!                             d2fdxdy, d2fdxdz, d2fdydz, d3fdxdydz)
    SUBROUTINE Tricub_Coef_r4(a, f, deriv)
!
        IMPLICIT NONE
!
        REAL(r4), DIMENSION(64), INTENT(OUT) :: a
        REAL(r4), DIMENSION(8),  INTENT(IN)  :: f
!       REAL(r4), DIMENSION(8),   INTENT(IN)  :: dfdx, dfdy, dfdz
!       REAL(r4), DIMENSION(8),   INTENT(IN)  :: d2fdxdy, d2fdxdz, d2fdydz
!       REAL(r4), DIMENSION(8),   INTENT(IN)  :: d3fdxdydz
        REAL(r4), DIMENSION(8,7), INTENT(IN)  :: deriv
!
        REAL(r4), DIMENSION(64) :: x
        INTEGER(i4) :: i
!!!
!
        DO i = 1, 8
            x(   i) = f(i)
            x( 8+i) = deriv(i,1) !dfdx(i)
            x(16+i) = deriv(i,2) !dfdy(i)
            x(24+i) = deriv(i,3) !dfdz(i)
            x(32+i) = deriv(i,4) !d2fdxdy(i)
            x(40+i) = deriv(i,5) !d2fdxdz(i)
            x(48+i) = deriv(i,6) !d2fdydz(i)
            x(56+i) = deriv(i,7) !d3fdxdydz(i)
        END DO
!
        CALL Tricub_Coef_Stacked_r4(a, x)
!
        RETURN
    END SUBROUTINE Tricub_Coef_r4
!
!-------------------------------------------------------------------------------
!
!   Tricub_Eval_Smpl_r8
!
!   Fortran implementation of the function
!   double tricubic_eval(double a[64], double x, double y, double z)
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_Eval_Smpl_r8(a, x, y, z) RESULT(ret)
!
!   This is the short version of tricubic_eval. It is used to compute
!   the value of the function at a given point (x,y,z). To compute partial
!   derivatives of f, use the full version with the extra args.

        IMPLICIT NONE
!
        REAL(r8), DIMENSION(64), INTENT(IN) :: a
        REAL(r8),                INTENT(IN) :: x, y, z
!
        REAL(r8) :: ret
        INTEGER(i4) :: i, j, k
!!!
!
        ret = 0.0_r8
        DO i = 0, 3
            DO j = 0, 3
                DO k = 0, 3
                    ret = ret + a(Tricub_IJK2N(i,j,k)+1)*(x**i)*(y**j)*(z**k)
                END DO
            END DO
        END DO
!
        RETURN
    END FUNCTION Tricub_Eval_Smpl_r8
!
!-------------------------------------------------------------------------------
!
!   Tricub_Eval_Smpl_r4
!
!   Single precision version of tricubic_eval_smpl_dbl
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_Eval_Smpl_r4(a, x, y, z) RESULT(ret)
!
!   This is the short version of tricubic_eval. It is used to compute
!   the value of the function at a given point (x,y,z). To compute partial
!   derivatives of f, use the full version with the extra args.

        IMPLICIT NONE
!
        REAL(r4), DIMENSION(64), INTENT(IN) :: a
        REAL(r4),                INTENT(IN) :: x, y, z
!
        REAL(r4) :: ret
        INTEGER(i4) :: i, j, k
!!!
!
        ret = 0.0_r4
        DO i = 0, 3
            DO j = 0, 3
                DO k = 0, 3
                    ret = ret + a(Tricub_IJK2N(i,j,k)+1)*(x**i)*(y**j)*(z**k)
                END DO
            END DO
        END DO

        RETURN
    END FUNCTION Tricub_Eval_Smpl_r4
!
!-------------------------------------------------------------------------------
!
!   Tricub_Eval_Full_r8
!
!   Fortran implementation of the function
!   double tricubic_eval(double a[64], double x, double y, double z, 
!   int derx, int dery, int derz)
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_Eval_Full_r8(a, x, y, z, derx, dery, derz) RESULT(ret)
!
!   The full version takes 3 extra integers args that allows to evaluate any
!   partial derivative of f at the point
!
!   derx = dery = derz = 0    => f
!   derx = 2, dery = derz = 0 => d2f/dx2
!   derx = dery = derz = 1    => d3f/dxdydz
!
!   NOTICE that (derx>3)||(dery>3)||(derz>3) => returns 0.0
!   This computes \frac{\partial ^{derx+dery+derz} d}
!   {\partial x^{derx} \partial y^{dery} \partial z^{derz}}
!
        IMPLICIT NONE
!
        REAL(r8),    DIMENSION(64), INTENT(IN) :: a
        REAL(r8),                   INTENT(IN) :: x, y, z
        INTEGER(i4),                INTENT(IN) :: derx, dery, derz
!
        REAL(r8) :: ret, cont
        INTEGER(i4) :: i, j, k, w
!!!
!
        ret = 0.0_r8
        DO i = derx, 3
            DO j = dery, 3
                DO k = derz, 3
                    cont = a(Tricub_IJK2N(i,j,k)+1)*  &
                           (x**(i-derx))*(y**(j-dery))*(z**(k-derz))
                    DO w = 0, derx-1
                        cont = cont * (i - w)
                    END DO
                    DO w = 0, dery-1
                        cont = cont * (j - w)
                    END DO
                    DO w = 0, derz-1
                        cont = cont * (k - w)
                    END DO
                    ret = ret + cont
                END DO
            END DO
        END DO
!
        RETURN
    END FUNCTION Tricub_Eval_Full_r8
!
!-------------------------------------------------------------------------------
!
!   Tricub_Eval_Full_r4
!
!   Single precision version of tricubic_eval_full_dbl
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_Eval_Full_r4(a, x, y, z, derx, dery, derz) RESULT(ret)
!
!   The full version takes 3 extra integers args that allows to evaluate any
!   partial derivative of f at the point
!
!   derx = dery = derz = 0    => f
!   derx = 2, dery = derz = 0 => d2f/dx2
!   derx = dery = derz = 1    => d3f/dxdydz
!
!   NOTICE that (derx>3)||(dery>3)||(derz>3) => returns 0.0
!   This computes \frac{\partial ^{derx+dery+derz} d}
!   {\partial x^{derx} \partial y^{dery} \partial z^{derz}}
!
        IMPLICIT NONE
!
        REAL(r4),    DIMENSION(64), INTENT(IN) :: a
        REAL(r4),                   INTENT(IN) :: x, y, z
        INTEGER(i4),                INTENT(IN) :: derx, dery, derz
!
        REAL(r4) :: ret, cont
        INTEGER(i4) :: i, j, k, w
!!!
!
        ret = 0.0_r4
        DO i = derx, 3
            DO j = dery, 3
                DO k = derz, 3
                    cont = a(Tricub_IJK2N(i,j,k)+1)*  &
                           (x**(i-derx))*(y**(j-dery))*(z**(k-derz))
                    DO w = 0, derx-1
                        cont = cont * (i - w)
                    END DO
                    DO w = 0, dery-1
                        cont = cont * (j - w)
                    END DO
                    DO w = 0, derz-1
                        cont = cont * (k - w)
                    END DO
                    ret = ret + cont
                END DO
            END DO
        END DO
!
        RETURN
    END FUNCTION Tricub_Eval_Full_r4
!
!-------------------------------------------------------------------------------
!  
!   FUNCTION Tricub_IJK2N
!
!   Fortran implementation of the function
!   int ijk2n(int i, int j, int k) in ltricubic_utils.cpp
!
!-------------------------------------------------------------------------------
!
    FUNCTION Tricub_IJK2N(i, j, k) RESULT(ret)
!
        IMPLICIT NONE
! 
        INTEGER(i4), INTENT(IN) :: i, j, k
        INTEGER(i4) :: ret
!!!
!
        ret = i + 4*j + 16*k
!
        RETURN
    END FUNCTION Tricub_IJK2N
!
!-------------------------------------------------------------------------------
!
!   SUBROUTINE Tricub_POINT2XYZ
!
!   Fortran implementation of the function
!   void point2xyz(int p, int *x, int *y, int *z) in ltricubic_utils.cpp
!
!-------------------------------------------------------------------------------
!
    SUBROUTINE Tricub_POINT2XYZ(p, x, y, z)
!
        IMPLICIT NONE
!
        INTEGER(i4), INTENT(IN)    :: p
        INTEGER(i4), INTENT(INOUT) :: x, y, z
!!!
!
        SELECT CASE(p)
            CASE (0); x = 0; y = 0; z = 0
            CASE (1); x = 1; y = 0; z = 0
            CASE (2); x = 0; y = 1; z = 0
            CASE (3); x = 1; y = 1; z = 0
            CASE (4); x = 0; y = 0; z = 1
            CASE (5); x = 1; y = 0; z = 1
            CASE (6); x = 0; y = 1; z = 1
            CASE (7); x = 1; y = 1; z = 1
            CASE DEFAULT; x = 0; y = 0; z = 0
        END SELECT

        RETURN
    END SUBROUTINE Tricub_POINT2XYZ
!
END MODULE Tricub

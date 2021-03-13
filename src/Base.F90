!-------------------------------------------------------------------------------
!
!  MODULE Base
!
!> @brief
!> - Declare type parameters
!
!  REVISION HISTORY
!
!  16AUG2015 : In-Sun Song : First crack
!  19OCT2019 : In-Sun Song : Reviewed and comments added
!  22FEB2021 : In-Sun Song : Fortran 2008 standard is applied
!
!-------------------------------------------------------------------------------
!
MODULE Base
!
   USE, INTRINSIC ::  &
       ISO_FORTRAN_ENV, ONLY: I4 => int32, I8 => int64,  &
                              R4 => real32, R8 => real64,  & 
                              R16 => real128
!
   IMPLICIT NONE
!
   PRIVATE
!
   PUBLIC :: I4, I8, R4, R8, R16
!
!  INTEGER, PARAMETER, PUBLIC :: R8 = SELECTED_REAL_KIND(12)
!  INTEGER, PARAMETER, PUBLIC :: R4 = SELECTED_REAL_KIND( 6)
!  INTEGER, PARAMETER, PUBLIC :: RN = KIND(1.0)   ! native real
!  INTEGER, PARAMETER, PUBLIC :: I8 = SELECTED_INT_KIND(13)
!  INTEGER, PARAMETER, PUBLIC :: I4 = SELECTED_INT_KIND( 6)
!  INTEGER, PARAMETER, PUBLIC :: IN = KIND(1)     ! native integer
!
!  16-byte real is necessary for computation of Gaussian quadrature
!  in the spherical harmonic dynamical core.
!
!  INTEGER, PARAMETER, PUBLIC :: R16 = 8
!  INTEGER, PARAMETER, PUBLIC :: R16 = SELECTED_REAL_KIND(12)
!  INTEGER, PARAMETER, PUBLIC :: R16 = SELECTED_REAL_KIND(24)
!
END MODULE Base

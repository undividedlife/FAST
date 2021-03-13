!-------------------------------------------------------------------------------
!  MODULE PhysConst
!-------------------------------------------------------------------------------
!
MODULE PhysConst
!
   USE Base,     ONLY: i4, r8, r16
!!!
!!! DO NOT USE PhysGrid for any reasons.
!!! Using PhysGrid leads to circular dependency between DynGrid and PhyGrid
!!!
!
   IMPLICIT NONE
!
   PRIVATE
!
   SAVE
!
   REAL(r8), PARAMETER, PUBLIC :: PI          = 3.14159265358979323846_r8 
!
   REAL(r8), PARAMETER, PUBLIC :: AVOGAD      = 6.02214e26_r8 
   REAL(r8), PARAMETER, PUBLIC :: BOLTZ       = 1.38065e-23_r8 
   REAL(r8), PARAMETER, PUBLIC :: CDAY        = 86400.0_r8 
   REAL(r8), PARAMETER, PUBLIC :: CPAIR       = 1.00464e3_r8 
   REAL(r8), PARAMETER, PUBLIC :: CPLIQ       = 4.188e3_r8 
   REAL(r8), PARAMETER, PUBLIC :: KARMAN      = 0.4_r8 
   REAL(r8), PARAMETER, PUBLIC :: LATICE      = 3.337e5_r8 
   REAL(r8), PARAMETER, PUBLIC :: LATVAP      = 2.501e6_r8
   REAL(r8), PARAMETER, PUBLIC :: PSTD        = 101325.0_r8 
   REAL(r8), PARAMETER, PUBLIC :: R_UNIVERSAL = AVOGAD*BOLTZ
   REAL(r8), PARAMETER, PUBLIC :: RHOH2O      = 1.000e3_r8 
   REAL(r8), PARAMETER, PUBLIC :: SPVAL       = 1.e30_r8  ! shr_const_spval 
   REAL(r8), PARAMETER, PUBLIC :: STEBOL      = 5.67e-8_r8 
   REAL(r8), PARAMETER, PUBLIC :: H2OTRIP     = 273.16_r8 
!
   REAL(r8), PARAMETER, PUBLIC :: C0          = 2.99792458e8_r8
   REAL(r8), PARAMETER, PUBLIC :: PLANCK      = 6.6260755e-34_r8
   REAL(r8), PARAMETER, PUBLIC :: GRAVIT      = 9.80616_r8 
   REAL(r8), PARAMETER, PUBLIC :: SDAY        = 86164.0_r8 
   REAL(r8), PARAMETER, PUBLIC :: MWH2O       = 18.016_r8 
   REAL(r8), PARAMETER, PUBLIC :: CPWV        = 1.810e3_r8 
   REAL(r8), PARAMETER, PUBLIC :: MWDRY       = 28.966_r8
   REAL(r8), PARAMETER, PUBLIC :: REARTH      = 6.37122e6_r8
   REAL(r8), PARAMETER, PUBLIC :: TMELT       = 273.15_r8
!
   REAL(r8), PARAMETER, PUBLIC :: RGA         = 1._r8/GRAVIT
   REAL(r8), PARAMETER, PUBLIC :: RA          = 1._r8/REARTH
   REAL(r8), PARAMETER, PUBLIC :: OMEGA       = 2.0_r8*PI/SDAY 
   REAL(r8), PARAMETER, PUBLIC :: RH2O        = R_UNIVERSAL/MWH2O
   REAL(r8), PARAMETER, PUBLIC :: RAIR        = R_UNIVERSAL/MWDRY
   REAL(r8), PARAMETER, PUBLIC :: EPSILO      = MWH2O/MWDRY
   REAL(r8), PARAMETER, PUBLIC :: ZVIR        = RH2O/RAIR - 1.0_r8 
   REAL(r8), PARAMETER, PUBLIC :: CPVIR       = CPWV/CPAIR - 1.0_r8 
   REAL(r8), PARAMETER, PUBLIC :: RHODAIR     = PSTD/(RAIR*TMELT) 
   REAL(r8),            PUBLIC :: ez   ! Coriolis expansion coeff
   REAL(r8), PARAMETER, PUBLIC :: CAPPA       = (R_UNIVERSAL/MWDRY)/CPAIR
   REAL(r8), PARAMETER, PUBLIC :: CPD_ON_CPV  = CPAIR/CPWV
!
!  Physical constants for dynamical core
!
   REAL(r8),  PARAMETER, PUBLIC :: DD_PI         = 3.141592653589793238462643383279_r8
   REAL(r16), PARAMETER, PUBLIC :: QQ_PI         = 3.141592653589793238462643383279_r16
   REAL(r8),  PARAMETER, PUBLIC :: G             = GRAVIT
   REAL(r8),  PARAMETER, PUBLIC :: RGAS          = RAIR
   REAL(r8),  PARAMETER, PUBLIC :: CP            = CPAIR
   REAL(r8),  PARAMETER, PUBLIC :: P0            = PSTD/100.0_r8
   REAL(r8),  PARAMETER, PUBLIC :: MWDAIR        = MWDRY
   REAL(r8),  PARAMETER, PUBLIC :: RWATER_VAPOR  = RH2O
   REAL(r8),  PARAMETER, PUBLIC :: CPWATER_VAPOR = CPWV
   REAL(r8),  PARAMETER, PUBLIC :: KAPPA         = CAPPA
   REAL(r8),  PARAMETER, PUBLIC :: RD_ON_RV      = EPSILO
   REAL(r8),  PARAMETER, PUBLIC :: RREARTH       = RA
!
   PUBLIC :: PhysConstInit
!
CONTAINS
!
!-------------------------------------------------------------------------------
!  SUBROUTINE PhysConstInit
!-------------------------------------------------------------------------------
!
   SUBROUTINE PhysConstInit
!
      IMPLICIT NONE
!
      INTEGER(i4) :: ierr
!!!
!   
      ez = OMEGA / SQRT(0.375_r8)
!
      RETURN
   END SUBROUTINE PhysConstInit
!
END MODULE

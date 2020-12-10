!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : UTILITIES
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief A module for general subroutines that will be used in several places
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE UTILS

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
  
  ! Define the subroutines here as public so that they may be called from
  ! other places.
  PUBLIC :: NETCDF_ERR_CHECK
  PUBLIC :: RAND_INIT
  PUBLIC :: RAND_NORMAL

  ! Define interface blocks for subrotuines in other files

  !> Error checking for netCDF subroutines
  INTERFACE NETCDF_ERR_CHECK

     SUBROUTINE NETCDF_ERR_CHECK_SBR(status)
       USE NETCDF
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       INTEGER(qb), INTENT(IN) :: status
     END SUBROUTINE NETCDF_ERR_CHECK_SBR

  END INTERFACE NETCDF_ERR_CHECK

  !> Initializes RANDOM_NUMBER with a given seed
  INTERFACE RAND_INIT

     SUBROUTINE RAND_INIT_SBR(seed)
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       INTEGER(qb), INTENT(IN) :: seed
     END SUBROUTINE RAND_INIT_SBR

  END INTERFACE RAND_INIT

  !> Generates a standard Gaussian random variable
  INTERFACE RAND_NORMAL
     
     SUBROUTINE RAND_NORMAL_SP(randNum) !< Single-precision real
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(sp), INTENT(INOUT) :: randNum
     END SUBROUTINE RAND_NORMAL_SP

     SUBROUTINE RAND_NORMAL_DP(randNum) !< Double-precision real
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(dp), INTENT(INOUT) :: randNum
     END SUBROUTINE RAND_NORMAL_DP

     SUBROUTINE RAND_NORMAL_C_SP(randNum) !< Single-precision complex
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       COMPLEX(sp), INTENT(INOUT) :: randNum
     END SUBROUTINE RAND_NORMAL_C_SP

     SUBROUTINE RAND_NORMAL_C_DP(randNum) !< Double-precision complex
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       COMPLEX(dp), INTENT(INOUT) :: randNum
     END SUBROUTINE RAND_NORMAL_C_DP

  END INTERFACE RAND_NORMAL
  
END MODULE UTILS

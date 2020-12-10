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
!> @brief Contains general random-number generation subroutines that will be
!! used in various places
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Initializes RANDOM_NUMBER with the given seed. Call this for repeatability.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE RAND_INIT_SBR(seed)

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  INTEGER(qb), INTENT(IN) :: seed !< The seed to use with the RNG
  INTEGER(qb) :: i !< Counter for DO loops
  INTEGER(qb), DIMENSION(33) :: seedArr !< The array generated for putting
  !! into RANDOM_SEED

  seedArr = 0_qb
  seedArr(1) = 104729_qb + 5_qb * seed ! This is pretty aribtrary

  DO i = 2, 33
     seedArr(i) = MOD(420_qb * seedArr(i-1) + 69_qb, 4294967_qb) ! This is
     !! also pretty arbitrary
  END DO
  
  CALL RANDOM_SEED(put = seedArr)

END SUBROUTINE RAND_INIT_SBR

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Generates a single-precision real number from a standard Gaussian
!! distribution.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  

SUBROUTINE RAND_NORMAL_SP(randNum)

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
  
  REAL(sp), INTENT(INOUT) :: randNum
  REAL(sp) :: r1 !< Real number to be used on Box-Mueller transform
  REAL(sp) :: r2 !< Real-number to be used in Box-Mueller transform
  REAL(sp) :: pi_sp !< Pi

  pi_sp = 4.0_sp * ATAN(1.0_sp)

  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)

  randNum = SQRT(-2.0_sp * LOG(r1)) * COS(2.0_sp * pi_sp * r2)
  
END SUBROUTINE RAND_NORMAL_SP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Generates a double-precision real number from a standard Gaussian
!! distribution.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE RAND_NORMAL_DP(randNum)
  
  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
  
  REAL(dp), INTENT(INOUT) :: randNum
  REAL(dp) :: r1 !< Real number to be used on Box-Mueller transform
  REAL(dp) :: r2 !< Real-number to be used in Box-Mueller transform
  REAL(dp) :: pi_dp !< Pi

  pi_dp = 4.0_dp * ATAN(1.0_dp)

  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)

  randNum = SQRT(-2.0_dp * LOG(r1)) * COS(2.0_dp * pi_dp * r2)
  
END SUBROUTINE RAND_NORMAL_DP


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Generates a single-precision complex number from a standard Gaussian
!! distribution, indepndent real and complex parts.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE RAND_NORMAL_C_SP(randNum) !< Single-precision complex
  
  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
  
  COMPLEX(sp), INTENT(INOUT) :: randNum
  REAL(sp) :: r1 !< Real number to be used on Box-Mueller transform
  REAL(sp) :: r2 !< Real-number to be used in Box-Mueller transform
  REAL(sp) :: pi_sp !< Pi

  pi_sp = 4.0_sp * ATAN(1.0_sp)

  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)

  randNum = CMPLX(SQRT(-2.0_sp * LOG(r1)) * COS(2.0_sp * pi_sp * r2), &
       & SQRT(-2.0_sp * LOG(r1)) * SIN(2.0_sp * pi_sp * r2), sp)
  
END SUBROUTINE RAND_NORMAL_C_SP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Generates a double-precision complex number from a standard Gaussian
!! distribution, indepndent real and complex parts.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE RAND_NORMAL_C_DP(randNum) !< Double-precision complex

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
  
  COMPLEX(dp), INTENT(INOUT) :: randNum
  REAL(dp) :: r1 !< Real number to be used on Box-Mueller transform
  REAL(dp) :: r2 !< Real-number to be used in Box-Mueller transform
  REAL(dp) :: pi_dp !< Pi

  pi_dp = 4.0_dp * ATAN(1.0_dp)

  CALL RANDOM_NUMBER(r1)
  CALL RANDOM_NUMBER(r2)

  randNum = CMPLX(SQRT(-2.0_dp * LOG(r1)) * COS(2.0_dp * pi_dp * r2), &
       & SQRT(-2.0_dp * LOG(r1)) * SIN(2.0_dp * pi_dp * r2), dp)
  
END SUBROUTINE RAND_NORMAL_C_DP

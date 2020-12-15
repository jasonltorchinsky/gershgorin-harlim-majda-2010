!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : TIME_STEPPER
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief An adaptive Forward Euler method for solving the Gershgorin-Majda
!! 2010 system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Uses an adaptive Forward Euler scheme to progress that state variable,
!! additive bias, and multiplicative bias to maximize the time-step size during
!! exponential decay periods (when the method is absolutely stable).
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ADAPTIVE_FORWARD_EULER_SBR(currTime, defaultTimeStepSize)

  USE INITIALIZE
  USE DET_FORCING
  USE UTILS
  
  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation.
  REAL(dp), INTENT(IN) :: defaultTimeStepSize !< The original time-step size
  REAL(dp) :: realNoise !< The noise to be added to the multiplicative bias
  COMPLEX(dp) :: cmplxNoise !< The noise to be added to the additive bias
  !! and state variable
  COMPLEX(dp) :: detForce !< The deterministic forcing for the state variable

  ! Determine if we are changing the time-step size
  IF (multBias .GT. 0) THEN ! Change time-step size
     timeStepSize = multBias / ((multBias * multBias + oscFreqU * oscFreqU))
  ELSE IF (multBias .LE. 0) THEN ! Choose the default time-step size
     timeStepSize = defaultTimeStepSize
  END IF
  
  
  ! Update deterministic forcing
  CALL DET_FORCING_DEFAULT(currTime, detForce)
  
  ! Update state variable
  CALL RAND_NORMAL(cmplxNoise) ! Get noise for state variable
  stateVar = ((1.0_dp + (-1.0_dp * multBias + (0.0_dp, 1.0_dp) * oscFreqU) &
       & * timeStepSize) * stateVar) + ((addBias + detForce) * timeStepSize) &
       & + (noiseStrU * SQRT(timeStepSize) * cmplxNoise)
  
  ! Update multiplicative bias
  CALL RAND_NORMAL(realNoise) ! Get the noise for the multiplicative bias
  multBias = ((1.0_dp - dampGamma * timeStepSize) * multBias) &
       & + (dampGamma * meanGamma * timeStepSize) &
       & + (noiseStrGamma * SQRT(timeStepSize) * realNoise)
  
  ! Update additive bias
  CALL RAND_NORMAL(cmplxNoise) ! Get noise for additive bias
  addBias = ((1.0_dp + (-1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB) &
       & * timeStepSize) * addBias) &
       & - ((dampB + (0.0_dp, 1.0_dp) * oscFreqB) * meanB * timeStepSize) &
       & + (noiseStrB * SQRT(timeStepSize) * cmplxNoise)

  

END SUBROUTINE ADAPTIVE_FORWARD_EULER_SBR

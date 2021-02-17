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
!> @brief An adaptive Backward Euler method for solving the Gershgorin-Majda
!! 2010 system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
!> @brief
!> Uses an adaptive Backward Euler scheme to progress that state variable,
!! additive bias, and multiplicative bias to ensure absolute stability during
!! periods of exponential growth.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ADAPTIVE_BACKWARD_EULER_SBR(currTime, defaultTimeStepSize)

  USE INITIALIZE
  USE DET_FORCING
  USE UTILS
  
  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation.
  REAL(dp), INTENT(IN) :: defaultTimeStepSize !< The original time-step size
  REAL(dp) :: realNoise !< The noise to be added to the multiplicative bias
  REAL(dp) :: prevMultBias !< The multiplicative bias of the previous step, in
  !! case we must adjust the step size
  REAL(dp) :: timeStepSizeBound !< Loweer bound for the time-step size
  COMPLEX(dp) :: cmplxNoise !< The noise to be added to the additive bias
  !! and state variable
  COMPLEX(dp) :: detForce !< The deterministic forcing for the state variable
  COMPLEX(dp) :: lambdaB !< The lambda_b parameter in the equations
  COMPLEX(dp) :: lambda !< Holds the value -gamma_{k+1} + i omegaU

  ! Always try the default time-step size first
  timeStepSize = defaultTimeStepSize
  
  ! Update multiplicative bias
  prevMultBias = multBias
  CALL RAND_NORMAL(realNoise) ! Get the noise for the multiplicative bias
  multBias = (1.0_dp / (1.0_dp + dampGamma * timeStepSize)) &
       & * (multBias + (dampGamma * meanGamma * timeStepSize) &
       & + (noiseStrGamma * SQRT(timeStepSize) * realNoise))
  

  ! If the multiplicative bias is negative, we may need to change the time-step
  ! size and redo it
  IF (multBias .LT. 0.0_dp) THEN
     timeStepSizeBound = -2.0_dp * multBias / (multBias * multBias &
          & + oscFreqU * oscFreqU) ! Find the lower bound for the time-step size
     DO WHILE (timeStepSize .LT. timeStepSizeBound) ! If our time-step size is
        ! too small, we need to try to find the multiplicative bias again
        timeStepSize = timeStepSizeBound
        CALL RAND_NORMAL(realNoise) ! Get the noise for the multiplicative bias
        multBias = (1.0_dp / (1.0_dp + dampGamma * timeStepSize)) &
             & * (prevMultBias + (dampGamma * meanGamma * timeStepSize) &
             & + (noiseStrGamma * SQRT(timeStepSize) * realNoise))
        timeStepSizeBound = -2.0_dp * multBias / (multBias * multBias &
             & + oscFreqU * oscFreqU)
     END DO
  END IF
  
  ! Update additive bias
  CALL RAND_NORMAL(cmplxNoise) ! Get noise for additive bias
  lambdaB = -1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB
  addBias = (1.0_dp / (1.0_dp - (lambdaB * timeStepSize))) &
       & * (addBias - (lambdaB * meanB  * timeStepSize) &
       & + (noiseStrB * SQRT(timeStepSize) * cmplxNoise))

  ! Update deterministic forcing
  CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
  
  ! Update state variable
  CALL RAND_NORMAL(cmplxNoise) ! Get noise for state variable
  lambda = -1.0_dp * multBias + (0.0_dp, 1.0_dp) * oscFreqU
  stateVar = (1.0_dp / (1.0_dp - lambda * timeStepSize)) &
       & * (stateVar + ((addBias + detForce) * timeStepSize) &
       & * (noiseStrU * SQRT(timeStepSize) * cmplxNoise))

  

END SUBROUTINE ADAPTIVE_BACKWARD_EULER_SBR

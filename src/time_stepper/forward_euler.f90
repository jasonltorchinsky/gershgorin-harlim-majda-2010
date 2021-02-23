!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-HARLIM-MAJDA_10 Project of the the multi-model communication
! research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE       : Multi-Model Communication
! PROJECT     : GERSHGORIN-HARLIM-MAJDA_10
! MODULE      : TIME_STEPPER
! URL         : https://github.com/jasonltorchinsky/gershgorin-harlim-majda-2010
! AFFILIATION : University of Wisconsin-Madison
! DATE        : Winter 2020
! REVISION    : 1.00
!
!> @author
!> Jason Torchinsky
!
!> @brief The Forward Euler method for solving the Gershgorin-Majda 2010 system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Torchinsky, University of Wisconsin-Madison
!> @brief
!> Uses the Forward Euler scheme to progress that state variable, additive
!! bias, and multiplicative bias.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE FORWARD_EULER_SBR(currTime)

  USE INITIALIZE
  USE DET_FORCING
  USE UTILS
  
  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation.
  REAL(dp) :: realNoise !< The noise to be added to the multiplicative bias
  COMPLEX(dp) :: cmplxNoise !< The noise to be added to the additive bias
  !! and state variable
  COMPLEX(dp) :: detForce !< The deterministic forcing for the state variable

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

  

END SUBROUTINE FORWARD_EULER_SBR

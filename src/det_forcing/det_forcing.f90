!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-HARLIM-MAJDA_10 Project of the the multi-model communication
! research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE       : Multi-Model Communication
! PROJECT     : GERSHGORIN-HARLIM-MAJDA_10
! MODULE      : DET_FORCING
! URL         : https://github.com/jasonltorchinsky/gershgorin-harlim-majda-2010
! AFFILIATION : University of Wisconsin-Madison
! DATE        : Winter 2020
! REVISION    : 1.00
!
!> @author
!> Jason Torchinsky
!
!> @brief The determinstic forcing the Gershgorin-Harlim-Majda 2010 system.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE DET_FORCING

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  !> Declare the subroutines here as public so that they may be called from
  !! elsewhere as necessary.
  PUBLIC :: DET_FORCING_DEFAULT

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> The standard deterministic forcing for the system, complex exponential.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE DET_FORCING_DEFAULT(time, detForce)

    USE INITIALIZE

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: time !< Time to evaluate the determinsitic
    !! forcing at
    COMPLEX(dp), INTENT(INOUT) :: detForce !< The deterministic forcing at the
    !! given time

    detForce = strForce * EXP( (0.0_dp, 1.0_dp) * oscFreqForce * time)


  END SUBROUTINE DET_FORCING_DEFAULT

END MODULE DET_FORCING

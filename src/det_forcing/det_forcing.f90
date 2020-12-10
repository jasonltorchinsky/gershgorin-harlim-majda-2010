!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : DET_FORCING
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The determinstic forcing the Gershgorin-Majda 2010 system.
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
  !> @author Jason Turner, University of Wisconsin-Madison
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

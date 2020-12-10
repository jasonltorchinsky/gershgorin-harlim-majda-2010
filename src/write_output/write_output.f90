!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : OUTPUT
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The output module for the Gershgorin-Majda 2010 system solver,
!! contains variable declarations and subroutines needed for output.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE WRITE_OUTPUT

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  ! Declare the subroutines here as public so that they may be called from the
  ! main driver.
  PUBLIC :: OUTPUT_STATE
  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Writes the current state to the netCDF file
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE OUTPUT_STATE(stepNum, time)

    USE INITIALIZE
    USE NETCDF
    USE UTILS

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: stepNum
    REAL(dp), INTENT(IN) :: time

    ! outputFreq = 0 means no output
    IF (outputFreq .NE. 0) THEN
       CALL NETCDF_ERR_CHECK( NF90_PUT_VAR(outFileID, runID, &
            & (/ REAL(stepNum, dp), time, REAL(stateVar, dp), AIMAG(stateVar), &
            &   REAL(addBias, dp), AIMAG(addBias), multBias /), &
            & (/ 1, stepNum/outputFreq + 1 /), (/ 7, 1 /)) ) ! Note the integer
            ! division
    END IF
    
  END SUBROUTINE OUTPUT_STATE
  
  
END MODULE WRITE_OUTPUT

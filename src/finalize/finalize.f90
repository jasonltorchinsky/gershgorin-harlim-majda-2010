!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : FINALIZE
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The finalization module for the Gershgorin-Majda 2010 system solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE FINALIZE

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer


  ! Declare the subroutines here as public so that they may be called from 
  ! elsewhere.
  PUBLIC :: FINALIZE_OUTPUT_FILE

  
CONTAINS


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @Author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Finalizes the output file for the simulation.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE FINALIZE_OUTPUT_FILE

    USE INITIALIZE
    USE NETCDF
    USE UTILS

    IMPLICIT NONE

    ! outputFreq = 0 means no output
    IF (outputFreq .NE. 0_qb) THEN 
       ! Close the file.
       CALL NETCDF_ERR_CHECK( NF90_CLOSE(outFileID) )
    END IF
    
  END SUBROUTINE FINALIZE_OUTPUT_FILE
  
  
END MODULE FINALIZE

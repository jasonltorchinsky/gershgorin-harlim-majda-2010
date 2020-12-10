!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : MAIN
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The main module for the Gershgorin-Majda 2010 system solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PROGRAM GERSHGORIN_MAJDA_10

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines kinds of real, complex, and integer
  INCLUDE 'ver.h' ! Defines program version information

  REAL(dp) :: startTime
  REAL(dp) :: endTime

  WRITE(*,"(A,I1,A,I2.2)") "Starting Gershgorin-Majda 2010 solver, version ", &
       & GERSHGORIN_MAJDA_10_VERSION_MAJOR, ".", &
       & GERSHGORIN_MAJDA_10_VERSION_MINOR

  CALL CPU_TIME(startTime)
  CALL MAIN
  CALL CPU_TIME(endTime)

  WRITE(*,"(A,A,ES12.6,A)") "Simulation completed successfully, " , &
       "total execution time: ", endTime - startTime, " seconds."

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, Unviersity of Wisconsin-Madison
  !> @brief
  !> The main program driver, calls all necessary stages to run the simulation.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE MAIN

    USE INITIALIZE
    USE WRITE_OUTPUT
    USE TIME_STEPPER
    USE FINALIZE

    IMPLICIT NONE

    INTEGER(qb) :: runNum !< Current simulation number.

    ! Initialize the system parameters, the system variables, the simulation
    ! parameters, and the output file.
    CALL INITIALIZE_SYSTEM_PARAMETERS
    CALL INITIALIZE_SYSTEM_CONDITIONS
    CALL INITIALIZE_SIMULATION_PARAMETERS
    CALL INITIALIZE_OUTPUT_FILE

    ! Run the simulation.
    DO runNum = 1, runCount
       CALL INITIALIZE_SIMULATION(runNum)
       CALL RUN_SIMULATION
    END DO
    
    ! Finalize the simulation by closing the output file.
    CALL FINALIZE_OUTPUT_FILE

    RETURN

  END SUBROUTINE MAIN

END PROGRAM GERSHGORIN_MAJDA_10

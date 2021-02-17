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
    USE UTILS

    USE MPI

    IMPLICIT NONE

    INTEGER(qb) :: runNum !< Current simulation number.
    INTEGER(qb) :: ierror !< Additional flag for MPI commands.

    ! Initialize MPI.
    CALL MPI_INIT(ierror)

    ! Get processor ID for writing messages
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, procID, ierror)
    IF (procID .EQ. 0_qb) THEN
       WRITE(*,"(A,I1,A,I2.2)") "Starting Gershgorin-Majda 2010 solver, version ", &
            & GERSHGORIN_MAJDA_10_VERSION_MAJOR, ".", &
            & GERSHGORIN_MAJDA_10_VERSION_MINOR
    END IF

    ! Initialize the system parameters, the system variables, the simulation
    ! parameters, and the output file.
    CALL INITIALIZE_SYSTEM_PARAMETERS
    CALL INITIALIZE_SIMULATION_PARAMETERS
    CALL INITIALIZE_OUTPUT_FILE
    
    ! Seed RANDOM_NUMBER, we seed it here so that each run will use
    ! different random numbers.
    CALL RAND_INIT(procID*seed)
    ! Run the simulation.
    DO runNum = 1, runCount
       IF ((MOD(runNum, runCount/10_qb) .EQ. 0_qb) &
            & .AND. (procID .EQ. 0_qb)) THEN
          PRINT *, "Starting run ", runNum, "of ", runCount, "."
       END IF
       CALL INITIALIZE_SYSTEM_CONDITIONS
       CALL INITIALIZE_SIMULATION(runNum)
       CALL RUN_SIMULATION
    END DO
    
    ! Finalize the simulation by closing the output file.
    CALL FINALIZE_OUTPUT_FILE
    CALL MPI_FINALIZE(ierror)

    ! Stop all but process 0, which will write final messages.
    IF (procID .NE. 0_qb) THEN
       STOP
    END IF

    RETURN

  END SUBROUTINE MAIN

END PROGRAM GERSHGORIN_MAJDA_10

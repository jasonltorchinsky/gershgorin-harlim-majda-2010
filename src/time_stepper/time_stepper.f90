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
!> @brief The time-stepper module for the Gershgorin-Majda 2010 system solver,
!! contains the DO loop for performing the simualtion and outputting the
!! results.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE TIME_STEPPER

  USE INITIALIZE

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  REAL(dp) :: defaultTimeStepSize !< Store the initial time-step size so
  !! that we may use it in other places in the module.

  !> Declare the subroutines here as public so that they may be called from
  !! elsewhere as necessary.
  PUBLIC :: RUN_SIMULATION
  PUBLIC :: FORWARD_EULER

  ! Declare the interfaces for subroutines outside of this module

  !> The Forward Euler time-stepping scheme.
  INTERFACE FORWARD_EULER
     
     SUBROUTINE FORWARD_EULER_SBR(currTime)
       USE INITIALIZE
       USE DET_FORCING
       USE UTILS
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation
     END SUBROUTINE FORWARD_EULER_SBR
     
  END INTERFACE FORWARD_EULER

  !> The Backward Euler time-stepping scheme.
  INTERFACE BACKWARD_EULER
     
     SUBROUTINE BACKWARD_EULER_SBR(currTime)
       USE INITIALIZE
       USE DET_FORCING
       USE UTILS
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation
     END SUBROUTINE BACKWARD_EULER_SBR
     
  END INTERFACE BACKWARD_EULER

  !> An adaptive Forward Euler time-stepping scheme, based on the absolute
  !! stability criteria.

  INTERFACE ADAPTIVE_FORWARD_EULER

     SUBROUTINE ADAPTIVE_FORWARD_EULER_SBR(currTime, defaultTimeStepSize)
       USE INITIALIZE
       USE DET_FORCING
       USE UTILS
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation
       REAL(dp), INTENT(IN) :: defaultTimeStepSize !< The original time-step
       !! size
     END SUBROUTINE ADAPTIVE_FORWARD_EULER_SBR

  END INTERFACE ADAPTIVE_FORWARD_EULER

  !> An adaptive Backward Euler time-stepping scheme, based on the absolute
  !! stability criteria.

  INTERFACE ADAPTIVE_BACKWARD_EULER

     SUBROUTINE ADAPTIVE_BACKWARD_EULER_SBR(currTime, defaultTimeStepSize)
       USE INITIALIZE
       USE DET_FORCING
       USE UTILS
       INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer
       REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation
       REAL(dp), INTENT(IN) :: defaultTimeStepSize !< The original time-step
       !! size
     END SUBROUTINE ADAPTIVE_BACKWARD_EULER_SBR

  END INTERFACE ADAPTIVE_BACKWARD_EULER

  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Performs the DO loop for progressing the simulation.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE RUN_SIMULATION

    USE INITIALIZE
    USE WRITE_OUTPUT
    USE UTILS

    IMPLICIT NONE

    INTEGER(qb) :: stepNum !< Current step number of the simulation
    REAL(dp) :: time !< Time of the current step

    ! Set initial step number and time
    stepNum = 0_qb
    time = 0.0_dp

    ! Set the default timeStepSize as the one the user input
    defaultTimeStepSize = timeStepSize

    ! If we are writing output, write the initial condition to file
    IF (outputFreq .NE. 0_qb) THEN
       CALL OUTPUT_STATE(stepNum, time)
    END IF
    
    ! Seed RANDOM_NUMBER
    CALL RAND_INIT(seed)

    ! Main time-stepping loop
    DO WHILE (stepNum .LT. timeStepCount)

       SELECT CASE (timeStepScheme) ! Choose the time-step-scheme

       CASE (0) ! Forward Euler
          CALL FORWARD_EULER(time)
          
       CASE(1) ! Backward Euler
          CALL BACKWARD_EULER(time)

       CASE(2) ! Adaptive Forward Euler
          CALL ADAPTIVE_FORWARD_EULER(time, defaultTimeStepSize)

       CASE(3) ! Adaptive Backward Euler
          CALL ADAPTIVE_BACKWARD_EULER(time, defaultTimeStepSize)
          
       CASE DEFAULT ! Choose Forward Euler by default
          ERROR STOP "Invalid time-stepper selected."
          
       END SELECT

       ! Update the step number and time
       stepNum = stepNum + 1
       time = time + timeStepSize

       ! Output the updated state
       IF (outputFreq .NE. 0_qb) THEN
          IF (MOD(stepNum, outputFreq) .EQ. 0) THEN
             CALL OUTPUT_STATE(stepNum, time)
          END IF
       END IF
       
    END DO

    
  END SUBROUTINE RUN_SIMULATION
  
END MODULE TIME_STEPPER

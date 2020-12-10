!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : INITIALIZE
! URL              : N/A, will be posted to GitHub.
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2020
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The initalization module for the Gershgorin-Majda 2010 system solver,
!! contains variable declarations and subroutines needed for initialization.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE INITIALIZE

  IMPLICIT NONE

  PRIVATE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  ! System variables (u, b, gamma)
  COMPLEX(dp), PUBLIC :: stateVar !< State variable in the GS10 system (u)
  COMPLEX(dp), PUBLIC :: addBias !< Additive bias correction in the GM10
  !! system (b)
  REAL(dp), PUBLIC :: multBias !< Multiplicative bais in GM10 system (gamma)

  ! System parameters
  REAL(dp), PUBLIC :: oscFreqU !< Oscillation frequency of the state variable
  !! in the GM10 system (omega)
  REAL(dp), PUBLIC :: noiseStrU !< Strength of white noise for the state
  !! variable (u) in the GM10 system (sigma)
  REAL(dp), PUBLIC :: dampB !< Strength of damping for the additive bias (b)
  !! in the GM10 system (gamma_b)
  REAL(dp), PUBLIC :: oscFreqB !< Oscillation frequency of the additive bias (b)
  !! in the GM10 system (omega_b)
  COMPLEX(dp), PUBLIC :: meanB !< Stationary mean value for the additive
  !! bias (b) in the GM10 system (widehat{b})
  REAL(dp), PUBLIC :: noiseStrB !< Strength of white noise for the additive
  !! bias (b) in the GM10 system (sigma_b)
  REAL(dp), PUBLIC :: dampGamma !< Strength of damping for the multiplicative
  !! bias (gamma) for the GM10 system (d_gamma)
  REAL(dp), PUBLIC :: meanGamma !< Stationary mean value for the multiplicative
  !! bias (gamma) in the GM10 system (widehat{gamma})
  REAL(dp), PUBLIC :: noiseStrGamma !< Strength of white noise for the
  !! multiplicative bias (gamma) in the GM10 system (sigma_gamma)
  REAL(dp), PUBLIC :: strForce !< Strength of the deterministic forcing (f) in
  !! the GM10 system (A_f, assuming the oscillatory determinstic forcing)
  REAL(dp), PUBLIC :: oscFreqForce !< Oscillation frequency of the determinstic
  !! forcing (f) in the GM10 system (omega_f, assuming oscillatory determinstic
  !! frequency)

  ! Simulation parameters
  REAL(dp), PUBLIC :: timeStepSize !< Size of time-steps to take
  INTEGER(qb), PUBLIC :: timeStepCount !< Number of time-steps to take
  INTEGER(qb), PUBLIC :: seed !< Seed to use for RNG
  INTEGER(qb), PUBLIC :: timeStepScheme !< Which time-stepping scheme to use
  !! (0 = Forward Euler)
  INTEGER(qb), PUBLIC :: runCount !< Number of runs to perform
  INTEGER(qb), PUBLIC :: outputFreq !< Frequency of output to file (every n steps)

  ! Output file parameters
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: outFileName = "out.nc" !< Output file name
  INTEGER(qb), PUBLIC :: outFileID !< NetCDF ID for output file
  INTEGER(qb), PUBLIC :: stepCountID !< NetCDF ID for step number dimension
  INTEGER(qb), PUBLIC :: varCountID !< NetCDF ID for the vars dimension (seven
  !! vars: time-step number, time, state variable value (real and complex parts),
  !! additive bias value (real and complex parts), and multiplicative bias value)
  INTEGER(qb), PUBLIC :: runID !< NetCDF ID for the output for the current run

  ! Declare the subroutines here as public so that they may be called from the
  ! main driver.
  PUBLIC :: INITIALIZE_SYSTEM_PARAMETERS
  PUBLIC :: INITIALIZE_SYSTEM_CONDITIONS
  PUBLIC :: INITIALIZE_SIMULATION_PARAMETERS
  PUBLIC :: INITIALIZE_OUTPUT_FILE
  PUBLIC :: INITIALIZE_SIMULATION
  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Reads the NAMELIST for input parameters.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_SYSTEM_PARAMETERS

    IMPLICIT NONE

    NAMELIST /GM10PARAMS/ oscFreqU, noiseStrU, dampB, oscFreqB, meanB, &
         & noiseStrB, dampGamma, meanGamma, noiseStrGamma, strForce, &
         & oscFreqForce

    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = GM10PARAMS)
    CLOSE(1000)
    
  END SUBROUTINE INITIALIZE_SYSTEM_PARAMETERS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @Author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Assigns the initial value to the state variable (u), additive bias (b), and
  !! multiplicative bias (gamma) in the GM10 system.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_SYSTEM_CONDITIONS

    IMPLICIT NONE

    COMPLEX(dp) :: initState !< Initial value of the state variable (u) in the
    !! GM10 system
    COMPLEX(dp) :: initAddBias !< Initial value of the additive bias (b) in the
    !! GM10 system
    REAL(dp) :: initMultBias !< Initial value of the multiplicative bias (gamma)
    !! in the GM10 system
    
    NAMELIST /GM10INITCOND/ initState, initAddBias, initMultBias

    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = GM10INITCOND)
    CLOSE(1000)

    stateVar = initState
    addBias = initAddBias
    multBias = initMultBias

  END SUBROUTINE INITIALIZE_SYSTEM_CONDITIONS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @Author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Reads the simulation parameters from file.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_SIMULATION_PARAMETERS

    IMPLICIT NONE
    
    NAMELIST /SIMPARAMS/ timeStepSize, timeStepCount, seed, timeStepScheme, &
         & runCount, outputFreq

    OPEN(1000, file = "NAMELIST")
    READ(1000, nml = SIMPARAMS)
    CLOSE(1000)

  END SUBROUTINE INITIALIZE_SIMULATION_PARAMETERS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @Author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Initializes the output file for the simulation.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_OUTPUT_FILE

    USE NETCDF
    USE UTILS

    IMPLICIT NONE

    ! outputFreq = 0 means no output
    IF (outputFreq .NE. 0_qb) THEN 
       ! Create the netCDF file. NF90_CLOBBER tells netCDF to overwrite the file
       ! if it already exists.
       CALL NETCDF_ERR_CHECK( NF90_CREATE(outFileName, NF90_CLOBBER, &
            & outFileID) )
       
       ! Define the dimensions for the output of the simulations. We want to
       ! have each variable represent a run of the simulation, so each variable
       ! will have time-step number, time, state variable value (real and complex
       ! parts), additive bias value (real and complex parts), and multiplicative
       ! bias value.
       CALL NETCDF_ERR_CHECK( NF90_DEF_DIM(outFileID, "step_count", &
            & timeStepCount/outputFreq + 1, stepCountID) ) ! Add one more step
            ! for the initial condition.
       CALL NETCDF_ERR_CHECK( NF90_DEF_DIM(outFileID, "var_count", 7, &
            & varCountID) )

       ! End define mode, this tells netCDF that we are done defining metadata.
       CALL NETCDF_ERR_CHECK( NF90_ENDDEF(outFileID) )
       
    END IF
    
  END SUBROUTINE INITIALIZE_OUTPUT_FILE

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @Author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Initializes the simulation, e.g., sets up the output.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  SUBROUTINE INITIALIZE_SIMULATION(runNum)

    USE NETCDF
    USE UTILS

    IMPLICIT NONE

    INTEGER(qb), INTENT(IN) :: runNum !< Run number
    CHARACTER(LEN=11) :: runName !< Name for the run netCDF variable

    ! outputFreq = 0 means no output
    IF (outputFreq .NE. 0_qb) THEN 
       WRITE(runName, "(A,I0.8)") "run", runNum

       ! Go back into define mode.
       CALL NETCDF_ERR_CHECK( NF90_REDEF(outFileID) )
       
       ! Create the variable to store all of the information from the
       ! simulation.
       CALL NETCDF_ERR_CHECK( NF90_DEF_VAR(outFileID, runName, NF90_DOUBLE, &
            & (/ varCountID, stepCountID /), runID ) )

       ! Put an attribute containing information about the run output
       CALL NETCDF_ERR_CHECK( NF90_PUT_ATT(outFileID, runID, "Info", &
            & "Contains the information for the run (each variable &
            &is contiguous): step number, time, state variable (real part) &
            &, state variable (complex part), additive bias (real part) &
            &, additive bias (complex part), and multiplicative bias.") )

       ! Set to not automatically fill with numbers, good for HPC applications
       CALL NETCDF_ERR_CHECK( NF90_DEF_VAR_FILL(outFileID, runID, 1, &
            & 0) )

       ! End define mode, this tells netCDF that we are done defining metadata.
       CALL NETCDF_ERR_CHECK( NF90_ENDDEF(outFileID) )
       

    END IF
    
  END SUBROUTINE INITIALIZE_SIMULATION
  
  
END MODULE INITIALIZE

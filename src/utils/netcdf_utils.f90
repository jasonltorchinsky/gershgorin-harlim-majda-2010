!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-HARLIM-MAJDA_10 Project of the the multi-model communication
! research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE       : Multi-Model Communication
! PROJECT     : GERSHGORIN-HARLIM-MAJDA_10
! MODULE      : UTILITIES
! URL         : https://github.com/jasonltorchinsky/gershgorin-harlim-majda-2010
! AFFILIATION : University of Wisconsin-Madison
! DATE        : Winter 2020
! REVISION    : 1.00
!
!> @author
!> Jason Torchinsky
!
!> @brief Contains general netCDF subroutines that will be used in various
!! places
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Unknown, NCAR
!> @brief
!> Checks the success/failure of a netCDF subroutine
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE NETCDF_ERR_CHECK_SBR(status)

  USE NETCDF

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  INTEGER(qb), INTENT(IN) :: status

  IF (status .NE. NF90_NOERR) THEN
     PRINT *, TRIM(NF90_STRERROR(status))
     STOP "Stopped due to netCDF error."
  END IF

END SUBROUTINE NETCDF_ERR_CHECK_SBR

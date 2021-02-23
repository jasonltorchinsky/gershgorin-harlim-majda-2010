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
! DATE        : Winter 2021
! REVISION    : 1.00
!
!> @author
!> Jason Torchinsky
!
!> @brief The analytic statistical solution for the Gershgorin-Majda 2010
!! system, using the trapezoidal rule for approximating the necessary integrals.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Torchinsky, University of Wisconsin-Madison
!> @brief
!> Uses the statistical analytic equations to progress that state variable,
!! additive bias, and multiplicative bias, with the trapezoidal rule for
!! approximating integrals. Although it is less readable, we will be using
!! the letter names for the variables given in the mathematical equations.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE TRAPEZOIDAL_ANALYTIC_SBR(currTime)

  USE INITIALIZE
  USE DET_FORCING
  USE UTILS

  IMPLICIT NONE

  INCLUDE 'kinds.h' ! Defines the kinds of real, complex, and integer

  REAL(dp), INTENT(IN) :: currTime !< Current time in the simulation.
  COMPLEX(dp) :: meanBNext !< Mean value of the additive bias (b) at the next
  !! time-step.
  REAL(dp) :: meanGammaNext !< Mean value of the multiplicative bias (gamma)
  !! at the next time-step.
  COMPLEX(dp) :: meanUNext !< Mean value of the state variable (u) at the
  !! next time step.
  REAL(dp) :: varAlpha !< Variance of the real part of b.
  REAL(dp) :: varBeta !< Variance of the imaginary part of b.
  REAL(dp) :: varGamma !< Variance of gamma.
  COMPLEX(dp) :: varU !< Variance of u, to be used in calculating
  !! other variances.
  COMPLEX(dp) :: covUUConjg !< Covariance of uConjg(u), to be used in calculating
  !! other variances.
  REAL(dp) :: varMu !< Variance of the real part of u.
  REAL(dp) :: varNu !< Varianve of the imaginary part of u.
  REAL(dp) :: covMuNu !< Covariance of the real and imaginary part of u.
  COMPLEX(dp) :: covUGamma !< Covariance of u and gamma.
  COMPLEX(dp) :: covUB !< Covariance of u and b.
  COMPLEX(dp) :: covUBConjg !< Covariance of u and CONJG(b).
  REAL(dp) :: covMuGamma !< Covariance of the real part of u and gamma.
  REAL(dp) :: covNuGamma !< Covariance of the imaginary part of u and gamma.
  REAL(dp) :: covMuAlpha !< Covariance of the real parts of u and b.
  REAL(dp) :: covNuAlpha !< Covariance of the imaginary part of u and the real
  !! part of b.
  REAL(dp) :: covMuBeta !< Covariance of the real part of u and the imaginary
  !! part of b.
  REAL(dp) :: covNuBeta !< Covariance of the imaginary parts of u and b.
  INTEGER(qb) :: i !< Counter for DO loops.
  REAL(dp) :: nextState(5) !< Vector to hold the value the real and imaginary
  !! parts of the additive bias, multiplicative bias, and state variable at the
  !! next time-step.
  REAL(dp) :: covMtx(5,5) !< Matrix of covariances of the real and imaginary
  !! parts of b, gamma, and u.
  REAL(dp) :: eValsCovMtx(5) !< Eigenvalues of the covaraince matrix.
  REAL(dp) :: work(170) !< Work array for the eigenvalue decomposition
  !! subroutine.
  INTEGER(qb) :: info !< Output integer for eigenvalue decompisiotn subroutine.


  ! Calculate the mean value of the additive bias (b), multiplicative bias (gamma)
  ! and state variable (u) at the next time-step
  meanBNext =  MEANBS_FUNC(currTime, currTime + timeStepSize)
  meanGammaNext = MEANGAMMAS_FUNC(currTime, currTime + timeStepSize)
  meanUNext = MEANUS_FUNC(currTime, currTime + timeStepSize)

  ! Calculate the variances and covariances of the real and imaginary parts of
  ! b, gamma, and u.
  varAlpha = (noiseStrB**2_qb / (4.0_dp * dampB)) &
       & * (1.0_dp - EXP(-2.0_dp * dampB * timeStepSize))
  varBeta = varAlpha
  varGamma = (noiseStrGamma**2_qb / (2.0_dp * dampGamma)) &
       & * (1.0_dp - EXP(-2.0_dp * dampGamma * timeStepSize))

  varU = VARU_FUNC(currTime, meanUNext)
  covUUConjg = COVUUCONJG_FUNC(currTime, meanUNext)
  varMu = 0.5_dp * REAL(varU + covUUConjg, dp)
  varNu = 0.5_dp * REAL(varU - covUUConjg, dp)
  covMuNu = 0.5_dp * REAL((0.0_dp, 1.0_dp) * (covUUConjg - varU))

  covUGamma = COVUGAMMA_FUNC(currTime, meanUNext)
  covMuGamma = REAL(covUGamma, dp)
  covNuGamma = -1.0_dp * REAL((0.0_dp, 1.0_dp) * covUGamma, dp)

  covUB = COVUB_FUNC(currTime, meanUNext)
  covUBConjg = COVUBCONJG_FUNC(currTime, meanUNext)
  covMuAlpha = 0.5_dp * REAL(covUB + covUBConjg, dp)
  covNuAlpha = -0.5_dp * REAL((0.0_dp, 1.0_dp) * (covUB + covUBConjg), dp)
  covMuBeta = 0.5_dp * REAL((0.0_dp, 1.0_dp) * (covUB - covUBConjg), dp)
  CovNuBeta = 0.5_dp * REAL(covUB - covUBConjg, dp)

  ! Generate five random Gaussian numbers.
  DO i = 1, 5
     CALL RAND_NORMAL_DP(nextState(i))
  END DO
  ! Set the covariance matrix
  covMtx = 0.0_dp
  covMtx(1,1) = varAlpha
  covMtx(4,1) = covMuAlpha
  covMtx(5,1) = covNuAlpha
  covMtx(2,2) = varBeta
  covMtx(4,2) = covMuBeta
  covMtx(5,2) = covNuBeta
  covMtx(3,3) = varGamma
  covMtx(4,3) = covMuGamma
  covMtx(5,3) = covNuGamma
  covMtx(1,4) = covMuAlpha
  covMtx(2,4) = covMuBeta
  covMtx(3,4) = covMuGamma
  covMtx(4,4) = varMu
  covMtx(5,4) = covMuNu
  covMtx(1,5) = covNuAlpha
  covMtx(2,5) = covNuBeta
  covMtx(3,5) = covNuGamma
  covMtx(4,5) = covMuNu
  covMtx(5,5) = varNu

  
  ! Obtain the eigen-decomposition of the covariance matrix.
  CALL DSYEV('V', 'U', 5, covMtx, 5, eValsCovMtx, work, 170, info)
  ! NOTE: covMtx is now the eigenvector matrix of covMtx.

  ! If the eigenvalue decomposition fails, print some failure information.
  IF (info .NE. 0_qb) THEN
     PRINT *, "LAPACK subroutine DSYEV has failed, info = ", info
     PRINT *, "Time-step: ", INT(currTime / timeStepSize, qb)
     ! Reset the covariance matrix
     covMtx = 0.0_dp
     covMtx(1,1) = varAlpha
     covMtx(4,1) = covMuAlpha
     covMtx(5,1) = covNuAlpha
     covMtx(2,2) = varBeta
     covMtx(4,2) = covMuBeta
     covMtx(5,2) = covNuBeta
     covMtx(3,3) = varGamma
     covMtx(4,3) = covMuGamma
     covMtx(5,3) = covNuGamma
     covMtx(1,4) = covMuAlpha
     covMtx(2,4) = covMuBeta
     covMtx(3,4) = covMuGamma
     covMtx(4,4) = varMu
     covMtx(5,4) = covMuNu
     covMtx(1,5) = covNuAlpha
     covMtx(2,5) = covNuBeta
     covMtx(3,5) = covNuGamma
     covMtx(4,5) = covMuNu
     covMtx(5,5) = varNu
     PRINT *, "Covariance matrix: "
     PRINT *, covMtx
     PRINT *, "Means for next time-step: "
     PRINT *, "b: ", meanBNext
     PRINT *, "gamma: ", meanGammaNext
     PRINT *, "u: ", meanUNext
     PRINT *, "Previous state: "
     PRINT *, "b: "
     PRINT *, addBias
     PRINT *, "gamma: "
     PRINT *, multBias
     PRINT *, "u: "
     PRINT *, stateVar
     ERROR STOP &
          & "Eigenvalue decomposition for the trapezoidal-analytic method failed."
  END IF

  ! Use the eigen-decomposition to adjust the random Gaussian numbers, and update
  ! b, gamma, and u.
  nextState =  MATMUL(covMtx, eValsCovMtx * nextState)
  addBias = meanBNext + nextState(1) + (0.0_dp, 1.0_dp) * nextState(2)
  multBias = meanGammaNext + nextState(3)
  stateVar = meanUNext + nextState(4) + (0.0_dp, 1.0_dp) * nextState(5)
  
CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of the additive bias (b) at the time s.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANBS_FUNC(currTime, s) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s !< Input for b(s).

    output = meanB + (addBias - meanB) &
         & * EXP((-1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB) &
         &       * (s - currTime))

  END FUNCTION MEANBS_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of the multiplicative bias (gamma) at time s.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(dp) FUNCTION MEANGAMMAS_FUNC(currTime, s) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s !< Input for gamma(s).

    output = meanGamma + (multBias - meanGamma) &
         & * EXP(-1.0_dp * dampGamma * (s - currTime))

  END FUNCTION MEANGAMMAS_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean value of J(s,t_{k+1}).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(dp) FUNCTION MEANJST_FUNC(currTime, s) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s !< Lower bound for the input of J(s,t_{k+1}).

    output = ((multBias - meanGamma) / (dampGamma)) &
         & * (EXP(-1.0_dp * dampGamma * (s - currTime)) &
         &    - EXP(-1.0_dp * dampGamma * timeStepSize))

  END FUNCTION MEANJST_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of J(s,t_{k+1}) and J(r,t_{k+1}).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(dp) FUNCTION COVJSTJRT_FUNC(currTime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s !< Lower bound for the input of J(s,t_{k+1}).
    REAL(dp), INTENT(IN) :: r !< Lower bound for the input of J(r,t_{k+1}).
    REAL(dp) :: nextTime !< Time of next time-step in simulation.

    nextTime = currTime + timeStepSize

    IF (r .LE. s) THEN

       output = (1.0_dp + EXP(-1.0_dp * dampGamma * (s - r))) &
            &    - 2.0_dp * dampGamma * (nextTime - s) &
            &    - EXP(-1.0_dp * dampGamma * (s + timeStepSize - currTime)) &
            &      * ((1.0_dp + EXP(dampGamma * (s - r))) &
            &         + (EXP(2.0_dp * dampGamma * (s - currTime)) &
            &            + EXP(dampGamma * (s + r - 2.0_dp * currTime))) &
            &         - (EXP(dampGamma * (nextTime - s)) &
            &            + EXP(-1.0_dp * dampGamma * (nextTime - r))))

    ELSE IF (s .LT. r) THEN

       output = (1.0_dp + EXP(-1.0_dp * dampGamma * (r - s))) &
            &    - 2.0_dp * dampGamma * (nextTime - r) &
            &    - EXP(-1.0_dp * dampGamma * (r + timeStepSize - currTime)) &
            &      * ((1.0_dp + EXP(dampGamma * (r - s))) &
            &         + (EXP(2.0_dp * dampGamma * (r - currTime)) &
            &            + EXP(dampGamma * (r + s - 2.0_dp * currTime))) &
            &         - (EXP(dampGamma * (nextTime - r)) &
            &            + EXP(-1.0_dp * dampGamma * (nextTime - s))))

    END IF

    output = -1.0_dp * ((noiseStrGamma ** 2_qb) &
            &                 / (2.0_dp * (dampGamma ** 3_qb))) * output


  END FUNCTION COVJSTJRT_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of the state variable (u) at time s.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANUS_FUNC(currTime, s) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s !< Input for u(s).
    INTEGER(qb) :: i !< Counters for DO loops.
    REAL(dp) :: integrateSize !< Size of the intervals used for the trapezoidal
    !! rule.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: meanJst !< <J(s,t_{k+1})>
    REAL(dp) :: covJstJrt !< Cos(J(r,t_{k+1}),J(r,t_{k+1}))
    COMPLEX(dp) :: meanBs !< <b(s)>
    COMPLEX(dp) :: detForce !< Deterministic forcing term.

    ! Set up trapezoidal integraion.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-integral term
    meanJst = MEANJST_FUNC(currTime, currTime)
    covJstJrt =  COVJSTJRT_FUNC(currTime, currTime, currTime)
    output = stateVar * EXP(-1.0_dp * meanJst + 0.5_dp * covJstJrt &
         &                  + (-1.0_dp * meanGamma &
         &                     + (0.0_dp, 1.0_dp) * oscFreqU) * timeStepSize)
    
    ! Add the integral term. Since we are doing a single trapezoidal rule,
    ! we have two terms outside of the summation. The first term (t_k) and last
    ! term (t_{k+1}).
    meanBs =  MEANBS_FUNC(currTime, currTime)
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    meanJst = MEANJST_FUNC(currTime, currTime)
    covJstJrt = COVJSTJRT_FUNC(currTime, currTime, currTime)
    output = output + 0.5_dp * integrateSize * (meanBs + detForce) &
         &            * EXP(-1.0_dp * meanJst + 0.5_dp * covJstJrt &
         &                  + (-1.0_dp * meanGamma &
         &                     + (0.0_dp, 1.0_dp) * oscFreqU) * timeStepSize)
    meanBs =  MEANBS_FUNC(currTime, currTime + timeStepSize)
    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    meanJst = MEANJST_FUNC(currTime, currTime + timeStepSize)
    covJstJrt = COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         & currTime + timeStepSize)
    output = output + 0.5_dp * integrateSize * (meanBs + detForce) &
         &            * EXP(-1.0_dp * meanJst + 0.5_dp * covJstJrt)
    ! The middle terms
    DO i = 1, integrateCount - 1
       meanBs =  MEANBS_FUNC(currTime, currTime + i * integrateSize)
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       meanJst = MEANJST_FUNC(currTime, currTime + i * integrateSize)
       covJstJrt = COVJSTJRT_FUNC(currTime, currTime + i * integrateSize, &
            & currTime + i * integrateSize)
       output = output + integrateSize * (meanBs + detForce) &
            &            * EXP(-1.0_dp * meanJst + 0.5_dp * covJstJrt &
            &                  + (-1.0_dp * meanGamma &
            &                     + (0.0_dp, 1.0_dp) * oscFreqU) &
            &                    * (timeStepSize - (i * integrateSize)))
    END DO

  END FUNCTION MEANUS_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of |A|^2 (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANABSA2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.

    output = ABS(stateVar)**2_qb * EXP(2.0_dp &
         &                             * (COVJSTJRT_FUNC(currTime, currTime, &
         &                                               currTime) &
         &                                - MEANJST_FUNC(currTime, currTime) &
         &                                - meanGamma * timeStepSize))

  END FUNCTION MEANABSA2_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of b(s) and b(r) for s, r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVBSBR_FUNC(currtime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s, r !< Times for b(s), b(r)

    output = ((noiseStrB**2_qb) / (2.0_dp * dampB)) &
         & * EXP(-1.0_dp * dampB * (s + r - 2.0_dp * currTime)) &
         & * (EXP(2.0_dp * dampB * (MIN(s,r) - currTime)) - 1.0_dp)

  END FUNCTION COVBSBR_FUNC


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the value of bvar(s,r) for given s, r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION BVAR_FUNC(currTime, s, r) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s, r !< Inputs for bvar(s,r)
    COMPLEX(dp) :: detForceS, detForceR !< Deterministic forcing at times s and r

    CALL DET_FORCING_DEFAULT(s, detForceS)
    CALL DET_FORCING_DEFAULT(r, detForceR)

    output = EXP(-1.0_dp * MEANJST_FUNC(currTime, s) - MEANJST_FUNC(currTime, r) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, s, s) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, r, r) &
         &       + COVJSTJRT_FUNC(currTime, s, r) &
         &       - meanGamma * (2.0_dp * (currTime + timeStepSize) - s - r) &
         &       + (0.0_dp, 1.0_dp) * oscFreqU * (s - r)) &
         & * (COVBSBR_FUNC(currTime, s, r) &
         &    + (MEANBS_FUNC(currTime, s) + detForceS) &
         &       * (CONJG(MEANBS_FUNC(currTime, r)) + CONJG(detForceR)))

!!$    IF (ISNAN(REAL(output, dp))) THEN
!!$       PRINT *, "s, r: ", s, r
!!$       PRINT *, -1.0_dp * MEANJST_FUNC(currTime, s) - MEANJST_FUNC(currTime, r)
!!$       PRINT *,  0.5_dp * COVJSTJRT_FUNC(currTime, s, s) 
!!$       PRINT *,  0.5_dp * COVJSTJRT_FUNC(currTime, r, r) 
!!$       PRINT *,  COVJSTJRT_FUNC(currTime, s, r) 
!!$       PRINT *,  meanGamma * (2.0_dp * (currTime + timeStepSize) - s - r) &
!!$            &       + (0.0_dp, 1.0_dp) * oscFreqU * (s - r)
!!$    END IF
    
  END FUNCTION BVAR_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of |B|^2 (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANABSB2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: i, j !< Counters for DO loops
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule, same for both integral directions.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the
    !! trapezoidal rule, same for both directions.

    ! Set up the trapezoidals integration.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the terms that are outside of the summation terms
    output = 0.25_dp * integrateSize**2_qb &
         & * (BVAR_FUNC(currTime, currTime, currTime) &
         &    + BVAR_FUNC(currTime, currTime, currTime + timeStepSize) &
         &    + BVAR_FUNC(currTime, currTime + timeStepSize, currTime) &
         &    + BVAR_FUNC(currTime, currTime + timeStepSize, &
         &                currTime + timeStepSize))

    ! Add the single summation term
    DO i = 1, integrateCount - 1
       output = output + 0.5_dp * integrateSize**2_qb &
            &            * (BVAR_FUNC(currTime, currTime, &
            &                         currTime + i * integrateSize) &
            &               + BVAR_FUNC(currTime, currTime + timeStepSize, &
            &                           currTime + i * integrateSize) &
            &               + BVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                           currTime) &
            &               + BVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                           currTime + timeStepSize))
    END DO

    ! Add the double summation term.
    DO j = 1, integrateCount - 1
       DO i = 1, integrateCount - 1
          output = output + integrateSize**2_qb &
               &            * BVAR_FUNC(currTime, currTime + i * integrateSize, &
               &                        currTime + j * integrateSize)
       END DO
    END DO


  END FUNCTION MEANABSB2_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of |C|^2 (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANABSC2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms of the trapezoidal rule.
    output = 0.5_dp * integrateSize &
         & * (EXP(2.0_dp &
         &        * (COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &           - MEANJST_FUNC(currTime, currTime) &
         &           - meanGamma * timeStepSize)) &
         &    + EXP(2.0_dp &
         &          * (COVJSTJRT_FUNC(currTime, &
         &                            currTime + timeStepSize, &
         &                            currTime + timeStepSize) &
         &             - MEANJST_FUNC(currTime, currTime + timeStepSize))))

    ! Add the summation term of the trapezoidal rule
    DO i = 1, integrateCount - 1
       output = output &
            & + integratesize &
            &   * EXP(2.0_dp &
            &         * (COVJSTJRT_FUNC(currTime, &
            &                           currTime + i * integrateSize, &
            &                           currTime + i * integrateSize) &
            &            - MEANJST_FUNC(currTime, &
            &                           currTime + i * integrateSize) &
            &            - meanGamma * (timeStepSize &
            &                                  - (i * integrateSize))))
    END DO

    ! Scale the by the noise strength of u
    output = output * noiseStrU**2_qb

  END FUNCTION MEANABSC2_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of ACONJG(B) (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANABCONJG_FUNC(currTime) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms.
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = 0.5_dp * stateVar &
         & * (CONJG(MEANBS_FUNC(currTime, currTime)) + CONJG(detForce)) &
         & * EXP(-2.0_dp * MEANJST_FUNC(currTime, currTime) &
         &       + 2.0_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       - 2.0_dp * meanGamma * timeStepSize)
    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output + 0.5_dp * stateVar &
         & * (CONJG(MEANBS_FUNC(currTime, currTime + timeStepSize)) + CONJG(detForce)) &
         & * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &       - MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                 currTime + timeStepSize) &
         &       + COVJSTJRT_FUNC(currTime, currTime, currTime + timeStepSize) &
         &       - meanGamma * timeStepSize &
         &       + (0.0_dp, 1.0_dp) * oscFreqU * timeStepSize)

    ! Add the summation term
    DO i = 1, integrateCount
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output + 0.5_dp * stateVar &
            & * (CONJG(MEANBS_FUNC(currTime, currTime + i * integrateSize)) &
            &    + CONJG(detForce)) &
            & * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
            &       - MEANJST_FUNC(currTime, currTime + i * integrateSize) &
            &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
            &       + 0.5_dp * COVJSTJRT_FUNC(currTime, &
            &                                 currTime + i * integrateSize, &
            &                                 currTime + i * integrateSize) &
            &       + COVJSTJRT_FUNC(currTime, currTime, &
            &                        currTime + i * integrateSize) &
            &       - meanGamma * (2.0_dp * timeStepSize - i * integrateSize) &
            &       + (0.0_dp, 1.0_dp) * oscFreqU * (i * integrateSize))       
    END DO

    ! Scale the output by the trapezoid size.
    output = output * integrateSize

  END FUNCTION MEANABCONJG_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the variance of u.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION VARU_FUNC(currTime, meanU) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean of u_{k+1}

    output = MEANABSA2_FUNC(currTime) + MEANABSB2_FUNC(currTime) &
         & + MEANABSC2_FUNC(currTime) &
         & + 2.0_dp * REAL(MEANABCONJG_FUNC(currTime), dp) &
         & - ABS(meanU)**2_qb

!!$    PRINT *, "<|A|^2>: ", MEANABSA2_FUNC(currTime)
!!$    PRINT *, "<|B|^2>: ", MEANABSB2_FUNC(currTime)
!!$    PRINT *, "<|C|^2>: ", MEANABSC2_FUNC(currTime)
!!$    PRINT *, 2.0_dp * REAL(MEANABCONJG_FUNC(currTime), dp)
!!$    PRINT *, ABS(meanU)**2_qb

  END FUNCTION VARU_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of A^2 (see Cov(u, CONJG(u))).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANA2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.

    output = (stateVar**2_qb) &
         & * EXP(2.0_dp * (COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &                 + (meanGamma + (0.0_dp, 1.0_dp) * oscFreqU) &
         &                   * timeStepSize &
         &                 - MEANJST_FUNC(currTime, currTime)))

  END FUNCTION MEANA2_FUNC


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of b(s) and CONJG(b(r)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVBSBRCONJG_FUNC(currTime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s, r !< Inputs to b(s) and CONJG(b(r)).
    COMPLEX(dp) :: lambdaB !< The parameter lamnbda_b in the equations.

    lambdaB = -1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB
    output = (noiseStrB**2_qb) / (2.0_dp * lambdaB) &
         & * EXP(lambdaB * (s + r - 2.0_dp * currTime))

  END FUNCTION COVBSBRCONJG_FUNC


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the value of bcovar(s, r) at specified s and r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION BCOVAR_FUNC(currTime, s, r) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s, r !< Inputs to bcovar(s, r)
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat in the equations.
    COMPLEX(dp) :: detForceS, detForceR !< The determinsitc forcing at s, r.

    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    CALL DET_FORCING_DEFAULT(s, detForceS)
    CALL DET_FORCING_DEFAULT(r, detForceR)
    output = (COVBSBRCONJG_FUNC(currTime, s, r) &
         &    + (MEANBS_FUNC(currTime, s) + detForceS) &
         &       * (MEANBS_FUNC(currTime, r) + detForceR)) &
         & * EXP(-1.0_dp * MEANJST_FUNC(currTime, s) &
         &       - MEANJST_FUNC(currTime, r) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, s, s) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, r, r) &
         &       + COVJSTJRT_FUNC(currTime, s, r) &
         &       + lambdaHat * (2.0_dp * (currTime + timeStepSize - s - r)))
    

  END FUNCTION BCOVAR_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of B^2 (see cov(u, CONJG(u))).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANB2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: i, j !< Counters for DO loops
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule, same for both integral directions.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the
    !! trapezoidal rule, same for both directions.

    ! Set up the trapezoidals integration.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

   

    ! Add the terms that are outside of the summation terms
    output = 0.25_dp * (BCOVAR_FUNC(currTime, currTime, currTime) &
         & + BCOVAR_FUNC(currTime, currTime, currTime + timeStepSize) &
         & + BCOVAR_FUNC(currTime, currTime + timeStepSize, currTime) &
         & + BCOVAR_FUNC(currTime, currTime + timeStepSize, &
         &             currTime + timeStepSize))

    ! Add the single summation term
    DO i = 1, integrateCount - 1
       output = output &
            & + 0.5_dp * (BCOVAR_FUNC(currTime, currTime, &
            &                       currTime + i * integrateSize) &
            &             + BCOVAR_FUNC(currTime, currTime + timeStepSize, &
            &                         currTime + i * integrateSize) &
            &             + BCOVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                         currTime) &
            &             + BCOVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                         currTime + timeStepSize))
    END DO


    ! Add the double summation term.
    DO j = 1, integrateCount - 1
       DO i = 1, integrateCount - 1
          output = output + BCOVAR_FUNC(currTime, currTime + i * integrateSize, &
               &                      currTime + j * integrateSize)
       END DO
    END DO


    ! Scale output by integrateSize
    output = output * integrateSize**2_qb

  END FUNCTION MEANB2_FUNC


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of C^2 (see cov(u, CONJG(u))).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANC2_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat in the equations.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms of the trapezoidal rule.
    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    output = 0.5_dp * (EXP(2.0_dp &
         &                * (COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &                   - MEANJST_FUNC(currTime, currTime) &
         &                   + lambdaHat * timeStepSize)) &
         &             + EXP(2.0_dp &
         &                   * (COVJSTJRT_FUNC(currTime, &
         &                      currTime + timeStepSize, &
         &                      currTime + timeStepSize) &
         &                   - MEANJST_FUNC(currTime, currTime + timeStepSize))))

    ! Add the summation term of the trapezoidal rule
    DO i = 1, integrateCount - 1
       output = output + EXP(2.0_dp &
            &                * (COVJSTJRT_FUNC(currTime, &
            &                   currTime + i * integrateSize, &
            &                   currTime + i * integrateSize) &
            &                   - MEANJST_FUNC(currTime, &
            &                                  currTime + i * integrateSize) &
            &                   + lambdaHat * (timeStepSize &
            &                                  - (i * integrateSize))))
    END DO

    ! Scale the numerical integral by the noise strength of u and the
    ! integration interval size
    output = output * noiseStrU**2_qb * integrateSize

  END FUNCTION MEANC2_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of AB (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANAB_FUNC(currTime) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: i !< Counter for DO loops.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat  in the equations.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms.
    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = 0.5_dp * stateVar &
         & * (MEANBS_FUNC(currTime, currTime) + detForce) &
         & * EXP(-2.0_dp * MEANJST_FUNC(currTime, currTime) &
         &       + 2.0_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 2.0_dp * lambdaHat * timeStepSize)
    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output + 0.5_dp * stateVar &
         & * (MEANBS_FUNC(currTime, currTime + timeStepSize) + detForce) &
         & * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &       - MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                 currTime + timeStepSize) &
         &       + COVJSTJRT_FUNC(currTime, currTime, currTime + timeStepSize) &
         &       + lambdaHat * timeStepSize)

    ! Add the summation term
    DO i = 1, integrateCount
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output + stateVar &
            & * (MEANBS_FUNC(currTime, currTime + i * integrateSize) &
            &    + detForce) &
            & * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
            &       - MEANJST_FUNC(currTime, currTime + i * integrateSize) &
            &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
            &       + 0.5_dp * COVJSTJRT_FUNC(currTime, &
            &                                 currTime + i * integrateSize, &
            &                                 currTime + i * integrateSize) &
            &       + COVJSTJRT_FUNC(currTime, currTime, &
            &                        currTime + i * integrateSize) &
            &       + lambdaHat * (2.0_dp * timeStepSize - i * integrateSize))       
    END DO

    ! Scale the output by the trapezoid size.
    output = output * integrateSize

  END FUNCTION MEANAB_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of uCONJG(u).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVUUCONJG_FUNC(currTime, meanU) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean of u_{k+1}

    output = MEANA2_FUNC(currTime) + MEANB2_FUNC(currTime) &
         & + MEANC2_FUNC(currTime) &
         & + 2.0_dp * MEANAB_FUNC(currTime) - meanU**2_qb

  END FUNCTION COVUUCONJG_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of ugamma (see Cov(u, gamma)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANUGAMMA_FUNC(currTime, meanU) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of u at the next time-step.
    INTEGER(qb) :: i !< Counter for DO loops.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat  in the equations.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU

    ! Add the non-integral terms.

    output = meanGamma * meanU &
         & + stateVar * (multBias - meanGamma) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat - dampGamma) * timeStepSize)

    ! Add the not-do loop terms of the integral
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (multBias - meanGamma) &
         &   * (MEANBS_FUNC(currTime, currTime) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat - dampGamma) * timeStepSize)

    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (multBias - meanGamma) &
         &   * (MEANBS_FUNC(currTime, currTime + timeStepSize) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                   currTime + timeStepSize) &
         &         - dampGamma * timeStepSize)

    DO i = 1, integrateCount - 1
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output &
            & + integrateSize * (multBias - meanGamma) &
            &   * (MEANBS_FUNC(currTime, currTime + i * integrateSize) + detForce) &
            &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime + i * integrateSize) &
            &         + 0.5_dp * COVJSTJRT_FUNC(currTime, &
            &                                   currTime + i * integrateSize, &
            &                                   currTime + i * integrateSize) &
            &         + lambdaHat * (timeStepSize - i * integrateSize) &
            &         - dampGamma * timeStepSize)
    END DO

  END FUNCTION MEANUGAMMA_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of u and gamma.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVUGAMMA_FUNC(currTime, meanU) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of U at the next time-step.

    output = MEANUGAMMA_FUNC(currTime, meanU) &
         & - meanU * MEANGAMMAS_FUNC(currTime, currTime + timeStepSize)

  END FUNCTION COVUGAMMA_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of uCONJG(b) (see Cov(u, b)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANUBCONJG_FUNC(currTime, meanU) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of u at the next time-step.
    INTEGER(qb) :: i !< Counter for DO loops.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat in the equations.
    COMPLEX(dp) :: lambdaB !< The parameter lambda_b in the equations.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    lambdaB = -1.0_dp * meanB + (0.0_dp, 1.0_dp) * oscFreqB

    ! Add the non-integral terms.

    output = CONJG(meanB) * meanU &
         & + stateVar * (CONJG(addBias) - CONJG(meanB)) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat + CONJG(lambdaB)) * timeStepSize)

    ! Add the not-do loop terms of the integral
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (CONJG(addBias) - CONJG(meanB)) &
         &   * (MEANBS_FUNC(currTime, currTime) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat + CONJG(lambdaB)) * timeStepSize)

    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (CONJG(addBias) - CONJG(meanB)) &
         &   * (MEANBS_FUNC(currTime, currTime + timeStepSize) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                   currTime + timeStepSize) &
         &         + CONJG(lambdaB) * timeStepSize)

    DO i = 1, integrateCount - 1
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output &
            & + integrateSize * (CONJG(addBias) - CONJG(meanB)) &
            &   * (MEANBS_FUNC(currTime, currTime + i * integrateSize) &
            &      + detForce) &
            &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, &
            &                                currTime + i * integrateSize) &
            &         + 0.5_dp * COVJSTJRT_FUNC(currTime, &
            &                                   currTime + i * integrateSize, &
            &                                   currTime + i * integrateSize) &
            &         + lambdaHat * (timeStepSize - i * integrateSize) &
            &         + CONJG(lambdaB) * timeStepSize)
    END DO

  END FUNCTION MEANUBCONJG_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of u and b.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVUB_FUNC(currTime, meanU) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of U at the next time-step.

    output = MEANUBCONJG_FUNC(currTime, meanU) &
         & - meanU * MEANBS_FUNC(currTime, currTime + timeStepSize)

  END FUNCTION COVUB_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of ub (see Cov(u, CONJG(b))).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANUB_FUNC(currTime, meanU) RESULT(output)

    USE DET_FORCING

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of u at the next time-step.
    INTEGER(qb) :: i !< Counter for DO loops.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat in the equations.
    COMPLEX(dp) :: lambdaB !< The parameter lambda_b in the equations.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 1000_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    lambdaB = -1.0_dp * meanB + (0.0_dp, 1.0_dp) * oscFreqB

    ! Add the non-integral terms.

    output = meanB * meanU &
         & + stateVar * (addBias - meanB) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat + lambdaB) * timeStepSize)

    ! Add the not-do loop terms of the integral
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (addBias - meanB) &
         &   * (MEANBS_FUNC(currTime, currTime) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &         + (lambdaHat + lambdaB) * timeStepSize)

    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output &
         & + 0.5_dp * integrateSize * (addBias - meanB) &
         &   * (MEANBS_FUNC(currTime, currTime + timeStepSize) + detForce) &
         &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &         + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                   currTime + timeStepSize) &
         &         + lambdaB * timeStepSize)

    DO i = 1, integrateCount - 1
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output &
            & + integrateSize * (addBias - meanB) &
            &   * (MEANBS_FUNC(currTime, currTime + i * integrateSize) + detForce) &
            &   * EXP(-1.0_dp * MEANJST_FUNC(currTime, currTime + i * integrateSize) &
            &         + 0.5_dp * COVJSTJRT_FUNC(currTime, &
            &                                   currTime + i * integrateSize, &
            &                                   currTime + i * integrateSize) &
            &         + lambdaHat * (timeStepSize - i * integrateSize) &
            &         + lambdaB * timeStepSize)
    END DO

  END FUNCTION MEANUB_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Torchinsky, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of u and CONJG(b).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVUBCONJG_FUNC(currTime, meanU) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    COMPLEX(dp), INTENT(IN) :: meanU !< Mean value of U at the next time-step.

    output = MEANUBCONJG_FUNC(currTime, meanU) &
         & - meanU * MEANBS_FUNC(currTime, currTime + timeStepSize)

  END FUNCTION COVUBCONJG_FUNC


END SUBROUTINE TRAPEZOIDAL_ANALYTIC_SBR

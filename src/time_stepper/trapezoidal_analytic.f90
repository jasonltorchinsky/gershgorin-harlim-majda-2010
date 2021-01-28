!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The GERSHGORIN-MAJDA_10 Project of the the multi-model communication research.
! To create the initial solver.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! TITLE            : Multi-Model Communication
! PROJECT          : GERSHGORIN-MAJDA_10
! MODULE           : TIME_STEPPER
! URL              : https://github.com/jasonlturner/gershgorin-majda-10
! AFFILIATION      : University of Wisconisn-Madison
! DATE             : Winter 2021
! REVISION         : 1.00
!
!> @author
!> Jason Turner
!
!> @brief The analytic statistical solution for the Gershgorin-Majda 2010
!! system, using the trapezoidal rule for approximating the necessary integrals.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!> @author Jason Turner, University of Wisconsin-Madison
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
  REAL(dp) :: varAlpha !< Variance of the real part of b
  REAL(dp) :: varBeta !< Variance of the imaginary part of b
  REAL(dp) :: varGamma !< Variance of gamma
  COMPLEX(dp) :: varU !< Variance of u, to be used in calculating
  !! other variances.
  COMPLEX(dp) :: covUUConjg !< Covariance of uConjg(u), to be used in calculating
  !! other variances.
  REAL(dp) :: varMu !< Variance of the real part of u
  REAL(dp) :: varNu !< Varianve of the imaginary part of u
  REAL(dp) :: covMuNu !< Covariance of the real and imaginary part of u.


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

  ! TO-DO: COVARIANCE OF COVARIANCE OF U WITH GAMMA,
  ! COVARIANCE OF U WITH B, COVARIANCE OF U WITH CONJG(B).
  

CONTAINS

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of the additive bias (b) at the time s.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANBS_FUNC(currTime, s) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s !< Input for b(s).

    output = meanB + (addBias - meanB) &
         & * EXP((-1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB) &
         &       * (s - currTime)

  END FUNCTION MEANBS_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
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
  !> @author Jason Turner, University of Wisconsin-Madison
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
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of J(s,t_{k+1}) and J(r,t_{k+1}).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  REAL(dp) FUNCTION COVJSTJRT_FUNC(currTime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s !< Lower bound for the input of J(s,t_{k+1}).
    REAL(dp), INTENT(IN) :: r !< Lower bound for the input of J(r,t_{k+1}).

    IF (s .LE. r) THEN

       output = -1.0_dp * ((noiseStrGamma ** 2_qb) &
            &                 / (2.0_dp * (dampGamma ** 3_qb))) &
            & * ((1.0_dp + EXP(-1.0_dp * dampGamma * (r - s))) &
            &    - 2.0_dp * dampGamma * (currTime + timeStepSize - s) &
            &    - EXP(-1.0_dp * dampGamma * (s + timeStepSize - currTime)) &
            &      * ((1.0_dp + EXP(dampGamma * (s - r))) &
            &         + (EXP(2.0_dp * dampGamma * (s - currTime)) &
            &            + EXP(dampGamma * (s + r - 2.0_dp * currTime))) &
            &         - (EXP(dampGamma * (currTime + timeStepSize - s)) &
            &            + EXP(-1.0_dp * dampGamma * (currTime &
            &                                         + timeStepSize - r)))))

    ELSE IF (r .LT. s) THEN

      output = -1.0_dp * ((noiseStrGamma ** 2_qb) &
            &                 / (2.0_dp * (dampGamma ** 3_qb))) &
            & * ((1.0_dp + EXP(-1.0_dp * dampGamma * (s - r))) &
            &    - 2.0_dp * dampGamma * (currTime + timeStepSize - r) &
            &    - EXP(-1.0_dp * dampGamma * (r + timeStepSize - currTime)) &
            &      * ((1.0_dp + EXP(dampGamma * (r - s))) &
            &         + (EXP(2.0_dp * dampGamma * (r - currTime)) &
            &            + EXP(dampGamma * (r + s - 2.0_dp * currTime))) &
            &         - (EXP(dampGamma * (currTime + timeStepSize - r)) &
            &            + EXP(-1.0_dp * dampGamma * (currTime &
            &                                         + timeStepSize - s)))))
       
    END IF

  END FUNCTION COVJSTJRT_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of the state variable (u) at time s.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANUS_FUNC(currTime, s) RESULT(output)

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
    integrateCount = 10_qb
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
  !> @author Jason Turner, University of Wisconsin-Madison
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
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of b(s) and b(r) for s, r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVBSBR_FUNC(currtime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation
    REAL(dp), INTENT(IN) :: s, r !< Times for b(s), b(r)

    output = ((noiseStrB**2_qb) / (2.0_dp * dampB)) &
         & * EXP(-1.0_dp & dampB * (s + r - 2.0_dp * currTime)) &
         & * (EXP(2.0_dp * dampB * (MIN(s,r) - currTime)) - 1.0_dp)
    
  END FUNCTION COVBSBR_FUNC

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the value of bvar(s,r) for given s, r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION BVAR_FUNC(currTime, s, r) RESULT(output)

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
    
  END FUNCTION BVAR_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
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
    integrateCount = 10_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the terms that are outside of the summation terms
    output = 0.25_dp * (BVAR_FUNC(currTime, currTime, currTime) &
         & + BVAR_FUNC(currTime, currTime, currTime + timeStepSize) &
         & + BVAR_FUNC(currTime, currTime + timeStepSize, currTime) &
         & + BVAR_FUNC(currTime, currTime + timeStepSize, &
         &             currTime + timeStepSize)

    ! Add the single summation term
    DO i = 1, integrateCount - 1
       output = output &
            & + 0.5_dp * (BVAR_FUNC(currTime, currTime, &
            &                       currTime + i * integrateSize) &
            &             + BVAR_FUNC(currTime, currTime + timeStepSize, &
            &                         currTime + i * integrateSize) &
            &             + BVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                         currTime) &
            &             + BVAR_FUNC(currTime, currTime + i * integrateSize, &
            &                         currTime + timeStepSize))
    END DO

    ! Add the double summation term.
    DO j = 1, integrateCount - 1
       DO i = 1, integrateCount - 1
          output = output + BVAR_FUNC(currTime, currTime + i * integrateSize, &
               &                      currTime + j * integrateSize)
       END DO
    END DO

    ! Scale output by integrateSize
    output = output * integrateSize**2_qb
    
  END FUNCTION MEANABSB2_FUNC
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
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
    integrateCount = 10_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms of the trapezoidal rule.
    output = 0.5_dp * (EXP(2.0_dp &
         &                * (COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &                   - MEANJST_FUNC(currTime, currTime) &
         &                   - meanGamma * timeStepSize)) &
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
            &                   - meanGamma * (timeStepSize &
            &                                  - (i * integrateSize))))
    END DO

    ! Scale the numerical integral by the noise strength of u and the
    ! integration interval size
    output = output * noiseStrU**2_qb * integrateSize

  END FUNCTION MEANABSC2_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of ACONJG(B) (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANABCONJG_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 10_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms.
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = 0.5_dp * stateVar &
         & * (CONJG(MEANBS_FUNC(currTime, currTime)) + CONJG(detForce)) &
         & * EXP(-2.0_dp * MEANJST(currTime, currTime) &
         &       + 2.0_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       - 2.0_dp * meanGamma * timeStepSize)
    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output + 0.5_dp * stateVar &
         & * (CONJG(MEANBS_FUNC(currTime, currTime + timeStepSize)) + CONJG(detForce)) &
         & * EXP(-1.0_dp * MEANJST(currTime, currTime) &
         &       - MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                 currTime + timeStepSize) &
         &       + COVJSTJRT_FUNC(currTime, currTime, currTime + timeStepSize) &
         &       - meanGamma * timeStepSize &
         &       + (0.0_dp, 1.0_dp) * oscFreqU * timeStepSize)

    ! Add the summation term
    DO i = 1, integrationCount
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output + 0.5_dp * stateVar &
            & * (CONJG(MEANBS_FUNC(currTime, currTime + i * integrateSize)) &
            &    + CONJG(detForce)) &
            & * EXP(-1.0_dp * MEANJST(currTime, currTime) &
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
  !> @author Jason Turner, University of Wisconsin-Madison
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

  END FUNCTION VARU_FUNC

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
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
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the covariance of b(s) and CONJG(b(r)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION COVBSCONJGBR_FUNC(currTime, s, r)) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s, r !< Inputs to b(s) and CONJG(b(r)).
    COMPLEX(dp) :: lambdaB !< The parameter lamnbda_b in the equations.

    lambdaB = -1.0_dp * dampB + (0.0_dp, 1.0_dp) * oscFreqB
    output = (noiseStrB**2_qb) / (2.0_dp * lambdaB) &
         & * EXP(lambdaB * (s + r - 2.0_dp * currTime))

  END FUNCTION COVBSCONJGBR_FUNC

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the value of bcovar(s, r) at specified s and r.
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION BCOVAR_FUNC(currTime, s, r) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    REAL(dp), INTENT(IN) :: s, r !< Inputs to bcovar(s, r)
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat in the equations.
    COMPLEX(dp) :: detForceS, detForceR !< The determinsitc forcing at s, r.

    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    CALL DET_FORCE_DEFAULT(s, detForceS)
    CALL DET_FORCE_DEFAULT(r, detForceR)
    output = (COVBSCONJGBR_FUNC(currTime, s, r) &
         &    + (MEANBS_FUNC(currTime, s) + detForceS) &
         &       * (MEANBS_FUNC(currTime, r) + detForceR)) &
         & * EXP(-1.0_dp * MEANJST_FUNC(currTime, s) &
         &       - MEANJST_FUNC(currTime, r) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, s, s) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, r, r) &
         &       + COVJSTJRT(currTime, s, r) &
         &       + lambdaHat * (2.0_dp * (currTime + timeStepSize - s - r)))

  END FUNCTION BCOVAR_FUNC

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !> @author Jason Turner, University of Wisconsin-Madison
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
    integrateCount = 10_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the terms that are outside of the summation terms
    output = 0.25_dp * (BCOVAR_FUNC(currTime, currTime, currTime) &
         & + BCOVAR_FUNC(currTime, currTime, currTime + timeStepSize) &
         & + BCOVAR_FUNC(currTime, currTime + timeStepSize, currTime) &
         & + BCOVAR_FUNC(currTime, currTime + timeStepSize, &
         &             currTime + timeStepSize)

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
  !> @author Jason Turner, University of Wisconsin-Madison
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
    integrateCount = 10_qb
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
  !> @author Jason Turner, University of Wisconsin-Madison
  !> @brief
  !> Calculates the mean of AB (see var(u)).
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  COMPLEX(dp) FUNCTION MEANAB_FUNC(currTime) RESULT(output)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: currTime !< Current time in simulation.
    INTEGER(qb) :: integrateCount !< Number of intervals used for the trapezoidal
    !! rule.
    REAL(dp) :: integrateSize !< Size of intervals used for the trapezoidal
    !! rule.
    COMPLEX(dp) :: lambdaHat !< The parameter lambda^hat  in the equations.
    COMPLEX(dp) :: detForce !< Deterministic forcing.

    ! Set up trapezoidal rule.
    integrateCount = 10_qb
    integrateSize = timeStepSize / REAL(integrateCount, dp)

    ! Add the non-summation terms.
    lambdaHat = -1.0_dp * meanGamma + (0.0_dp, 1.0_dp) * oscFreqU
    CALL DET_FORCING_DEFAULT(currTime, detForce)
    output = 0.5_dp * stateVar &
         & * (MEANBS_FUNC(currTime, currTime) + detForce) &
         & * EXP(-2.0_dp * MEANJST(currTime, currTime) &
         &       + 2.0_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 2.0_dp * lambdaHat * timeStepSize)
    CALL DET_FORCING_DEFAULT(currTime + timeStepSize, detForce)
    output = output + 0.5_dp * stateVar &
         & * (MEANBS_FUNC(currTime, currTime + timeStepSize) + detForce) &
         & * EXP(-1.0_dp * MEANJST(currTime, currTime) &
         &       - MEANJST_FUNC(currTime, currTime + timeStepSize) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime, currTime) &
         &       + 0.5_dp * COVJSTJRT_FUNC(currTime, currTime + timeStepSize, &
         &                                 currTime + timeStepSize) &
         &       + COVJSTJRT_FUNC(currTime, currTime, currTime + timeStepSize) &
         &       + lambdaHat * timeStepSize)

    ! Add the summation term
    DO i = 1, integrationCount
       CALL DET_FORCING_DEFAULT(currTime + i * integrateSize, detForce)
       output = output + 0.5_dp * stateVar &
            & * (MEANBS_FUNC(currTime, currTime + i * integrateSize) &
            &    + detForce) &
            & * EXP(-1.0_dp * MEANJST(currTime, currTime) &
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
  !> @author Jason Turner, University of Wisconsin-Madison
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
  
END SUBROUTINE TRAPEZOIDAL_ANALYTIC_SBR

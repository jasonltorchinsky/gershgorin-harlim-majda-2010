###############################################################################
# Import packages:

import numpy as np
from scipy import integrate

###############################################################################
# Helper functions:

def complex_quad(func, a, b, funcArgs):
    """
    Integrate a complex-valued function.
    """
    def real_func(x, *funcArgs):
        return np.real(func(x, *funcArgs))
    def imag_func(x, *funcArgs):
        return np.imag(func(x, *funcArgs))
    realIntegral = integrate.quad(real_func, a, b, args=funcArgs)
    imagIntegral = integrate.quad(imag_func, a, b, args=funcArgs)
    return (realIntegral[0] + 1.0j * imagIntegral[0], realIntegral[1],
            imagIntegral[1])

def complex_dblquad(func, a, b, gfun, hfun, funcArgs):
    """
    Double-integrate a complex-valued function.
    """
    def real_func(x, *funcArgs):
        return np.real(func(x, *funcArgs))
    def imag_func(x, *funcArgs):
        return np.imag(func(x, *funcArgs))
    realIntegral = integrate.dblquad(real_func, a, b, gfun, hfun, args=funcArgs)
    imagIntegral = integrate.dblquad(imag_func, a, b, gfun, hfun, args=funcArgs)
    return (realIntegral[0] + 1.0j * imagIntegral[0], realIntegral[1],
            imagIntegral[1])

###############################################################################
# Functions for calculating statistics later:

def fs(s, aF, omegaF):
    """
    f(s) from the equations.
    """

    output = aF * np.exp(1.0j * omegaF * s)

    return output

def mean_bs(s, bHat, b0, lambdaB, t0):
    """
    <b(s)> from the equations.
    """

    output = bHat + (b0 - bHat) * np.exp(lambdaB * (s - t0))

    return output

def mean_Jst(s, t, gamma0, gammaHat, dGamma, t0):
    """
    <J(s,t)> from the equations.
    """

    output = (((gamma0 - gammaHat) / dGamma)
              * (np.exp(-dGamma * (s - t0)) - np.exp(-dGamma * (t - t0)))
             )
    
    return output

def var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0):
    """
    Var(J(s,t)) from the equations.
    """

    firstTerm = ((varGamma0 / dGamma**2)
                 * (np.exp(-dGamma * (s - t0)) - np.exp(-dGamma * (t - t0)))**2
                )
    secondTerm = (((sigmaGamma**2) / (dGamma**3))
                  * (-1 + dGamma * (t - s)
                     + np.exp(-dGamma * (t + s - 2.0 * t0))
                       * (1 + np.exp(2.0 * dGamma * (s - t0))
                          - np.cosh(dGamma * (t - s))
                         )
                    )
                 )

    return firstTerm + secondTerm

def cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0):
    """
    Cov(u0, J(t0,t)) from the equations.
    """

    output = (covU0Gamma0 / dGamma) * (1.0 - np.exp(-dGamma * (t - t0)))

    return output

def cov_bs_Jst(s, t, covB0Gamma0, dGamma, lambdaB, t0):
    """
    Cov(b(s), J(s,t)) from the equations.
    """

    output = ((covB0Gamma0 / dGamma) * np.exp(lambdaB * (s - t0))
              * (np.exp(-dGamma * (s - t0)) - np.exp(-dGamma * (t - t0)))
             )

    return output

def integrand1_meanU(t, s, lambdaHat, bHat, b0, lambdaB, covB0Gamma0,
                     gamma0, varGamma0, gammaHat, dGamma, sigmaGamma, t0):
    """
    The first integrand of <u(t)>.
    """

    meanBs = mean_bs(s, bHat, b0, lambdaB, t0)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0,sigmaGamma, dGamma, t0)
    covBsJst = cov_bs_Jst(s, t, covB0Gamma0, dGamma, lambdaB, t0)
    
    frstMult = np.exp(lambdaHat * (t - s))
    scndMult = meanBs - covBsJst
    thrdMult = np.exp(-meanJst + 0.5 * varJst)
    
    return frstMult * scndMult * thrdMult

def integrand2_meanU(t, s, aF, omegaF, lambdaHat, gamma0, gammaHat,
                     dGamma, varGamma0, sigmaGamma, t0):
    """
    The second integrand of <u(t)>.
    """

    fsVal = fs(s, aF, omegaF)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    
    scndMult = np.exp(lambdaHat * (t - s))
    thrdMult = np.exp(-meanJst + 0.5 * varJst)
    
    return fsVal * scndMult * thrdMult

def cov_Jst_Jrt(s, r, t, varGamma0, dGamma, sigmaGamma, t0):
    """
    Cov(J(s,t),J(r,t)) from the equations.
    """

    if (r <= s):
        term1 = ((varGamma0 / dGamma**2)
                 * (np.exp(-dGamma * (t - t0)) - np.exp(-dGamma * (s - t0)))
                 * (np.exp(-dGamma * (t - t0)) - np.exp(-dGamma * (r - t0)))
                )
        term2 = (-1.0 - np.exp(-dGamma * (s - r)) + 2.0 * dGamma * (t - s)
                 + np.exp(-dGamma * (s + t - 2.0 * t0))
                   * (1.0 + np.exp(dGamma * (s - r))
                      + np.exp(2.0 * dGamma * (s - t0))
                      + np.exp(dGamma * (s + r - 2.0 * t0))
                      - np.exp(dGamma * (t - r))
                      - np.exp(-dGamma * (t - s)))
                )
        term2 = 0.5 * (sigmaGamma**2 / dGamma**3) * term2
    elif (s < r):
        term1 = ((varGamma0 / dGamma**2)
                 * (np.exp(-dGamma * (t - t0)) - np.exp(-dGamma * (r - t0)))
                 * (np.exp(-dGamma * (t - t0)) - np.exp(-dGamma * (s - t0)))
                )
        term2 = (-1.0 - np.exp(-dGamma * (r - s)) + 2.0 * dGamma * (t - r)
                  + np.exp(-dGamma * (r + t - 2.0 * t0))
                    * (1.0 + np.exp(dGamma * (r - s))
                       + np.exp(2.0 * dGamma * (r - t0))
                       + np.exp(dGamma * (r + s - 2.0* t0))
                       - np.exp(dGamma * (t - s))
                       - np.exp(-dGamma * (t - r)))
                )
        term2 = 0.5 * (sigmaGamma**2 / dGamma**3) * term2

    return term1 + term2

def cov_bs_Jrt(s, r, t, covB0Gamma0, dGamma, lambdaB, t0):
    """
    Cov(b(s),J(r,t)) from the equations.
    """

    output = ((covB0Gamma0 / dGamma) * np.exp(lambdaB * (s - t0))
              * (np.exp(-dGamma * (r - t0)) - np.exp(-dGamma * (t - t0)))
             )

    return output

def cov_conjbr_Jst(s, r, t, covB0Gamma0, dGamma, lambdaB, t0):
    """
    Cov(conj(b(r)),J(s,t)) from the equations.
    """

    output = ((np.conj(covB0Gamma0) / dGamma)
              * np.exp(np.conj(lambdaB) * (r - t0))
              * (np.exp(-dGamma * (s - t0)) - np.exp(-dGamma * (t - t0)))
             )

    return output

def cov_bs_br(s, r, gammaB, omegaB, varB0, sigmaB, t0):
    """
    Cov(b(s), b(r)) from the equations.
    """

    factor1 = np.exp(-gammaB * (s + r - 2.0 * t0) + 1.0j * omegaB * (s - r))
    factor2 = (varB0
               + (sigmaB**2 / (2.0 * gammaB))
                 * (np.exp(2.0 * gammaB * (min(s, r) - t0)) - 1)
              )

    return factor1 * factor2

def bvar(s, r, t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
         omegaU, b0, gammaB, omegaB, varB0, sigmaB, bHat, lambdaB, aF, omegaF,
         covB0Gamma0, t0):
    """
    bvar(s,r) from the equations.
    """

    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    meanJrt = mean_Jst(r, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    varJrt = var_Jst(r, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covJstJrt = cov_Jst_Jrt(s, r, t, varGamma0, dGamma, sigmaGamma, t0)
    covBsBr = cov_bs_br(s, r, gammaB, omegaB, varB0, sigmaB, t0)
    meanBs = mean_bs(s, bHat, b0, lambdaB, t0)
    valFs = fs(s, aF, omegaF)
    covBsJrt = cov_bs_Jrt(s, r, t, covB0Gamma0, dGamma, lambdaB, t0)
    covBsJst = cov_bs_Jrt(s, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    conjMeanBr = np.conj(mean_bs(r, bHat, b0, lambdaB, t0))
    conjValFr = np.conj(fs(r, aF, omegaF))
    covConjBrJrt = cov_conjbr_Jst(r, r, t, covB0Gamma0, dGamma, lambdaB, t0)
    covConjBrJst = cov_conjbr_Jst(s, r, t, covB0Gamma0, dGamma, lambdaB, t0)
    
    factor1 = np.exp(-meanJst - meanJrt + 0.5 * varJst + 0.5 * varJrt
                     + covJstJrt - gammaHat * (2.0 * t - s - r)
                     + 1.0j * omegaU * (s - r))
    factor2 = (covBsBr
               + (meanBs + valFs - covBsJrt - covBsJst)
                 * (conjMeanBr + conjValFr - covConjBrJrt - covConjBrJst)
              )

    return factor1 * factor2

def integrand_absc2(s, t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma, t0):
    """
    Integrand of |C|^2 from the equations.
    """

    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)

    output = np.exp(2.0 * (varJst - meanJst - gammaHat * (t - s)))

    return output

def cov_u0_Jst(s, t, covU0Gamma0, dGamma, t0):
    """
    Cov(u0, J(s,t)) from the equations.
    """

    output = ((covU0Gamma0 / dGamma)
              * (np.exp(-dGamma * (s - t0)) - np.exp(-dGamma * (t - t0)))
             )

    return output

def integrand_absaconjb(s, t, covU0B0, lambdaB, u0, covU0Gamma0, dGamma,
                         bHat, b0, aF, omegaF, covB0Gamma0, gamma0, gammaHat,
                         varGamma0, sigmaGamma, omegaU, t0):
    """
    Integrand of |Aconj(B)| from the equations.
    """

    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    covU0Jst = cov_u0_Jst(s, t, covU0Gamma0, dGamma, t0)
    conjMeanBs = np.conj(mean_bs(s, bHat, b0, lambdaB, t0))
    conjValFs = np.conj(fs(s, aF, omegaF))
    covConjBsJt0t = cov_conjbr_Jst(t0, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    covConjBsJst = cov_conjbr_Jst(s, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    meanJt0t = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    varJt0t = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covJt0tJst = cov_Jst_Jrt(t0, s, t, varGamma0, dGamma, sigmaGamma, t0)

    factor1Term1 = covU0B0 * np.exp(np.conj(lambdaB) * (s - t0))
    factor1Term2 = ((u0 - covU0Jt0t - covU0Jst)
                    * (conjMeanBs + conjValFs - covConjBsJt0t - covConjBsJst)
                   )
    factor2 = np.exp(-meanJt0t - meanJst + 0.5 * varJt0t + 0.5 * varJst
                     + covJt0tJst - gammaHat * (2.0 * t - s - t0)
                     + 1.0j * omegaU * (s - t0)
                    )

    return (factor1Term1 + factor1Term2) * factor2

def cov_bs_conjbr(s, r, covB0ConjB0, lambdaB, t0):
    """
    Cov(b(s),conj(b(r))) in equations.
    """

    return covB0ConjB0 * np.exp(lambdaB * (s + r - 2.0 * t0))

def bcovar(s, r, t, gamma0, dGamma, varGamma0, sigmaGamma, gammaHat,
           b0, gammaB, omegaB, covB0ConjB0, sigmaB, bHat, lambdaB, aF, omegaF,
           covB0Gamma0, lambdaHat, t0):
    """
    bcovar(s,r) from the equations.
    """

    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    meanJrt = mean_Jst(r, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    varJrt = var_Jst(r, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covJstJrt = cov_Jst_Jrt(s, r, t, varGamma0, dGamma, sigmaGamma, t0)
    covBsConjBr = cov_bs_conjbr(s, r, covB0ConjB0, lambdaB, t0)
    meanBs = mean_bs(s, bHat, b0, lambdaB, t0)
    valFs = fs(s, aF, omegaF)
    covBsJrt = cov_bs_Jrt(s, r, t, covB0Gamma0, dGamma, lambdaB, t0)
    covBsJst = cov_bs_Jrt(s, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    meanBr = mean_bs(r, bHat, b0, lambdaB, t0)
    valFr = fs(r, aF, omegaF)
    covBrJrt = cov_bs_Jrt(r, r, t, covB0Gamma0, dGamma, lambdaB, t0)
    covBrJst = cov_bs_Jrt(r, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    
    factor1 = np.exp(-meanJst - meanJrt + 0.5 * varJst + 0.5 * varJrt
                     + covJstJrt + lambdaHat * (2.0 * t - s - r))
    factor2 = (covBsConjBr
               + (meanBs + valFs - covBsJrt - covBsJst)
                 * (meanBr + valFr - covBrJrt - covBrJst)
              )

    return factor1 * factor2

def integrand_c2(s, t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
                 lambdaHat, t0):
    """
    Integrand of C^2 from the equations.
    """

    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)

    output = np.exp(2.0 * (varJst + lambdaHat * (t - s) - meanJst))

    return output

def integrand_ab(s, t, covU0B0, lambdaB, u0, covU0Gamma0, dGamma,
                 bHat, b0, aF, omegaF, covB0Gamma0, gamma0, gammaHat,
                 varGamma0, sigmaGamma, lambdaHat, t0):
    """
    Integrand of AB from the equations.
    """

    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    covU0Jst = cov_u0_Jst(s, t, covU0Gamma0, dGamma, t0)
    meanBs = np.conj(mean_bs(s, bHat, b0, lambdaB, t0))
    valFs = np.conj(fs(s, aF, omegaF))
    covBsJt0t = cov_bs_Jrt(s, t0, t, covB0Gamma0, dGamma, lambdaB, t0)
    covBsJst = cov_bs_Jrt(s, s, t, covB0Gamma0, dGamma, lambdaB, t0)
    meanJt0t = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    meanJst = mean_Jst(s, t, gamma0, gammaHat, dGamma, t0)
    varJt0t = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    varJst = var_Jst(s, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covJt0tJst = cov_Jst_Jrt(t0, s, t, varGamma0, dGamma, sigmaGamma, t0)

    factor1Term1 = covU0B0 * np.exp(lambdaB * (s - t0))
    factor1Term2 = ((u0 - covU0Jt0t - covU0Jst)
                    * (meanBs + valFs - covBsJt0t - covBsJst)
                   )
    factor2 = np.exp(-meanJt0t - meanJst + 0.5 * varJt0t + 0.5 * varJst
                     + covJt0tJst + lambdaHat * (2.0 * t - s - t0)
                    )

    return (factor1Term1 + factor1Term2) * factor2

###############################################################################
# Functions for calculating final statistics:

def mean_gamma_anal(t, gammaHat, gamma0, dGamma, t0):
    """
    Returns the mean of gamma at time t obtained via the analytic solution.
    """

    term2 = (gamma0 - gammaHat) * np.exp(-dGamma * (t - t0))

    return gammaHat + term2

def mean_b_anal(t, bHat, b0, lambdaB, t0):
    """
    Returns the mean of b at time t obtained via the analytic solution.
    """

    return mean_bs(t, bHat, b0, lambdaB, t0)

def var_gamma_anal(t, varGamma0, dGamma, sigmaGamma, t0):
    """
    Returns the variance of gamma at time t obtained via the analytic solution.
    """

    term1 = varGamma0 * np.exp(-2.0 * dGamma * (t - t0))
    term2 = ((sigmaGamma**2 / (2.0 * dGamma))
             * (1.0 - np.exp(-2.0 * dGamma * (t - t0)))
            )

    return term1 + term2

def var_b_anal(t, varB0, gammaB, sigmaB, t0):
    """
    Returns the variance of b at time t obtained via the analytic solution.
    """

    term1 = varB0 * np.exp(-2.0 * gammaB * (t - t0))
    term2 = ((sigmaB**2 / (2.0 * gammaB))
             * (1.0 - np.exp(-2.0 * gammaB * (t - t0)))
            )

    return term1 + term2

def cov_b_conjb_anal(t, covB0ConjB0, lambdaB, t0):
    """
    Returns the covariance of b and conj(b) at time t obtained via the analytic 
    solution.
    """

    return covB0ConjB0 * np.exp(2.0 * lambdaB * (t - t0))

def cov_b_gamma_anal(t, covB0Gamma0, lambdaB, dGamma, t0):
    """
    Returns the covariance of b and gamma at time t obtained via the analytic 
    solution.
    """

    return covB0Gamma0 * np.exp((lambdaB - dGamma) * (t - t0))

def mean_u_anal(t, u0, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
                covU0Gamma0, lambdaHat, b0, bHat, lambdaB, covB0Gamma0,
                aF, omegaF, t0):

    """
    Returns the mean of u at time t obtained via the analytic solution.
    """

    meanJst = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    term1 = (np.exp(lambdaHat * (t - t0))
             * (u0 - cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0))
             * np.exp(-meanJst + 0.5 * varJst)
            )

    args = (t0,lambdaHat,bHat,b0,lambdaB,covB0Gamma0,gamma0,varGamma0,
            gammaHat,dGamma,sigmaGamma,t0)
    term2 = complex_quad(integrand1_meanU, t0, t, funcArgs=args)[0]

    args = (t0,aF,omegaF,lambdaHat,gamma0,gammaHat,dGamma,varGamma0,
            sigmaGamma,t0)
    term3 = complex_quad(integrand2_meanU, t0, t, funcArgs=args)[0]

    return term1 + term2 + term3

def var_u_anal(t, u0, varU0, covU0Gamma0, covU0B0, dGamma, gamma0, gammaHat,
               varGamma0, sigmaGamma, omegaU, b0, gammaB, omegaB, varB0,
               sigmaB, bHat, lambdaB, aF, omegaF, covB0Gamma0, lambdaHat,
               sigmaU, t0):
    """
    Returns the variance of u at time t obtained via the analytic solution.
    """

    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    meanJt0t = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    varJt0t = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    meanAbsA2Factor1 = (u0**2 + varU0 - 4.0 * np.real(np.conj(u0) * covU0Jt0t)
                        + 4.0 * np.abs(covU0Jt0t)**2
                       )
    meanAbsA2Factor2 = np.exp(2.0 * (varJt0t - meanJt0t - gammaHat * (t - t0)))
    meanAbsA2 = meanAbsA2Factor1 * meanAbsA2Factor2

    args = (t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
            omegaU, b0, gammaB, omegaB, varB0, sigmaB, bHat, lambdaB,
            aF, omegaF, covB0Gamma0, t0)
    meanAbsB2 = complex_dblquad(bvar, t0, t, t0, t, args)[0]

    args = (t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma, t0)
    meanAbsC2 = sigmaU**2 * complex_quad(integrand_absc2, t0, t, args)[0]

    args = (t, covU0B0, lambdaB, u0, covU0Gamma0, dGamma,
            bHat, b0, aF, omegaF, covB0Gamma0, gamma0, gammaHat,
            varGamma0, sigmaGamma, omegaU, t0)
    meanAbsAConjB = complex_quad(integrand_absaconjb,t0, t, args)[0]

    meanU = mean_u_anal(t, u0, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
                        covU0Gamma0, lambdaHat, b0, bHat, lambdaB, covB0Gamma0,
                        aF, omegaF, t0)
    absMeanU2 = np.abs(meanU)**2

    output = (meanAbsA2 + meanAbsB2 + meanAbsC2
              + 2.0 * np.real(meanAbsAConjB) - absMeanU2
             )
    
    return output

def cov_u_conju_anal(t, u0, covU0ConjU0, covU0Gamma0, covU0B0, dGamma,
                     gamma0, gammaHat, varGamma0, sigmaGamma, omegaU, b0,
                     gammaB, omegaB, varB0, sigmaB, bHat, lambdaB, covB0ConjB0,
                     aF, omegaF, covB0Gamma0, lambdaHat, sigmaU, t0):
    """
    Returns the covariance of u and conj(u) at time t obtained via the analytic 
    solution.
    """

    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    meanJt0t = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    varJt0t = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    meanA2Factor1 = covU0ConjU0 + (u0 - 2.0 * covU0Jt0t)**2
    meanA2Factor2 = np.exp(2.0 * (varJt0t + lambdaHat * (t - t0) - meanJt0t))
    meanA2 = meanA2Factor1 * meanA2Factor2

    args = (t, gamma0, dGamma, varGamma0, sigmaGamma, gammaHat,
            b0, gammaB, omegaB, covB0ConjB0, sigmaB, bHat, lambdaB, aF, omegaF,
            covB0Gamma0, lambdaHat, t0)
    meanB2 = complex_dblquad(bcovar, t0, t, t0, t, args)[0]

    args = (t, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma, lambdaHat,
            t0)
    meanC2 = sigmaU**2 * complex_quad(integrand_c2, t0, t, args)[0]

    args = (t, covU0B0, lambdaB, u0, covU0Gamma0, dGamma,
                 bHat, b0, aF, omegaF, covB0Gamma0, gamma0, gammaHat,
                 varGamma0, sigmaGamma, lambdaHat, t0)
    meanAB = complex_quad(integrand_ab, t0, t, args)[0]

    meanU = mean_u_anal(t, u0, gamma0, gammaHat, dGamma, varGamma0, sigmaGamma,
                        covU0Gamma0, lambdaHat, b0, bHat, lambdaB, covB0Gamma0,
                        aF, omegaF, t0)

    
    return meanA2 + meanB2 + meanC2 + 2.0 * meanAB - meanU**2

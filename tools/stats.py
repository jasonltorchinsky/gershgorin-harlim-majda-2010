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

def cov_b_conjgb_anal(t, covB0ConjgB0, lambdaB, t0):
    """
    Returns the covariance of b and conjg(b) at time t obtained via the analytic 
    solution.
    """

    return covB0ConjgB0 * np.exp(2.0 * lambdaB * (t - t0))

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

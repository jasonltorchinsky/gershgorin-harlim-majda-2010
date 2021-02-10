"""
This script is based on the example script "NCL_color_1.py" provided by 
GeoCAT-Examples.

The original author and date of publication is unknown.

Modifying Author: Jason Torchinsky
Modifying Date: Winter 2021
"""

###############################################################################
# Import packages:

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
from scipy import integrate

import geocat.datafiles as gdf
from geocat.viz import cmaps as gvcmaps
import geocat.viz.util as gvutil

###############################################################################
# Plotting functions:

def stat_plot(num, anal, time, stat, varName1, varName2=""):
    """
    Creates a plot of the statistic of a variable, given an array of numerical 
    means, analytic statistic, the times of the data points, and the variable 
    name.
    """

    fig = plt.figure(figsize=(10, 8))

    # Generate axes 
    ax = plt.axes()

    ax.set_xscale("linear")
    dt = time[1] - time[0] # (! ASSUME CONSTANT TIME-STEP SIZE !)
    t0 = time[0]
    tf = time[len(time)-1]
    xTicks = np.linspace(t0,tf+dt,num=10)
    ax.set_xticks(xTicks)
    
    ax.set_yscale("linear")


    # Style for the line for the numerical mean line
    kwargsNum = dict(
                     color='k',
                     label = "Numerical"
                    )

    # Style for the line for the analytic mean line
    kwargsAnal = dict(
                      color='r',
                      label="Analytic"
                     )
    

    # Generate the plot.
    plt.plot(time, num, **kwargsNum)
    plt.plot(time, anal, **kwargsAnal)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.minorticks_off()

    # Set the title, labels, and legend
    if stat=="mean":
        titleStr = "Mean of " + varName1
        fileName = stat + varName1 + ".png"
    elif stat=="var":
        titleStr = "Variance of " + varName1
        fileName = stat + varName1 + ".png"
    elif stat=="cov":
        titleStr = "Covariance of " + varName1 + " and " + varName2
        fileName = stat + varName1 + varName2 + ".png"

    plt.title(titleStr,
              {'fontsize': 22},
              loc='center',
              pad=30.0)
    ax.set_xlabel("Time",
                  {'fontsize': 18})
    plt.legend(loc='lower right',
               fontsize=18)

    # Save the plot
    plt.savefig(fileName,dpi=1000)

###############################################################################
# Functions for calculating later statistics:

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
# Parameters for reading in data (! UPDATE FOR EACH SIMULATION !):

nProcs = 16
nRunsPerProc = 300
nRuns = nProcs * nRunsPerProc
mainFileStr = "../build/out"

###############################################################################
# Simulation parameters (! UPDATE FOR EACH SIMULATION !):

t0 = 0.0
dt = 1.0

gammaHat = 1.5
dGamma = 0.015
sigmaGamma = 0.7745
gamma0 = 0.0
varGamma0 = 0.0

bHat = 0.0
gammaB = 0.15
omegaB = 1.78
sigmaB = 0.7745
lambdaB = complex(-gammaB, omegaB)
b0 = 0.0 + 0.0j
varB0 = 0.0 + 0.0j
covB0ConjgB0 = 0.0 + 0.0j

covB0Gamma0 = 0.0 + 0.0j

omegaU = 1.78
sigmaU = 0.1549
lambdaHat = complex(-gammaHat, omegaU)
u0 = 0.0 + 0.0j

covU0Gamma0 = 0.0 + 0.0j

aF = 0.0
omegaF = 0.0


###############################################################################
# Read in data, using the first file to get the time of the time-steps:

print("Reading output files...")

# Open a netCDF data file using xarray default engine and load the data into
# xarray.
with xr.open_dataset(mainFileStr + "000.nc") as ds:
    time = np.asarray(ds.run00000001)[:,1]
    reU = np.zeros(shape=(nRuns, time.size))
    imU = np.zeros_like(reU)
    reB = np.zeros_like(reU)
    imB = np.zeros_like(reU)
    gamma = np.zeros_like(reU)
    
for proc in range(0,nProcs):
    procStr = '{0:03d}'.format(proc)
    fileName = mainFileStr + procStr + ".nc"
    with xr.open_dataset(fileName) as ds:
        runStrs = list(ds.keys())
        for runNum in range(0,nRunsPerProc):
            varIdx = proc * nRunsPerProc + runNum
            runStr = runStrs[runNum]
            reU[varIdx,:] = ds[runStr][:,2]
            imU[varIdx,:] = ds[runStr][:,3]
            reB[varIdx,:] = ds[runStr][:,4]
            imB[varIdx,:] = ds[runStr][:,5]
            gamma[varIdx,:] = ds[runStr][:,6]

###############################################################################
# Get statistics for the numerical results:

print("Calculating statistics for numerical results...")

# Mean of gamma:

meanGamma = np.mean(gamma, axis=0)

# Mean of b:

meanReB = np.mean(reB, axis=0)
meanImB = np.mean(imB, axis=0)

# Variance of gamma

varGamma = np.var(gamma, axis=0)

# Variance of b

varReB = np.var(reB, axis=0)
varImB = np.var(imB, axis=0)
# To get the covariance, we need to take the covariance at each time and put
# them into an array.
covReBImB = np.zeros_like(time)
for t in range(0,len(time)):
    covReBImB[t] = np.cov(reB[:,t], imB[:,t], rowvar=False)[0,1]

# Covariance of gamma and b

covReBGamma = np.zeros_like(time)
for t in range(0,len(time)):
    covReBGamma[t] = np.cov(reB[:,t], gamma[:,t], rowvar=False)[0,1]
covImBGamma = np.zeros_like(time)
for t in range(0,len(time)):
    covImBGamma[t] = np.cov(imB[:,t], gamma[:,t], rowvar=False)[0,1]

# Mean of u

meanReU = np.mean(reU, axis=0)
meanImU = np.mean(imU, axis=0)

###############################################################################
# Get statistics from analytic expressions:

print("Calculating analytic statistics...")

# Mean of gamma:

meanGammaAnal = (gammaHat
                 + (gamma0 - gammaHat) * np.exp(-dGamma * (time - t0))
                )

# Mean of b:

meanBAnal = (bHat
             + (b0 - bHat) * np.exp(lambdaB * (time - t0))
            )

meanReBAnal = np.real(meanBAnal)
meanImBAnal = np.imag(meanBAnal)

# Variance of gamma:

varGammaAnal = (varGamma0 * np.exp(-2.0 * dGamma * (time - t0))
                + (sigmaGamma**2 / (2.0 * dGamma))
                  * (1.0 - np.exp(-2.0 * dGamma * (time - t0)))
               )

# Variance of b:

varBAnal = (varB0 * np.exp(-2.0 * gammaB * (time - t0))
            + (sigmaB**2 / (2.0 * gammaB))
              * (1.0 - np.exp(-2.0 * gammaB * (time - t0)))
           )
covBConjgBAnal = covB0ConjgB0 * np.exp(2.0 * lambdaB * (time - t0))

varReBAnal = 0.5 * np.real(varBAnal + covBConjgBAnal)
varImBAnal = 0.5 * np.real(varBAnal - covBConjgBAnal)
covReBImBAnal = -0.5 * np.imag(varBAnal - covBConjgBAnal)

# Covariance of b and gamma

covBGammaAnal = covB0Gamma0 * np.exp((lambdaB - dGamma) * (time - t0))

covReBGammaAnal = np.real(covBGammaAnal)
covImBGammaAnal = np.imag(covBGammaAnal)

# Mean of u

meanUAnal1 = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    meanJst = mean_Jst(t0, t, gamma0, gammaHat, dGamma, t0)
    varJst = var_Jst(t0, t, gamma0, varGamma0, sigmaGamma, dGamma, t0)
    covU0Jt0t = cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0)
    meanUAnal1[tIdx] = (np.exp(lambdaHat * (t - t0))
                        * (u0 - cov_u0_Jt0t(t, covU0Gamma0, dGamma, t0))
                        * np.exp(-meanJst + 0.5 * varJst)
                       )

meanUAnal2 = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    args = (t0,lambdaHat,bHat,b0,lambdaB,covB0Gamma0,gamma0,varGamma0,
            gammaHat,dGamma,sigmaGamma,t0)
    meanUAnal2[tIdx] = complex_quad(integrand1_meanU, t0, t, funcArgs=args)[0]

meanUAnal3 = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    args = (t0,aF,omegaF,lambdaHat,gamma0,gammaHat,dGamma,varGamma0,
            sigmaGamma,t0)
    meanUAnal3[tIdx] = complex_quad(integrand2_meanU, t0, t, funcArgs=args)[0]

meanUAnal = meanUAnal1 + meanUAnal2 + meanUAnal3

meanReUAnal = np.real(meanUAnal)
meanImUAnal = np.imag(meanUAnal)

###############################################################################
# Plot numerical versus analytic statistics:

print("Generating plots...")

# Mean of gamma
stat_plot(meanGamma, meanGammaAnal, time, "mean", "\u03B3")

# Mean of b
stat_plot(meanReB, meanReBAnal, time, "mean", "Re(b)")
stat_plot(meanImB, meanImBAnal, time, "mean", "Im(b)")

# Variance of gamma
stat_plot(varGamma, varGammaAnal, time, "var", "\u03B3")

# Variance of b
stat_plot(varReB, varReBAnal, time, "var", "Re(b)")
stat_plot(varImB, varImBAnal, time,  "var","Im(b)")
stat_plot(covReBImB, covReBImBAnal, time, "cov", "Re(b)", "Im(b)")

# Covariance of b and gamma
stat_plot(covReBGamma, covReBGammaAnal, time, "cov", "Re(b)", "\u03B3")
stat_plot(covImBGamma, covImBGammaAnal, time, "cov", "Im(b)", "\u03B3")

# Mean of u
stat_plot(meanReU, meanReUAnal, time, "mean", "Re(u)")
stat_plot(meanImU, meanImUAnal, time, "mean", "Im(u)")

print("Done!")

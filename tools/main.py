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

from plot import stat_plot
from stats import *

    
###############################################################################
# Parameters for reading in data (! UPDATE FOR EACH SIMULATION !):

nProcs = 16
nRunsPerProc = 500
nRuns = nProcs * nRunsPerProc
mainFileStr = "../build/out"

###############################################################################
# Simulation parameters (! UPDATE FOR EACH SIMULATION !):

t0 = 0.0
dt = 0.001

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
covB0ConjB0 = 0.0 + 0.0j

covB0Gamma0 = 0.0 + 0.0j

omegaU = 1.78
sigmaU = 0.1549
lambdaHat = complex(-gammaHat, omegaU)
u0 = 0.0 + 0.0j
varU0 = 0.0 + 0.0j
covU0ConjU0 = 0.0 + 0.0j

covU0Gamma0 = 0.0 + 0.0j
covU0B0 = 0.0 + 0.0j

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
covReBImB = np.empty_like(time)
for tIdx in range(0,len(time)):
    covReBImB[tIdx] = np.cov(reB[:,tIdx], imB[:,tIdx], rowvar=False)[0,1]

# Covariance of gamma and b

covReBGamma = np.empty_like(time)
for tIdx in range(0,len(time)):
    covReBGamma[tIdx] = np.cov(reB[:,tIdx], gamma[:,tIdx], rowvar=False)[0,1]
covImBGamma = np.empty_like(time)
for tIdx in range(0,len(time)):
    covImBGamma[tIdx] = np.cov(imB[:,tIdx], gamma[:,tIdx], rowvar=False)[0,1]

# Mean of u

meanReU = np.mean(reU, axis=0)
meanImU = np.mean(imU, axis=0)

# Variance of u

varReU = np.var(reU, axis=0)
varImU = np.var(imU, axis=0)
covReUImU = np.empty_like(time)
for tIdx in range(0,len(time)):
    covReUImU[tIdx] = np.cov(reU[:,tIdx], imU[:,tIdx], rowvar=False)[0,1]

# Plot all u trajectories because why not
"""
fig = plt.figure(figsize=(10, 8))

# Generate axes 
ax = plt.axes()

ax.set_xscale("linear")
t0 = time[0]
tf = time[len(time)-1]
xTicks = np.linspace(t0,tf,num=10)
ax.set_xticks(xTicks)

ax.set_yscale("linear")

for runIdx in range(0, len(reU[:,0])):
    plt.plot(time, imU[runIdx,:])

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.minorticks_off()
# Save the plot
plt.savefig("plots/imUTraj.png",dpi=300)
"""

###############################################################################
# Get statistics from analytic expressions:

print("Calculating analytic statistics...")

# Mean of gamma:

meanGammaAnal = np.empty_like(time, dtype=np.float64)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    meanGammaAnal[tIdx] = mean_gamma_anal(t, gammaHat, gamma0, dGamma, t0)
    
print("  Calculated mean \u03B3.")

# Mean of b:

meanBAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    meanBAnal[tIdx] = mean_b_anal(t, bHat, b0, lambdaB, t0)

meanReBAnal = np.real(meanBAnal)
meanImBAnal = np.imag(meanBAnal)

print("  Calculated mean b.")

# Variance of gamma:

varGammaAnal = np.empty_like(time, dtype=np.float64)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    varGammaAnal[tIdx] = var_gamma_anal(t, varGamma0, dGamma, sigmaGamma, t0)
    
print("  Calculated variance \u03B3.")

# Variance of b:

varBAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    varBAnal[tIdx] = var_b_anal(t, varB0, gammaB, sigmaB, t0)
covBConjBAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    covBConjBAnal[tIdx] = cov_b_conjb_anal(t, covB0ConjB0, lambdaB, t0)

varReBAnal = 0.5 * np.real(varBAnal + covBConjBAnal)
varImBAnal = 0.5 * np.real(varBAnal - covBConjBAnal)
covReBImBAnal = -0.5 * np.imag(varBAnal - covBConjBAnal)

print("  Calculated variance b and covariance b, conj(b).")

# Covariance of b and gamma

covBGammaAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    covBGammaAnal[tIdx] = cov_b_gamma_anal(t, covB0Gamma0, lambdaB, dGamma, t0)

covReBGammaAnal = np.real(covBGammaAnal)
covImBGammaAnal = np.imag(covBGammaAnal)

print("  Calculated covariance b and \u03B3.")

# Mean of u

meanUAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    meanUAnal[tIdx] = mean_u_anal(t, u0, gamma0, gammaHat, dGamma, varGamma0,
                                  sigmaGamma, covU0Gamma0, lambdaHat, b0,
                                  bHat, lambdaB, covB0Gamma0, aF, omegaF, t0)

meanReUAnal = np.real(meanUAnal)
meanImUAnal = np.imag(meanUAnal)

print("  Calculated mean u.")

# Variance of u

varUAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    varUAnal[tIdx] = var_u_anal(t, u0, varU0, covU0Gamma0, covU0B0, dGamma, gamma0, gammaHat,
                                varGamma0, sigmaGamma, omegaU, b0, gammaB, omegaB, varB0,
                                sigmaB, bHat, lambdaB, aF, omegaF, covB0Gamma0, lambdaHat,
                                sigmaU, t0)
    
covUConjUAnal = np.empty_like(time, dtype=np.complex128)
for tIdx in range(0,len(time)):
    t = time[tIdx]
    covUConjUAnal[tIdx] = cov_u_conju_anal(t, u0, covU0ConjU0, covU0Gamma0, covU0B0, dGamma,
                                           gamma0, gammaHat, varGamma0, sigmaGamma, omegaU, b0,
                                           gammaB, omegaB, varB0, sigmaB, bHat, lambdaB, covB0ConjB0,
                                           aF, omegaF, covB0Gamma0, lambdaHat, sigmaU, t0)

varReUAnal = 0.5 * np.real(varUAnal + covUConjUAnal)
varImUAnal = 0.5 * np.real(varUAnal - covUConjUAnal)
covReUImUAnal = -0.5 * np.imag(varUAnal - covUConjUAnal)

print("  Calculated variance u and covariance u, conj(u).")

###############################################################################
# Plot numerical versus analytic statistics:

print("Generating plots...")

# Mean of gamma
stat_plot(meanGamma, meanGammaAnal, time, "mean", "\u03B3")
print(" Generated plot for mean \u03B3.")

# Mean of b
stat_plot(meanReB, meanReBAnal, time, "mean", "Re(b)")
stat_plot(meanImB, meanImBAnal, time, "mean", "Im(b)")
print(" Generated plots for mean b.")

# Variance of gamma
stat_plot(varGamma, varGammaAnal, time, "var", "\u03B3")
print(" Generated plot for variance of \u03B3.")

# Variance of b
stat_plot(varReB, varReBAnal, time, "var", "Re(b)")
stat_plot(varImB, varImBAnal, time, "var","Im(b)")
stat_plot(covReBImB, covReBImBAnal, time, "cov", "Re(b)", "Im(b)")
print(" Generated plots for variance of b.")

# Covariance of b and gamma
stat_plot(covReBGamma, covReBGammaAnal, time, "cov", "Re(b)", "\u03B3")
stat_plot(covImBGamma, covImBGammaAnal, time, "cov", "Im(b)", "\u03B3")
print(" Generated plots for covariance of b and \u03B3.")

# Mean of u
stat_plot(meanReU, meanReUAnal, time, "mean", "Re(u)")
stat_plot(meanImU, meanImUAnal, time, "mean", "Im(u)")
print(" Generated plots for mean u.")

# Variance of u
stat_plot(varReU, np.max(varReU) / np.max(varReUAnal) * varReUAnal, time, "var", "Re(u)")
stat_plot(varImU, np.max(varImU) / np.max(varImUAnal) *varImUAnal, time, "var","Im(u)")
stat_plot(covReUImU, np.max(covReUImU) / np.max(covReUImUAnal) * covReUImUAnal, time, "cov", "Re(u)", "Im(u)")
print(" Generated plots for variance of u.")

print("Done!")

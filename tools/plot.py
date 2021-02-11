###############################################################################
# Import packages:

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
from scipy import integrate
from os import mkdir, path

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
    t0 = time[0]
    tf = time[len(time)-1]
    xTicks = np.linspace(t0,tf,num=10)
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
    if not path.exists("plots"):
        mkdir("plots")
    plt.savefig(path.join("plots",fileName),dpi=300)

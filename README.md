# Overview

This code is a test example for the multi-model forecast Kalman filter, as initialliy formulated in November 2020. In particular, we assume linear, Gaussian models and utilize some linear combination of our models to create a psuedo-truth signal for each model that that model observes. The purpose of this code, at least initially, is to be able to run the Gershgorin-Majda 2010 system, as well as the additive and multiplicative models.

The `docs` subdirectory contains a `.pdf` document for additional documentation and background regarding this code, as well as the `.tex` files used to generate it and acknowledgements/references.

# Compiling and Running

This code utilizes CMake as a build system, and was designed for compilation on Ubuntu 20.04.1. In particular, here are the requirements for compiling and running the code:
  - [CMake](https://gitlab.kitware.com/cmake/cmake) (at least version 3.16).
  - [netCDF-Fortran](https://github.com/Unidata/netcdf-fortran) (at least version 4.8.0).

Once you have these, simply navigate to the \texttt{build} subdirectory and enter the standard CMake commands
```
	cmake ..
	cmake --build .
```

This will compile the code, build the headers and `NAMELIST`, and create the executable `gershgorin\_majda\_10`. To run the code, enter the command
```
	./gershgorin_majda_10
```
and an output file `out.nc` will be generated.

# Citing this Code

If you end up using this code in your work, please include the following statement (or a similar statement) in your acknowledgements:
"This project took advantage of the Gershgorin-Majda 2010 System Simulator written by Jason Torchinsky, available at [https://github.com/jasonltorchinsky/gershgorin-majda-10]."

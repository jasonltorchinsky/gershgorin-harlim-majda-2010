# Overview

This code is an implementation of the Gershgorin-Harlim-Majda 2010 system. It inlcudes a Forward Euler, Backward Euler, adaptive Forward Euler, adaptive Backward Euler, and an analytic statistics solver (using the trapezoidal method for required integrals). It was originally intended for use in the multi-model communication project (working title).

The `docs\build` subdirectory contains a `.pdf` document for additional documentation and background regarding this code, as well as the `.tex` files used to generate it and acknowledgements/references.

# Compiling and Running

This code utilizes CMake as a build system, and was written on Ubuntu 20.04.1. In particular, here are the requirements for compiling and running the code:
  - [CMake](https://gitlab.kitware.com/cmake/cmake) (at least version 3.16).
  - [LAPACK](https://github.com/Reference-LAPACK/lapack) (at least version 3.9.0).
  - [MPI](https://www.mpi-forum.org/docs/) (at least version 4.0).
  - [netCDF-Fortran](https://github.com/Unidata/netcdf-fortran) (at least version 4.8.0).

The minimum required versions are not "hard" minimums; older versions may work, but the code was developed using these versions. Once you have these, simply navigate to the `build` subdirectory and enter the standard CMake commands
```
	mkdir build
	cd build
	cmake ..
	cmake --build .
```

This will compile the code, build the headers and `NAMELIST`, and create the executable `gershgorin_majda_10`. To run the code, enter the command
```
	mpirun -np X gershgorin_harlim_majda_10
```
(where `X` is the number of processors) and several output files `outXXX.nc` will be generated.

# Citing this Code

If you end up using this code in your work, please include the following statement (or a similar statement) in your acknowledgements:
"This project took advantage of the Gershgorin-Harlim-Majda 2010 System Simulator written by Jason Torchinsky, available at [https://github.com/jasonltorchinsky/gershgorin-harlim-majda-2010]."

# Acknowledgements

This project took advantage of CMake software developed by Kitware, LAPACK software, MPI software, and netCDF software developed by UCAR/Unidata [http://doi.org/10.5065/D6H70CW6]. This work was funded by the Department of Energy's Computational Science Graduate Fellowship, under grant number DE-SC0020347. We would also like to acknowledge Samuel Stechmann for advising the development of this simulator, and Nelia Mann for insight into debugging some features.

Please refer to the `main.pdf` document in the `docs\build` subdirectory for further acknowledgements and citations.



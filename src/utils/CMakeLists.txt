# CMakeLists.txt for src/utils subdirectory.
# Author: Jason Turner
# Date: 2012-12-XX

# Defines netcdf_utils as a library
add_library(utils
	    utils.f90
	    netcdf_utils.f90
	    rng_utils.f90
           )

target_include_directories(utils
			   PUBLIC $ENV{INCLUDE_DIR}
			  )

# Find the netcdf library and link it to here
find_package(netCDF-Fortran REQUIRED)

target_link_libraries(utils
		      PRIVATE netCDF::netcdff
		     )

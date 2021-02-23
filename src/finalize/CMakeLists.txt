# CMakeLists.txt for src/finalize subdirectory.
# Author: Jason Torchinsky
# Date: 2012-12-XX

# Defines time_stepper as a library
add_library(finalize
	    finalize.f90
	   )

# Find the netCDF library
find_package(netCDF-Fortran REQUIRED)

# Add libraries to link to the module
target_link_libraries(finalize
		      PRIVATE netCDF::netcdff
		      PUBLIC initialize
		      PUBLIC utils
		     )


# Add directories the module needs to include
target_include_directories(finalize
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{INITIALIZE_DIR}
			   PUBLIC $ENV{UTILS_DIR}
			  )

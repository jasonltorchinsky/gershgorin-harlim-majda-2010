# CMakeLists.txt for src/initialize subdirectory.
# Author: Jason Turner
# Date: 2012-12-XX

# Defines initialize as a library
add_library(initialize
	    initialize.f90
	   )

# Find the netCDF library
find_package(netCDF-Fortran REQUIRED)

# Add libraries to link to the module
target_link_libraries(initialize
		      PRIVATE netCDF::netcdff
		      PUBLIC utils
		     )

# Add directories the module needs to include
target_include_directories(initialize
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{UTILS_DIR}
			  )

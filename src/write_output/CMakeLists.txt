# CMakeLists.txt for the src/write_output subdirectory
# Author: Jason Torchinsky
# Date: 2020-12-XX

# Defines output as a library
add_library(write_output
	    write_output.f90
	   )

# Find the netCDF library
find_package(netCDF-Fortran REQUIRED)

# Add libraries to link to the module
target_link_libraries(write_output
		      PRIVATE netCDF::netcdff
		      PUBLIC initialize
		      PUBLIC utils
		     )

# Add directories the module needs to include
target_include_directories(write_output
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{INITIALIZE_DIR}
			   PUBLIC $ENV{UTILS_DIR}
			  )
			  

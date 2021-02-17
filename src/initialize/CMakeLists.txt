# CMakeLists.txt for src/initialize subdirectory.
# Author: Jason Turner
# Date: 2012-12-XX

# Defines initialize as a library
add_library(initialize
	    initialize.f90
	   )

# Find the needed libraries.
find_package(netCDF-Fortran REQUIRED) # netCDF-Fortran
find_package(MPI REQUIRED) # MPI

# Add libraries to link to the module
target_link_libraries(initialize
		      PRIVATE netCDF::netcdff
		      PRIVATE ${MPI_Fortran_LIBRARIES}
		      PUBLIC utils
		     )

# Add directories the module needs to include
target_include_directories(initialize
			   PUBLIC ${MPI_Fortran_INCLUDE_PATH}
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{UTILS_DIR}
			  )

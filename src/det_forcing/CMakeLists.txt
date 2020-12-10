# CMakeLists.txt for src/det_forcing subdirectory.
# Author: Jason Turner
# Date: 2012-12-XX

# Defines time_stepper as a library
add_library(det_forcing
	    det_forcing.f90
	   )

# Add libraries to link to the module
target_link_libraries(det_forcing
		      PUBLIC initialize
		     )


# Add directories the module needs to include
target_include_directories(det_forcing
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{INITIALIZE_DIR}
			  )

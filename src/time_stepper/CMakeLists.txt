# CMakeLists.txt for src/time_stepper subdirectory.
# Author: Jason Torchinsky
# Date: 2012-12-XX

# Defines time_stepper as a library
add_library(time_stepper
	    time_stepper.f90
	    forward_euler.f90
	    backward_euler.f90
	    adaptive_forward_euler.f90
	    adaptive_backward_euler.f90
	    trapezoidal_analytic.f90
	   )

# Add libraries to link to the module
target_link_libraries(time_stepper
		      PUBLIC initialize
		      PUBLIC write_output
		      PUBLIC det_forcing
		      PUBLIC utils
		      PUBLIC lapack
		      PUBLIC blas
		     )


# Add directories the module needs to include
target_include_directories(time_stepper
			   PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{INITIALIZE_DIR}
			   PUBLIC $ENV{WRITE_OUTPUT_DIR}
			   PUBLIC $ENV{DET_FORCING_DIR}
			   PUBLIC $ENV{UTILS_DIR}
			  )

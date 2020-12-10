# CMakeLists.txt for src directory
# Author: Jason Turner
# Date: 2020-12-XX

# Define variables for directories for use throughout the project.
set(ENV{INCLUDE_DIR} ${CMAKE_BINARY_DIR}/inc)
set(ENV{SOURCE_DIR} ${CMAKE_BINARY_DIR}/src)
set(ENV{INITIALIZE_DIR} $ENV{SOURCE_DIR}/initialize)
set(ENV{WRITE_OUTPUT_DIR} $ENV{SOURCE_DIR}/write_output)
set(ENV{TIME_STEPPER_DIR} $ENV{SOURCE_DIR}/time_stepper)
set(ENV{DET_FORCING_DIR} $ENV{SOURCE_DIR}/det_forcing)
set(ENV{FINALIZE_DIR} $ENV{SOURCE_DIR}/finalize)
set(ENV{UTILS_DIR} $ENV{SOURCE_DIR}/utils)

# Add subdirectories containing libraries/modules
add_subdirectory(initialize)
add_subdirectory(input)
add_subdirectory(write_output)
add_subdirectory(time_stepper)
add_subdirectory(finalize)
add_subdirectory(det_forcing)
add_subdirectory(utils)

# Set the version information file and kinds file
configure_file(ver.h.in $ENV{INCLUDE_DIR}/ver.h)
configure_file(kinds.h.in $ENV{INCLUDE_DIR}/kinds.h)

# Define the project's executable
add_executable(gershgorin_majda_10)

# Add source files to the executable
target_sources(gershgorin_majda_10
               PUBLIC main.f90
	      )

# Add libraries to link to the executable
target_link_libraries(gershgorin_majda_10
                      PUBLIC initialize
		      PUBLIC write_output
		      PUBLIC time_stepper
		      PUBLIC finalize
		     )

# Add directories to include with the executable
target_include_directories(gershgorin_majda_10
                           PUBLIC $ENV{INCLUDE_DIR}
			   PUBLIC $ENV{INITIALIZE_DIR}
			   PUBLIC $ENV{WRITE_OUTPUT_DIR}
			   PUBLIC $ENV{TIME_STEPPER_DIR}
			   PUBLIC $ENV{FINALIZE_DIR}
			  )



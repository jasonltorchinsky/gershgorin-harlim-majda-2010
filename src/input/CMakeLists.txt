# CMakeLists.txt for src/in subdirectory
# Author: Jason Turner
# Date: 2020-12-XX

# Copy over the namelist file (including the default parameters)
configure_file(NAMELIST.in ${CMAKE_BINARY_DIR}/NAMELIST COPYONLY)

# Code to simulate the Gershgorin-Harlim-Majda 2010 system.
# Author: Jason Torchinsky
# Date: 2021-02-XX

cmake_minimum_required(VERSION 3.16)

# Set the project name and version
enable_language(Fortran)
project(gershgorin_harlim_majda_10
        VERSION 1.0
        DESCRIPTION "Simulates the Gershgorin-Harlim-Majda 2010 system, as well as the additive and multiplicative models."
       )
       
# Set the build type, typically will be debug for the developpment process
set(CMAKE_BUILD_TYPE Debug)

# Specify the executable directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})


# Add the subdirectories containing source code
add_subdirectory(src)

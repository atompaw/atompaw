# Cmake config file for ATOMPAW
# Authors: M. Torrent, P. Kestener

cmake_minimum_required(VERSION 3.18)

# Suppress developer warnings
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE INTERNAL "" FORCE)

# Retrieve ATOMPAW version
file(STRINGS "VERSION" FILE_VERSION)
string(REGEX MATCH " *Version *([0-9][0-9.]*)" match "${FILE_VERSION}")
if(CMAKE_MATCH_1)
    set(ATP_VERSION "${CMAKE_MATCH_1}")
else()
    set(ATP_VERSION "undefined")
endif()

# Define project
project(ATOMPAW
  LANGUAGES C Fortran
  VERSION ${ATP_VERSION}
  HOMEPAGE_URL https://github.com/atompaw/atompaw
  DESCRIPTION "Generator of projectors and basis functions for electronic structure calculations using the Projector Augmented-Wave (PAW) method")

# Add subdirectory with CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Set build type (Release, Debug, ...)
include(set_default_build_type)
#include(prevent_build_in_place)

# Gnu compatibility,
include(GNUInstallDirs)

# Export compiles commands to json (for LSP)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

# Look for PkgConfig
find_package(PkgConfig)

# Perform some operations related to compilers
include(compilers)

# Dependency 1: BLAS/LAPACK
include(get_linalg)

# Dependency 2: LIBXC
include(get_libxc)

# Generate config.h and make sure it is found (autotools compatibility)
include(generate_config_h)
add_compile_definitions(HAVE_CONFIG_H)
add_compile_options(-I${CMAKE_BINARY_DIR})

# Add source files
add_subdirectory(src)

# Installation of docs and scripts
install(FILES
  doc/atompaw-usersguide.pdf
  doc/ATOMPAW_Explore_Userguide.pdf
  DESTINATION share/doc/atompaw/doc)
install(FILES
  doc/exploresetup/ParameterExplore.cpp
  doc/exploresetup/ParameterExploreReadMe.txt
  DESTINATION share/doc/atompaw/doc/exploresetup)
install(FILES
  doc/scripts/explorelogderiv
  doc/scripts/explorewfn
  doc/scripts/plotdensity
  doc/scripts/plotlogderiv
  doc/scripts/plotpotential
  doc/scripts/plotwfn
  DESTINATION bin/scripts)

# Create dist target
include(create_dist)

# Final printing
message(STATUS "")
message(STATUS "=================================================")
message(STATUS " ${PROJECT_NAME} version ${PROJECT_VERSION}")
message(STATUS " Build configuration")
message(STATUS "=================================================")
message(STATUS "CMake version ${CMAKE_VERSION}")
message(STATUS "CMake install prefix   : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "CMake system processor : ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "CMake system name (OS) : ${CMAKE_SYSTEM_NAME}")
message(STATUS "==== BLAS/LAPACK config:")
if (BLAS_LAPACK_OK)
  message(STATUS "  Lapack libraries   : ${LAPACK_LIBRARIES}")
else()
  message(STATUS "  BLAS/LAPACK not working!")
endif()
if (USE_LIBXC)
  message(STATUS "==== LIBXC config:")
  if (LIBXC_OK)
    message(STATUS "  Libxc version      : ${Libxc_VERSION}")
    message(STATUS "  Libxc Libraries    : ${Libxc_LIBRARIES}")
    message(STATUS "  Libxc Includes     : ${Libxc_INCLUDE_DIRS}")
  else()
    message(STATUS "  Libxc not found!")
  endif()
endif()
message(STATUS "")

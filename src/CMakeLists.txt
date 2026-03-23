# Cmake file for ATOMPAW/src
# Authors: M. Torrent

# Get source files from Makefile.am
include(get_sources_from_makefile)
get_sources_from_makefile(${CMAKE_SOURCE_DIR}/src/Makefile.am atp_srcs SOURCES)

# Dynamically build pkginfo.f90
configure_file(${CMAKE_SOURCE_DIR}/cmake/pkginfo.f90.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/pkginfo.f90)

# Define a common library for all executables
add_library(atompaw_lib "")
target_sources(atompaw_lib PRIVATE ${SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/pkginfo.f90)
target_link_libraries(atompaw_lib PRIVATE BLAS::BLAS LAPACK::LAPACK)
if (LIBXC_OK) 
  target_link_libraries(atompaw_lib PRIVATE Libxc::xc)
endif()

# Define atompaw executable
add_executable(atompaw atompaw_prog.F90)
target_link_libraries(atompaw PRIVATE atompaw_lib)
install(TARGETS atompaw DESTINATION bin)

# Define atompaw_aeonly executable
add_executable(atompaw_aeonly atompaw_aeonly.F90)
target_link_libraries(atompaw_aeonly PRIVATE atompaw_lib)
install(TARGETS atompaw_aeonly DESTINATION bin)
